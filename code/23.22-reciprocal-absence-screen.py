#!/usr/bin/env python3
"""
23.22-reciprocal-absence-screen.py  (Phase 4, step 7, of code/23-next-phase-research-plan.md)

RECIPROCAL GENOMIC-ABSENCE SCREEN FOR THE BRAKER-NEW GENES. Follow-up to 23.19, which classes each
BRAKER gene against that ecotype's OWN Liftoff annotation (shared / liftoff_locus / new) and against
the SaNama reference proteins (ref_homolog / ref_paralog / no_ref_hit). The `new` + `no_ref_hit`
genes are reported as novel candidates (lean 1,137 hq, siscowet 1,262 hq).

What 23.19 has NOT asked, and this script does: **is the candidate actually absent from the other
ecotype's genome?**

A gene called in lean and missing from siscowet's ANNOTATION is not an ecotype-specific gene. It is
equally consistent with BRAKER simply failing to predict it over there -- a real risk when the two
assemblies were annotated in separate runs with separate repeat libraries. The only defensible claim
is GENOMIC absence, so each candidate protein is mapped back to the other ecotype's purged assembly
with miniprot (splice-aware protein-to-genome) and kept only when no credible locus exists there.

This is the per-gene check that 23.19's own closing caveat asks for before ecotype differences in
gene content are interpreted, and it applies the same >=2-independent-lines standard that Phase 3
pass 2 applied to the methylation calls.

Classification of each candidate (protein -> other ecotype's assembly):
  present_in_other   best hit identity >= MIN_IDENT and query coverage >= MIN_COV.
                     NOT ecotype-specific. Informative in its own right: the locus exists in the
                     other assembly but was not annotated there, i.e. an annotation gap, and the
                     count of these is a direct estimate of how leaky "BRAKER-new" is as a proxy
                     for ecotype-specific.
  weak_hit_in_other  a hit below those thresholds. Ambiguous -- a diverged copy, a pseudogene, or a
                     fragment. Reported separately rather than forced into either bucket.
  absent_from_other  no hit at all. The surviving ecotype-specific candidate.

Independent corroboration. Phase 2 (23.4) called ecotype-specific sequence from READ coverage:
analyses/23-reciprocal-pav/{eco}_specific.present_regions.on_{eco}.bed holds regions covered by this
ecotype's reads and not by the other ecotype's. A candidate whose locus overlaps one of those regions
is supported by two independent methods (protein-level absence here; read-level absence in Phase 2)
and is flagged pav_supported. That intersection, not the raw absent_from_other count, is the set
worth carrying into Phase 5.

LIMITS, which matter for how the output is read:
  * Absence of a miniprot hit is evidence of absence only up to assembly completeness. A gene lying
    in an unassembled or collapsed region of the other assembly looks absent and is not. The
    pav_supported flag is the guard against this, because read-level evidence does not depend on the
    other assembly being complete.
  * This screen tests PRESENCE/ABSENCE, not copy number. 23.19 flags residual haplotypes in siscowet
    (2.84 vs 2.58 Gb, higher BUSCO duplication); duplicate-derived models inflate the candidate list
    at the input, and a gene present twice in one ecotype and once in the other is invisible here.
  * Thresholds are deliberately permissive (60/60). The screen is built to be conservative about
    CALLING something ecotype-specific, so a marginal hit counts as presence.

Inputs (per ecotype, from 23.19 / 23.4):
  analyses/23-denovo-annotation/final/{eco}/braker_gene_classes.tsv
  analyses/23-denovo-annotation/final/{eco}/{eco}.braker_new.aa
  analyses/23-denovo-annotation/final/{eco}/{eco}.merged.genes.bed
  output/21.1-gene-anchored-synteny/{other}.purged.fa
  analyses/23-reciprocal-pav/{eco}_specific.present_regions.on_{eco}.bed   (optional)

Outputs (analyses/23-denovo-annotation/final/{eco}/):
  {eco}.reciprocal_screen.tsv      every candidate: call, best hit, identity, coverage, pav_supported
  {eco}.ecotype_specific.ids       absent_from_other candidates (one gene ID per line)
  {eco}.ecotype_specific.aa        their proteins
  reciprocal_screen_summary.txt

Usage (raven, after 23.19):
  /home/sr320/miniconda3/bin/python code/23.22-reciprocal-absence-screen.py [lean siscowet]

Author: Generated for project-lake-trout
Date: 2026-11-01
"""

import argparse
import subprocess
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parent.parent
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
BEDTOOLS = "/usr/bin/bedtools"
MINIPROT = "miniprot"

MIN_IDENT = 0.60
MIN_COV = 0.60
OTHER = {"lean": "siscowet", "siscowet": "lean"}


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", **kw)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def parse_paf(gff_path):
    """Best hit per query from miniprot's ##PAF comment lines.

    Fields after the ##PAF token are standard PAF, so with the token as field 0:
    1 qname, 2 qlen, 3 qstart, 4 qend, 5 strand, 6 tname, 7 tlen, 8 tstart, 9 tend,
    10 nmatch, 11 alnlen, 12 mapq.  identity = nmatch/alnlen, coverage = (qend-qstart)/qlen.
    """
    best = {}
    with open(gff_path) as fh:
        for line in fh:
            if not line.startswith("##PAF"):
                continue
            f = line.split("\t")
            if len(f) < 13:
                continue
            q, qlen, qs, qe = f[1], int(f[2]), int(f[3]), int(f[4])
            tname, nmatch, alnlen = f[6], int(f[10]), int(f[11])
            if qlen <= 0 or alnlen <= 0:
                continue
            ident, cov = nmatch / alnlen, (qe - qs) / qlen
            score = ident * cov
            if q not in best or score > best[q]["score"]:
                best[q] = {"hit_contig": tname, "ident": ident, "cov": cov, "score": score}
    return best


def screen(eco, p4, keep_lq):
    other = OTHER[eco]
    out = p4 / "final" / eco
    work = out / "reciprocal_screen"
    work.mkdir(parents=True, exist_ok=True)

    cl = pd.read_csv(out / "braker_gene_classes.tsv", sep="\t")
    sel = cl["class"].eq("new") & cl["ref_status"].eq("no_ref_hit")
    if not keep_lq:
        sel &= cl["quality"].eq("hq")
    cand = cl[sel].copy()
    if cand.empty:
        print(f"{eco}: no candidates")
        return None

    # 23.19 writes one protein per new gene, named by transcript; map back to gene via the table.
    aa = read_fasta(out / f"{eco}.braker_new.aa")
    tx_col = "braker_tx" if "braker_tx" in cand.columns else "tx"
    ids = set(cand[tx_col]) if tx_col in cand.columns else set(cand["braker_gene"])
    q_fa = work / "candidates.aa"
    with open(q_fa, "w") as f:
        for name, seq in aa.items():
            if name in ids:
                f.write(f">{name}\n{seq}\n")

    target = SYN / f"{other}.purged.fa"
    gff = work / f"vs_{other}.gff"
    sh(f"{MINIPROT} -t 16 --gff {target} {q_fa} > {gff}")
    best = parse_paf(gff)

    key = tx_col if tx_col in cand.columns else "braker_gene"
    cand["hit_contig"] = cand[key].map(lambda q: best.get(q, {}).get("hit_contig"))
    cand["hit_ident"] = cand[key].map(lambda q: best.get(q, {}).get("ident"))
    cand["hit_cov"] = cand[key].map(lambda q: best.get(q, {}).get("cov"))

    def call(r):
        if pd.isna(r["hit_ident"]):
            return "absent_from_other"
        if r["hit_ident"] >= MIN_IDENT and r["hit_cov"] >= MIN_COV:
            return "present_in_other"
        return "weak_hit_in_other"

    cand["screen_call"] = cand.apply(call, axis=1)

    # Independent read-level corroboration from Phase 2 (23.4), if present.
    pav = BASE / "analyses" / "23-reciprocal-pav" / f"{eco}_specific.present_regions.on_{eco}.bed"
    cand["pav_supported"] = False
    genes_bed = out / f"{eco}.merged.genes.bed"
    if pav.exists() and genes_bed.exists():
        gb = pd.read_csv(genes_bed, sep="\t", header=None,
                         names=["chrom", "start", "end", "gene_id", "provenance", "strand"])
        sub = gb[gb["gene_id"].isin(set(cand["braker_gene"]))]
        if not sub.empty:
            a = work / "cand.bed"
            sub.sort_values(["chrom", "start"]).to_csv(a, sep="\t", header=False, index=False)
            hit = work / "cand_in_pav.bed"
            sh(f"{BEDTOOLS} intersect -u -a {a} -b {pav} > {hit}")
            if hit.stat().st_size:
                supported = set(pd.read_csv(hit, sep="\t", header=None, usecols=[3])[3])
                cand["pav_supported"] = cand["braker_gene"].isin(supported)
    else:
        print(f"{eco}: WARNING no Phase-2 PAV regions at {pav} -- pav_supported left False")

    cand.to_csv(out / f"{eco}.reciprocal_screen.tsv", sep="\t", index=False)

    keep = cand[cand["screen_call"].eq("absent_from_other")]
    keep[["braker_gene"]].to_csv(out / f"{eco}.ecotype_specific.ids",
                                 sep="\t", header=False, index=False)
    keep_ids = set(keep[key])
    with open(out / f"{eco}.ecotype_specific.aa", "w") as f:
        for name, seq in aa.items():
            if name in keep_ids:
                f.write(f">{name} ecotype_specific={eco}\n{seq}\n")

    n = len(cand)
    counts = cand["screen_call"].value_counts()
    both = int((cand["screen_call"].eq("absent_from_other") & cand["pav_supported"]).sum())
    return {
        "ecotype": eco, "other": other, "candidates": n,
        "absent_from_other": int(counts.get("absent_from_other", 0)),
        "weak_hit_in_other": int(counts.get("weak_hit_in_other", 0)),
        "present_in_other": int(counts.get("present_in_other", 0)),
        "pav_supported": both,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ecotypes", nargs="*", default=["lean", "siscowet"])
    ap.add_argument("--p4", default="analyses/23-denovo-annotation",
                    help="Phase 4 output root (default: %(default)s)")
    ap.add_argument("--keep-lq", action="store_true",
                    help="also screen lq candidates (default: hq only, matching 23.19's headline)")
    a = ap.parse_args()

    p4 = (BASE / a.p4) if not Path(a.p4).is_absolute() else Path(a.p4)
    rows = [r for eco in a.ecotypes if (r := screen(eco, p4, a.keep_lq))]
    if not rows:
        return
    df = pd.DataFrame(rows)
    summary = p4 / "final" / "reciprocal_screen_summary.txt"
    with open(summary, "w") as f:
        f.write("Reciprocal genomic-absence screen of BRAKER-new candidates (23.22)\n")
        f.write(f"thresholds: present if identity >= {MIN_IDENT:.2f} and query coverage >= {MIN_COV:.2f}\n\n")
        f.write(df.to_string(index=False))
        f.write("\n\nabsent_from_other = no miniprot locus in the other assembly (ecotype-specific candidate).\n")
        f.write("present_in_other  = the locus exists there but was not annotated; an annotation gap,\n")
        f.write("                    and a direct measure of how leaky 'BRAKER-new' is as a proxy.\n")
        f.write("pav_supported     = also inside a Phase-2 (23.4) read-level ecotype-specific region;\n")
        f.write("                    two independent lines. Carry THIS set into Phase 5, not the raw count.\n")
    print(df.to_string(index=False))
    print(f"\nwrote {summary}")


if __name__ == "__main__":
    main()
