#!/usr/bin/env python3
"""
23.20-functional-annotation.py  (Phase 4, step 5, of code/23-next-phase-research-plan.md)

FUNCTIONAL ANNOTATION OF THE NEW BRAKER3 GENES and the per-ecotype gene function table. Liftoff
genes inherit the reference annotation (analyses/18-annotation/gene_function_table.tsv) through
their reference gene ID; only the BRAKER-unique (`new`) genes from 23.19 need new annotation.
This follows plan 20 Step 5, with Swiss-Prot + Pfam standing in for eggNOG-mapper / InterProScan,
whose databases (~50 GB / ~15 GB) are not on raven.

New genes (final/{eco}/{eco}.braker_new.aa):
  1. Homology: DIAMOND blastp --sensitive vs UniProt Swiss-Prot. Best hit with e <= 1e-5 and query
     cover >= 50 % is kept; it is NAMED when identity >= NAME_PID and e <= 1e-10 (symbol =
     Swiss-Prot primary gene name, product = Swiss-Prot protein name + "-like").
  2. GO: terms of the named Swiss-Prot hit (UniProt reviewed-entry GO table), by aspect.
  3. Domains: hmmsearch --cut_ga vs Pfam-A; per gene the domains in sequence order.
  4. TE flag: a transposable-element Pfam domain (RVT_*, rve, DDE_*, Transposase*, MULE, HAT*,
     Exo_endo_phos*, gag/RNase H families ...) or a TE-like hit name, or CDS repeat fraction
     >= 0.5 (23.19); a DNA_pol_B domain without a near-identical RefSeq match (Polinton / Maverick
     self-synthesizing transposons; real POLA/D/E genes match the reference); or a protein sequence
     identical in >= MULTICOPY_MIN new genes (high-copy element).
Output columns match the reference table (gene_id, symbol, geneid, biotype, chrom, start, end,
strand, is_named, product, n_go, go_bp, go_mf, go_cc, go_terms) plus provenance / evidence columns;
start is 0-based (BED convention, as in final/{eco}/{eco}.merged.genes.bed) for every row.

Outputs (analyses/23-denovo-annotation/final/{eco}/):
  {eco}.braker_new.function.tsv   new genes: annotation + evidence (quality, ref_status, hits, Pfam, TE)
  gene_function_table.tsv         all merged genes on native coordinates: Liftoff rows from the
                                  reference table (gene_id = native Liftoff ID, ref_gene_id =
                                  canonical), BRAKER rows from above
  function_summary.txt

Databases (downloaded once to analyses/23-denovo-annotation/evidence/function/, approved
2026-10-17): uniprot_sprot.fasta.gz (UniProt FTP), sprot_go.tsv.gz (UniProt REST, reviewed entries,
accession + gene + protein name + GO), Pfam-A.hmm (EBI FTP). Release info is written to
db_versions.txt.

Run on raven after 23.19: /home/sr320/miniconda3/bin/python code/23.20-functional-annotation.py
  smoke test: ... lean --p4 analyses/23-denovo-annotation/smoke

Author: Generated for project-lake-trout
Date: 2026-10-17
"""

import argparse
import gzip
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
DBDIR = BASE / "analyses" / "23-denovo-annotation" / "evidence" / "function"
REF_TABLE = BASE / "analyses" / "18-annotation" / "gene_function_table.tsv"
DIAMOND = "/home/shared/diamond-2.1.8"
HMMSEARCH = "/home/shared/hmmer-3.4/src/hmmsearch"
HMMPRESS = "/home/shared/hmmer-3.4/src/hmmpress"
THREADS = 16

URL_SPROT = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
URL_SPROT_REL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt"
URL_GO = ("https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cgene_primary"
          "%2Cprotein_name%2Cgo_p%2Cgo_f%2Cgo_c&format=tsv&query=%28reviewed%3Atrue%29")
URL_PFAM = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
URL_PFAM_REL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt"

HIT_EVALUE, HIT_QCOV = 1e-5, 50
MULTICOPY_MIN = 3
NAME_PID, NAME_EVALUE = 30.0, 1e-10
TE_PFAM = re.compile(r"^(RVT_\d|RVT_thumb|rve|rve_\d|DDE_|Transposase|MULE|HAT|hAT|Dimer_Tnp_hAT|"
                     r"Exo_endo_phos|DUF4371|Tnp_|THAP|Retrotrans_gag|gag|zf-CCHC|RNase_H|"
                     r"Transpos|Phage_integrase|Chromo|Helitron|Tc1|CENP-B_dimeris|Plant_tran|"
                     r"LINE|Mariner|Piggy|DUF1758|Peptidase_A17)", re.I)
TE_NAME = re.compile(r"transpos|retrovir|reverse transcriptase|gag-pol|\bpol polyprotein|"
                     r"integrase|LINE-1|retroelement|piggyBac|Tigger|pogo|Jerky|harbinger", re.I)


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", **kw)


# =============================================================================
def fetch_databases():
    DBDIR.mkdir(parents=True, exist_ok=True)
    sprot = DBDIR / "uniprot_sprot.fasta.gz"
    go = DBDIR / "sprot_go.tsv.gz"
    pfam = DBDIR / "Pfam-A.hmm"
    if not sprot.exists():
        sh(f"curl -sSfL -o '{sprot}.tmp' '{URL_SPROT}' && mv '{sprot}.tmp' '{sprot}'")
    if not (DBDIR / "sprot.dmnd").exists():
        sh(f"'{DIAMOND}' makedb --in '{sprot}' -d '{DBDIR / 'sprot'}' --quiet")
    if not go.exists():
        sh(f"curl -sSfL -o '{go}.tmp' '{URL_GO}' && mv '{go}.tmp' '{go}'")
    if not pfam.exists():
        sh(f"curl -sSfL '{URL_PFAM}' | gunzip > '{pfam}.tmp' && mv '{pfam}.tmp' '{pfam}'")
    if not Path(f"{pfam}.h3m").exists():
        sh(f"'{HMMPRESS}' '{pfam}' > /dev/null")
    ver = DBDIR / "db_versions.txt"
    if not ver.exists():
        rel = sh(f"curl -sSfL '{URL_SPROT_REL}'; echo; curl -sSfL '{URL_PFAM_REL}' | head -5",
                 capture_output=True, text=True).stdout
        ver.write_text(rel)
    return sprot, go, pfam


def load_go(path):
    df = pd.read_csv(path, sep="\t", dtype=str, compression="gzip").fillna("")
    df.columns = ["acc", "gene", "protein", "go_p", "go_f", "go_c"]
    ids = lambda s: re.findall(r"GO:\d{7}", s)
    return {r.acc: {"gene": r.gene, "protein": re.split(r" \(", r.protein)[0],
                    "P": ids(r.go_p), "F": ids(r.go_f), "C": ids(r.go_c)}
            for r in df.itertuples(index=False)}


def read_fasta_headers(path):
    out = {}
    for line in open(path):
        if line.startswith(">"):
            name, *rest = line[1:].split()
            out[name] = dict(kv.split("=", 1) for kv in rest if "=" in kv)
    return out


# =============================================================================
def annotate_new(eco, final, sprot_db, go, pfam):
    aa = final / f"{eco}.braker_new.aa"
    classes = pd.read_csv(final / "braker_gene_classes.tsv", sep="\t")
    new = classes[classes["class"] == "new"].copy()
    new["gene_id"] = f"{eco}_braker_" + new["braker_gene"]
    seqs = {}
    name = None
    for line in open(aa):
        if line.startswith(">"):
            name = line[1:].split()[0].rsplit(".", 1)[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip().rstrip("*")
    copies = pd.Series(list(seqs.values())).value_counts()
    multicopy = {g for g, sq in seqs.items() if sq and copies.get(sq, 0) >= MULTICOPY_MIN}

    # 1. Swiss-Prot
    hits = final / "_new_vs_sprot.tsv"
    sh(f"'{DIAMOND}' blastp -q '{aa}' -d '{sprot_db}' -o '{hits}' --sensitive --evalue {HIT_EVALUE} "
       f"--query-cover {HIT_QCOV} --max-target-seqs 1 --threads {THREADS} --quiet "
       f"--outfmt 6 qseqid sseqid pident qcovhsp evalue bitscore stitle")
    h = pd.read_csv(hits, sep="\t", header=None, quoting=3,
                    names=["q", "sseqid", "sp_pident", "sp_qcov", "sp_evalue", "sp_bits", "stitle"]) \
        if hits.stat().st_size else pd.DataFrame(columns=["q", "sseqid", "sp_pident", "sp_qcov",
                                                          "sp_evalue", "sp_bits", "stitle"])
    h = h.sort_values("sp_bits", ascending=False).drop_duplicates("q")
    h["sp_acc"] = h["sseqid"].str.split("|").str[1]
    h["sp_desc"] = h["stitle"].str.replace(r"^\S+ ", "", regex=True).str.replace(r" OS=.*", "", regex=True)
    h["sp_species"] = h["stitle"].str.extract(r"OS=(.+?) OX=", expand=False)
    h["gene_id"] = h["q"].map(lambda q: q.rsplit(".", 1)[0])

    # 3. Pfam
    dom = final / "_new_pfam.domtbl"
    sh(f"'{HMMSEARCH}' --cut_ga --cpu {THREADS} --noali -o /dev/null --domtblout '{dom}' '{pfam}' '{aa}'")
    pf = defaultdict(list)
    for line in open(dom):
        if line.startswith("#"):
            continue
        p = line.split()
        pf[p[0].rsplit(".", 1)[0]].append((int(p[17]), p[3]))     # (ali_from, domain name)
    pfam_str = {g: ",".join(dict.fromkeys(d for _, d in sorted(v))) for g, v in pf.items()}

    n = new.merge(h[["gene_id", "sp_acc", "sp_desc", "sp_species", "sp_pident", "sp_qcov", "sp_evalue"]],
                  on="gene_id", how="left")
    n["pfam_domains"] = n["gene_id"].map(pfam_str).fillna("")
    named = (n["sp_pident"] >= NAME_PID) & (n["sp_evalue"] <= NAME_EVALUE)
    rows = []
    for r, is_named in zip(n.itertuples(index=False), named):
        g = go.get(r.sp_acc, {}) if is_named else {}
        sym = (g.get("gene") or "").lower()
        prod = g.get("protein") or (r.sp_desc if isinstance(r.sp_desc, str) else "")
        doms = [d for d in r.pfam_domains.split(",") if d]
        te = bool(any(TE_PFAM.match(d) for d in doms)
                  or (isinstance(r.sp_desc, str) and TE_NAME.search(r.sp_desc))
                  or (r.repeat_frac >= 0.5)
                  or ("DNA_pol_B" in doms and r.ref_status != "ref_homolog")
                  or (r.gene_id in multicopy))
        gos = {k: sorted(set(g.get(k, []))) for k in ("P", "F", "C")}
        rows.append({
            "gene_id": r.gene_id,
            "symbol": f"{sym}-like" if is_named and sym else r.gene_id,
            "geneid": "", "biotype": "protein_coding",
            "chrom": r.seqid, "start": int(r.start), "end": int(r.end), "strand": r.strand,
            "is_named": int(bool(is_named and sym)),
            "product": (f"{prod}-like" if is_named and prod else
                        ("uncharacterized protein" + (f" with {r.pfam_domains} domain(s)" if r.pfam_domains else ""))),
            "n_go": sum(len(v) for v in gos.values()),
            "go_bp": len(gos["P"]), "go_mf": len(gos["F"]), "go_cc": len(gos["C"]),
            "go_terms": ";".join(gos["P"] + gos["F"] + gos["C"]),
            "provenance": "braker", "quality": r.quality, "ref_status": r.ref_status,
            "ref_hit": r.ref_hit if isinstance(r.ref_hit, str) else "",
            "ref_pident": r.ref_pident,
            "sp_acc": r.sp_acc if isinstance(r.sp_acc, str) else "",
            "sp_desc": r.sp_desc if isinstance(r.sp_desc, str) else "",
            "sp_species": r.sp_species if isinstance(r.sp_species, str) else "",
            "sp_pident": r.sp_pident, "sp_qcov": r.sp_qcov, "sp_evalue": r.sp_evalue,
            "pfam_domains": r.pfam_domains, "te_like": int(te),
            "multicopy_identical": int(r.gene_id in multicopy),
            "aa_len": r.aa_len, "n_introns": r.n_introns, "introns_rna": r.introns_rna,
            "prothint_cov": r.prothint_cov, "repeat_frac": r.repeat_frac,
        })
    out = pd.DataFrame(rows)
    out.to_csv(final / f"{eco}.braker_new.function.tsv", sep="\t", index=False, float_format="%.4g")
    for f in (hits, dom):
        f.unlink(missing_ok=True)
    return out


def ecotype_table(eco, final, newfun):
    """Native-coordinate function table: Liftoff genes from the reference table + new BRAKER genes."""
    ref = pd.read_csv(REF_TABLE, sep="\t", dtype=str).set_index("gene_id")
    bed = pd.read_csv(final / f"{eco}.merged.genes.bed", sep="\t", header=None,
                      names=["chrom", "start", "end", "gene_id", "provenance", "strand"])
    lift = bed[bed["provenance"] == "liftoff"].copy()
    lift["ref_gene_id"] = lift["gene_id"].str.replace(r"_\d+$", "", regex=True)
    cols = ["symbol", "geneid", "biotype", "is_named", "product", "n_go", "go_bp", "go_mf", "go_cc", "go_terms"]
    lift = lift.join(ref[cols], on="ref_gene_id")
    base_cols = ["gene_id", "symbol", "geneid", "biotype", "chrom", "start", "end", "strand", "is_named",
                 "product", "n_go", "go_bp", "go_mf", "go_cc", "go_terms"]
    lift = lift[base_cols + ["provenance", "ref_gene_id"]]
    nf = newfun.copy()
    nf["ref_gene_id"] = ""
    nf["start"] = nf["start"] - 1                # GFF 1-based -> BED-style 0-based, as for Liftoff rows
    extra = [c for c in nf.columns if c not in base_cols + ["provenance", "ref_gene_id"]]
    table = pd.concat([lift, nf[base_cols + ["provenance", "ref_gene_id"] + extra]], ignore_index=True)
    table = table.sort_values(["chrom", "start"])
    table.to_csv(final / "gene_function_table.tsv", sep="\t", index=False, float_format="%.4g")
    return table


def summarize(eco, final, nf, table):
    lines = [f"{eco}: merged gene function table {len(table):,} genes "
             f"(liftoff {int((table['provenance'] == 'liftoff').sum()):,}, braker {len(nf):,})"]
    for q in ("hq", "lq"):
        x = nf[nf["quality"] == q]
        if not len(x):
            continue
        lines.append(f"  new {q} {len(x):,}: Swiss-Prot hit {(x['sp_acc'] != '').mean():.0%}, named "
                     f"{x['is_named'].mean():.0%}, Pfam {(x['pfam_domains'] != '').mean():.0%}, GO "
                     f"{(x['n_go'] > 0).mean():.0%}, TE-like {x['te_like'].mean():.0%}")
        for st in ("no_ref_hit", "ref_paralog", "ref_homolog"):
            y = x[x["ref_status"] == st]
            if len(y):
                lines.append(f"    {st:11s} {len(y):5,}: TE-like {y['te_like'].mean():.0%}, named "
                             f"{y['is_named'].mean():.0%}")
    cand = nf[(nf["quality"] == "hq") & (nf["ref_status"] == "no_ref_hit") & (nf["te_like"] == 0)]
    lines.append(f"  novel candidates (hq, no SaNama RefSeq hit, not TE-like): {len(cand):,}; "
                 f"of these named via Swiss-Prot {int(cand['is_named'].sum()):,}, Pfam-only "
                 f"{int(((cand['is_named'] == 0) & (cand['pfam_domains'] != '')).sum()):,}, unannotated "
                 f"{int(((cand['is_named'] == 0) & (cand['pfam_domains'] == '')).sum()):,}")
    for r in cand.sort_values("aa_len", ascending=False).head(25).itertuples():
        lines.append(f"    {r.gene_id:28s} {r.chrom}:{r.start}-{r.end} aa={r.aa_len:<5} "
                     f"{r.symbol[:16]:16s} {str(r.product)[:50]:50s} {r.pfam_domains[:30]}")
    text = "\n".join(lines) + "\n"
    (final / "function_summary.txt").write_text(text)
    print(text, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ecotypes", nargs="*", default=["lean", "siscowet"])
    ap.add_argument("--p4", type=Path, default=BASE / "analyses" / "23-denovo-annotation")
    a = ap.parse_args()
    p4 = a.p4 if a.p4.is_absolute() else BASE / a.p4
    sprot, go_path, pfam = fetch_databases()
    go = load_go(go_path)
    for eco in a.ecotypes:
        final = p4 / "final" / eco
        if not (final / "braker_gene_classes.tsv").exists():
            print(f"{eco}: no 23.19 output in {final}; skipped")
            continue
        nf = annotate_new(eco, final, DBDIR / "sprot", go, pfam)
        table = ecotype_table(eco, final, nf)
        summarize(eco, final, nf, table)


if __name__ == "__main__":
    main()
