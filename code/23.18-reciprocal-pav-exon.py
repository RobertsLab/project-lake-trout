#!/usr/bin/env python3
"""
23.18-reciprocal-pav-exon.py  (follow-up to 23.17, code/23-next-phase-research-plan.md)

WHERE IN THE GENE ARE THE RECIPROCAL-PAV SEGMENTS? 23.17 showed that a reciprocal-PAV "gene" (23.4) is
usually a large gene containing a short native presence/absence segment (median 3 % of the gene). This
asks whether those segments remove coding sequence or sit in introns.

Per ecotype (segments are on the native assembly that carries them, 23.4 present_regions BEDs):
  1. Gene models from the native Liftoff GFF (output/21.1-gene-anchored-synteny/{eco}.purged.liftoff.gff3):
     gene -> transcripts (mRNA and ncRNA) -> exon / CDS, all copies kept by their Liftoff ID.
  2. For each reciprocal-PAV gene: union of CDS, union of exons, intron = gene span - exons.
     Absent bp = present_region ∩ gene, split into CDS / non-coding exon (UTR, ncRNA) / intron.
  3. Gene classes: cds_hit (>= 1 bp CDS absent), exon_hit (non-coding exon only), intron_only.
     Also: fraction of the gene's CDS absent, number of exons touched.
  4. Expectation: for each gene, the segments placed at random inside the same gene would put
     absent_bp x (CDS bp / gene bp) into CDS. Observed / expected CDS bp summed over genes, plus a
     per-gene permutation (N_PERM circular shifts of the gene's segments inside the gene span).

Inputs (raven): analyses/23-reciprocal-pav/{eco}_specific.{genes.tsv,present_regions.on_{eco}.bed};
  output/21.1-gene-anchored-synteny/{eco}.purged.liftoff.gff3; gene_evidence_refined.tsv (--genes) for
  tier / product.

Outputs (analyses/23-lean-baseline/):
  reciprocal_pav_exon.tsv          one row per reciprocal-PAV gene
  reciprocal_pav_exon_summary.tsv  per ecotype x class, and observed / expected CDS bp

Author: Generated for project-lake-trout
Date: 2026-10-16
"""

import argparse
import collections
import csv
import random
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================
BASE = Path(__file__).resolve().parent.parent
RPAV = BASE / "analyses" / "23-reciprocal-pav"
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
GENES = BASE / "analyses" / "23-integration" / "refined" / "gene_evidence_refined.tsv"
OUT = BASE / "analyses" / "23-lean-baseline"
ECOTYPES = ["lean", "siscowet"]
N_PERM = 1000
SEED = 1


# =============================================================================
# Interval helpers (half-open, 0-based)
# =============================================================================
def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def overlap(a, b):
    """bp of overlap between two merged interval lists."""
    i = j = n = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            n += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return n


def length(iv):
    return sum(e - s for s, e in iv)


def attrs(field):
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


# =============================================================================
def load_models(gff, wanted):
    """gene_id -> {chrom, start, end, exons[], cds[], exon_list[]} for wanted genes only."""
    genes, tx2gene = {}, {}
    rows = []
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            a = attrs(p[8])
            ftype, s, e = p[2], int(p[3]) - 1, int(p[4])
            if ftype == "gene" and a.get("ID") in wanted:
                genes[a["ID"]] = {"chrom": p[0], "start": s, "end": e, "exons": [], "cds": [],
                                  "biotype": a.get("gene_biotype", "")}
            elif ftype in ("exon", "CDS"):
                rows.append((ftype, a.get("Parent", ""), s, e))
            elif "Parent" in a and "ID" in a:
                tx2gene[a["ID"]] = a["Parent"]
    for ftype, parent, s, e in rows:
        g = tx2gene.get(parent, parent)
        if g in genes:
            genes[g]["exons" if ftype == "exon" else "cds"].append((s, e))
    for g in genes.values():
        g["exon_list"] = sorted(set(g["exons"]))
        g["exons"], g["cds"] = merge(g["exons"]), merge(g["cds"])
    return genes


def load_bed(path):
    d = collections.defaultdict(list)
    for line in open(path):
        c, s, e = line.split()[:3]
        d[c].append((int(s), int(e)))
    return {c: merge(v) for c, v in d.items()}


def clip(iv, s, e):
    return merge([(max(a, s), min(b, e)) for a, b in iv if a < e and b > s])


def shifted_cds(seg, g, rng):
    """CDS bp hit when the gene's segments are circularly shifted by a random offset."""
    span = g["end"] - g["start"]
    k = rng.randrange(span)
    moved = []
    for s, e in seg:
        s2, e2 = (s - g["start"] + k) % span, (e - g["start"] + k) % span
        if s2 < e2:
            moved.append((g["start"] + s2, g["start"] + e2))
        else:  # wraps
            moved += [(g["start"] + s2, g["end"]), (g["start"], g["start"] + e2)]
    return overlap(merge(moved), g["cds"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default=str(GENES))
    args = ap.parse_args()
    info = {r["gene_id"]: r for r in csv.DictReader(open(args.genes), delimiter="\t")}
    rng = random.Random(SEED)
    rows = []
    for eco in ECOTYPES:
        ids = [l.split("\t")[3] for l in open(RPAV / f"{eco}_specific.genes.tsv").read().split("\n")[1:] if l]
        models = load_models(SYN / f"{eco}.purged.liftoff.gff3", set(ids))
        pres = load_bed(RPAV / f"{eco}_specific.present_regions.on_{eco}.bed")
        for gid in ids:
            g = models.get(gid)
            if g is None:
                print(f"  missing model: {gid}")
                continue
            seg = clip(pres.get(g["chrom"], []), g["start"], g["end"])
            gene_bp = g["end"] - g["start"]
            absent = length(seg)
            exon_abs = overlap(seg, g["exons"])
            cds_abs = overlap(seg, g["cds"])
            cds_bp, exon_bp = length(g["cds"]), length(g["exons"])
            exons_hit = sum(1 for x in g["exon_list"] if overlap(seg, [x]))
            if cds_abs:
                cls = "cds_hit"
            elif exon_abs:
                cls = "noncoding_exon_hit"
            else:
                cls = "intron_only"
            perm = [shifted_cds(seg, g, rng) for _ in range(N_PERM)] if cds_bp else [0]
            base = gid.split("_")[0] if gid.rsplit("_", 1)[-1].isdigit() else gid
            r = info.get(gid) or info.get(base) or {}
            rows.append({
                "gene_id": gid, "present_in": eco, "tier": r.get("tier", ""),
                "product": r.get("product", ""), "biotype": g["biotype"],
                "extra_copy": int(gid != base),
                "gene_bp": gene_bp, "exon_bp": exon_bp, "cds_bp": cds_bp, "absent_bp": absent,
                "absent_frac": round(absent / gene_bp, 4),
                "absent_cds_bp": cds_abs, "absent_noncoding_exon_bp": exon_abs - cds_abs,
                "absent_intron_bp": absent - exon_abs,
                "cds_frac_absent": round(cds_abs / cds_bp, 4) if cds_bp else "",
                "exons_hit": exons_hit, "exons_total": len(g["exon_list"]),
                "expected_cds_bp": round(absent * cds_bp / gene_bp, 1),
                "perm_p_cds_ge_obs": round(sum(p >= cds_abs for p in perm) / len(perm), 4),
                "class": cls,
            })
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "reciprocal_pav_exon.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summ = []
    for eco in ECOTYPES + ["ALL"]:
        sub = [r for r in rows if eco == "ALL" or r["present_in"] == eco]
        for subset, rs in (("all", sub), ("tier_A", [r for r in sub if r["tier"] == "A"]),
                           ("protein_coding", [r for r in sub if r["cds_bp"]])):
            if not rs:
                continue
            c = collections.Counter(r["class"] for r in rs)
            obs = sum(r["absent_cds_bp"] for r in rs)
            exp = sum(r["expected_cds_bp"] for r in rs)
            summ.append({
                "present_in": eco, "subset": subset, "genes": len(rs),
                "cds_hit": c["cds_hit"], "noncoding_exon_hit": c["noncoding_exon_hit"],
                "intron_only": c["intron_only"],
                "cds_frac_absent_ge_0_5": sum(1 for r in rs if r["cds_frac_absent"] != ""
                                              and r["cds_frac_absent"] >= 0.5),
                "obs_cds_bp": obs, "exp_cds_bp": round(exp, 1),
                "obs_over_exp": round(obs / exp, 3) if exp else "",
            })
    with open(OUT / "reciprocal_pav_exon_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0]), delimiter="\t")
        w.writeheader()
        w.writerows(summ)
    for s in summ:
        print("\t".join(str(v) for v in s.values()))


if __name__ == "__main__":
    main()
