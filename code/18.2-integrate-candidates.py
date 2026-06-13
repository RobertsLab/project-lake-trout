#!/usr/bin/env python3
"""
18.2-integrate-candidates.py  (Step 3 of code/18-diff-annotation-phenotype-plan.md)

Integrate the methylation (Step 1) and PAV (Step 2) gene assignments into one ranked
candidate table, add liver-RNAseq expression evidence, and flag CONVERGENT genes hit by
both differential methylation and a stringent siscowet-specific deletion.

Join key throughout: `gene_id` (= `gene-XXX`), shared by the annotation table, both
assignment tables, the burden table, and data/whole_tx_table.csv.

Master set = genes with >=1 DMR, DMC, or stringent-PAV assignment within 5 kb. Lenient PAV
burden is attached as context (burden-only genes live in pav_gene_burden.tsv, not here).

EXPRESSION CAVEAT: data/whole_tx_table.csv is liver RNAseq from a *separate* parasite study
(12 lean + 12 siscowet, different individuals than the PacBio methylation/PAV fish). It is
orthogonal support, one tissue only — not confirmation. Carried as a column note.
"""

import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
ANN = BASE / "analyses" / "18-annotation"
RNA = BASE / "data" / "whole_tx_table.csv"
OUT = ANN / "integrated_candidate_genes.tsv"

CLASS_RANK = {"promoter": 4, "exon": 3, "intron": 2,
              "upstream<=5kb": 1, "downstream<=5kb": 1}
RANK2CLASS = {v: k for k, v in
              {"promoter": 4, "exon": 3, "intron": 2, "flank<=5kb": 1}.items()}


def load_genes():
    g = {}
    with open(ANN / "gene_function_table.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g[row["gene_id"]] = row
    return g


def best_class(classes):
    if not classes:
        return ""
    return max(classes, key=lambda c: CLASS_RANK.get(c, 0))


def load_dmr():
    d = defaultdict(lambda: {"n": 0, "classes": [], "dirs": set(),
                             "maxdiff": 0.0, "promoter": 0})
    with open(ANN / "dmr_gene_assignments.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            x = d[r["gene_id"]]
            x["n"] += 1
            x["classes"].append(r["location_class"])
            x["dirs"].add(r["direction"])
            x["maxdiff"] = max(x["maxdiff"], abs(float(r["meth_diff"])))
            x["promoter"] = max(x["promoter"], int(r["in_promoter"]))
    return d


def load_dmc():
    d = defaultdict(lambda: {"n": 0, "promoter_n": 0})
    with open(ANN / "dmc_gene_assignments.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            x = d[r["gene_id"]]
            x["n"] += 1
            x["promoter_n"] += int(r["in_promoter"])
    return d


def load_pav_strin():
    d = defaultdict(lambda: {"n": 0, "classes": [], "exonic": 0, "maxsize": 0})
    with open(ANN / "pav_gene_assignments.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            x = d[r["gene_id"]]
            x["n"] += 1
            x["classes"].append(r["location_class"])
            x["exonic"] = max(x["exonic"], int(r["exon_overlap"]))
            x["maxsize"] = max(x["maxsize"], int(r["size_bp"]))
    return d


def load_burden():
    d = {}
    with open(ANN / "pav_gene_burden.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["gene_id"]] = r
    return d


def load_expression():
    """gene_id -> (lean_mean_fpkm, sisco_mean_fpkm) summing transcript FPKM per gene."""
    with open(RNA) as fh:
        hdr = next(csv.reader(fh))
    fpkm_cols = [h for h in hdr if h.startswith("FPKM.")]
    lean_c = [h for h in fpkm_cols if "LL" in h]
    sis_c = [h for h in fpkm_cols if "SL" in h]
    gene_lean = defaultdict(lambda: [0.0] * len(lean_c))
    gene_sis = defaultdict(lambda: [0.0] * len(sis_c))
    with open(RNA) as fh:
        for r in csv.DictReader(fh):
            gid = r["gene_id"]
            for i, c in enumerate(lean_c):
                try:
                    gene_lean[gid][i] += float(r[c])
                except (ValueError, KeyError):
                    pass
            for i, c in enumerate(sis_c):
                try:
                    gene_sis[gid][i] += float(r[c])
                except (ValueError, KeyError):
                    pass
    expr = {}
    keys = set(gene_lean) | set(gene_sis)
    for gid in keys:
        lv = gene_lean.get(gid, [0.0])
        sv = gene_sis.get(gid, [0.0])
        lm = sum(lv) / len(lv) if lv else 0.0
        sm = sum(sv) / len(sv) if sv else 0.0
        expr[gid] = (lm, sm)
    return expr


def concordance(dmr, log2fc, expressed):
    """promoter hyper_siscowet -> expect siscowet DOWN (log2fc<0); hypo -> UP."""
    if not dmr or not dmr["promoter"] or not expressed or log2fc is None:
        return "NA"
    if abs(log2fc) < 0.5:
        return "flat"
    dirs = dmr["dirs"]
    has_hyper = "hyper_siscowet" in dirs
    has_hypo = "hypo_siscowet" in dirs
    if has_hyper and has_hypo:
        return "mixed"
    if has_hyper:
        return "concordant" if log2fc < 0 else "discordant"
    if has_hypo:
        return "concordant" if log2fc > 0 else "discordant"
    return "NA"


def main():
    genes = load_genes()
    dmr = load_dmr()
    dmc = load_dmc()
    pav = load_pav_strin()
    burden = load_burden()
    expr = load_expression()

    candidates = set(dmr) | set(dmc) | set(pav)
    print(f"candidate genes (DMR/DMC/stringent-PAV within 5kb): {len(candidates)}",
          file=sys.stderr)

    cols = [
        "gene_id", "symbol", "product", "biotype", "n_go", "chrom", "start", "end",
        "strand",
        "dmr_n", "dmr_best_class", "dmr_promoter", "dmr_directions", "dmr_max_absdiff",
        "dmc_n", "dmc_promoter_n",
        "pav_strin_n", "pav_strin_best_class", "pav_strin_exonic", "pav_strin_maxsize",
        "lean_del_bp", "sisco_del_bp", "lean_ins_n", "sisco_ins_n",
        "lean_fpkm", "sisco_fpkm", "log2fc_sisco_lean", "expressed",
        "convergent", "meth_expr_concordant", "caution", "rank_score",
    ]

    # repetitive / multicopy ncRNA biotypes where promoter "DMRs" are usually rDNA-array
    # methylation artifacts, not gene regulation -> deprioritize in the ranked view.
    REPETITIVE = {"rRNA", "snRNA", "snoRNA", "tRNA"}

    rows = []
    n_conv = 0
    for gid in candidates:
        g = genes.get(gid, {})
        dm = dmr.get(gid)
        dc = dmc.get(gid)
        pv = pav.get(gid)
        bd = burden.get(gid, {})

        lm, sm = expr.get(gid, (0.0, 0.0))
        log2fc = math.log2((sm + 0.1) / (lm + 0.1))
        expressed = (lm >= 1.0 or sm >= 1.0)

        convergent = int(bool(dm) and bool(pv))
        if convergent:
            n_conv += 1
        conc = concordance(dm, log2fc if expressed else None, expressed)
        biotype = g.get("biotype", "")
        caution = "repetitive_ncRNA" if biotype in REPETITIVE else ""

        score = 0
        if convergent:
            score += 3
        if dm and dm["promoter"]:
            score += 2
        if pv and pv["exonic"]:
            score += 2
        if dm and best_class(dm["classes"]) == "exon":
            score += 1
        if dc and dc["promoter_n"] > 0:
            score += 1
        if expressed:
            score += 1
        if conc == "concordant":
            score += 2

        rows.append({
            "gene_id": gid, "symbol": g.get("symbol", ""),
            "product": g.get("product", ""), "biotype": g.get("biotype", ""),
            "n_go": g.get("n_go", ""), "chrom": g.get("chrom", ""),
            "start": g.get("start", ""), "end": g.get("end", ""),
            "strand": g.get("strand", ""),
            "dmr_n": dm["n"] if dm else 0,
            "dmr_best_class": best_class(dm["classes"]) if dm else "",
            "dmr_promoter": dm["promoter"] if dm else 0,
            "dmr_directions": ",".join(sorted(dm["dirs"])) if dm else "",
            "dmr_max_absdiff": f"{dm['maxdiff']:.2f}" if dm else "",
            "dmc_n": dc["n"] if dc else 0,
            "dmc_promoter_n": dc["promoter_n"] if dc else 0,
            "pav_strin_n": pv["n"] if pv else 0,
            "pav_strin_best_class": best_class(pv["classes"]) if pv else "",
            "pav_strin_exonic": pv["exonic"] if pv else 0,
            "pav_strin_maxsize": pv["maxsize"] if pv else 0,
            "lean_del_bp": bd.get("lean_del_bp", 0),
            "sisco_del_bp": bd.get("sisco_del_bp", 0),
            "lean_ins_n": bd.get("lean_ins_n", 0),
            "sisco_ins_n": bd.get("sisco_ins_n", 0),
            "lean_fpkm": f"{lm:.2f}", "sisco_fpkm": f"{sm:.2f}",
            "log2fc_sisco_lean": f"{log2fc:.2f}" if expressed else "NA",
            "expressed": int(expressed),
            "convergent": convergent, "meth_expr_concordant": conc,
            "caution": caution, "rank_score": score,
        })

    # primary: score; then protein-coding/non-repetitive first; then effect sizes
    rows.sort(key=lambda r: (r["rank_score"],
                             0 if r["caution"] else 1,
                             1 if r["biotype"] == "protein_coding" else 0,
                             float(r["dmr_max_absdiff"] or 0),
                             r["pav_strin_maxsize"]), reverse=True)

    with open(OUT, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"convergent (DMR + stringent PAV) genes: {n_conv}", file=sys.stderr)
    print(f"wrote {OUT} ({len(rows)} genes)", file=sys.stderr)


if __name__ == "__main__":
    main()
