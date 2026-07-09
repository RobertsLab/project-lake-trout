#!/usr/bin/env python3
"""
23.9-native-meth-dmg.py  (Phase 3, part B, of code/23-next-phase-research-plan.md)

Gene-anchored differential methylation on ECOTYPE-NATIVE coordinates. Phase-3 part A
(23.8-native-meth-pipeline.py) called 5mC for each ecotype on its OWN purged assembly, so lean and
siscowet methylation live on different contig coordinate systems and cannot be merged by chrom:pos
(the constraint flagged in plan 13.3 and the merge used by 14-diff-meth.py).

We bridge the two coordinate systems with the shared **Liftoff gene ID** (both assemblies were lifted
from the same reference, so a gene carries the same ID in each). For every sample we assign each
covered CpG to a gene using that ecotype's native Liftoff annotation and take the coverage-weighted
gene methylation (sum Nmod / sum Nvalid over the gene's CpGs). A gene then has one methylation value
per sample, keyed by gene ID -- directly comparable across ecotypes without any base-level alignment.

We then run the 14-diff-meth statistical recipe per gene (Mann-Whitney U, BH-FDR, 15% difference)
to call **differentially methylated genes (DMGs)** on native coordinates, and compare them back to
the reference-based DMRs (analyses/18-annotation/dmr_gene_assignments.tsv) to see which reference
signals are confirmed vs. reference artifacts.

Inputs (from 23.8): analyses/23-native-methylation/meth/{sample}.cpg_counts.bed  (chrom,start,end,Nmod,Nvalid)
Gene models (native): output/21.1-gene-anchored-synteny/{lean,siscowet}.purged.liftoff.genes.bed

Outputs (analyses/23-native-methylation/):
  native_methylation_by_gene.tsv   gene x sample coverage-weighted methylation (+ per-ecotype n_cpgs)
  dmg_results.tsv                  per-gene test (all genes tested)
  dmg_significant.tsv              significant DMGs
  reference_dmr_confirmation.tsv   reference-DMR genes: confirmed native DMG?  concordant direction?
  native_dmg_report.txt

Author: Generated for project-lake-trout
Date: 2026-07-07
"""

import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parent.parent
METH = BASE / "analyses" / "23-native-methylation" / "meth"
OUT = BASE / "analyses" / "23-native-methylation"
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
ANN = BASE / "analyses" / "18-annotation"
BEDTOOLS = "/usr/bin/bedtools"

SAMPLES_LEAN = ["bc2041", "bc2068", "bc2069", "bc2070"]
SAMPLES_SISCOWET = ["bc2071", "bc2072", "bc2073", "bc2096"]
ECO_GENES = {
    "lean": SYN / "lean.purged.liftoff.genes.bed",
    "siscowet": SYN / "siscowet.purged.liftoff.genes.bed",
}
ECO_SAMPLES = {"lean": SAMPLES_LEAN, "siscowet": SAMPLES_SISCOWET}

# statistical choices mirror 14-diff-meth.py
MIN_CPGS_PER_GENE = 3     # covered CpGs required to score a gene in a sample
MIN_SAMPLES = 3           # samples per ecotype group with a gene value
DIFF_THRESHOLD = 15.0     # percent methylation difference (15%)
P_VALUE_THRESHOLD = 0.05
Q_VALUE_THRESHOLD = 0.10


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


def nonempty(p):
    p = Path(p)
    return p.exists() and p.stat().st_size > 0


def sorted_genes(eco):
    """Native gene BED sorted, unique gene per row (chrom,start,end,gene_id)."""
    src = ECO_GENES[eco]
    dst = OUT / f"_{eco}.genes.sorted.bed"
    sh(f"cut -f1-4 '{src}' | sort -k1,1 -k2,2n > '{dst}'")
    return dst


def gene_methylation(eco, sample, genes_bed):
    """coverage-weighted gene methylation for one sample via bedtools map.
    Returns dict gene_id -> (pct_meth, n_cpgs) for genes with >= MIN_CPGS_PER_GENE covered CpGs."""
    counts = METH / f"{sample}.cpg_counts.bed"
    if not nonempty(counts):
        print(f"  WARNING missing {counts.name}; skipping {sample}")
        return {}
    counts_sorted = OUT / f"_{sample}.counts.sorted.bed"
    sh(f"sort -k1,1 -k2,2n '{counts}' > '{counts_sorted}'")
    mapped = OUT / f"_{sample}.gene_meth.bed"
    # per gene: sum Nmod (col4), sum Nvalid (col5), count CpGs
    sh(f"'{BEDTOOLS}' map -a '{genes_bed}' -b '{counts_sorted}' -c 4,5,4 -o sum,sum,count "
       f"-null 0 > '{mapped}'")
    out = {}
    with open(mapped) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            gene = p[3]
            try:
                nmod, nvalid, ncpg = float(p[4]), float(p[5]), int(p[6])
            except (ValueError, IndexError):
                continue
            if ncpg >= MIN_CPGS_PER_GENE and nvalid > 0:
                out[gene] = (100.0 * nmod / nvalid, ncpg)
    counts_sorted.unlink(missing_ok=True)
    mapped.unlink(missing_ok=True)
    return out


def main():
    if not nonempty(METH / f"{SAMPLES_LEAN[0]}.cpg_counts.bed"):
        sys.exit("FATAL: native methylation counts not found. Run 23.8 first.")
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- per-sample gene methylation (each ecotype on its own genome) ----
    gene_meth = {}     # sample -> {gene: pct}
    gene_ncpg = defaultdict(dict)  # eco -> {gene: representative n_cpgs}
    for eco, samples in ECO_SAMPLES.items():
        gbed = sorted_genes(eco)
        for s in samples:
            print(f"  [{eco}] {s}: aggregating CpGs to genes ...", flush=True)
            gm = gene_methylation(eco, s, gbed)
            gene_meth[s] = {g: v[0] for g, v in gm.items()}
            for g, v in gm.items():
                gene_ncpg[eco][g] = max(gene_ncpg[eco].get(g, 0), v[1])
        gbed.unlink(missing_ok=True)

    all_samples = SAMPLES_LEAN + SAMPLES_SISCOWET
    # ---- gene x sample matrix (genes seen in either ecotype) ----
    genes = set()
    for s in all_samples:
        genes |= set(gene_meth.get(s, {}))
    rows = []
    for g in sorted(genes):
        row = {"gene_id": g}
        for s in all_samples:
            row[s] = gene_meth.get(s, {}).get(g, np.nan)
        row["n_cpgs_lean"] = gene_ncpg["lean"].get(g, 0)
        row["n_cpgs_siscowet"] = gene_ncpg["siscowet"].get(g, 0)
        rows.append(row)
    mat = pd.DataFrame(rows).set_index("gene_id")
    mat.to_csv(OUT / "native_methylation_by_gene.tsv", sep="\t")
    print(f"  genes with methylation in >=1 sample: {len(mat):,}")

    # ---- gene annotation ----
    fn = {}
    with open(ANN / "gene_function_table.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fn[r["gene_id"]] = (r.get("symbol", ""), r.get("product", ""))

    # ---- per-gene Mann-Whitney (lean vs siscowet), on the shared gene-ID axis ----
    results = []
    for g, r in mat.iterrows():
        lean_vals = pd.to_numeric(r[SAMPLES_LEAN], errors="coerce").dropna()
        sisc_vals = pd.to_numeric(r[SAMPLES_SISCOWET], errors="coerce").dropna()
        if len(lean_vals) < MIN_SAMPLES or len(sisc_vals) < MIN_SAMPLES:
            continue
        lean_mean, sisc_mean = lean_vals.mean(), sisc_vals.mean()
        diff = sisc_mean - lean_mean            # + = hyper in siscowet
        try:
            stat, pval = stats.mannwhitneyu(lean_vals, sisc_vals, alternative="two-sided")
        except ValueError:
            continue
        sym, prod = fn.get(g, ("", ""))
        results.append({
            "gene_id": g, "symbol": sym, "product": prod,
            "lean_mean": round(lean_mean, 2), "siscowet_mean": round(sisc_mean, 2),
            "meth_diff": round(diff, 2), "n_lean": len(lean_vals), "n_siscowet": len(sisc_vals),
            "n_cpgs_lean": int(r["n_cpgs_lean"]), "n_cpgs_siscowet": int(r["n_cpgs_siscowet"]),
            "pvalue": pval,
        })
    if not results:
        sys.exit("No genes had enough samples in both ecotypes to test.")
    df = pd.DataFrame(results)
    df["qvalue"] = multipletests(df["pvalue"], method="fdr_bh")[1]
    df["significant"] = (df["pvalue"] < P_VALUE_THRESHOLD) & \
                        (df["meth_diff"].abs() >= DIFF_THRESHOLD)
    df["significant_fdr"] = (df["qvalue"] < Q_VALUE_THRESHOLD) & \
                            (df["meth_diff"].abs() >= DIFF_THRESHOLD)
    df["direction"] = np.where(df["meth_diff"] > 0, "hyper_siscowet", "hypo_siscowet")
    df = df.sort_values(["significant_fdr", "significant", "pvalue"],
                        ascending=[False, False, True])
    df.to_csv(OUT / "dmg_results.tsv", sep="\t", index=False)
    sig = df[df["significant"]]
    sig.to_csv(OUT / "dmg_significant.tsv", sep="\t", index=False)
    print(f"  genes tested: {len(df):,}  DMG (p<0.05, |diff|>=15%): {len(sig):,}  "
          f"(FDR<0.10: {int(df['significant_fdr'].sum()):,})")

    # ---- compare to reference-based DMRs (which reference signals survive natively) ----
    ref_dmr_genes = set()
    dmr_dir = {}
    dpath = ANN / "dmr_gene_assignments.tsv"
    if nonempty(dpath):
        with open(dpath) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                g = r.get("gene_id")
                if g:
                    ref_dmr_genes.add(g)
                    dmr_dir.setdefault(g, r.get("direction", ""))
    tested = set(df["gene_id"])
    sig_set = set(sig["gene_id"])
    native_dir = dict(zip(df["gene_id"], df["direction"]))
    conf_rows = []
    for g in sorted(ref_dmr_genes):
        in_tested = g in tested
        is_dmg = g in sig_set
        nd = native_dir.get(g, "")
        conf_rows.append({
            "gene_id": g, "symbol": fn.get(g, ("", ""))[0],
            "ref_dmr_direction": dmr_dir.get(g, ""),
            "native_testable": int(in_tested),
            "native_dmg": int(is_dmg),
            "native_direction": nd if in_tested else "",
            "direction_concordant": int(bool(is_dmg) and nd and
                                         nd == dmr_dir.get(g, "")),
        })
    conf = pd.DataFrame(conf_rows)
    conf.to_csv(OUT / "reference_dmr_confirmation.tsv", sep="\t", index=False)
    n_ref = len(ref_dmr_genes)
    n_testable = int(conf["native_testable"].sum()) if n_ref else 0
    n_conf = int(conf["native_dmg"].sum()) if n_ref else 0
    n_conc = int(conf["direction_concordant"].sum()) if n_ref else 0

    # ---- report ----
    with open(OUT / "native_dmg_report.txt", "w") as f:
        f.write("Ecotype-native differential methylation (gene-anchored)\n")
        f.write("=" * 56 + "\n\n")
        f.write(f"Genes tested (>= {MIN_SAMPLES} samples/ecotype, "
                f">= {MIN_CPGS_PER_GENE} CpGs/gene): {len(df):,}\n")
        f.write(f"DMGs (p<{P_VALUE_THRESHOLD}, |diff|>={DIFF_THRESHOLD:.0f}%): {len(sig):,}\n")
        f.write(f"  hyper in siscowet: {(sig['direction']=='hyper_siscowet').sum():,}\n")
        f.write(f"  hypo  in siscowet: {(sig['direction']=='hypo_siscowet').sum():,}\n")
        f.write(f"DMGs at FDR<{Q_VALUE_THRESHOLD}: {int(df['significant_fdr'].sum()):,}\n\n")
        f.write("Reference-DMR confirmation on native coordinates\n")
        f.write("-" * 56 + "\n")
        f.write(f"Reference-DMR genes: {n_ref:,}\n")
        f.write(f"  testable natively (enough coverage both ecotypes): {n_testable:,}\n")
        f.write(f"  confirmed as native DMG: {n_conf:,}\n")
        f.write(f"  direction-concordant: {n_conc:,}\n")
        if n_testable:
            f.write(f"  confirmation rate among testable: {n_conf/n_testable:.1%}\n")
        f.write("\nTop 20 native DMGs by |meth_diff|:\n")
        for _, r in sig.reindex(sig["meth_diff"].abs().sort_values(ascending=False).index).head(20).iterrows():
            f.write(f"  {r['gene_id']} {r['symbol']}  diff={r['meth_diff']:+.1f}%  "
                    f"p={r['pvalue']:.1e}  [{r['product'][:45]}]\n")
    print("\n" + (OUT / "native_dmg_report.txt").read_text())
    print(f"Outputs -> {OUT}")


if __name__ == "__main__":
    main()
