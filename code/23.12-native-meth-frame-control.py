#!/usr/bin/env python3
"""
23.12-native-meth-frame-control.py  (Phase 3, part D, of code/23-next-phase-research-plan.md)

FRAME CONTROL for the native methylation differences of 23.11. In 23.11 each ecotype is read on its
OWN assembly and features are paired by Liftoff gene ID, so ecotype is perfectly confounded with
assembly: a gene whose Liftoff placement or upstream sequence differs between the two assemblies can
look "differentially methylated" without any methylation difference. Label permutation cannot catch
this (every permuted split mixes assemblies).

The Phase-2 cross alignments (23.4: lean reads -> siscowet.purged.fa, siscowet reads -> lean.purged.fa)
carry the MM/ML tags, so each assembly can serve as a single shared FRAME holding both ecotypes'
reads on identical coordinates:

    lean frame      lean reads on lean (23.8)       vs  siscowet reads on lean (23.4 + modkit)
    siscowet frame  siscowet reads on siscowet (23.8) vs lean reads on siscowet (23.4 + modkit)

No gene-ID pairing is needed inside a frame. The frame tests have the opposite weakness (the
visiting ecotype's reads carry divergence / mapping bias), so the two designs fail in different ways.
A feature is called ROBUST when it is FDR-significant in 23.11 and, in BOTH frames, has the same sign,
p < 0.05 and |diff| >= 10 points. CpG content is also compared between the two assemblies
(|log2 ratio| of median covered CpGs) as a homology flag.

Robust promoter / body genes then go to GO over-representation (18.3 machinery) against the genes
whose features were testable, with 23.7 tandem-cluster collapsing, and the lipid / calcium
hypotheses are read off as in 23.7.

Inputs (raven): 23.11 outputs; analyses/23-native-methylation/meth_cross/{s}.on_{frame}.CpG.bedmethyl
  (modkit 0.5.0 pileup --cpg --combine-strands on analyses/23-reciprocal-pav/align/*.bam; see
  meth_cross/run_modkit_cross.sh)

Outputs (analyses/23-native-methylation/):
  dmf_frame_{lean,siscowet}.tsv.gz       per-feature tests inside each frame
  dmf_robust.tsv                         23.11 features with both frame tests + robust flag
  ref_dmr_frame_tests.tsv                reference-DMR loci tested inside each frame
  go_robust_{promoter,body,all}.tsv      GO ORA (collapsed) of robust genes
  native_meth_frame_report.txt

Author: Generated for project-lake-trout
Date: 2026-10-14
"""

import importlib.util
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "analyses" / "23-native-methylation"
CROSS = OUT / "meth_cross"


def load_module(fname, modname):
    spec = importlib.util.spec_from_file_location(modname, BASE / "code" / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m11 = load_module("23.11-native-meth-counts.py", "native_meth_counts_2311")

FRAME_OWN = {"lean": m11.SAMPLES_LEAN, "siscowet": m11.SAMPLES_SISCOWET}
ROBUST_P = 0.05
CPG_RATIO_MAX = 1.0            # |log2| of median covered CpGs, siscowet vs lean


def cross_counts(sample, frame):
    """All-coverage 5mC counts for a visiting sample on the other ecotype's assembly."""
    dst = m11.COUNTS / f"{sample}.on_{frame}.cpg_counts_all.sorted.bed"
    if nonempty(dst):
        return dst
    src = CROSS / f"{sample}.on_{frame}.CpG.bedmethyl"
    if not nonempty(src):
        sys.exit(f"FATAL: {src} missing -- run meth_cross/run_modkit_cross.sh first")
    m11.sh(f"awk -v OFS='\\t' '$4==\"m\" && $10>=1{{print $1,$2,$3,$12,$10}}' '{src}' "
           f"| LC_ALL=C sort -S 8G --parallel=8 -k1,1 -k2,2n > '{dst}.tmp' && mv '{dst}.tmp' '{dst}'")
    return dst


def nonempty(p):
    p = Path(p)
    return p.exists() and p.stat().st_size > 0


def frame_test(frame):
    """Both ecotypes' reads on one assembly; features from that assembly's 23.11 feature BED."""
    fbed = m11.COUNTS / f"{frame}.features.bed"
    per_sample = {}
    for s in m11.SAMPLES:
        if s in FRAME_OWN[frame]:
            counts = m11.COUNTS / f"{s}.cpg_counts_all.sorted.bed"
            tsv = m11.COUNTS / f"{s}.features.tsv"
        else:
            counts = cross_counts(s, frame)
            tsv = m11.COUNTS / f"{s}.on_{frame}.features.tsv"
        per_sample[s] = m11.map_counts(fbed, counts, tsv).set_index("name")
    names = sorted(set.intersection(*(set(d.index) for d in per_sample.values())))
    arr = {c: np.column_stack([per_sample[s].loc[names, c].to_numpy(float) for s in m11.SAMPLES])
           for c in ("nmod", "nvalid", "ncpg")}
    res, testable = m11.test_matrix(arr["nmod"], arr["nvalid"], arr["ncpg"])
    res.insert(0, "feature", names)
    res["region"] = [n.split("|")[1] for n in names]
    res = res[testable].reset_index(drop=True)
    for region in ("body", "promoter"):
        k = res["region"] == region
        res.loc[k, "qvalue"] = multipletests(res.loc[k, "pvalue"], method="fdr_bh")[1]
    res.to_csv(OUT / f"dmf_frame_{frame}.tsv.gz", sep="\t", index=False, compression="gzip",
               float_format="%.6g")
    print(f"  [{frame} frame] tested {len(res):,}  FDR<0.1 & |diff|>=10: "
          f"{int(((res['qvalue'] < 0.1) & (res['meth_diff'].abs() >= 10)).sum()):,}", flush=True)
    return res.set_index("feature")


def refdmr_frame_test(frame):
    """The 23.11 reference-DMR loci placed on one assembly, tested with both ecotypes' reads."""
    bed = m11.COUNTS / f"refdmr.on_{frame}.bed"
    per_sample = {}
    for s in m11.SAMPLES:
        if s in FRAME_OWN[frame]:
            counts, tsv = m11.COUNTS / f"{s}.cpg_counts_all.sorted.bed", m11.COUNTS / f"refdmr.{s}.tsv"
        else:
            counts, tsv = cross_counts(s, frame), m11.COUNTS / f"refdmr.{s}.on_{frame}.tsv"
        per_sample[s] = m11.map_counts(bed, counts, tsv).drop_duplicates("name").set_index("name")
    names = sorted(per_sample[m11.SAMPLES[0]].index)
    arr = {c: np.column_stack([per_sample[s][c].reindex(names).fillna(0).to_numpy(float)
                               for s in m11.SAMPLES]) for c in ("nmod", "nvalid", "ncpg")}
    res, testable = m11.test_matrix(arr["nmod"], arr["nvalid"], arr["ncpg"], min_cpgs=2)
    res.insert(0, "dmr_id", names)
    return res[testable].set_index("dmr_id")


def go_robust(robust, tested_genes):
    if len(robust) < 10:
        return [f"  {len(robust)} robust features -- too few for GO; skipped"]
    m7 = load_module("23.7-refine-evidence-model.py", "refine_evidence_237")
    go = load_module("18.3-go-enrichment.py", "go_enrichment_183")
    fn = {r["gene_id"]: r for r in m7.read_tsv(m7.ANN / "gene_function_table.tsv")}
    rep_of, _ = m7.build_clusters(fn)
    name, ns, parents, alt = go.parse_obo(go.OBO)
    _, anc_fn = go.build_ancestors(parents)
    gene_go = go.load_gene_go(name, alt, anc_fn)
    # background = genes with a testable feature, collapsed to tandem-cluster representatives
    rep_go = defaultdict(set)
    for g in tested_genes:
        if g in gene_go:
            rep_go[rep_of.get(g, g)] |= gene_go[g]
    t2g = defaultdict(set)
    for g, terms in rep_go.items():
        for t in terms:
            t2g[t].add(g)
    universe = set(rep_go)
    lines, pheno = [], []
    for key in ("promoter", "body", "all"):
        genes = set(robust["gene_id"] if key == "all" else
                    robust.loc[robust["region"] == key, "gene_id"])
        study = {rep_of.get(g, g) for g in genes} & universe
        rows, n, _ = go.ora(study, t2g, universe, name, ns)
        sig = [r for r in rows if r["fdr"] < 0.1]
        lines.append(f"  {key:8s} genes {len(genes):4d} -> clusters in background {n:4d}; "
                     f"GO terms FDR<0.1: {len(sig)}")
        for r in sig[:8]:
            lines.append(f"      {r['name'][:55]:55s} k={r['study_k']:3d} fold={r['fold_enrichment']:5.1f} "
                         f"FDR={r['fdr']:.1e}")
        for r in rows:
            grp = ("lipid" if m7.PHENO_LIPID.search(r["name"]) else
                   "calcium" if m7.PHENO_CALCIUM.search(r["name"]) else None)
            if grp and r["fdr"] < 0.25:
                pheno.append(f"      [{key}] {grp}: {r['name']} k={r['study_k']} FDR={r['fdr']:.2g}")
        go.write(rows, OUT / f"go_robust_{key}.tsv")      # formats fdr as text in place
    lines.append("  lipid / calcium terms at FDR<0.25:")
    lines += pheno or ["      none"]
    return lines


def main():
    for f in ("dmf_results.tsv.gz", "native_feature_counts.tsv.gz"):
        if not nonempty(OUT / f):
            sys.exit(f"FATAL: {f} missing -- run 23.11 first")
    harm = pd.read_csv(OUT / "dmf_results.tsv.gz", sep="\t")
    counts = pd.read_csv(OUT / "native_feature_counts.tsv.gz", sep="\t").set_index("feature")

    lean_cpg = counts[[f"ncpg_{s}" for s in m11.SAMPLES_LEAN]].median(axis=1)
    sisc_cpg = counts[[f"ncpg_{s}" for s in m11.SAMPLES_SISCOWET]].median(axis=1)
    harm["log2_cpg_ratio"] = np.log2((sisc_cpg.reindex(harm["feature"]).to_numpy() + 1) /
                                     (lean_cpg.reindex(harm["feature"]).to_numpy() + 1))

    frames = {f: frame_test(f) for f in ("lean", "siscowet")}
    for f, r in frames.items():
        for c in ("meth_diff", "pvalue", "qvalue"):
            harm[f"{f}_frame_{c}"] = r[c].reindex(harm["feature"]).to_numpy()

    def frame_ok(f):
        d, p = harm[f"{f}_frame_meth_diff"], harm[f"{f}_frame_pvalue"]
        return (np.sign(d) == np.sign(harm["meth_diff"])) & (p < ROBUST_P) & (d.abs() >= m11.DIFF_THRESHOLD)

    harm["lean_frame_ok"] = frame_ok("lean")
    harm["siscowet_frame_ok"] = frame_ok("siscowet")
    harm["cpg_homology_ok"] = harm["log2_cpg_ratio"].abs() <= CPG_RATIO_MAX
    harm["robust"] = harm["significant_fdr"] & harm["lean_frame_ok"] & harm["siscowet_frame_ok"]
    sig = harm[harm["significant_fdr"]].copy()
    sig.to_csv(OUT / "dmf_robust.tsv", sep="\t", index=False, float_format="%.6g")
    robust = sig[sig["robust"]]

    lines = ["Native methylation frame control (23.12)", "=" * 42, ""]
    for region in ("promoter", "body"):
        s = sig[sig["region"] == region]
        lines += [f"  {region}: 23.11 FDR hits {len(s):,}",
                  f"    testable in lean frame {s['lean_frame_pvalue'].notna().sum():,}, "
                  f"siscowet frame {s['siscowet_frame_pvalue'].notna().sum():,}",
                  f"    same sign & p<0.05 & |diff|>=10 -- lean frame {int(s['lean_frame_ok'].sum()):,}, "
                  f"siscowet frame {int(s['siscowet_frame_ok'].sum()):,}, both (ROBUST) "
                  f"{int(s['robust'].sum()):,}",
                  f"    robust with CpG content within 2x: "
                  f"{int((s['robust'] & s['cpg_homology_ok']).sum()):,}  "
                  f"(hyper {int((s['robust'] & (s['meth_diff'] > 0)).sum())} / "
                  f"hypo {int((s['robust'] & (s['meth_diff'] < 0)).sum())} in siscowet)", ""]
        tested = s.dropna(subset=["lean_frame_meth_diff", "siscowet_frame_meth_diff"])
        if len(tested) > 2:
            lines.append(f"    Spearman 23.11 diff vs lean frame {tested['meth_diff'].corr(tested['lean_frame_meth_diff'], method='spearman'):.2f}, "
                         f"vs siscowet frame {tested['meth_diff'].corr(tested['siscowet_frame_meth_diff'], method='spearman'):.2f}")
            lines.append("")
    for f, r in frames.items():
        pi0 = min(1.0, (r["pvalue"] > 0.5).mean() / 0.5)
        lines.append(f"  {f} frame: {len(r):,} features tested; p<0.05 in "
                     f"{(r['pvalue'] < 0.05).mean():.1%} (5 % under the null); pi0 estimate {pi0:.2f}; "
                     f"min q {r['qvalue'].min():.2f}; |diff| median {r['meth_diff'].abs().median():.1f}, "
                     f"99th pct {r['meth_diff'].abs().quantile(0.99):.1f}")
    lines.append("")

    loci = pd.read_csv(OUT / "ref_dmr_native_loci.tsv", sep="\t").set_index("dmr_id")
    lines.append("Reference DMRs (14-diff-meth) inside each frame")
    rows = []
    for f in ("lean", "siscowet"):
        r = refdmr_frame_test(f)
        j = loci[["ref_meth_diff"]].join(r[["meth_diff", "pvalue"]], how="inner")
        j["frame"] = f
        rows.append(j.reset_index())
        agree = (np.sign(j["meth_diff"]) == np.sign(j["ref_meth_diff"]))
        conf = agree & (j["pvalue"] < ROBUST_P) & (j["meth_diff"].abs() >= m11.DIFF_THRESHOLD)
        rho = j["ref_meth_diff"].corr(j["meth_diff"], method="spearman") if len(j) > 2 else np.nan
        lines.append(f"  {f} frame: testable {len(j)}; sign agrees {int(agree.sum())}/{len(j)}; "
                     f"Spearman rho {rho:.2f}; same sign & p<0.05 & |diff|>=10: {int(conf.sum())}")
    pd.concat(rows).to_csv(OUT / "ref_dmr_frame_tests.tsv", sep="\t", index=False, float_format="%.6g")
    lines.append("")
    tested_genes = set(harm["gene_id"])
    lines += ["GO over-representation of robust genes (background: genes with a testable feature, "
              "tandem clusters collapsed)"]
    lines += go_robust(robust, tested_genes)
    lines += ["", "Top 30 robust features by 23.11 p:"]
    for r in robust.sort_values("pvalue").head(30).itertuples():
        lines.append(f"  {r.gene_id:24s} {r.region:8s} diff={r.meth_diff:+6.1f} "
                     f"(lean frame {r.lean_frame_meth_diff:+6.1f}, sisc frame {r.siscowet_frame_meth_diff:+6.1f}) "
                     f"{str(r.symbol)[:12]:12s} {str(r.product)[:40]}")
    text = "\n".join(lines) + "\n"
    (OUT / "native_meth_frame_report.txt").write_text(text)
    print("\n" + text)


if __name__ == "__main__":
    main()
