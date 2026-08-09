#!/usr/bin/env python3
"""
14.1-within-vs-across-meth.py

Direct quantification of DNA methylation variation WITHIN ecotype vs ACROSS ecotype.

Input: per-sample CpG methylation (percent) files from 10-mCG-call:
       analyses/10-mCG-call/{sample}.cov10.bedGraph  (chrom, start, end, pct_meth)
       already coverage-filtered (>=10x).

Approach:
  1. Build a CpG x sample matrix of percent-methylation, keeping only CpG sites
     covered in ALL 8 samples (complete cases) so every pairwise comparison uses
     the identical set of sites.
  2. For every pair of samples compute:
       - mean absolute methylation difference (|%meth_i - %meth_j|)  -> "distance"
       - Pearson & Spearman correlation
       - Euclidean distance
  3. Split the 28 pairs into WITHIN-Lean, WITHIN-Siscowet, and ACROSS (Lean-Siscowet).
  4. Report within vs across summaries + a within/across ratio, and hierarchical
     clustering / PCA to see whether individuals group by ecotype.

Outputs (analyses/14-diff-meth/):
  within_vs_across_summary.txt
  pairwise_sample_distances.csv
  sample_correlation_matrix.csv
"""

import sys
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
MDIR = BASE / "analyses" / "10-mCG-call"
ODIR = BASE / "analyses" / "14-diff-meth"

LEAN = ["bc2041", "bc2069", "bc2070", "bc2068"]
SISCOWET = ["bc2071", "bc2073", "bc2072", "bc2096"]
SAMPLES = LEAN + SISCOWET
GROUP = {s: "Lean" for s in LEAN}
GROUP.update({s: "Siscowet" for s in SISCOWET})


def load_sample(sample):
    """Return a Series indexed by 'chrom:pos' -> pct_meth for one sample."""
    fp = MDIR / f"{sample}.cov10.bedGraph"
    df = pd.read_csv(
        fp, sep="\t", header=None,
        names=["chrom", "start", "end", "pct"],
        usecols=[0, 1, 2, 3],
        dtype={"chrom": str, "start": np.int64, "pct": np.float64},
        comment="t",  # skip any 'track' header lines
    )
    idx = df["chrom"] + ":" + df["start"].astype(str)
    s = pd.Series(df["pct"].values, index=idx, name=sample)
    # collapse any accidental duplicate positions
    s = s[~s.index.duplicated(keep="first")]
    return s


def main():
    print("Loading per-sample methylation ...")
    series = {}
    for s in SAMPLES:
        series[s] = load_sample(s)
        print(f"  {s}: {len(series[s]):,} CpGs")

    print("Building complete-case matrix (CpGs covered in all 8 samples) ...")
    mat = pd.DataFrame(series)
    n_before = len(mat)
    mat = mat.dropna(how="any")
    print(f"  {n_before:,} union CpGs -> {len(mat):,} shared by all 8 samples")

    if len(mat) < 1000:
        sys.exit("Too few shared CpGs; check inputs.")

    X = mat[SAMPLES].values  # sites x samples
    n_sites = X.shape[0]

    # ---- pairwise metrics ----
    rows = []
    corr = pd.DataFrame(np.eye(len(SAMPLES)), index=SAMPLES, columns=SAMPLES)
    for a, b in combinations(SAMPLES, 2):
        va, vb = mat[a].values, mat[b].values
        mad = np.mean(np.abs(va - vb))          # mean abs % difference
        eucl = np.sqrt(np.mean((va - vb) ** 2))  # RMS % difference
        pear = np.corrcoef(va, vb)[0, 1]
        ga, gb = GROUP[a], GROUP[b]
        if ga == gb:
            cls = f"within_{ga}"
        else:
            cls = "across"
        rows.append(dict(sample_a=a, sample_b=b, group_a=ga, group_b=gb,
                         comparison=cls, mean_abs_diff=mad, rms_diff=eucl,
                         pearson_r=pear))
        corr.loc[a, b] = corr.loc[b, a] = pear

    pw = pd.DataFrame(rows)
    pw.to_csv(ODIR / "pairwise_sample_distances.csv", index=False)
    corr.to_csv(ODIR / "sample_correlation_matrix.csv")

    within = pw[pw.comparison.str.startswith("within")]
    across = pw[pw.comparison == "across"]

    def summ(d, col):
        return d[col].mean(), d[col].std(ddof=1), d[col].min(), d[col].max()

    lines = []
    def out(txt=""):
        print(txt)
        lines.append(txt)

    out("=" * 64)
    out("WITHIN-ECOTYPE vs ACROSS-ECOTYPE DNA METHYLATION VARIATION")
    out("=" * 64)
    out(f"Shared CpG sites (covered >=10x in all 8 samples): {n_sites:,}")
    out(f"Samples: Lean={LEAN}  Siscowet={SISCOWET}")
    out("")
    out(f"Pair counts: within-Lean={ (pw.comparison=='within_Lean').sum() }, "
        f"within-Siscowet={ (pw.comparison=='within_Siscowet').sum() }, "
        f"across={ len(across) }")
    out("")

    for col, label, unit in [
        ("mean_abs_diff", "Mean absolute methylation difference", "% points"),
        ("pearson_r", "Genome-wide Pearson correlation", "r"),
    ]:
        wm, ws, wlo, whi = summ(within, col)
        am, a_s, alo, ahi = summ(across, col)
        out(f"{label} ({unit}):")
        out(f"    WITHIN ecotype : mean={wm:.3f}  sd={ws:.3f}  range[{wlo:.3f},{whi:.3f}]")
        out(f"    ACROSS ecotype : mean={am:.3f}  sd={a_s:.3f}  range[{alo:.3f},{ahi:.3f}]")
        if col == "mean_abs_diff":
            out(f"    ACROSS/WITHIN ratio (distance): {am/wm:.3f}x")
        out("")

    # per-group within breakdown
    out("Within-group detail (mean_abs_diff, % points):")
    for g in ["within_Lean", "within_Siscowet"]:
        d = pw[pw.comparison == g]
        out(f"    {g}: mean={d.mean_abs_diff.mean():.3f}  "
            f"range[{d.mean_abs_diff.min():.3f},{d.mean_abs_diff.max():.3f}]")
    out("")

    # Welch t-test across vs within (distance)
    try:
        from scipy import stats
        t, p = stats.ttest_ind(across.mean_abs_diff, within.mean_abs_diff,
                               equal_var=False)
        out(f"Welch t-test (across vs within mean_abs_diff): t={t:.3f}, p={p:.4g}")
    except Exception as e:
        out(f"(scipy unavailable for t-test: {e})")
    out("")

    # ---- clustering: does each sample's nearest neighbour share its ecotype? ----
    out("Nearest-neighbour ecotype check (by mean_abs_diff):")
    concord = 0
    for s in SAMPLES:
        sub = pw[(pw.sample_a == s) | (pw.sample_b == s)].copy()
        sub["other"] = np.where(sub.sample_a == s, sub.sample_b, sub.sample_a)
        nn = sub.loc[sub.mean_abs_diff.idxmin()]
        same = GROUP[nn["other"]] == GROUP[s]
        concord += same
        out(f"    {s} ({GROUP[s]}): nearest = {nn['other']} ({GROUP[nn['other']]}) "
            f"{'[same]' if same else '[DIFFERENT]'}")
    out(f"  {concord}/{len(SAMPLES)} samples nearest-neighbour a same-ecotype fish")
    out("")

    # ---- PCA (top 2 PCs on samples) ----
    try:
        # center sites, PCA over samples
        Xc = X - X.mean(axis=1, keepdims=True)
        # sample x sample covariance via SVD on sites x samples
        U, S, Vt = np.linalg.svd(Xc.T @ Xc)  # samples x samples
        pcs = (Xc.T @ Xc @ U)  # not needed; use eigen of gram
        evals = S
        # coordinates: eigenvectors of gram matrix scaled
        coords = U[:, :2] * np.sqrt(S[:2])
        var_expl = S / S.sum()
        out(f"PCA on samples (PC1 {var_expl[0]*100:.1f}%, PC2 {var_expl[1]*100:.1f}% of variance):")
        for i, s in enumerate(SAMPLES):
            out(f"    {s} ({GROUP[s]}): PC1={coords[i,0]:+.1f}  PC2={coords[i,1]:+.1f}")
    except Exception as e:
        out(f"(PCA skipped: {e})")
    out("")

    out("INTERPRETATION")
    out("-" * 64)
    ratio = across.mean_abs_diff.mean() / within.mean_abs_diff.mean()
    if ratio < 1.15:
        out(f"Across/within distance ratio = {ratio:.2f}x: across-ecotype methylation")
        out("differences are comparable to variation among individuals of the SAME")
        out("ecotype -- consistent with 0 CpGs surviving FDR in the group contrast.")
    else:
        out(f"Across/within distance ratio = {ratio:.2f}x: across-ecotype differences")
        out("exceed within-ecotype variation, though FDR-level per-site signal is weak.")

    (ODIR / "within_vs_across_summary.txt").write_text("\n".join(lines) + "\n")
    print(f"\nWrote {ODIR/'within_vs_across_summary.txt'}")
    print(f"Wrote {ODIR/'pairwise_sample_distances.csv'}")
    print(f"Wrote {ODIR/'sample_correlation_matrix.csv'}")


if __name__ == "__main__":
    main()
