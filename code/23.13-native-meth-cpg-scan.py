#!/usr/bin/env python3
"""
23.13-native-meth-cpg-scan.py  (Phase 3, part E, of code/23-next-phase-research-plan.md)

CpG-LEVEL DMR SCAN INSIDE EACH SHARED FRAME. 23.12 showed that the only valid cross-ecotype
comparison puts both ecotypes' reads on one assembly (a "frame") and repeats it in both directions.
At the promoter / gene-body level nothing differed, but 2 kb features can dilute a localized
difference, and the reference DMRs that place uniquely did replicate. This script repeats the
14-diff-meth design (DMCs merged into DMRs) inside each frame, with a calibrated error rate.

Per frame (lean assembly, siscowet assembly):
  1. CpG matrix: the eight samples' all-coverage 5mC counts on the frame's coordinates (own
     ecotype from 23.11 counts/, visiting ecotype from 23.12 cross counts). A CpG is tested when
     >= MIN_SAMPLES_COV of 8 samples have Nvalid >= MIN_COV. The rule ignores ecotype labels, so
     permuted runs test exactly the same CpGs.
  2. Per-CpG quasi-binomial test (23.11 `quasibinomial_two_group`) and the difference in mean
     per-sample % methylation (samples with Nvalid >= MIN_COV).
  3. DMRs, as in 14-diff-meth: DMCs (p < p_thr, |diff| >= diff_thr) on the same contig with the same
     direction and <= MERGE_GAP bp apart are chained; chains with >= min_cpgs DMCs are DMRs. Region
     statistics (pooled counts of every tested CpG in the DMR) are descriptive only, because the
     region was selected on the same data.
  4. Error rate by permutation: the whole calling step is rerun on the 18 fully balanced relabelings
     (each new group = 2 lean + 2 siscowet fish). Empirical FDR = mean permuted DMR count / observed
     count, over a small grid of (p_thr, diff_thr, min_cpgs). The primary setting is the one with
     the most observed DMRs at empirical FDR <= EFDR_MAX in BOTH frames (14-diff-meth recipe is
     reported regardless).
  5. Cross-frame replication: each primary DMR (+/- 500 bp) is placed on the other assembly
     (23.11 `place_intervals`) and re-tested there with pooled counts from all eight fish. A DMR
     replicates when the placed locus has the same sign, p < 0.05 and |diff| >= 10 points.
  6. Annotation: nearest Liftoff gene on the frame assembly (promoter = within 2 kb upstream of the
     TSS), overlap with the placed reference DMRs (23.11 refdmr.on_{frame}.bed), and GO
     over-representation of replicated-DMR genes against genes with >= 3 tested CpGs (23.12
     `go_robust`, tandem clusters collapsed).

Inputs (raven): analyses/23-native-methylation/counts/{s}.cpg_counts_all.sorted.bed and
  {s}.on_{frame}.cpg_counts_all.sorted.bed (from 23.11 / 23.12), purged assemblies + Liftoff BEDs.

Outputs (analyses/23-native-methylation/cpg_scan/):
  dmr_scan_grid.tsv          per frame x setting: observed DMRs, permuted mean / max, empirical FDR
  dmrs_{lean,siscowet}.tsv   primary-setting DMRs with region stats, gene, ref-DMR overlap, and the
                             test of their placed locus in the other frame
  dmrs_replicated.tsv        DMRs replicated in the other frame
  go_cpgscan_{promoter,body,all}.tsv
  cpg_scan_report.txt
  (counts/frame_{frame}.cpg_matrix.npz caches the filtered matrices; gitignored)

Run on raven with /home/sr320/miniconda3/bin/python. Memory: ~30-40 M CpGs x 8 samples per frame,
tens of GB peak; 19 calling passes per frame.

Author: Generated for project-lake-trout
Date: 2026-10-14
"""

import importlib.util
import itertools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent


def load_module(fname, modname):
    spec = importlib.util.spec_from_file_location(modname, BASE / "code" / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m11 = load_module("23.11-native-meth-counts.py", "native_meth_counts_2311")
m12 = load_module("23.12-native-meth-frame-control.py", "native_meth_frame_2312")

OUT = m11.OUT / "cpg_scan"
COUNTS = m11.COUNTS
SAMPLES = m11.SAMPLES
IS_SISC = np.array([s in m11.SAMPLES_SISCOWET for s in SAMPLES])
FRAMES = ("lean", "siscowet")

MIN_COV = 5               # reads per sample for a CpG value
MIN_SAMPLES_COV = 7       # of 8 samples with >= MIN_COV (label-free filter)
MERGE_GAP = 250           # bp between consecutive DMCs (14-diff-meth)
GRID = [                  # (p_thr, diff_thr, min_cpgs); first row = 14-diff-meth recipe
    (0.05, 15, 3), (0.05, 15, 5), (0.01, 15, 3), (0.01, 25, 3), (0.01, 25, 5), (0.001, 25, 3),
]
EFDR_MAX = 0.10
REPL_P, REPL_DIFF = 0.05, 10.0
CHUNK = 4_000_000


# =============================================================================
# CpG matrix per frame
# =============================================================================
def sample_counts(sample, frame):
    own = sample in m11.ECO_SAMPLES[frame]
    return COUNTS / (f"{sample}.cpg_counts_all.sorted.bed" if own
                     else f"{sample}.on_{frame}.cpg_counts_all.sorted.bed")


def load_frame(frame):
    """Filtered CpG x sample matrices for one frame, sorted by (contig index, position)."""
    cache = COUNTS / f"frame_{frame}.cpg_matrix.npz"
    contigs = pd.read_csv(f"{m11.ECO_FA[frame]}.fai", sep="\t", header=None, usecols=[0])[0].tolist()
    if cache.exists():
        z = np.load(cache)
        return contigs, z["key"], z["nmod"], z["nvalid"]
    cidx = {c: i for i, c in enumerate(contigs)}
    per = []
    for s in SAMPLES:
        path = sample_counts(s, frame)
        if not m11.nonempty(path):
            sys.exit(f"FATAL: {path} missing -- run 23.11 and 23.12 first")
        print(f"  [{frame}] reading {path.name}", flush=True)
        df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 3, 4],
                         names=["chrom", "pos", "nmod", "nvalid"],
                         dtype={"chrom": "category", "pos": np.int64,
                                "nmod": np.int32, "nvalid": np.int32})
        cat_idx = pd.Series(df["chrom"].cat.categories).map(cidx).fillna(-1).astype(np.int64)
        codes = cat_idx.to_numpy()[df["chrom"].cat.codes.to_numpy()]
        key = (codes.astype(np.int64) << 32) | df["pos"].to_numpy()
        per.append((key, df["nmod"].to_numpy(), df["nvalid"].to_numpy()))
    union = np.unique(np.concatenate([k for k, _, _ in per]))
    nmod = np.zeros((len(union), len(SAMPLES)), dtype=np.int32)
    nvalid = np.zeros_like(nmod)
    for j, (k, m, v) in enumerate(per):
        idx = np.searchsorted(union, k)
        nmod[idx, j], nvalid[idx, j] = m, v
    keep = (nvalid >= MIN_COV).sum(1) >= MIN_SAMPLES_COV
    key, nmod, nvalid = union[keep], nmod[keep], nvalid[keep]
    print(f"  [{frame}] CpGs in union {len(union):,}; tested (>= {MIN_COV} reads in >= "
          f"{MIN_SAMPLES_COV}/8 fish) {len(key):,}", flush=True)
    np.savez(cache, key=key, nmod=nmod, nvalid=nvalid)
    return contigs, key, nmod, nvalid


# =============================================================================
# Tests and DMR calling
# =============================================================================
def cpg_tests(nmod, nvalid, is_sisc):
    """Per-CpG quasi-binomial p and mean % difference (siscowet - lean), chunked."""
    n = len(nmod)
    p = np.empty(n)
    diff = np.empty(n, dtype=np.float32)
    for a in range(0, n, CHUNK):
        b = min(n, a + CHUNK)
        ok = nvalid[a:b] >= MIN_COV
        nm = np.where(ok, nmod[a:b], np.nan)
        nv = np.where(ok, nvalid[a:b], np.nan)
        _, _, _, pp, _, _ = m11.quasibinomial_two_group(nm, nv, is_sisc)
        pct = 100.0 * nm / nv
        with np.errstate(all="ignore"):
            diff[a:b] = np.nanmean(pct[:, is_sisc], 1) - np.nanmean(pct[:, ~is_sisc], 1)
        p[a:b] = pp
    return p, diff


def call_dmrs(key, p, diff, p_thr, diff_thr, min_cpgs):
    """Chain DMCs (same contig, same direction, gap <= MERGE_GAP). Returns arrays of
    (first DMC row, last DMC row, n DMCs) into the tested-CpG arrays."""
    dmc = np.flatnonzero((p < p_thr) & (np.abs(diff) >= diff_thr))
    if len(dmc) == 0:
        return np.empty((0, 3), dtype=np.int64)
    k = key[dmc]
    chrom, pos = k >> 32, k & 0xFFFFFFFF
    sign = np.sign(diff[dmc])
    brk = np.ones(len(dmc), dtype=bool)
    brk[1:] = (chrom[1:] != chrom[:-1]) | (sign[1:] != sign[:-1]) | (pos[1:] - pos[:-1] > MERGE_GAP)
    starts = np.flatnonzero(brk)
    ends = np.append(starts[1:], len(dmc)) - 1
    sizes = ends - starts + 1
    good = sizes >= min_cpgs
    return np.column_stack([dmc[starts[good]], dmc[ends[good]], sizes[good]])


def balanced_permutations():
    """The 18 relabelings with exactly two fish of each ecotype in each group."""
    lean = np.flatnonzero(~IS_SISC)
    sisc = np.flatnonzero(IS_SISC)
    perms, seen = [], set()
    for a in itertools.combinations(lean, 2):
        for b in itertools.combinations(sisc, 2):
            g = np.zeros(len(SAMPLES), dtype=bool)
            g[list(a) + list(b)] = True
            canon = tuple(g) if g[0] else tuple(~g)
            if canon not in seen:
                seen.add(canon)
                perms.append(g)
    return perms


def region_stats(key, nmod, nvalid, intervals, contigs, min_cpgs=1):
    """Pool counts of every tested CpG inside each (contig, start, end) interval; test.
    Returns the 23.11 test_matrix frame (+ n_cpgs_tested) aligned to the intervals."""
    cidx = {c: i for i, c in enumerate(contigs)}
    lo = np.array([(cidx.get(c, -1) << 32) | s for c, s, _ in intervals], dtype=np.int64)
    hi = np.array([(cidx.get(c, -1) << 32) | e for c, _, e in intervals], dtype=np.int64)
    i0, i1 = np.searchsorted(key, lo, "left"), np.searchsorted(key, hi, "left")
    cm = np.vstack([np.zeros((1, nmod.shape[1])), np.cumsum(nmod, 0, dtype=np.float64)])
    cv = np.vstack([np.zeros((1, nvalid.shape[1])), np.cumsum(nvalid, 0, dtype=np.float64)])
    cc = np.vstack([np.zeros((1, nvalid.shape[1])),
                    np.cumsum(nvalid >= MIN_COV, 0, dtype=np.float64)])
    rm, rv, rc = cm[i1] - cm[i0], cv[i1] - cv[i0], cc[i1] - cc[i0]
    saved = m11.MIN_NVALID
    m11.MIN_NVALID = MIN_COV                       # region totals, not 2 kb features
    res, testable = m11.test_matrix(rm, rv, rc, min_cpgs=min_cpgs)
    m11.MIN_NVALID = saved
    res["n_cpgs_tested"] = i1 - i0
    res["testable"] = testable
    return res


# =============================================================================
# Annotation helpers
# =============================================================================
def nearest_gene(frame, dmrs):
    """bedtools closest against the frame's Liftoff genes; promoter / body / distance."""
    tmp_a = OUT / f"_{frame}.dmrs.bed"
    tmp_g = OUT / f"_{frame}.genes.bed"
    dmrs[["contig", "start", "end", "dmr_id"]].sort_values(["contig", "start"]).to_csv(
        tmp_a, sep="\t", header=False, index=False)
    m11.sh(f"LC_ALL=C sort -k1,1 -k2,2n '{m11.ECO_GENES[frame]}' > '{tmp_g}'")
    out = OUT / f"_{frame}.closest.tsv"
    m11.sh(f"LC_ALL=C sort -k1,1 -k2,2n '{tmp_a}' | '{m11.BEDTOOLS}' closest -a - -b '{tmp_g}' "
           f"-D b -t first > '{out}'")
    cl = pd.read_csv(out, sep="\t", header=None,
                     names=["c", "s", "e", "dmr_id", "gc", "gs", "ge", "gene_id_raw", "gsc",
                            "strand", "distance"])
    cl["gene_id"] = cl["gene_id_raw"].str.replace(m11.COPY_SUFFIX, "", regex=True)
    # -D b: negative = upstream of the gene
    cl["region"] = np.where(cl["distance"] == 0, "body",
                            np.where((cl["distance"] < 0) & (cl["distance"] >= -m11.PROMOTER_BP),
                                     "promoter", "intergenic"))
    for f in (tmp_a, tmp_g, out):
        f.unlink(missing_ok=True)
    return cl.drop_duplicates("dmr_id").set_index("dmr_id")[["gene_id", "distance", "region"]]


def tested_gene_background(frame, key, contigs):
    """Genes (canonical IDs) with >= 3 tested CpGs in promoter + body on the frame assembly."""
    genes = pd.read_csv(m11.ECO_GENES[frame], sep="\t", header=None,
                        names=["chrom", "start", "end", "gene_id", "score", "strand"])
    s = np.where(genes["strand"] == "-", genes["start"], np.maximum(0, genes["start"] - m11.PROMOTER_BP))
    e = np.where(genes["strand"] == "-", genes["end"] + m11.PROMOTER_BP, genes["end"])
    cidx = {c: i for i, c in enumerate(contigs)}
    ci = genes["chrom"].map(cidx).fillna(-1).astype(np.int64).to_numpy()
    n = np.searchsorted(key, (ci << 32) | e) - np.searchsorted(key, (ci << 32) | s)
    ids = genes.loc[n >= 3, "gene_id"].str.replace(m11.COPY_SUFFIX, "", regex=True)
    return set(ids)


def overlaps_ref_dmr(frame, dmrs):
    bed = COUNTS / f"refdmr.on_{frame}.bed"
    if not m11.nonempty(bed):
        return pd.Series("", index=dmrs["dmr_id"])
    ref = pd.read_csv(bed, sep="\t", header=None, names=["c", "s", "e", "ref_dmr"])
    hits = {}
    for r in ref.itertuples(index=False):
        m = (dmrs["contig"] == r.c) & (dmrs["start"] < r.e) & (dmrs["end"] > r.s)
        for d in dmrs.loc[m, "dmr_id"]:
            hits.setdefault(d, []).append(r.ref_dmr)
    return pd.Series({d: ",".join(hits.get(d, [])) for d in dmrs["dmr_id"]})


# =============================================================================
# Main
# =============================================================================
def scan_frame(frame):
    contigs, key, nmod, nvalid = load_frame(frame)
    p, diff = cpg_tests(nmod, nvalid, IS_SISC)
    observed = {g: call_dmrs(key, p, diff, *g) for g in GRID}
    perm_counts = {g: [] for g in GRID}
    for i, g_lab in enumerate(balanced_permutations()):
        pp, dd = cpg_tests(nmod, nvalid, g_lab)
        for g in GRID:
            perm_counts[g].append(len(call_dmrs(key, pp, dd, *g)))
        print(f"  [{frame}] permutation {i + 1}/18 done", flush=True)
    grid = []
    for g in GRID:
        obs, pc = len(observed[g]), np.array(perm_counts[g])
        grid.append({"frame": frame, "p_thr": g[0], "diff_thr": g[1], "min_cpgs": g[2],
                     "observed": obs, "perm_mean": pc.mean(), "perm_max": int(pc.max()),
                     "efdr": pc.mean() / obs if obs else np.nan})
    n_tested = len(key)
    return {"contigs": contigs, "key": key, "nmod": nmod, "nvalid": nvalid,
            "observed": observed, "grid": pd.DataFrame(grid), "n_tested": n_tested,
            "p_frac": float((p < 0.05).mean())}


def dmr_table(frame, fr, setting):
    rows = fr["observed"][setting]
    if len(rows) == 0:
        return pd.DataFrame(columns=["dmr_id", "contig", "start", "end", "n_dmcs", "n_cpgs_tested",
                                     "lean_mean", "siscowet_mean", "meth_diff", "direction"])
    key, contigs = fr["key"], fr["contigs"]
    ci = (key[rows[:, 0]] >> 32).astype(int)
    df = pd.DataFrame({
        "dmr_id": [f"{frame[0].upper()}DMR_{i + 1}" for i in range(len(rows))],
        "contig": [contigs[c] for c in ci],
        "start": (key[rows[:, 0]] & 0xFFFFFFFF).astype(np.int64),
        "end": (key[rows[:, 1]] & 0xFFFFFFFF).astype(np.int64) + 1,
        "n_dmcs": rows[:, 2],
    })
    st = region_stats(key, fr["nmod"], fr["nvalid"],
                      list(zip(df["contig"], df["start"], df["end"])), contigs)
    df["n_cpgs_tested"] = st["n_cpgs_tested"].to_numpy()
    df["lean_mean"] = st["lean_mean"].round(2).to_numpy()
    df["siscowet_mean"] = st["siscowet_mean"].round(2).to_numpy()
    df["meth_diff"] = st["meth_diff"].round(2).to_numpy()
    df["direction"] = np.where(df["meth_diff"] > 0, "hyper_siscowet", "hypo_siscowet")
    return df


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    frames = {}
    for f in FRAMES:
        print("=" * 70 + f"\n{f} frame\n" + "=" * 70, flush=True)
        frames[f] = scan_frame(f)
    grid = pd.concat([frames[f]["grid"] for f in FRAMES])
    grid.to_csv(OUT / "dmr_scan_grid.tsv", sep="\t", index=False, float_format="%.4g")

    # primary setting: most DMRs (summed over frames) with eFDR <= EFDR_MAX in both frames
    wide = grid.pivot_table(index=["p_thr", "diff_thr", "min_cpgs"], columns="frame",
                            values=["observed", "efdr"])
    ok = wide[(wide[("efdr", "lean")] <= EFDR_MAX) & (wide[("efdr", "siscowet")] <= EFDR_MAX)]
    if len(ok):
        primary = (ok[("observed", "lean")] + ok[("observed", "siscowet")]).idxmax()
        controlled = True
    else:
        primary = (wide[("efdr", "lean")] + wide[("efdr", "siscowet")]).idxmin()
        controlled = False
    print(f"  primary setting {primary} (eFDR <= {EFDR_MAX} in both frames: {controlled})")

    tables = {f: dmr_table(f, frames[f], primary) for f in FRAMES}
    # cross-frame replication
    for f, other in (("lean", "siscowet"), ("siscowet", "lean")):
        df = tables[f]
        placed = m11.place_intervals(
            m11.ECO_FA[f], list(zip(df["dmr_id"], df["contig"], df["start"], df["end"])),
            {other: m11.ECO_FA[other]}, f"cpgscan_{f}") if len(df) else {}
        loci = [placed.get((d, other)) for d in df["dmr_id"]]
        df[f"{other}_locus"] = [f"{x[0]}:{x[1]}-{x[2]}" if x else "" for x in loci]
        has = [x is not None for x in loci]
        o = frames[other]
        st = region_stats(o["key"], o["nmod"], o["nvalid"],
                          [(x[0], x[1], x[2]) for x in loci if x], o["contigs"], min_cpgs=2) \
            if any(has) else None
        for c in ("meth_diff", "pvalue", "n_cpgs_tested", "testable"):
            col = np.full(len(df), np.nan, dtype=object)
            if st is not None:
                col[np.flatnonzero(has)] = st[c].to_numpy()
            df[f"other_frame_{c}"] = col
        d_o = pd.to_numeric(df["other_frame_meth_diff"], errors="coerce")
        p_o = pd.to_numeric(df["other_frame_pvalue"], errors="coerce")
        df["replicated"] = (df["other_frame_testable"] == True) & \
                           (np.sign(d_o) == np.sign(df["meth_diff"])) & \
                           (p_o < REPL_P) & (d_o.abs() >= REPL_DIFF)
        if len(df):
            ann = nearest_gene(f, df)
            df = df.join(ann, on="dmr_id")
            df["ref_dmr_overlap"] = overlaps_ref_dmr(f, df).reindex(df["dmr_id"]).to_numpy()
            fn = pd.read_csv(m11.ANN / "gene_function_table.tsv", sep="\t", dtype=str,
                             usecols=["gene_id", "symbol", "product"]).set_index("gene_id")
            df["symbol"] = df["gene_id"].map(fn["symbol"]).fillna("")
            df["product"] = df["gene_id"].map(fn["product"]).fillna("")
        tables[f] = df
        df.to_csv(OUT / f"dmrs_{f}.tsv", sep="\t", index=False, float_format="%.6g")

    rep = pd.concat([tables[f].assign(frame=f) for f in FRAMES], ignore_index=True)
    rep = rep[rep["replicated"] == True]
    rep.to_csv(OUT / "dmrs_replicated.tsv", sep="\t", index=False, float_format="%.6g")

    # ---- report ----
    lines = ["Native CpG-level DMR scan inside each frame (23.13)", "=" * 52, "",
             f"CpG filter: >= {MIN_COV} reads in >= {MIN_SAMPLES_COV}/8 fish (label-free); "
             f"DMCs chained within {MERGE_GAP} bp, same direction", ""]
    for f in FRAMES:
        lines.append(f"  {f} frame: {frames[f]['n_tested']:,} CpGs tested; "
                     f"p<0.05 at {frames[f]['p_frac']:.1%} of CpGs")
    lines += ["", "Grid (18 balanced permutations; eFDR = permuted mean / observed):",
              "  frame     p_thr  diff  minCpG  observed  perm_mean  perm_max  eFDR"]
    for r in grid.itertuples():
        lines.append(f"  {r.frame:9s} {r.p_thr:<6g} {r.diff_thr:<5g} {r.min_cpgs:<7d} "
                     f"{r.observed:<9d} {r.perm_mean:<10.1f} {r.perm_max:<9d} {r.efdr:.2f}")
    lines += ["", f"Primary setting: p<{primary[0]}, |diff|>={primary[1]}, >={primary[2]} DMCs "
                  f"({'empirical FDR <= %.2f in both frames' % EFDR_MAX if controlled else 'NOT FDR-controlled: lowest eFDR in grid'})"]
    for f in FRAMES:
        df = tables[f]
        if not len(df):
            lines.append(f"  {f}: 0 DMRs")
            continue
        t = df[df["other_frame_testable"] == True]
        lines += [f"  {f}: {len(df)} DMRs (hyper {int((df['meth_diff'] > 0).sum())} / hypo "
                  f"{int((df['meth_diff'] < 0).sum())} in siscowet); median {df['n_dmcs'].median():.0f} "
                  f"DMCs, {(df['end'] - df['start']).median():.0f} bp",
                  f"    placed and testable in the other frame: {len(t)}; replicated "
                  f"(same sign, p<{REPL_P}, |diff|>={REPL_DIFF:.0f}): {int(df['replicated'].sum())}",
                  f"    region: " + ", ".join(f"{k} {v}" for k, v in df["region"].value_counts().items()),
                  f"    overlapping a placed reference DMR: {int((df['ref_dmr_overlap'] != '').sum())}"]
    lines += ["", "GO over-representation, genes of replicated DMRs (background: genes with >= 3 "
                  "tested CpGs in either frame; tandem clusters collapsed)"]
    background = set()
    for f in FRAMES:
        background |= tested_gene_background(f, frames[f]["key"], frames[f]["contigs"])
    if "gene_id" in rep and len(rep):
        go_in = rep.dropna(subset=["gene_id"]).drop_duplicates("gene_id")
        lines += m12.go_robust(go_in, background, prefix="go_cpgscan", out=OUT)
    else:
        lines.append("  no replicated DMRs -- GO skipped")
    lines += ["", "Replicated DMRs (by |diff|):"]
    for r in rep.reindex(rep["meth_diff"].abs().sort_values(ascending=False).index).head(40).itertuples():
        lines.append(f"  {r.frame:8s} {r.dmr_id:10s} {r.contig}:{r.start}-{r.end} n={r.n_dmcs:<3d} "
                     f"diff={r.meth_diff:+6.1f} other={float(r.other_frame_meth_diff):+6.1f}  "
                     f"{r.region:10s} {str(r.symbol)[:12]:12s} {str(r.product)[:38]}")
    text = "\n".join(lines) + "\n"
    (OUT / "cpg_scan_report.txt").write_text(text)
    print("\n" + text)


if __name__ == "__main__":
    main()
