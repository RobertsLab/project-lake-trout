#!/usr/bin/env python3
"""
23.15-dmr-fixed-variants.py  (Phase 3 follow-up, code/23-next-phase-research-plan.md)

ARE THE NATIVE DMRs NEXT TO FIXED SEQUENCE DIFFERENCES? The 23.13 DMRs separate all four lean
from all four siscowet fish, which fits a stable epigenetic ecotype difference but equally fits
cis-genetic control (a fixed variant nearby sets the methylation state). This script asks whether
sequence differences that are fixed between the ecotypes in these eight fish are enriched near the
DMRs relative to matched control windows.

Per frame (lean assembly, siscowet assembly), using the same eight BAMs whose methylation defined
the frame (own ecotype: 23.8 align/{s}.to_own.sorted.bam; visiting ecotype: 23.4
align/{s}.to_{frame}.sorted.bam):

  1. Windows. Every primary 23.13 DMR, plus N_CONTROLS control windows per DMR: same length, a
     tested-CpG count within +/-30 % (min +/-1) of the DMR's (CpG-dense sequence mutates faster),
     on the same frame, >= EXCLUDE_BP from any DMR, with room for the 5 kb flanks. Control centres
     are drawn from the tested CpGs of the 23.13 matrix (counts/frame_{frame}.cpg_matrix.npz).
  2. Variants. bcftools 1.14 mpileup -X pacbio-ccs -q 20 -a AD,DP over window +/- FLANK_BP, all 8
     BAMs (--ignore-RG, one sample per BAM), then call -mv; sites with QUAL >= 20 kept.
  3. Callable bases: samtools depth -a -Q 20, depth >= MIN_DP in all 8 fish.
  4. A site is FIXED when every fish has >= MIN_DP reads (AD sum) and all four fish of one ecotype
     have alt fraction <= 0.2 while all four of the other have >= 0.8. Each variant is binned by its
     distance to the window core: core (inside), near (<= 1 kb), far (1-5 kb).
     Two less strict classes are also counted, SNVs only: DIFF75 / DIFF50 = |mean alt fraction
     siscowet - mean alt fraction lean| >= 0.75 / 0.5 (all 8 fish >= MIN_DP).
  5. Statistics per frame, class and bin. Primary: the fraction of DMR windows with >= 1 site of
     the class vs matched controls, with a matched resampling p (10,000 draws of one control per
     DMR). Window-level, so one divergent haplotype block cannot dominate. The pooled rate per
     callable kb is reported alongside for context. Repeated for replicated DMRs and by direction.
  QC: the share of covered variant sites where all 8 fish look heterozygous (collapsed WGD
     paralogs in the purged assemblies); such sites sit near AF 0.5 and cannot be called fixed.
  6. Copy-number check. If the ecotypes differ in how many paralog copies collapse onto a locus,
     both allele fractions and methylation shift. Per window: mean Nvalid over its tested CpGs per
     fish (23.13 matrix), normalised by that fish's median, then log2(siscowet mean / lean mean).
     |log2| >= DEPTH_DEV flags a depth difference; DMR vs control by matched resampling, and
     within DMRs with vs without a differentiated SNV (|dAF| >= 0.5) within 1 kb.

Outputs (analyses/23-native-methylation/fixed_variants/):
  windows_{frame}.tsv            DMR + control windows (window_id, type, dmr_id, contig, start, end)
  window_counts_{frame}.tsv      per window x bin: callable bp, fixed SNV, fixed indel, all variants
  fixed_sites_{frame}.tsv        every fixed site with per-fish alt fractions and its window/bin
  window_depth_{frame}.tsv       per window: siscowet / lean depth ratio (log2), diff50 hit
  fixed_variant_report.txt
  (vcf/ and callable/ intermediates are gitignored)

Run on raven: /home/sr320/miniconda3/bin/python code/23.15-dmr-fixed-variants.py

Author: Generated for project-lake-trout
Date: 2026-10-15
"""

import importlib.util
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
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

OUT = m11.OUT / "fixed_variants"
SCAN = m11.OUT / "cpg_scan"
VCF = OUT / "vcf"
CALL = OUT / "callable"
BCFTOOLS = "/home/shared/bcftools-1.14/bcftools"
SAMTOOLS = m11.SAMTOOLS
BEDTOOLS = m11.BEDTOOLS
SAMPLES = m11.SAMPLES
IS_SISC = np.array([s in m11.SAMPLES_SISCOWET for s in SAMPLES])
FRAMES = ("lean", "siscowet")

N_CONTROLS = 10
EXCLUDE_BP = 20_000
FLANK_BP = 5_000
NEAR_BP = 1_000
MIN_DP = 4
MIN_QUAL = 20
MIN_MAPQ = 20
LO_AF, HI_AF = 0.2, 0.8
N_CHUNKS = 24
N_RESAMPLE = 10_000
SEED = 2315
DEPTH_DEV = 0.5           # |log2 siscowet / lean depth|, ~1.4-fold


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def bams(frame):
    out = []
    for s in SAMPLES:
        if s in m11.ECO_SAMPLES[frame]:
            out.append(m11.OUT / "align" / f"{s}.to_own.sorted.bam")
        else:
            out.append(BASE / "analyses" / "23-reciprocal-pav" / "align" / f"{s}.to_{frame}.sorted.bam")
    for b in out:
        if not b.exists():
            sys.exit(f"FATAL: {b} missing")
    return out


# =============================================================================
# 1. Windows
# =============================================================================
def build_windows(frame, rng):
    path = OUT / f"windows_{frame}.tsv"
    if path.exists():
        return pd.read_csv(path, sep="\t")
    dmrs = pd.read_csv(SCAN / f"dmrs_{frame}.tsv", sep="\t")
    dmrs["status"] = np.where(dmrs["replicated"].astype(str) == "True", "replicated",
                              np.where(dmrs["other_frame_testable"].astype(str) == "True",
                                       "not_replicated", "single_frame"))
    z = np.load(m11.COUNTS / f"frame_{frame}.cpg_matrix.npz")
    key = z["key"]
    fai = pd.read_csv(f"{m11.ECO_FA[frame]}.fai", sep="\t", header=None, usecols=[0, 1])
    contigs, sizes = fai[0].tolist(), fai[1].to_numpy()
    cidx = {c: i for i, c in enumerate(contigs)}

    # DMR exclusion zones, as key intervals
    dk = np.array([(cidx[c] << 32) | s for c, s in zip(dmrs["contig"], dmrs["start"])])
    dk_sorted = np.sort(dk)

    # candidate control starts: random tested CpGs
    cand = np.sort(rng.choice(len(key), size=min(len(key), 3_000_000), replace=False))
    ck = key[cand]
    cc, cp = (ck >> 32).astype(np.int64), (ck & 0xFFFFFFFF).astype(np.int64)
    # distance to nearest DMR start on the same contig
    j = np.searchsorted(dk_sorted, ck)
    near = np.full(len(ck), np.iinfo(np.int64).max)
    for jj in (j - 1, j):
        ok = (jj >= 0) & (jj < len(dk_sorted))
        d = np.where(ok, dk_sorted[np.clip(jj, 0, len(dk_sorted) - 1)], -1)
        same = ok & ((d >> 32) == cc)
        near = np.where(same, np.minimum(near, np.abs((d & 0xFFFFFFFF) - cp)), near)
    pool_ok = near >= EXCLUDE_BP

    rows, used = [], set()
    for r in dmrs.itertuples(index=False):
        L = int(r.end - r.start)
        n = int(r.n_cpgs_tested)
        tol = max(1, int(round(0.3 * n)))
        cnt = np.searchsorted(key, ck + L) - np.searchsorted(key, ck)
        fits = (cp - FLANK_BP >= 0) & (cp + L + FLANK_BP <= sizes[cc])
        ok = np.flatnonzero(pool_ok & fits & (np.abs(cnt - n) <= tol))
        if used:
            ok = np.setdiff1d(ok, np.fromiter(used, dtype=np.int64), assume_unique=True)
        pick = rng.choice(ok, size=min(N_CONTROLS, len(ok)), replace=False) if len(ok) else []
        rows.append({"window_id": r.dmr_id, "type": "dmr", "dmr_id": r.dmr_id,
                     "contig": r.contig, "start": int(r.start), "end": int(r.end),
                     "n_cpgs_tested": n, "status": r.status, "direction": r.direction})
        for k, i in enumerate(pick):
            used.add(int(i))
            rows.append({"window_id": f"{r.dmr_id}_c{k + 1}", "type": "control", "dmr_id": r.dmr_id,
                         "contig": contigs[cc[i]], "start": int(cp[i]), "end": int(cp[i] + L),
                         "n_cpgs_tested": int(cnt[i]), "status": r.status,
                         "direction": r.direction})
    w = pd.DataFrame(rows)
    w.to_csv(path, sep="\t", index=False)
    return w


# =============================================================================
# 2-3. Variants and callable bases, chunked in parallel
# =============================================================================
def flanked_bed(frame, w):
    sizes = dict(pd.read_csv(f"{m11.ECO_FA[frame]}.fai", sep="\t", header=None, usecols=[0, 1]).values)
    b = pd.DataFrame({"c": w["contig"],
                      "s": (w["start"] - FLANK_BP).clip(lower=0),
                      "e": [min(sizes[c], e + FLANK_BP) for c, e in zip(w["contig"], w["end"])]})
    path = OUT / f"_{frame}.regions.bed"
    b.sort_values(["c", "s"]).to_csv(path, sep="\t", header=False, index=False)
    merged = OUT / f"_{frame}.regions.merged.bed"
    sh(f"LC_ALL=C sort -k1,1 -k2,2n '{path}' | '{BEDTOOLS}' merge -i - > '{merged}'")
    return merged


def run_chunk(frame, i, chunk_bed, bam_list):
    vcf = VCF / f"{frame}.chunk{i:02d}.bcf"
    if not m11.nonempty(vcf):
        sh(f"'{BCFTOOLS}' mpileup -X pacbio-ccs -q {MIN_MAPQ} -a FORMAT/AD,FORMAT/DP --ignore-RG "
           f"-f '{m11.ECO_FA[frame]}' -R '{chunk_bed}' -b '{bam_list}' -Ou 2>/dev/null "
           f"| '{BCFTOOLS}' call -mv -Ob -o '{vcf}.tmp' 2>/dev/null && mv '{vcf}.tmp' '{vcf}'")
    cal = CALL / f"{frame}.chunk{i:02d}.callable.bed"
    if not m11.nonempty(cal):
        sh(f"'{SAMTOOLS}' depth -a -Q {MIN_MAPQ} -b '{chunk_bed}' -f '{bam_list}' "
           f"| awk -v m={MIN_DP} 'BEGIN{{OFS=\"\\t\"}} {{ok=1; for(i=3;i<=NF;i++) if($i<m){{ok=0;break}} "
           f"if(ok) print $1,$2-1,$2}}' | '{BEDTOOLS}' merge -i - > '{cal}.tmp' && mv '{cal}.tmp' '{cal}'")
    return vcf, cal


def variants_and_callable(frame, w):
    merged = flanked_bed(frame, w)
    bam_list = OUT / f"_{frame}.bams.txt"
    bam_list.write_text("\n".join(str(b) for b in bams(frame)) + "\n")
    regions = pd.read_csv(merged, sep="\t", header=None)
    chunks = np.array_split(np.arange(len(regions)), N_CHUNKS)
    beds = []
    for i, idx in enumerate(chunks):
        p = OUT / f"_{frame}.chunk{i:02d}.bed"
        regions.iloc[idx].to_csv(p, sep="\t", header=False, index=False)
        beds.append(p)
    with ThreadPoolExecutor(max_workers=N_CHUNKS) as ex:
        res = list(ex.map(lambda a: run_chunk(frame, *a, bam_list), enumerate(beds)))
    print(f"  [{frame}] {len(regions):,} merged regions, {(regions[2] - regions[1]).sum() / 1e6:.1f} Mb "
          f"genotyped in {N_CHUNKS} chunks", flush=True)
    return [r[0] for r in res], [r[1] for r in res]


def read_sites(vcfs):
    rows = []
    for v in vcfs:
        out = subprocess.run(
            f"'{BCFTOOLS}' query -i 'QUAL>={MIN_QUAL}' -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n' '{v}'",
            shell=True, check=True, capture_output=True, text=True, executable="/bin/bash").stdout
        for line in out.splitlines():
            p = line.split("\t")
            ref, alts = p[2], p[3].split(",")
            ad = [[int(x) if x != "." else 0 for x in s.split(",")] for s in p[4:]]
            # major alt allele across the 8 fish
            alt_tot = [sum(a[k] for a in ad if len(a) > k) for k in range(1, len(alts) + 1)]
            k = int(np.argmax(alt_tot)) + 1
            refc = np.array([a[0] for a in ad], dtype=float)
            altc = np.array([a[k] if len(a) > k else 0 for a in ad], dtype=float)
            dp = refc + altc
            is_indel = len(ref) != len(alts[k - 1])
            rows.append((p[0], int(p[1]) - 1, is_indel, *dp, *(altc / np.where(dp > 0, dp, np.nan))))
    cols = ["contig", "pos", "indel"] + [f"dp_{s}" for s in SAMPLES] + [f"af_{s}" for s in SAMPLES]
    s = pd.DataFrame(rows, columns=cols)
    dp = s[[f"dp_{x}" for x in SAMPLES]].to_numpy()
    af = s[[f"af_{x}" for x in SAMPLES]].to_numpy()
    full = (dp >= MIN_DP).all(1)
    lean_lo, lean_hi = (af[:, ~IS_SISC] <= LO_AF).all(1), (af[:, ~IS_SISC] >= HI_AF).all(1)
    sisc_lo, sisc_hi = (af[:, IS_SISC] <= LO_AF).all(1), (af[:, IS_SISC] >= HI_AF).all(1)
    s["covered"] = full
    s["fixed"] = full & ((lean_lo & sisc_hi) | (lean_hi & sisc_lo))
    s["alt_in"] = np.where(s["fixed"] & lean_hi, "lean",
                           np.where(s["fixed"] & sisc_hi, "siscowet", ""))
    d = np.nanmean(af[:, IS_SISC], 1) - np.nanmean(af[:, ~IS_SISC], 1)
    s["diff75"] = full & ~s["indel"] & (np.abs(d) >= 0.75)
    s["diff50"] = full & ~s["indel"] & (np.abs(d) >= 0.5)
    s["all_het"] = full & ((af >= 0.25) & (af <= 0.75)).all(1)
    return s


# =============================================================================
# 4. Per-window bins
# =============================================================================
def bin_intervals(w):
    """Window x bin intervals: core; near = core +/- 1 kb minus core; far = +/- 5 kb minus near."""
    rows = []
    for r in w.itertuples(index=False):
        s, e = r.start, r.end
        rows.append((r.contig, s, e, r.window_id, "core"))
        rows.append((r.contig, max(0, s - NEAR_BP), s, r.window_id, "near"))
        rows.append((r.contig, e, e + NEAR_BP, r.window_id, "near"))
        rows.append((r.contig, max(0, s - FLANK_BP), max(0, s - NEAR_BP), r.window_id, "far"))
        rows.append((r.contig, e + NEAR_BP, e + FLANK_BP, r.window_id, "far"))
    b = pd.DataFrame(rows, columns=["contig", "start", "end", "window_id", "bin"])
    return b[b["end"] > b["start"]]


def window_counts(frame, w, sites, cal_beds):
    bins = bin_intervals(w)
    bpath = OUT / f"_{frame}.bins.bed"
    bins.sort_values(["contig", "start"]).to_csv(bpath, sep="\t", header=False, index=False)
    cal = OUT / f"_{frame}.callable.bed"
    sh(f"cat {' '.join(repr(str(c)) for c in cal_beds)} | LC_ALL=C sort -k1,1 -k2,2n "
       f"| '{BEDTOOLS}' merge -i - > '{cal}'")
    ov = OUT / f"_{frame}.bins_callable.tsv"
    sh(f"LC_ALL=C sort -k1,1 -k2,2n '{bpath}' | '{BEDTOOLS}' intersect -a - -b '{cal}' -wao > '{ov}'")
    o = pd.read_csv(ov, sep="\t", header=None)
    callable_bp = o.groupby([3, 4])[o.columns[-1]].sum()

    # assign callable variant sites to bins (a site can sit in bins of several windows)
    cs = sites[sites["covered"]].copy()
    cpath = OUT / f"_{frame}.sites.bed"
    pd.DataFrame({"c": cs["contig"], "s": cs["pos"], "e": cs["pos"] + 1,
                  "i": np.arange(len(cs))}).sort_values(["c", "s"]).to_csv(
        cpath, sep="\t", header=False, index=False)
    hp = OUT / f"_{frame}.sites_bins.tsv"
    sh(f"LC_ALL=C sort -k1,1 -k2,2n '{bpath}' | '{BEDTOOLS}' intersect -a '{cpath}' -b - -wa -wb > '{hp}'")
    h = pd.read_csv(hp, sep="\t", header=None, names=["c", "s", "e", "i", "bc", "bs", "be", "window_id", "bin"])
    h = h.merge(cs.reset_index(drop=True)[["fixed", "indel", "alt_in", "diff75", "diff50"]],
                left_on="i", right_index=True)
    h["fixed_snv"] = h["fixed"] & ~h["indel"]
    h["fixed_indel"] = h["fixed"] & h["indel"]
    agg = h.groupby(["window_id", "bin"]).agg(
        variants=("i", "size"), fixed_snv=("fixed_snv", "sum"), fixed_indel=("fixed_indel", "sum"),
        diff75=("diff75", "sum"), diff50=("diff50", "sum"))
    idx = pd.MultiIndex.from_product([w["window_id"], ["core", "near", "far"]],
                                     names=["window_id", "bin"])
    wc = pd.DataFrame(index=idx)
    wc["callable_bp"] = callable_bp.reindex(idx).fillna(0).astype(int).to_numpy()
    wc = wc.join(agg).fillna(0).reset_index()
    wc = wc.merge(w[["window_id", "type", "dmr_id", "status", "direction"]], on="window_id")
    wc.to_csv(OUT / f"window_counts_{frame}.tsv", sep="\t", index=False)

    fixed = h[h["fixed"]].merge(cs.reset_index(drop=True), left_on="i", right_index=True,
                                suffixes=("", "_s"))
    keep = ["window_id", "bin", "contig", "pos", "indel", "alt_in"] + \
           [f"af_{s}" for s in SAMPLES] + [f"dp_{s}" for s in SAMPLES]
    fixed[keep].to_csv(OUT / f"fixed_sites_{frame}.tsv", sep="\t", index=False, float_format="%.3g")
    return wc


# =============================================================================
# 5. Statistics
# =============================================================================
CLASSES = (("fixed", "fixed SNV/indel"), ("diff75", "SNV |dAF|>=0.75"), ("diff50", "SNV |dAF|>=0.5"))


def compare(wc, rng, subset=None):
    """DMR vs matched controls per class and bin: fraction of windows with >= 1 site (primary,
    resampling p with one control per DMR) and pooled rate per callable kb (context)."""
    if subset is not None:
        wc = wc[wc["dmr_id"].isin(subset)]
    wc = wc.assign(fixed=wc["fixed_snv"] + wc["fixed_indel"])
    piv = wc.pivot_table(index=["window_id", "type", "dmr_id"], columns="bin",
                         values=["fixed", "diff75", "diff50", "callable_bp"], aggfunc="sum").fillna(0)
    piv.columns = [f"{a}_{b}" for a, b in piv.columns]
    piv = piv.reset_index()
    for cls, _ in CLASSES:
        piv[f"{cls}_core1kb"] = piv[f"{cls}_core"] + piv[f"{cls}_near"]
    d = piv[piv["type"] == "dmr"].set_index("dmr_id")
    c = piv[piv["type"] == "control"]
    groups = {k: g for k, g in c.groupby("dmr_id")}
    ids = [x for x in d.index if x in groups]
    if not ids:
        return ["    (no DMRs with controls)"]
    d = d.loc[ids]
    lines = [f"    {'class':17s} {'bin':10s} {'DMR windows hit':>17s} {'control':>8s} {'p':>7s}   "
             f"{'rate/kb DMR':>11s} {'control':>8s}"]
    for cls, label in CLASSES:
        for b in ("core1kb", "far"):
            col = f"{cls}_{b}"
            hit_d = (d[col] > 0).to_numpy()
            ctrl = [(groups[x][col] > 0).to_numpy() for x in ids]
            picks = [rng.integers(0, len(x), N_RESAMPLE) for x in ctrl]
            null = np.mean([h[p] for h, p in zip(ctrl, picks)], axis=0)
            p = (1 + (null >= hit_d.mean()).sum()) / (1 + N_RESAMPLE)
            bins = ("core", "near") if b == "core1kb" else ("far",)
            cb_d = sum(d[f"callable_bp_{x}"].sum() for x in bins)
            cb_c = sum(c[f"callable_bp_{x}"].sum() for x in bins)
            rd = 1000 * d[col].sum() / max(1, cb_d)
            rc = 1000 * c[col].sum() / max(1, cb_c)
            lines.append(f"    {label:17s} {b:10s} {int(hit_d.sum()):5d}/{len(d):<4d} ({hit_d.mean():5.1%}) "
                         f"{np.mean(np.concatenate(ctrl)):7.1%} {p:7.4f}   {rd:11.3f} {rc:8.3f}")
    return lines


def depth_check(frame, w, wc, rng):
    z = np.load(m11.COUNTS / f"frame_{frame}.cpg_matrix.npz")
    key, nv = z["key"], z["nvalid"].astype(float)
    med = np.median(nv, 0)
    contigs = pd.read_csv(f"{m11.ECO_FA[frame]}.fai", sep="\t", header=None, usecols=[0])[0].tolist()
    ci = {c: i for i, c in enumerate(contigs)}
    k0 = np.array([(ci[c] << 32) | s for c, s in zip(w["contig"], w["start"])], dtype=np.int64)
    k1 = np.array([(ci[c] << 32) | e for c, e in zip(w["contig"], w["end"])], dtype=np.int64)
    i0, i1 = np.searchsorted(key, k0), np.searchsorted(key, k1)
    cs = np.vstack([np.zeros((1, nv.shape[1])), np.cumsum(nv, 0)])
    dep = (cs[i1] - cs[i0]) / np.maximum(i1 - i0, 1)[:, None] / med
    out = w[["window_id", "type", "dmr_id", "status", "direction"]].copy()
    out["log2_depth_ratio"] = np.log2(dep[:, IS_SISC].mean(1) / dep[:, ~IS_SISC].mean(1))
    near = wc[wc["bin"].isin(["core", "near"])].groupby("window_id")["diff50"].sum() > 0
    out["diff50_core1kb"] = out["window_id"].map(near).fillna(False).astype(bool)
    out["depth_dev"] = out["log2_depth_ratio"].abs() >= DEPTH_DEV
    out.to_csv(OUT / f"window_depth_{frame}.tsv", sep="\t", index=False, float_format="%.4g")

    d = out[out["type"] == "dmr"].set_index("dmr_id")
    groups = {k: g["depth_dev"].to_numpy() for k, g in out[out["type"] == "control"].groupby("dmr_id")}
    ids = [x for x in d.index if x in groups]
    obs = d.loc[ids, "depth_dev"].mean()
    picks = [rng.integers(0, len(groups[x]), N_RESAMPLE) for x in ids]
    null = np.mean([groups[x][p] for x, p in zip(ids, picks)], axis=0)
    p = (1 + (null >= obs).sum()) / (1 + N_RESAMPLE)
    ctrl = np.mean(np.concatenate([groups[x] for x in ids]))
    dh, dn = d[d["diff50_core1kb"]], d[~d["diff50_core1kb"]]
    hypo = d[d["direction"] == "hypo_siscowet"]
    return [f"  copy-number check (|log2 siscowet/lean depth| >= {DEPTH_DEV}):",
            f"    DMR windows {obs:.1%} vs matched controls {ctrl:.1%} (one-sided p {p:.4f})",
            f"    DMRs with a |dAF|>=0.5 SNV within 1 kb {dh['depth_dev'].mean():.1%} (n={len(dh)}) vs "
            f"without {dn['depth_dev'].mean():.1%} (n={len(dn)})",
            f"    hypo-in-siscowet DMRs: siscowet depth lower {int((hypo['log2_depth_ratio'] <= -DEPTH_DEV).sum())}, "
            f"higher {int((hypo['log2_depth_ratio'] >= DEPTH_DEV).sum())}"]


def main():
    for d in (OUT, VCF, CALL):
        d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    report = ["Fixed sequence differences near native DMRs (23.15)", "=" * 52, "",
              f"Fixed = >= {MIN_DP} reads in all 8 fish; all 4 of one ecotype alt fraction <= {LO_AF} "
              f"and all 4 of the other >= {HI_AF}. bcftools mpileup -X pacbio-ccs -q {MIN_MAPQ}, "
              f"QUAL >= {MIN_QUAL}. Controls: {N_CONTROLS} per DMR, same length, tested-CpG count "
              f"within 30 %, >= {EXCLUDE_BP // 1000} kb from any DMR.", ""]
    for frame in FRAMES:
        print("=" * 70 + f"\n{frame} frame\n" + "=" * 70, flush=True)
        w = build_windows(frame, rng)
        n_ctrl = w[w["type"] == "control"].groupby("dmr_id").size()
        vcfs, cals = variants_and_callable(frame, w)
        sites = read_sites(vcfs)
        wc = window_counts(frame, w, sites, cals)
        d = w[w["type"] == "dmr"]
        report += [f"{frame} frame: {len(d)} DMRs, {len(w) - len(d)} controls "
                   f"(median {n_ctrl.median():.0f} per DMR; {int((n_ctrl < N_CONTROLS).sum())} DMRs "
                   f"with fewer than {N_CONTROLS})",
                   f"  variant sites with >= {MIN_DP} reads in all fish: {int(sites['covered'].sum()):,}; "
                   f"fixed: {int(sites['fixed'].sum()):,} "
                   f"(alt in lean {int((sites['alt_in'] == 'lean').sum())}, "
                   f"alt in siscowet {int((sites['alt_in'] == 'siscowet').sum())})",
                   f"  QC: covered sites where all 8 fish look heterozygous (collapsed paralogs): "
                   f"{int(sites['all_het'].sum()):,} ({sites['all_het'].sum() / max(1, sites['covered'].sum()):.0%}); "
                   f"differentiated SNVs |dAF|>=0.75: {int(sites['diff75'].sum()):,}, >=0.5: "
                   f"{int(sites['diff50'].sum()):,}",
                   "  all DMRs:"]
        report += compare(wc, rng)
        for label, ids in (("replicated DMRs", d.loc[d["status"] == "replicated", "dmr_id"]),
                           ("hypomethylated in siscowet", d.loc[d["direction"] == "hypo_siscowet", "dmr_id"]),
                           ("hypermethylated in siscowet", d.loc[d["direction"] == "hyper_siscowet", "dmr_id"])):
            report.append(f"  {label} (n={len(ids)}):")
            report += compare(wc, rng, set(ids))
        report += depth_check(frame, w, wc, rng)
        hits = wc[(wc["type"] == "dmr") & ((wc["fixed_snv"] + wc["fixed_indel"]) > 0)]
        report.append("  DMR windows with a fixed difference within 5 kb (window: bin=count):")
        for wid, g in hits.groupby("window_id"):
            dm = d.set_index("window_id").loc[wid]
            report.append(f"    {wid:9s} {dm['status']:13s} {dm['direction']:15s} " + ", ".join(
                f"{r.bin}={int(r.fixed_snv + r.fixed_indel)}" for r in g.itertuples()))
        report.append("")
    text = "\n".join(report) + "\n"
    (OUT / "fixed_variant_report.txt").write_text(text)
    print("\n" + text)


if __name__ == "__main__":
    main()
