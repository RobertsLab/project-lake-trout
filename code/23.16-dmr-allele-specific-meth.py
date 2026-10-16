#!/usr/bin/env python3
"""
23.16-dmr-allele-specific-meth.py  (Phase 3 follow-up, code/23-next-phase-research-plan.md)

ALLELE-SPECIFIC METHYLATION AT DIFFERENTIATED NATIVE DMRs. 23.15 found differentiated SNVs
(|dAF| >= 0.5 between ecotypes) enriched near the 23.13 DMRs. That fits cis control (the haplotype
sets its own methylation) but also linkage with a trans or environmental difference. HiFi reads carry
both the SNV base and the 5mC call (MM/ML), so in a fish heterozygous for a differentiated SNV the two
haplotypes can be read separately:

    ASM = % methylation on the siscowet-associated allele - % on the lean-associated allele
          (at the DMR core CpGs, in the same fish)

  cis:   ASM ~ the DMR difference (siscowet - lean), in lean AND siscowet heterozygotes
  trans / environment: ASM ~ 0 (both haplotypes follow the fish's ecotype)

Per frame (both ecotypes' reads on one purged assembly, as in 23.13 / 23.15):
  1. SNVs from the 23.15 bcftools calls with >= 4 reads in all 8 fish. Differentiated = |dAF| >= 0.5;
     the siscowet-associated allele is the one more frequent in siscowet. (A differentiated SNV cannot
     look heterozygous in all 8 fish, so the 23.15 paralog flag does not apply; collapsed-paralog copy
     number is instead handled by re-running the tests without DMRs that show an ecotype depth
     difference, 23.15 window_depth.)
  2. Windows: 23.13 DMRs and the 23.15 matched controls (10 per DMR) with a differentiated SNV within
     5 kb of the core. Controls give the null spread of ASM for the same SNV class at sequence with no
     ecotype methylation difference.
  3. For each heterozygous fish (alt fraction 0.2-0.8, >= 3 reads per allele): primary MAPQ >= 20
     reads spanning the SNV are assigned to an allele by their base; 5mC calls (ML >= 128 =
     methylated) at the window's tested core CpGs (23.13 matrix) are pooled per allele. Forward reads
     report the CpG C, reverse reads the G at +1 (checked on these BAMs).
     An allele needs >= MIN_READS phased reads with core calls and >= MIN_CALLS calls.
  4. Per (window, fish) the SNV with most phased reads is kept; per window, ASM is averaged over fish.
  5. Tests on DMR windows: slope of ASM on the DMR difference through the origin (cis = 1, trans = 0)
     with a bootstrap CI over DMRs; sign agreement; separately in lean and siscowet heterozygotes.
     Each DMR is classed cis-like (same sign, |ASM| >= half the DMR difference and above the 95th
     percentile of control |ASM|) or trans-like (|ASM| below the control median while |diff| >= 20).

Inputs (raven): 23.13 cpg_scan/dmrs_{frame}.tsv; 23.15 fixed_variants/{windows,window_depth}_{frame}.tsv
  and vcf/{frame}.chunk*.bcf; counts/frame_{frame}.cpg_matrix.npz; the 16 frame BAMs.

Outputs (analyses/23-native-methylation/asm/):
  asm_alleles_{frame}.tsv   window x SNV x fish x allele: reads, calls, % methylated
  asm_windows_{frame}.tsv   window-level ASM, DMR difference, class
  asm_report.txt

Run on raven with the pysam env: /home/sr320/miniconda3/envs/syri_env/bin/python
(kept free of the 23.11 module, which needs statsmodels).

Author: Generated for project-lake-trout
Date: 2026-10-16
"""

import glob
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy import stats

BASE = Path(__file__).resolve().parent.parent
NM = BASE / "analyses" / "23-native-methylation"
SCAN, FV, COUNTS = NM / "cpg_scan", NM / "fixed_variants", NM / "counts"
OUT = NM / "asm"
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
BCFTOOLS = "/home/shared/bcftools-1.14/bcftools"

SAMPLES_LEAN = ["bc2041", "bc2068", "bc2069", "bc2070"]
SAMPLES_SISCOWET = ["bc2071", "bc2072", "bc2073", "bc2096"]
SAMPLES = SAMPLES_LEAN + SAMPLES_SISCOWET
IS_SISC = np.array([s in SAMPLES_SISCOWET for s in SAMPLES])
ECO = {"lean": SAMPLES_LEAN, "siscowet": SAMPLES_SISCOWET}
FA = {"lean": SYN / "lean.purged.fa", "siscowet": SYN / "siscowet.purged.fa"}
FRAMES = ("lean", "siscowet")

MIN_DP, MIN_QUAL, MIN_MAPQ = 4, 20, 20
DIFF_AF = 0.5
SNV_WINDOW = 5000
HET_LO, HET_HI, MIN_ALLELE_AD = 0.2, 0.8, 3
MIN_READS, MIN_CALLS = 3, 10
METH_PROB = 128
N_BOOT = 10_000
SEED = 2316


def bam_path(sample, frame):
    if sample in ECO[frame]:
        return NM / "align" / f"{sample}.to_own.sorted.bam"
    return BASE / "analyses" / "23-reciprocal-pav" / "align" / f"{sample}.to_{frame}.sorted.bam"


# =============================================================================
# 1. SNVs
# =============================================================================
def load_snvs(frame):
    rows = []
    for v in sorted(glob.glob(str(FV / "vcf" / f"{frame}.chunk*.bcf"))):
        out = subprocess.run(
            f"'{BCFTOOLS}' query -i 'QUAL>={MIN_QUAL} && TYPE=\"snp\"' "
            f"-f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%AD]\\n' '{v}'",
            shell=True, check=True, capture_output=True, text=True).stdout
        for line in out.splitlines():
            p = line.split("\t")
            alts = p[3].split(",")
            ad = [[int(x) if x != "." else 0 for x in s.split(",")] for s in p[4:]]
            tot = [sum(a[k] for a in ad if len(a) > k) for k in range(1, len(alts) + 1)]
            k = int(np.argmax(tot)) + 1
            if len(p[2]) != 1 or len(alts[k - 1]) != 1:
                continue
            refc = [a[0] for a in ad]
            altc = [a[k] if len(a) > k else 0 for a in ad]
            rows.append((p[0], int(p[1]) - 1, p[2].upper(), alts[k - 1].upper(), *refc, *altc))
    cols = ["contig", "pos", "ref", "alt"] + [f"ref_{s}" for s in SAMPLES] + [f"alt_{s}" for s in SAMPLES]
    s = pd.DataFrame(rows, columns=cols).drop_duplicates(["contig", "pos"])
    rc = s[[f"ref_{x}" for x in SAMPLES]].to_numpy(float)
    ac = s[[f"alt_{x}" for x in SAMPLES]].to_numpy(float)
    dp = rc + ac
    af = ac / np.where(dp > 0, dp, np.nan)
    full = (dp >= MIN_DP).all(1)
    d = np.nanmean(af[:, IS_SISC], 1) - np.nanmean(af[:, ~IS_SISC], 1)
    s["dAF"] = d
    s["diff"] = full & (np.abs(d) >= DIFF_AF)
    s["all_het"] = full & ((af >= 0.25) & (af <= 0.75)).all(1)
    s["sisc_allele"] = np.where(d > 0, s["alt"], s["ref"])
    s["lean_allele"] = np.where(d > 0, s["ref"], s["alt"])
    for i, x in enumerate(SAMPLES):
        s[f"het_{x}"] = (af[:, i] >= HET_LO) & (af[:, i] <= HET_HI) & \
                        (rc[:, i] >= MIN_ALLELE_AD) & (ac[:, i] >= MIN_ALLELE_AD)
    return s[s["diff"]].reset_index(drop=True)


# =============================================================================
# 2-3. Phasing and per-allele methylation
# =============================================================================
def core_cpgs(frame, windows):
    z = np.load(COUNTS / f"frame_{frame}.cpg_matrix.npz")
    key = z["key"]
    contigs = pd.read_csv(f"{FA[frame]}.fai", sep="\t", header=None, usecols=[0])[0].tolist()
    ci = {c: i for i, c in enumerate(contigs)}
    out = {}
    for r in windows.itertuples(index=False):
        lo = np.searchsorted(key, (ci[r.contig] << 32) | r.start)
        hi = np.searchsorted(key, (ci[r.contig] << 32) | r.end)
        out[r.window_id] = set((key[lo:hi] & 0xFFFFFFFF).astype(int).tolist())
    return out


def allele_meth(bam, contig, snv_pos, alleles, cpgs):
    """Per allele: [reads, methylated calls, calls] at the given CpG positions."""
    res = {a: [0, 0, 0] for a in alleles}
    lo, hi = min(cpgs | {snv_pos}), max(cpgs | {snv_pos}) + 2
    for r in bam.fetch(contig, snv_pos, snv_pos + 1):
        if r.is_secondary or r.is_supplementary or r.mapping_quality < MIN_MAPQ:
            continue
        if r.reference_start > lo or r.reference_end < hi:
            continue                       # read must span the SNV and the whole core
        pairs = r.get_aligned_pairs(matches_only=True)
        q_of = {rp: q for q, rp in pairs}
        q = q_of.get(snv_pos)
        if q is None:
            continue
        base = r.query_sequence[q].upper()
        if base not in res:
            continue
        mb = r.modified_bases
        if not mb:
            continue
        ref_of = dict(pairs)
        m = n = 0
        shift = 1 if r.is_reverse else 0
        for key, lst in mb.items():
            if key[2] != "m":
                continue
            for qq, prob in lst:
                rp = ref_of.get(qq)
                if rp is None:
                    continue
                if (rp - shift) in cpgs:
                    n += 1
                    m += prob >= METH_PROB
        if n:
            res[base][0] += 1
            res[base][1] += m
            res[base][2] += n
    return res


def phase_frame(frame, windows, snvs, cpgs):
    bams = {s: pysam.AlignmentFile(str(bam_path(s, frame))) for s in SAMPLES}
    by_contig = {c: g for c, g in snvs.groupby("contig")}
    rows = []
    for w in windows.itertuples(index=False):
        g = by_contig.get(w.contig)
        if g is None or not cpgs.get(w.window_id):
            continue
        g = g[(g["pos"] >= w.start - SNV_WINDOW) & (g["pos"] < w.end + SNV_WINDOW)]
        for snv in g.itertuples(index=False):
            dist = 0 if w.start <= snv.pos < w.end else min(abs(snv.pos - w.start), abs(snv.pos - w.end + 1))
            for s in SAMPLES:
                if not getattr(snv, f"het_{s}"):
                    continue
                res = allele_meth(bams[s], w.contig, snv.pos, (snv.sisc_allele, snv.lean_allele),
                                  cpgs[w.window_id])
                for allele, side in ((snv.sisc_allele, "sisc"), (snv.lean_allele, "lean")):
                    nr, nm, nc = res[allele]
                    rows.append({"window_id": w.window_id, "type": w.type, "dmr_id": w.dmr_id,
                                 "snv": f"{w.contig}:{snv.pos + 1}", "snv_dist": dist,
                                 "dAF": round(snv.dAF, 3), "all_het": snv.all_het, "fish": s,
                                 "fish_ecotype": "siscowet" if s in SAMPLES_SISCOWET else "lean",
                                 "allele_side": side, "allele": allele, "reads": nr,
                                 "meth_calls": nm, "calls": nc,
                                 "pct": 100.0 * nm / nc if nc else np.nan})
    for b in bams.values():
        b.close()
    return pd.DataFrame(rows)


# =============================================================================
# 4-5. Window ASM and tests
# =============================================================================
def window_asm(al):
    piv = al.pivot_table(index=["window_id", "type", "dmr_id", "snv", "snv_dist", "fish", "fish_ecotype"],
                         columns="allele_side", values=["reads", "calls", "pct"]).dropna()
    piv.columns = [f"{a}_{b}" for a, b in piv.columns]
    piv = piv[(piv["reads_sisc"] >= MIN_READS) & (piv["reads_lean"] >= MIN_READS) &
              (piv["calls_sisc"] >= MIN_CALLS) & (piv["calls_lean"] >= MIN_CALLS)].reset_index()
    piv["asm"] = piv["pct_sisc"] - piv["pct_lean"]
    piv["phased_reads"] = piv["reads_sisc"] + piv["reads_lean"]
    best = piv.sort_values(["phased_reads", "snv_dist"], ascending=[False, True]) \
              .drop_duplicates(["window_id", "fish"])
    return best


def summarize(best, dmrs, depth):
    w = best.groupby(["window_id", "type", "dmr_id"]).agg(
        asm=("asm", "mean"), n_fish=("fish", "size"),
        n_lean_het=("fish_ecotype", lambda x: int((x == "lean").sum())),
        n_sisc_het=("fish_ecotype", lambda x: int((x == "siscowet").sum())),
        asm_in_lean=("asm", lambda x: x[best.loc[x.index, "fish_ecotype"] == "lean"].mean()),
        asm_in_sisc=("asm", lambda x: x[best.loc[x.index, "fish_ecotype"] == "siscowet"].mean()),
        min_snv_dist=("snv_dist", "min")).reset_index()
    w = w.merge(dmrs[["dmr_id", "meth_diff", "status", "direction", "symbol", "product"]],
                on="dmr_id", how="left")
    w = w.merge(depth[["window_id", "log2_depth_ratio", "depth_dev"]], on="window_id", how="left")
    return w


def slope_through_origin(x, y):
    return float(np.sum(x * y) / np.sum(x * x)) if len(x) else np.nan


def classify_and_test(w, rng, label):
    d = w[w["type"] == "dmr"].copy()
    c = w[w["type"] == "control"]
    lines = [f"  {label}: DMR windows with phased heterozygotes {len(d)}; control windows {len(c)}"]
    if len(d) < 3 or len(c) < 10:
        lines.append("    too few windows for tests")
        return lines, d
    ca = c["asm"].abs()
    q95, q50 = ca.quantile(0.95), ca.median()
    x, y = d["meth_diff"].to_numpy(), d["asm"].to_numpy()
    beta = slope_through_origin(x, y)
    boots = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(d), len(d))
        boots.append(slope_through_origin(x[i], y[i]))
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rho, rp = stats.spearmanr(x, y)
    agree = np.sign(x) == np.sign(y)
    sp = stats.binomtest(int(agree.sum()), len(d), 0.5, alternative="greater").pvalue
    d["class"] = np.where(agree & (np.abs(y) >= 0.5 * np.abs(x)) & (np.abs(y) > q95), "cis_like",
                          np.where((np.abs(y) < q50) & (np.abs(x) >= 20), "trans_like", "intermediate"))
    lines += [f"    control |ASM|: median {q50:.1f}, 95th pct {q95:.1f} points; DMR |ASM| median "
              f"{d['asm'].abs().median():.1f} (Mann-Whitney p "
              f"{stats.mannwhitneyu(d['asm'].abs(), ca, alternative='greater').pvalue:.2g})",
              f"    slope ASM ~ DMR difference: {beta:.2f} (95 % CI {lo:.2f} to {hi:.2f}; cis = 1, "
              f"trans = 0); Spearman rho {rho:.2f} (p {rp:.2g})",
              f"    sign agreement {int(agree.sum())}/{len(d)} (binomial p {sp:.2g})"]
    for eco in ("lean", "sisc"):
        col = f"asm_in_{eco}"
        k = d[col].notna()
        if k.sum() >= 3:
            lines.append(f"    in {'lean' if eco == 'lean' else 'siscowet'} heterozygotes (n={int(k.sum())} DMRs): "
                         f"slope {slope_through_origin(d.loc[k, 'meth_diff'].to_numpy(), d.loc[k, col].to_numpy()):.2f}")
    lines.append("    classes: " + ", ".join(f"{k} {v}" for k, v in d["class"].value_counts().items()))
    for cls in ("cis_like", "trans_like"):
        k = d["class"] == cls
        if k.any():
            lines.append(f"      {cls}: replicated {int((d.loc[k, 'status'] == 'replicated').sum())}/"
                         f"{int(k.sum())}, depth difference {int(d.loc[k, 'depth_dev'].fillna(False).sum())}/"
                         f"{int(k.sum())}")
    return lines, d


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    report = ["Allele-specific methylation at differentiated native DMRs (23.16)", "=" * 64, "",
              f"Differentiated SNV |dAF| >= {DIFF_AF} within {SNV_WINDOW // 1000} kb of the window core; "
              f"heterozygote = alt fraction {HET_LO}-{HET_HI} with >= {MIN_ALLELE_AD} reads per allele; "
              f">= {MIN_READS} phased reads and >= {MIN_CALLS} core CpG calls per allele; reads must span "
              f"SNV and core; ML >= {METH_PROB}/255 = methylated. ASM = siscowet-allele % - lean-allele %.",
              ""]
    for frame in FRAMES:
        print(f"== {frame} frame", flush=True)
        dmrs = pd.read_csv(SCAN / f"dmrs_{frame}.tsv", sep="\t")
        dmrs["status"] = np.where(dmrs["replicated"].astype(str) == "True", "replicated",
                                  np.where(dmrs["other_frame_testable"].astype(str) == "True",
                                           "not_replicated", "single_frame"))
        windows = pd.read_csv(FV / f"windows_{frame}.tsv", sep="\t")
        depth = pd.read_csv(FV / f"window_depth_{frame}.tsv", sep="\t")
        snvs = load_snvs(frame)
        cpgs = core_cpgs(frame, windows)
        cache = OUT / f"asm_alleles_{frame}.tsv"
        if cache.exists():
            al = pd.read_csv(cache, sep="\t")
        else:
            al = phase_frame(frame, windows, snvs, cpgs)
            al.to_csv(cache, sep="\t", index=False, float_format="%.4g")
        n_dmr_snv = al.loc[al["type"] == "dmr", "window_id"].nunique()
        report += [f"{frame} frame: {len(snvs):,} differentiated SNVs; DMR windows with a heterozygote "
                   f"at one within 5 kb: {n_dmr_snv}"]
        best = window_asm(al)
        w = summarize(best, dmrs, depth)
        lines, d = classify_and_test(w, rng, "all windows")
        report += lines
        keep = ~w["depth_dev"].fillna(False).astype(bool)
        lines_nd, _ = classify_and_test(w[keep], rng, "excluding windows with an ecotype depth difference")
        report += lines_nd
        rep = w[(w["type"] == "control") | (w["status"] == "replicated")]
        lines_rep, _ = classify_and_test(rep, rng, "replicated DMRs only")
        report += lines_rep
        w = w.merge(d[["window_id", "class"]], on="window_id", how="left")
        w.to_csv(OUT / f"asm_windows_{frame}.tsv", sep="\t", index=False, float_format="%.4g")
        top = d.reindex(d["meth_diff"].abs().sort_values(ascending=False).index).head(25)
        report.append("    DMRs (largest differences first): DMR diff | ASM (lean hets / siscowet hets) | class")
        for _, r in top.iterrows():
            report.append(f"      {r['dmr_id']:9s} {r['meth_diff']:+6.1f} | {r['asm']:+6.1f} "
                          f"({r['asm_in_lean']:+.0f} / {r['asm_in_sisc']:+.0f}) n={r['n_fish']} | "
                          f"{r['class']:12s} {r['status']:13s} {str(r['symbol'])[:12]:12s} "
                          f"{str(r['product'])[:34]}")
        report.append("")
    text = re.sub(r"[+-]?\bnan\b", "–", "\n".join(report)) + "\n"
    (OUT / "asm_report.txt").write_text(text)
    print("\n" + text)


if __name__ == "__main__":
    main()
