#!/usr/bin/env python3
"""
23.11-native-meth-counts.py  (Phase 3, part C, of code/23-next-phase-research-plan.md)

COUNT-BASED native differential methylation. The first native pass (23.9) was power-limited by two
choices it inherited from 14-diff-meth: (i) it only used CpGs with >=10 reads in a sample, although
native per-sample depth is ~8-20x, so bc2069 kept 0.68 M CpGs and only 10,838 genes were testable;
(ii) it tested with Mann-Whitney on n=4 vs 4, whose smallest possible p is 0.029, so nothing can
pass FDR. Here:

  A. Feature-level test on ALL covered CpGs. Every 5mC CpG call (Nvalid >= 1) from the 23.8 modkit
     bedMethyl is pooled into native features -- gene body and promoter (2 kb upstream of the TSS,
     strand-aware) from each ecotype's Liftoff annotation. Features are joined across ecotypes by
     gene ID, restricted to genes that are single-copy in BOTH assemblies (no Liftoff `_N` copies),
     so the pairing is 1:1. Pooling reads across CpGs is what rescues the low-depth samples.

     Test: quasi-binomial GLM, methylated ~ ecotype, weights = Nvalid. For a two-group design the
     MLE is closed form (pooled group proportions), so it is vectorised over all features:
        b      = logit(p_sisc) - logit(p_lean)
        var(b) = 1/(N_lean p_lean q_lean) + 1/(N_sisc p_sisc q_sisc)
        phi    = max(1, Pearson X^2 / (n - 2))   (fish-to-fish overdispersion)
        t      = b / sqrt(phi var(b)),  p from t with n - 2 df,  BH-FDR.
     Unlike Mann-Whitney the p-value is continuous, so FDR is attainable with n = 4 vs 4, and the
     dispersion term keeps a coverage-heavy feature from looking significant when fish disagree.
     The 14-diff-meth Mann-Whitney on per-sample percentages is kept alongside for comparability.

  B. Locus-level check of the 302 reference DMRs (analyses/14-diff-meth/dmrs.bed). Each DMR plus
     500 bp flanks is pulled from the reference and placed on each purged assembly with
     minimap2 -x asm5 (unique primary hit, MAPQ >= 20, >= 80 % of the query aligned). The DMR core is
     projected onto native coordinates, counts are pooled per sample, and the same two tests are run.
     The key quantity is whether native methylation differences reproduce the reference-based ones
     (sign agreement, correlation), not just how many pass a threshold.
     Placement control: 1,000 random 1.2 kb reference windows (bedtools random, seed 23, <5 % N) are
     placed the same way, so the MAPQ profile of the DMR loci can be compared with a genome baseline
     (MAPQ 0 = the locus exists in >= 2 near-identical copies in the native assembly).

Inputs (from 23.8, on raven):
  analyses/23-native-methylation/meth/{sample}.CpG.bedmethyl
  output/21.1-gene-anchored-synteny/{lean,siscowet}.purged.{fa,liftoff.genes.bed}
  data/GCF_016432855.1_SaNama_1.0_genomic.fa, analyses/14-diff-meth/dmrs.bed
  analyses/18-annotation/gene_function_table.tsv

Outputs (analyses/23-native-methylation/):
  counts/{sample}.cpg_counts_all.sorted.bed     chrom,start,end,Nmod,Nvalid (Nvalid>=1; gitignored)
  counts/{sample}.features.tsv, counts/refdmr.*  intermediates (gitignored)
  native_feature_counts.tsv.gz                  feature x sample Nmod / Nvalid / n_cpg
  dmf_results.tsv.gz                            per-feature tests (all tested)
  dmf_significant.tsv                           FDR < 0.10 and |diff| >= 10 points
  ref_dmr_native_loci.tsv                       302 reference DMRs on native coordinates + tests
  ref_dmr_placement_mapq.tsv                    MAPQ classes, reference DMRs vs random windows
  native_meth_counts_report.txt

Compute: one awk/sort pass per 3 GB bedMethyl, eight bedtools maps, two minimap2 runs of 302 short
queries against ~2.7 Gb assemblies. Tens of minutes on raven; run with
/home/sr320/miniconda3/bin/python (statsmodels).

Author: Generated for project-lake-trout
Date: 2026-10-14
"""

import csv
import gzip
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# =============================================================================
# Configuration
# =============================================================================
BASE = Path(__file__).resolve().parent.parent
OUT = BASE / "analyses" / "23-native-methylation"
METH = OUT / "meth"
COUNTS = OUT / "counts"
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
ANN = BASE / "analyses" / "18-annotation"
REF_FA = BASE / "data" / "GCF_016432855.1_SaNama_1.0_genomic.fa"
REF_DMRS = BASE / "analyses" / "14-diff-meth" / "dmrs.bed"

SYRI_BIN = Path("/home/sr320/miniconda3/envs/syri_env/bin")
MINIMAP2 = str(SYRI_BIN / "minimap2")
SAMTOOLS = "/usr/bin/samtools"
BEDTOOLS = "/usr/bin/bedtools"
THREADS = 32

SAMPLES_LEAN = ["bc2041", "bc2068", "bc2069", "bc2070"]
SAMPLES_SISCOWET = ["bc2071", "bc2072", "bc2073", "bc2096"]
SAMPLES = SAMPLES_LEAN + SAMPLES_SISCOWET
ECO_SAMPLES = {"lean": SAMPLES_LEAN, "siscowet": SAMPLES_SISCOWET}
ECO_FA = {"lean": SYN / "lean.purged.fa", "siscowet": SYN / "siscowet.purged.fa"}
ECO_GENES = {"lean": SYN / "lean.purged.liftoff.genes.bed",
             "siscowet": SYN / "siscowet.purged.liftoff.genes.bed"}

PROMOTER_BP = 2000        # upstream of TSS
MIN_NVALID = 20           # pooled 5mC calls per sample per feature
MIN_CPGS = 3              # covered CpGs per sample per feature
MIN_SAMPLES = 3           # samples per ecotype passing the two filters above
DIFF_THRESHOLD = 10.0     # percentage points; gene-level means are diluted vs single CpGs
Q_VALUE_THRESHOLD = 0.10
P_VALUE_THRESHOLD = 0.05

DMR_FLANK = 500
N_RANDOM = 1000
RANDOM_LEN = 1200
MIN_MAPQ = 20
MIN_QUERY_FRAC = 0.80

COPY_SUFFIX = re.compile(r"_\d+$")


def sh(cmd):
    print("  $", cmd[:160], flush=True)
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def nonempty(p):
    p = Path(p)
    return p.exists() and p.stat().st_size > 0


# =============================================================================
# Step 1: all-coverage CpG counts per sample
# =============================================================================
def all_counts(sample):
    """chrom,start,end,Nmod,Nvalid for every 5mC CpG with Nvalid >= 1, sorted for bedtools.
    bedMethyl columns: 4 = mod code, 10 = Nvalid_cov, 12 = Nmod (as used in 23.8)."""
    dst = COUNTS / f"{sample}.cpg_counts_all.sorted.bed"
    if nonempty(dst):
        return dst
    src = METH / f"{sample}.CpG.bedmethyl"
    if not nonempty(src):
        sys.exit(f"FATAL: {src} missing -- run 23.8 first")
    sh(f"awk -v OFS='\\t' '$4==\"m\" && $10>=1{{print $1,$2,$3,$12,$10}}' '{src}' "
       f"| LC_ALL=C sort -S 8G --parallel=8 -k1,1 -k2,2n > '{dst}.tmp' && mv '{dst}.tmp' '{dst}'")
    return dst


# =============================================================================
# Step 2: native features on the shared single-copy gene axis
# =============================================================================
def shared_single_copy():
    """Gene IDs present exactly once (no Liftoff `_N` copies) in both assemblies."""
    per_eco = {}
    for eco, bed in ECO_GENES.items():
        ids = pd.read_csv(bed, sep="\t", header=None, usecols=[3])[3]
        base = ids.str.replace(COPY_SUFFIX, "", regex=True)
        counts = base.value_counts()
        per_eco[eco] = set(counts[counts == 1].index) & set(ids)
    return per_eco["lean"] & per_eco["siscowet"]


def feature_bed(eco, keep):
    """Gene body + strand-aware promoter features for one ecotype, sorted, chrom-length clipped."""
    dst = COUNTS / f"{eco}.features.bed"
    if nonempty(dst):
        return dst
    sizes = dict(pd.read_csv(f"{ECO_FA[eco]}.fai", sep="\t", header=None, usecols=[0, 1]).values)
    genes = pd.read_csv(ECO_GENES[eco], sep="\t", header=None,
                        names=["chrom", "start", "end", "gene_id", "score", "strand"])
    genes = genes[genes["gene_id"].isin(keep)]
    rows = []
    for g in genes.itertuples(index=False):
        rows.append((g.chrom, g.start, g.end, f"{g.gene_id}|body"))
        if g.strand == "-":
            ps, pe = g.end, min(g.end + PROMOTER_BP, sizes[g.chrom])
        else:
            ps, pe = max(0, g.start - PROMOTER_BP), g.start
        if pe > ps:
            rows.append((g.chrom, ps, pe, f"{g.gene_id}|promoter"))
    df = pd.DataFrame(rows, columns=["chrom", "start", "end", "name"])
    df.sort_values(["chrom", "start", "end"]).to_csv(dst, sep="\t", header=False, index=False)
    return dst


def map_counts(features, counts, dst):
    """Per feature: sum Nmod, sum Nvalid, count CpGs."""
    if not nonempty(dst):
        sh(f"LC_ALL=C sort -k1,1 -k2,2n '{features}' | '{BEDTOOLS}' map -a - -b '{counts}' "
           f"-c 4,5,4 -o sum,sum,count -null 0 > '{dst}'")
    df = pd.read_csv(dst, sep="\t", header=None,
                     names=["chrom", "start", "end", "name", "nmod", "nvalid", "ncpg"])
    for c in ("nmod", "nvalid", "ncpg"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)
    return df


# =============================================================================
# Statistics
# =============================================================================
def quasibinomial_two_group(nmod, nvalid, is_sisc):
    """Vectorised quasi-binomial Wald test for methylated ~ group.

    nmod, nvalid: (features x samples) arrays, NaN where a sample fails the filters.
    is_sisc: boolean (samples,) group indicator. Returns logOR, se, t, p, df, phi.
    """
    ok = ~np.isnan(nvalid)
    nm = np.where(ok, nmod, 0.0)
    nv = np.where(ok, nvalid, 0.0)
    L, S = ~is_sisc, is_sisc
    # pooled group proportions (+0.5 keeps logits finite at 0 % / 100 %)
    mL, vL = nm[:, L].sum(1) + 0.5, nv[:, L].sum(1) + 1.0
    mS, vS = nm[:, S].sum(1) + 0.5, nv[:, S].sum(1) + 1.0
    pL, pS = mL / vL, mS / vS
    b = np.log(pS / (1 - pS)) - np.log(pL / (1 - pL))
    var = 1.0 / (vL * pL * (1 - pL)) + 1.0 / (vS * pS * (1 - pS))
    # Pearson dispersion over samples with data
    pg = np.where(is_sisc[None, :], pS[:, None], pL[:, None])
    resid2 = np.where(ok, (nm - nv * pg) ** 2 / (nv * pg * (1 - pg) + 1e-12), 0.0)
    n = ok.sum(1)
    df = np.maximum(n - 2, 1)
    phi = np.maximum(1.0, resid2.sum(1) / df)
    se = np.sqrt(phi * var)
    t = b / se
    p = 2 * stats.t.sf(np.abs(t), df)
    return b, se, t, p, df, phi


def mwu_rows(pct, is_sisc):
    """14-diff-meth Mann-Whitney on per-sample percentages, row-wise."""
    out = np.full(pct.shape[0], np.nan)
    for i, row in enumerate(pct):
        a, c = row[~is_sisc], row[is_sisc]
        a, c = a[~np.isnan(a)], c[~np.isnan(c)]
        if len(a) >= MIN_SAMPLES and len(c) >= MIN_SAMPLES:
            out[i] = stats.mannwhitneyu(a, c, alternative="two-sided").pvalue
    return out


def test_matrix(nmod, nvalid, ncpg, min_cpgs=MIN_CPGS):
    """Apply per-sample filters, then both tests. Returns a DataFrame aligned to the rows."""
    is_sisc = np.array([s in SAMPLES_SISCOWET for s in SAMPLES])
    passed = (nvalid >= MIN_NVALID) & (ncpg >= min_cpgs)
    nm = np.where(passed, nmod, np.nan)
    nv = np.where(passed, nvalid, np.nan)
    pct = 100.0 * nm / nv
    n_lean = passed[:, ~is_sisc].sum(1)
    n_sisc = passed[:, is_sisc].sum(1)
    b, se, t, p, df, phi = quasibinomial_two_group(nm, nv, is_sisc)
    with np.errstate(all="ignore"):
        lean_mean = np.nanmean(pct[:, ~is_sisc], axis=1)
        sisc_mean = np.nanmean(pct[:, is_sisc], axis=1)
    res = pd.DataFrame({
        "n_lean": n_lean, "n_siscowet": n_sisc,
        "lean_mean": lean_mean, "siscowet_mean": sisc_mean,
        "meth_diff": sisc_mean - lean_mean,           # + = hyper in siscowet
        "log_odds_ratio": b, "dispersion": phi, "t": t, "df": df, "pvalue": p,
        "pvalue_mwu": mwu_rows(pct, is_sisc),
    })
    for i, s in enumerate(SAMPLES):
        res[f"pct_{s}"] = pct[:, i]
    testable = (n_lean >= MIN_SAMPLES) & (n_sisc >= MIN_SAMPLES)
    return res, testable


# =============================================================================
# Part A: feature-level native test
# =============================================================================
def part_a(counts, fn):
    print("\n" + "=" * 70 + "\nPart A: feature-level count-based test\n" + "=" * 70, flush=True)
    keep = shared_single_copy()
    print(f"  genes single-copy in both assemblies: {len(keep):,}")

    per_sample = {}
    for eco, samples in ECO_SAMPLES.items():
        fbed = feature_bed(eco, keep)
        for s in samples:
            df = map_counts(fbed, counts[s], COUNTS / f"{s}.features.tsv")
            per_sample[s] = df.set_index("name")[["nmod", "nvalid", "ncpg"]]

    names = sorted(set.intersection(*(set(d.index) for d in per_sample.values())))
    nmod = np.column_stack([per_sample[s].loc[names, "nmod"].to_numpy(float) for s in SAMPLES])
    nvalid = np.column_stack([per_sample[s].loc[names, "nvalid"].to_numpy(float) for s in SAMPLES])
    ncpg = np.column_stack([per_sample[s].loc[names, "ncpg"].to_numpy(float) for s in SAMPLES])

    long = pd.DataFrame({"feature": names})
    for i, s in enumerate(SAMPLES):
        long[f"nmod_{s}"], long[f"nvalid_{s}"], long[f"ncpg_{s}"] = nmod[:, i], nvalid[:, i], ncpg[:, i]
    long.to_csv(OUT / "native_feature_counts.tsv.gz", sep="\t", index=False,
                compression="gzip", float_format="%.0f")

    res, testable = test_matrix(nmod, nvalid, ncpg)
    res.insert(0, "feature", names)
    res.insert(1, "gene_id", [n.split("|")[0] for n in names])
    res.insert(2, "region", [n.split("|")[1] for n in names])
    res = res[testable].reset_index(drop=True)
    for region in ("body", "promoter"):         # FDR within each feature class
        m = res["region"] == region
        res.loc[m, "qvalue"] = multipletests(res.loc[m, "pvalue"], method="fdr_bh")[1]
    res["symbol"] = res["gene_id"].map(lambda g: fn.get(g, ("", ""))[0])
    res["product"] = res["gene_id"].map(lambda g: fn.get(g, ("", ""))[1])
    res["direction"] = np.where(res["meth_diff"] > 0, "hyper_siscowet", "hypo_siscowet")
    res["significant_fdr"] = (res["qvalue"] < Q_VALUE_THRESHOLD) & \
                             (res["meth_diff"].abs() >= DIFF_THRESHOLD)
    res = res.sort_values("pvalue")
    res.to_csv(OUT / "dmf_results.tsv.gz", sep="\t", index=False, compression="gzip",
               float_format="%.6g")
    res[res["significant_fdr"]].to_csv(OUT / "dmf_significant.tsv", sep="\t", index=False,
                                       float_format="%.6g")
    return keep, res


# =============================================================================
# Part B: reference DMRs placed on native coordinates
# =============================================================================
def place_ref_dmrs():
    """Map each reference DMR (+flanks) to both purged assemblies; project the core."""
    dmrs = pd.read_csv(REF_DMRS, sep="\t")
    dmrs.columns = [c.lstrip("#") for c in dmrs.columns]
    sizes = dict(pd.read_csv(f"{REF_FA}.fai", sep="\t", header=None, usecols=[0, 1]).values)
    qfa = COUNTS / "refdmr.flanked.fa"
    regions = COUNTS / "refdmr.regions.txt"
    meta = {}
    with open(regions, "w") as f:
        for r in dmrs.itertuples(index=False):
            fs = max(0, r.start - DMR_FLANK)
            fe = min(sizes[r.chrom], r.end + DMR_FLANK)
            f.write(f"{r.chrom}:{fs + 1}-{fe}\n")
            # core offsets inside the flanked query (0-based, half-open)
            meta[f"{r.chrom}:{fs + 1}-{fe}"] = (r.name, r.start - fs, r.end - fs, fe - fs)
    if not nonempty(qfa):
        sh(f"'{SAMTOOLS}' faidx '{REF_FA}' -r '{regions}' > '{qfa}'")

    placed = {}
    for eco, fa in ECO_FA.items():
        paf = COUNTS / f"refdmr.to_{eco}.paf"
        if not nonempty(paf):
            sh(f"'{MINIMAP2}' -x asm5 -t {THREADS} --secondary=no '{fa}' '{qfa}' > '{paf}'")
        hits = {}
        with open(paf) as f:
            for line in f:
                p = line.split("\t")
                q, qlen, qs, qe, strand = p[0], int(p[1]), int(p[2]), int(p[3]), p[4]
                t, ts, te, mapq = p[5], int(p[7]), int(p[8]), int(p[11])
                if "tp:A:P" not in line:
                    continue
                hits.setdefault(q, []).append((q, qlen, qs, qe, strand, t, ts, te, mapq))
        for q, hs in hits.items():
            if len(hs) != 1:                    # split / ambiguous placement -> skip
                continue
            q, qlen, qs, qe, strand, t, ts, te, mapq = hs[0]
            dmr_id, cs, ce, _ = meta[q]
            if mapq < MIN_MAPQ or (qe - qs) / qlen < MIN_QUERY_FRAC or cs < qs or ce > qe:
                continue
            if strand == "+":
                ns, ne = ts + (cs - qs), te - (qe - ce)
            else:
                ns, ne = ts + (qe - ce), te - (cs - qs)
            if ne > ns:
                placed[(dmr_id, eco)] = (t, ns, ne, mapq)
    return dmrs, placed


def mapq_class(q):
    return "0" if q == 0 else "1-4" if q < 5 else "5-19" if q < 20 else "20-59" if q < 60 else "60"


def paf_mapq_profile(paf, label, eco):
    """MAPQ class counts of primary hits covering >= MIN_QUERY_FRAC of the query."""
    best = {}
    with open(paf) as f:
        for line in f:
            p = line.split("\t")
            if "tp:A:P" not in line or (int(p[3]) - int(p[2])) / int(p[1]) < MIN_QUERY_FRAC:
                continue
            best[p[0]] = max(best.get(p[0], -1), int(p[11]))
    prof = pd.Series([mapq_class(q) for q in best.values()]).value_counts()
    return [{"set": label, "assembly": eco, "mapq_class": k, "n": int(prof.get(k, 0)),
             "frac": prof.get(k, 0) / max(1, len(best))}
            for k in ("0", "1-4", "5-19", "20-59", "60")]


def placement_control():
    """Random reference windows placed like the DMRs: the genome-wide MAPQ baseline."""
    sizes = COUNTS / "refchroms.sizes"
    rfa = COUNTS / "random.clean.fa"
    if not nonempty(rfa):
        sh(f"grep '^NC_0523' '{REF_FA}.fai' | cut -f1,2 > '{sizes}'")
        sh(f"'{BEDTOOLS}' random -l {RANDOM_LEN} -n {int(N_RANDOM * 1.2)} -seed 23 -g '{sizes}' "
           f"| sort -k1,1 -k2,2n | awk '{{print $1\":\"$2+1\"-\"$3}}' > '{COUNTS}/random.regions.txt'")
        sh(f"'{SAMTOOLS}' faidx '{REF_FA}' -r '{COUNTS}/random.regions.txt' > '{COUNTS}/random.fa'")
        recs = open(COUNTS / "random.fa").read().split(">")[1:]
        keep = [r for r in recs
                if "".join(r.split("\n")[1:]).upper().count("N") < 0.05 * RANDOM_LEN][:N_RANDOM]
        rfa.write_text("".join(">" + r for r in keep))
    rows = []
    for eco, fa in ECO_FA.items():
        paf = COUNTS / f"random.to_{eco}.paf"
        if not nonempty(paf):
            sh(f"'{MINIMAP2}' -x asm5 -t {THREADS} --secondary=no '{fa}' '{rfa}' > '{paf}'")
        rows += paf_mapq_profile(COUNTS / f"refdmr.to_{eco}.paf", "reference_dmr", eco)
        rows += paf_mapq_profile(paf, "random_window", eco)
    prof = pd.DataFrame(rows)
    prof.to_csv(OUT / "ref_dmr_placement_mapq.tsv", sep="\t", index=False, float_format="%.4f")
    return prof


def part_b(counts):
    print("\n" + "=" * 70 + "\nPart B: reference DMRs on native coordinates\n" + "=" * 70,
          flush=True)
    dmrs, placed = place_ref_dmrs()
    for eco in ECO_SAMPLES:
        bed = COUNTS / f"refdmr.on_{eco}.bed"
        with open(bed, "w") as f:
            for (dmr_id, e), (t, ns, ne, _) in placed.items():
                if e == eco:
                    f.write(f"{t}\t{ns}\t{ne}\t{dmr_id}\n")
    both = [d for d in dmrs["name"] if (d, "lean") in placed and (d, "siscowet") in placed]
    print(f"  reference DMRs: {len(dmrs)}  placed uniquely on lean: "
          f"{sum(1 for k in placed if k[1] == 'lean')}  siscowet: "
          f"{sum(1 for k in placed if k[1] == 'siscowet')}  both: {len(both)}")

    per_sample = {}
    for eco, samples in ECO_SAMPLES.items():
        for s in samples:
            df = map_counts(COUNTS / f"refdmr.on_{eco}.bed", counts[s],
                            COUNTS / f"refdmr.{s}.tsv")
            per_sample[s] = df.drop_duplicates("name").set_index("name")
    nmod = np.column_stack([per_sample[s]["nmod"].reindex(both).to_numpy(float) for s in SAMPLES])
    nvalid = np.column_stack([per_sample[s]["nvalid"].reindex(both).to_numpy(float) for s in SAMPLES])
    ncpg = np.column_stack([per_sample[s]["ncpg"].reindex(both).to_numpy(float) for s in SAMPLES])
    # DMRs were called on >=3 CpGs; require >=2 covered natively (a SNP can remove a CpG)
    res, testable = test_matrix(np.nan_to_num(nmod), np.nan_to_num(nvalid), np.nan_to_num(ncpg),
                                min_cpgs=2)
    res.insert(0, "dmr_id", both)
    res["native_testable"] = testable
    for eco in ECO_SAMPLES:
        res[f"{eco}_locus"] = [f"{placed[(d, eco)][0]}:{placed[(d, eco)][1]}-{placed[(d, eco)][2]}"
                               for d in both]
    ref = dmrs.rename(columns={"name": "dmr_id", "meth_diff": "ref_meth_diff",
                               "direction": "ref_direction", "n_cpgs": "ref_n_cpgs"})
    ref = ref[["dmr_id", "chrom", "start", "end", "ref_n_cpgs", "ref_meth_diff", "ref_direction"]]
    out = ref.merge(res, on="dmr_id", how="left")
    out["placed_both"] = out["dmr_id"].isin(both)
    out["native_testable"] = out["native_testable"].astype("boolean").fillna(False).astype(bool)
    out["sign_agrees"] = out["native_testable"] & (np.sign(out["meth_diff"]) ==
                                                   np.sign(out["ref_meth_diff"]))
    out["native_confirmed"] = out["sign_agrees"] & (out["pvalue"] < P_VALUE_THRESHOLD) & \
                              (out["meth_diff"].abs() >= DIFF_THRESHOLD)
    out.to_csv(OUT / "ref_dmr_native_loci.tsv", sep="\t", index=False, float_format="%.6g")
    return out, placement_control()


# =============================================================================
# Report
# =============================================================================
def report(n_keep, feat, loci, prof):
    lines = ["Ecotype-native differential methylation, count-based (23.11)",
             "=" * 62, "",
             f"Part A -- features on {n_keep:,} genes single-copy in both assemblies",
             f"  per-sample filter: >= {MIN_NVALID} pooled calls and >= {MIN_CPGS} CpGs; "
             f">= {MIN_SAMPLES} samples per ecotype", ""]
    for region in ("body", "promoter"):
        r = feat[feat["region"] == region]
        sig = r[r["significant_fdr"]]
        lines += [f"  {region}: tested {len(r):,}",
                  f"    quasi-binomial p<0.05: {(r['pvalue'] < 0.05).sum():,}   "
                  f"FDR<{Q_VALUE_THRESHOLD}: {(r['qvalue'] < Q_VALUE_THRESHOLD).sum():,}   "
                  f"FDR<{Q_VALUE_THRESHOLD} & |diff|>={DIFF_THRESHOLD:.0f}: {len(sig):,} "
                  f"(hyper {int((sig['direction'] == 'hyper_siscowet').sum())} / "
                  f"hypo {int((sig['direction'] == 'hypo_siscowet').sum())})",
                  f"    Mann-Whitney p<0.05 & |diff|>=15 (23.9 rule): "
                  f"{((r['pvalue_mwu'] < 0.05) & (r['meth_diff'].abs() >= 15)).sum():,}",
                  f"    median dispersion: {r['dispersion'].median():.2f}", ""]
    t = loci[loci["native_testable"]]
    lines += ["Part B -- the 302 reference DMRs (14-diff-meth) on native coordinates",
              f"  placed uniquely on both assemblies: {int(loci['placed_both'].sum())}",
              f"  testable natively: {len(t)}"]
    if len(t):
        rho, rp = stats.spearmanr(t["ref_meth_diff"], t["meth_diff"])
        lines += [f"  sign agreement with reference diff: {int(t['sign_agrees'].sum())}/{len(t)} "
                  f"({t['sign_agrees'].mean():.0%}; 50 % expected by chance)",
                  f"  Spearman(ref diff, native diff): rho={rho:.2f} (p={rp:.1e})",
                  f"  confirmed (same sign, p<0.05, |diff|>={DIFF_THRESHOLD:.0f}): "
                  f"{int(t['native_confirmed'].sum())}",
                  f"  native |diff| median {t['meth_diff'].abs().median():.1f} vs reference "
                  f"|diff| median {t['ref_meth_diff'].abs().median():.1f}"]
    lines += ["  placement MAPQ (primary hit >= 80 % of query): fraction MAPQ 0 / MAPQ 60"]
    for eco in ECO_SAMPLES:
        f = prof[prof["assembly"] == eco].set_index(["set", "mapq_class"])["frac"]
        lines.append(f"    {eco:8s} reference DMRs {f[('reference_dmr', '0')]:.0%} / "
                     f"{f[('reference_dmr', '60')]:.0%}   random windows "
                     f"{f[('random_window', '0')]:.0%} / {f[('random_window', '60')]:.0%}")
    lines += ["", "Top 25 differentially methylated features (FDR, by p):"]
    for r in feat[feat["significant_fdr"]].head(25).itertuples():
        lines.append(f"  {r.gene_id:24s} {r.region:8s} diff={r.meth_diff:+6.1f}  "
                     f"q={r.qvalue:.1e}  {str(r.symbol)[:14]:14s} {str(r.product)[:40]}")
    text = "\n".join(lines) + "\n"
    (OUT / "native_meth_counts_report.txt").write_text(text)
    print("\n" + text)


def main():
    COUNTS.mkdir(parents=True, exist_ok=True)
    fn = {}
    with open(ANN / "gene_function_table.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fn[r["gene_id"]] = (r.get("symbol", ""), r.get("product", ""))

    print("=" * 70 + "\nStep 1: all-coverage CpG counts\n" + "=" * 70, flush=True)
    counts = {s: all_counts(s) for s in SAMPLES}

    keep, feat = part_a(counts, fn)
    loci, prof = part_b(counts)
    report(len(keep), feat, loci, prof)
    print(f"Outputs -> {OUT}")


if __name__ == "__main__":
    main()
