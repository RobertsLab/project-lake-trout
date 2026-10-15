#!/usr/bin/env python3
"""
23.8-native-meth-pipeline.py  (Phase 3, part A, of code/23-next-phase-research-plan.md)

Ecotype-NATIVE methylation calling. Every prior methylation result (analyses/14-diff-meth) mapped
all reads to the lean-background reference, so cross-ecotype methylation was reference-anchored. Here
each ecotype's HiFi reads are aligned to its OWN purged assembly and 5mC is called there, following
the 10-mCG-call.Rmd `modkit pileup` pattern.

    lean reads    (bc2041,bc2068,bc2069,bc2070)  ->  lean.purged.fa
    siscowet reads(bc2071,bc2072,bc2073,bc2096)  ->  siscowet.purged.fa

pbmm2 preserves the MM/ML base-modification tags from the CCS BAMs, so modkit can read 5mC directly.
The purged pbmm2 CCS indexes built in Phase 2 (analyses/23-reciprocal-pav/align/*.mmi) are reused.

Per sample:
  1. pbmm2 align --preset CCS --sort   (MM/ML preserved)
  2. modkit pileup --cpg --combine-strands --ref <own purged.fa>   -> bedMethyl
  3. cov10 %-methylation bedGraph  (col10 Nvalid_cov >= 10; col11 = percent_modified)
  4. per-CpG counts  (chrom,start,end,Nmod=col12,Nvalid=col10) for coverage-weighted gene
     aggregation in 23.9

Outputs (analyses/23-native-methylation/):
  align/{sample}.to_own.sorted.bam(.bai), .flagstat.txt
  meth/{sample}.CpG.bedmethyl            raw modkit bedMethyl
  meth/{sample}.cov10.bedGraph           chrom,start,end,percent (Nvalid_cov>=10)
  meth/{sample}.cpg_counts.bed           chrom,start,end,Nmod,Nvalid  (Nvalid>=10)
  native_meth_summary.csv

Compute: 8 pbmm2 alignments of multi-GB HiFi sets to ~2.5-2.8 Gb genomes -- multi-hour. Run
detached (setsid nohup) and resumable (skips a stage whose output already exists).

Author: Generated for project-lake-trout
Date: 2026-07-07
"""

import csv
import glob
import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
READS_DIR = BASE / "data" / "pacbio-reads"
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
OUT = BASE / "analyses" / "23-native-methylation"
ALIGN = OUT / "align"
METH = OUT / "meth"
LOGS = OUT / "logs"
# reuse the purged CCS indexes built in Phase 2 (23.4); fall back to building if absent
PHASE2_ALIGN = BASE / "analyses" / "23-reciprocal-pav" / "align"

PBMM2 = BASE / "code" / "04-pacbio" / "tools" / "pbmm2"
SYRI_BIN = Path("/home/sr320/miniconda3/envs/syri_env/bin")
SAMTOOLS = "/usr/bin/samtools"
MODKIT = str(BASE / "code" / "04-pacbio" / "tools" / "modkit")

THREADS = "46"
MODKIT_THREADS = "38"
MIN_COV = 10  # Nvalid_cov filter (matches 14-diff-meth cov10)

LEAN_FA = SYN / "lean.purged.fa"
SISCO_FA = SYN / "siscowet.purged.fa"
ECO = {
    "lean":     {"samples": ["bc2041", "bc2068", "bc2069", "bc2070"], "fa": LEAN_FA,
                 "mmi": PHASE2_ALIGN / "lean.purged.mmi"},
    "siscowet": {"samples": ["bc2071", "bc2072", "bc2073", "bc2096"], "fa": SISCO_FA,
                 "mmi": PHASE2_ALIGN / "siscowet.purged.mmi"},
}


def sh(cmd, log=None):
    print("  $", cmd, flush=True)
    if log:
        with open(log, "a") as lf:
            subprocess.run(cmd, shell=True, check=True, stdout=lf, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, shell=True, check=True)


def nonempty(p):
    p = Path(p)
    return p.exists() and p.stat().st_size > 0


def reads_bam(sample):
    hits = sorted(glob.glob(str(READS_DIR / f"*{sample}*.bam")))
    if not hits:
        sys.exit(f"FATAL: no read BAM for {sample} in {READS_DIR}")
    return hits[0]


def get_index(eco):
    """Reuse the Phase-2 purged CCS index; build one if it is missing."""
    mmi = ECO[eco]["mmi"]
    if nonempty(mmi):
        return mmi
    local = ALIGN / (Path(ECO[eco]["fa"]).stem + ".mmi")
    if not nonempty(local):
        print(f"[index] building pbmm2 CCS index for {eco}", flush=True)
        sh(f"'{PBMM2}' index --preset CCS '{ECO[eco]['fa']}' '{local}'",
           log=LOGS / "pbmm2_index.log")
    return local


def align_sample(sample, mmi):
    bam = ALIGN / f"{sample}.to_own.sorted.bam"
    if nonempty(bam):
        print(f"[align] {sample} exists, skipping", flush=True)
        return bam
    reads = reads_bam(sample)
    print(f"[align] {sample} -> own assembly ({Path(reads).name})", flush=True)
    sh(f"'{PBMM2}' align --preset CCS --sort -j {THREADS} '{mmi}' '{reads}' '{bam}'",
       log=LOGS / f"{sample}.pbmm2.log")
    sh(f"'{SAMTOOLS}' index '{bam}'")
    sh(f"'{SAMTOOLS}' flagstat '{bam}' > '{ALIGN}/{sample}.flagstat.txt'")
    return bam


def call_methylation(sample, bam, fa):
    """modkit pileup -> bedMethyl -> cov10 bedGraph + per-CpG counts.

    modkit bedMethyl columns used: 10 = Nvalid_cov, 11 = percent_modified, 12 = Nmod.
    (Verified against a sample's output before the aggregation step in 23.9.)
    """
    bedmethyl = METH / f"{sample}.CpG.bedmethyl"
    if not nonempty(bedmethyl):
        sh(f"'{MODKIT}' pileup '{bam}' '{bedmethyl}' --cpg --combine-strands "
           f"--ref '{fa}' --threads {MODKIT_THREADS}",
           log=LOGS / f"{sample}.modkit.log")

    cov10 = METH / f"{sample}.cov10.bedGraph"
    counts = METH / f"{sample}.cpg_counts.bed"
    # col4 == "m" keeps only 5mC (reads also carry 6mA 'a', though --cpg already excludes it);
    # col10 Nvalid_cov >= MIN_COV, col11 percent_modified, col12 Nmod.
    if not nonempty(cov10):
        sh(f"awk -v OFS='\\t' -v m={MIN_COV} '$4==\"m\" && $10>=m{{print $1,$2,$3,$11}}' "
           f"'{bedmethyl}' > '{cov10}'")
    if not nonempty(counts):
        # chrom,start,end,Nmod,Nvalid  (coverage-weighted gene aggregation in 23.9)
        sh(f"awk -v OFS='\\t' -v m={MIN_COV} '$4==\"m\" && $10>=m{{print $1,$2,$3,$12,$10}}' "
           f"'{bedmethyl}' > '{counts}'")
    return cov10, counts


def n_lines(p):
    return sum(1 for _ in open(p)) if nonempty(p) else 0


def main():
    for d in (ALIGN, METH, LOGS):
        d.mkdir(parents=True, exist_ok=True)
    if not nonempty(MODKIT):
        sys.exit(f"FATAL: modkit not found at {MODKIT}. Install it first "
                 f"(conda install -n syri_env -c bioconda modkit).")

    # PASS 1 -- align all samples first (the expensive step) so a downstream calling issue can
    # never strand a pending alignment.
    bams = {}
    for eco, cfg in ECO.items():
        print("=" * 70)
        print(f"ALIGN {eco}: reads -> own purged assembly")
        print("=" * 70)
        if not nonempty(f"{cfg['fa']}.fai"):
            sh(f"'{SAMTOOLS}' faidx '{cfg['fa']}'")
        mmi = get_index(eco)
        for s in cfg["samples"]:
            bams[s] = align_sample(s, mmi)

    # PASS 2 -- call methylation per sample (per-sample failure is non-fatal).
    summary = []
    for eco, cfg in ECO.items():
        for s in cfg["samples"]:
            try:
                cov10, counts = call_methylation(s, bams[s], cfg["fa"])
                n_cpg = n_lines(counts)
            except subprocess.CalledProcessError as e:
                print(f"  [{s}] modkit/calling FAILED: {e}", flush=True)
                n_cpg = -1
            print(f"  [{s}] CpGs with cov>={MIN_COV}: {n_cpg:,}", flush=True)
            summary.append({"ecotype": eco, "sample": s,
                            "cpgs_cov10": n_cpg, "bam": bams[s].name})

    with open(OUT / "native_meth_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    print(f"\nDone. Summary -> {OUT/'native_meth_summary.csv'}")
    print("Next: python3 code/23.9-native-meth-dmg.py")


if __name__ == "__main__":
    main()
