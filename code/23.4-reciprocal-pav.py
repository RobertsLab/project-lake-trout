#!/usr/bin/env python3
"""
23.4-reciprocal-pav.py  (Phase 2, Step 1 of code/23-next-phase-research-plan.md)

RECIPROCAL PAV. The reference-based PAV (code/15-diff-pav.py) mapped every ecotype's reads to the
lean-background reference, so siscowet-specific signal was divergence-inflated. Here we run the same
coverage/CIGAR logic as code/11-pav.Rmd BOTH DIRECTIONS, each ecotype's HiFi reads against the
OTHER ecotype's purged assembly, on native coordinates:

    lean reads    (bc2041,bc2068,bc2069,bc2070)  ->  siscowet.purged.fa
    siscowet reads(bc2071,bc2072,bc2073,bc2096)  ->  lean.purged.fa

Directional interpretation (BEDs on the TARGET genome's native contig coordinates):
  * A region of the TARGET genome with ZERO coverage in ALL FOUR querying-ecotype samples, and
    FLANKED by covered windows (interior -- not a contig end / unmappable edge), is sequence that is
    present in the TARGET ecotype but absent from the QUERYING ecotype. Because it is measured
    against the target ecotype's OWN assembled genome (definitely present) rather than a divergent
    reference, it is a confident presence/absence call:
        lean reads -> siscowet  absent-consensus  =>  SISCOWET-specific present (siscowet coords)
        siscowet   -> lean      absent-consensus  =>  LEAN-specific present     (lean coords)
  * CIGAR insertions (>=50 bp) shared by >= MIN_INS_SAMPLES querying-ecotype samples are sequence in
    the querying ecotype's reads that is NOVEL relative to the target genome:
        lean reads -> siscowet  insertions  =>  LEAN-specific novel sequence (positioned on siscowet)
        siscowet   -> lean      insertions  =>  SISCOWET-specific novel sequence (positioned on lean)

Each specific-present set is intersected with the matching native Liftoff gene BED
(output/21.1-gene-anchored-synteny/*.purged.liftoff.genes.bed) to name affected genes. SV / reference
cross-validation is done at the GENE level in 23.6 (avoids native<->anchored coordinate bridging).

Outputs (analyses/23-reciprocal-pav/):
  align/{sample}.to_{target}.sorted.bam(.bai), .flagstat.txt
  cov/{sample}.regions.bed.gz, calls/{sample}.{absent,present}_merged.bed, calls/{sample}.insertions.bed
  siscowet_specific.present_regions.on_siscowet.bed   lean_specific.present_regions.on_lean.bed
  lean_specific.novel_seq.on_siscowet.bed             siscowet_specific.novel_seq.on_lean.bed
  {ecotype}_specific.genes.tsv
  reciprocal_pav_summary.csv

Compute: 8 pbmm2 alignments of multi-GB HiFi read sets to ~2.5-2.8 Gb genomes -- multi-hour.
Sequential, resumable (skips a stage whose output already exists).

Author: Generated for project-lake-trout
Date: 2026-07-04
"""

import glob
import gzip
import subprocess
import sys
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================
BASE = Path(__file__).resolve().parent.parent
READS_DIR = BASE / "data" / "pacbio-reads"
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
OUT = BASE / "analyses" / "23-reciprocal-pav"
ALIGN = OUT / "align"
COV = OUT / "cov"
CALLS = OUT / "calls"
LOGS = OUT / "logs"

# Tools
PBMM2 = BASE / "code" / "04-pacbio" / "tools" / "pbmm2"
SYRI_BIN = Path("/home/sr320/miniconda3/envs/syri_env/bin")
SAMTOOLS = "/usr/bin/samtools"
BEDTOOLS = "/usr/bin/bedtools"
MOSDEPTH = str(SYRI_BIN / "mosdepth")
GAWK = "/usr/bin/gawk"

THREADS = "46"
MOSDEPTH_THREADS = "8"

# ecotype -> (samples, target assembly, target label, native gene bed)
LEAN_SAMPLES = ["bc2041", "bc2068", "bc2069", "bc2070"]
SISCO_SAMPLES = ["bc2071", "bc2072", "bc2073", "bc2096"]
LEAN_FA = SYN / "lean.purged.fa"
SISCO_FA = SYN / "siscowet.purged.fa"
LEAN_GENES = SYN / "lean.purged.liftoff.genes.bed"
SISCO_GENES = SYN / "siscowet.purged.liftoff.genes.bed"

# Two directions. querying ecotype's reads -> target ecotype's genome.
DIRECTIONS = [
    {"query_eco": "lean", "samples": LEAN_SAMPLES,
     "target_eco": "siscowet", "target_fa": SISCO_FA, "target_genes": SISCO_GENES,
     "tlabel": "siscowet"},
    {"query_eco": "siscowet", "samples": SISCO_SAMPLES,
     "target_eco": "lean", "target_fa": LEAN_FA, "target_genes": LEAN_GENES,
     "tlabel": "lean"},
]

WINDOW = 1000          # mosdepth window (matches 11-pav)
MIN_INS_SIZE = 50      # CIGAR insertion/deletion min size (matches 11-pav)
MIN_ABSENT_BP = 1000   # minimum consensus-absent region to report as present/absent event
FLANK_BP = 2000        # a consensus-absent region must have covered flanks within this distance
MIN_INS_SAMPLES = 3    # querying-ecotype samples that must share a novel-sequence insertion cluster


# =============================================================================
# Helpers
# =============================================================================
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


def chrom_sizes(target_fa, dest):
    if not nonempty(dest):
        sh(f"cut -f1,2 '{target_fa}.fai' > '{dest}'")
    return dest


# =============================================================================
# Stage 1: align + coverage + per-sample calls
# =============================================================================
def build_index(target_fa):
    mmi = ALIGN / (Path(target_fa).stem + ".mmi")
    if not nonempty(mmi):
        print(f"[index] building pbmm2 CCS index for {Path(target_fa).name}", flush=True)
        sh(f"'{PBMM2}' index --preset CCS '{target_fa}' '{mmi}'",
           log=LOGS / "pbmm2_index.log")
    return mmi


def align_sample(sample, mmi, tlabel):
    bam = ALIGN / f"{sample}.to_{tlabel}.sorted.bam"
    if nonempty(bam):
        print(f"[align] {sample}->{tlabel} exists, skipping", flush=True)
        return bam
    reads = reads_bam(sample)
    log = LOGS / f"{sample}.to_{tlabel}.pbmm2.log"
    print(f"[align] {sample} -> {tlabel}  ({Path(reads).name})", flush=True)
    sh(f"'{PBMM2}' align --preset CCS --sort -j {THREADS} "
       f"'{mmi}' '{reads}' '{bam}'", log=log)
    sh(f"'{SAMTOOLS}' index '{bam}'")
    sh(f"'{SAMTOOLS}' flagstat '{bam}' > '{ALIGN}/{sample}.to_{tlabel}.flagstat.txt'")
    return bam


def coverage_and_calls(sample, bam, tlabel):
    """mosdepth windows -> absent/present merged BEDs; CIGAR -> insertions/deletions (per 11-pav)."""
    regions = COV / f"{sample}.to_{tlabel}.regions.bed.gz"
    if not nonempty(regions):
        prefix = COV / f"{sample}.to_{tlabel}"
        sh(f"'{MOSDEPTH}' --threads {MOSDEPTH_THREADS} --no-per-base --by {WINDOW} "
           f"'{prefix}' '{bam}'", log=LOGS / f"{sample}.to_{tlabel}.mosdepth.log")

    absent = CALLS / f"{sample}.to_{tlabel}.absent_merged.bed"
    present = CALLS / f"{sample}.to_{tlabel}.present_merged.bed"
    if not nonempty(absent):
        sh(f"zcat '{regions}' | awk -v OFS='\\t' '$4==0{{print $1,$2,$3}}' | "
           f"'{BEDTOOLS}' merge -i - > '{absent}'")
    if not nonempty(present):
        sh(f"zcat '{regions}' | awk -v OFS='\\t' '$4>0{{print $1,$2,$3}}' | "
           f"'{BEDTOOLS}' merge -i - > '{present}'")

    ins = CALLS / f"{sample}.to_{tlabel}.insertions.bed"
    dele = CALLS / f"{sample}.to_{tlabel}.deletions.bed"
    if not nonempty(ins) or not nonempty(dele):
        # Parse CIGAR once; emit I>=min and D>=min. gawk 3-arg match (as in 11-pav).
        awk_prog = (
            r'{chrom=$3; if(chrom=="*")next; ref=$4; cig=$6;'
            r' while(match(cig,/([0-9]+)([MIDNSHP=X])/,a)){len=a[1]+0; op=a[2];'
            f'  if(op=="I"&&len>={MIN_INS_SIZE}) print chrom"\\t"ref"\\t"ref+1"\\tI\\t"len > "/dev/stdout";'
            f'  if(op=="D"&&len>={MIN_INS_SIZE}) print chrom"\\t"ref"\\t"ref+len"\\tD\\t"len > "/dev/stderr";'
            r'  if(op=="M"||op=="D"||op=="N"||op=="="||op=="X") ref+=len;'
            r'  cig=substr(cig,RSTART+RLENGTH)}}'
        )
        raw_i = CALLS / f"{sample}.to_{tlabel}.ins_raw.bed"
        raw_d = CALLS / f"{sample}.to_{tlabel}.del_raw.bed"
        sh(f"'{SAMTOOLS}' view '{bam}' | '{GAWK}' '{awk_prog}' "
           f"1> '{raw_i}' 2> '{raw_d}'")
        for raw, merged in [(raw_i, ins), (raw_d, dele)]:
            sh(f"sort -k1,1 -k2,2n '{raw}' | '{BEDTOOLS}' merge -i - -c 4 -o count "
               f"> '{merged}'")
            Path(raw).unlink(missing_ok=True)
    return absent, present, ins


# =============================================================================
# Stage 2: consensus per direction
# =============================================================================
def multiinter_consensus(beds, out_bed, need):
    """Regions covered by >= `need` of the input BEDs (bedtools multiinter num column)."""
    sorted_beds = []
    for b in beds:
        sb = Path(str(b) + ".srt")
        sh(f"sort -k1,1 -k2,2n '{b}' > '{sb}'")
        sorted_beds.append(str(sb))
    inputs = " ".join(f"'{b}'" for b in sorted_beds)
    sh(f"'{BEDTOOLS}' multiinter -i {inputs} | awk -v OFS='\\t' '$4>={need}{{print $1,$2,$3}}' | "
       f"'{BEDTOOLS}' merge -i - > '{out_bed}'")
    for sb in sorted_beds:
        Path(sb).unlink(missing_ok=True)


def interior_absent(absent_all, present_all, sizes, out_bed):
    """Keep consensus-absent regions >= MIN_ABSENT_BP that have consensus-present sequence within
    FLANK_BP on BOTH sides (interior -- excludes contig ends / unmappable edges)."""
    named = Path(str(out_bed) + ".named")
    sh(f"awk -v OFS='\\t' '($3-$2)>={MIN_ABSENT_BP}{{print $1,$2,$3,\"a\"NR}}' "
       f"'{absent_all}' > '{named}'")
    lflank = Path(str(out_bed) + ".lf")
    rflank = Path(str(out_bed) + ".rf")
    sh(f"'{BEDTOOLS}' flank -i '{named}' -g '{sizes}' -l {FLANK_BP} -r 0 > '{lflank}'")
    sh(f"'{BEDTOOLS}' flank -i '{named}' -g '{sizes}' -l 0 -r {FLANK_BP} > '{rflank}'")
    lsup = Path(str(out_bed) + ".ls")
    rsup = Path(str(out_bed) + ".rs")
    # names whose left/right flank overlaps a consensus-present region
    sh(f"'{BEDTOOLS}' intersect -a '{lflank}' -b '{present_all}' -u | cut -f4 | sort -u > '{lsup}'")
    sh(f"'{BEDTOOLS}' intersect -a '{rflank}' -b '{present_all}' -u | cut -f4 | sort -u > '{rsup}'")
    both = Path(str(out_bed) + ".both")
    sh(f"comm -12 '{lsup}' '{rsup}' > '{both}'")
    # keep named regions whose id is supported on both sides
    sh(f"awk 'NR==FNR{{k[$1]=1;next}} ($4 in k)' '{both}' '{named}' | "
       f"cut -f1-3 | sort -k1,1 -k2,2n > '{out_bed}'")
    for f in [named, lflank, rflank, lsup, rsup, both]:
        Path(f).unlink(missing_ok=True)


def consensus_insertions(ins_beds, out_bed):
    """Insertion clusters shared by >= MIN_INS_SAMPLES querying-ecotype samples.
    Tag each sample's insertions with the sample id, merge nearby (within 100 bp), and keep clusters
    hitting enough distinct samples."""
    tagged = Path(str(out_bed) + ".tagged")
    with open(tagged, "w") as o:
        for b in ins_beds:
            sid = Path(b).name.split(".")[0]
            with open(b) as f:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 3:
                        o.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{sid}\n")
    sh(f"sort -k1,1 -k2,2n '{tagged}' | "
       f"'{BEDTOOLS}' merge -i - -d 100 -c 4 -o count_distinct,distinct | "
       f"awk -v OFS='\\t' '$4>={MIN_INS_SAMPLES}{{print $1,$2,$3,$5}}' > '{out_bed}'")
    Path(tagged).unlink(missing_ok=True)


def annotate_genes(regions_bed, genes_bed, out_tsv):
    """Genes (native coords) overlapping the specific-present regions."""
    sh(f"'{BEDTOOLS}' intersect -a '{genes_bed}' -b '{regions_bed}' -u | "
       f"sort -k1,1 -k2,2n > '{out_tsv}.tmp'")
    n = 0
    with open(f"{out_tsv}.tmp") as f, open(out_tsv, "w") as o:
        o.write("chrom\tstart\tend\tgene_id\tscore\tstrand\n")
        for line in f:
            o.write(line)
            n += 1
    Path(f"{out_tsv}.tmp").unlink(missing_ok=True)
    return n


def count_lines(p):
    return sum(1 for _ in open(p)) if nonempty(p) else 0


def total_bp(p):
    if not nonempty(p):
        return 0
    s = 0
    with open(p) as f:
        for line in f:
            c = line.split("\t")
            if len(c) >= 3:
                s += int(c[2]) - int(c[1])
    return s


# =============================================================================
# Main
# =============================================================================
def main():
    for d in (ALIGN, COV, CALLS, LOGS):
        d.mkdir(parents=True, exist_ok=True)

    summary = []
    for D in DIRECTIONS:
        qeco, teco, tlabel = D["query_eco"], D["target_eco"], D["tlabel"]
        print("=" * 70)
        print(f"DIRECTION  {qeco} reads -> {teco} assembly")
        print("=" * 70)

        sizes = chrom_sizes(D["target_fa"], OUT / f"{tlabel}.chrom.sizes")
        mmi = build_index(D["target_fa"])

        absents, presents, inserts = [], [], []
        for s in D["samples"]:
            bam = align_sample(s, mmi, tlabel)
            a, p, i = coverage_and_calls(s, bam, tlabel)
            absents.append(a)
            presents.append(p)
            inserts.append(i)

        # consensus (all 4 querying samples)
        absent_all = CALLS / f"{qeco}_reads_on_{tlabel}.absent_all4.bed"
        present_all = CALLS / f"{qeco}_reads_on_{tlabel}.present_all4.bed"
        multiinter_consensus(absents, absent_all, need=len(D["samples"]))
        multiinter_consensus(presents, present_all, need=len(D["samples"]))

        # TARGET-ecotype-specific present regions (interior consensus-absent on target coords)
        present_regions = OUT / f"{teco}_specific.present_regions.on_{tlabel}.bed"
        interior_absent(absent_all, present_all, sizes, present_regions)

        # QUERYING-ecotype-specific novel sequence (consensus insertions, on target coords)
        novel = OUT / f"{qeco}_specific.novel_seq.on_{tlabel}.bed"
        consensus_insertions(inserts, novel)

        # genes overlapping the target-ecotype-specific present regions (native genes)
        genes_tsv = OUT / f"{teco}_specific.genes.tsv"
        n_genes = annotate_genes(present_regions, D["target_genes"], genes_tsv)

        summary.append({
            "direction": f"{qeco}_reads_to_{teco}_asm",
            "target_specific_ecotype": teco,
            "present_regions": count_lines(present_regions),
            "present_bp": total_bp(present_regions),
            "present_genes": n_genes,
            "querying_specific_ecotype": qeco,
            "novel_seq_clusters": count_lines(novel),
            "novel_seq_bp": total_bp(novel),
        })
        print(f"  [{teco}-specific present] regions={summary[-1]['present_regions']:,} "
              f"bp={summary[-1]['present_bp']:,} genes={n_genes:,}")
        print(f"  [{qeco}-specific novel]   clusters={summary[-1]['novel_seq_clusters']:,} "
              f"bp={summary[-1]['novel_seq_bp']:,}")

    import csv
    with open(OUT / "reciprocal_pav_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        w.writerows(summary)
    print(f"\nDone. Summary -> {OUT/'reciprocal_pav_summary.csv'}")


if __name__ == "__main__":
    main()
