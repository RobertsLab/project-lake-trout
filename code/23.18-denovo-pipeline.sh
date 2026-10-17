#!/bin/bash
# 23.18-denovo-pipeline.sh  (Phase 4 of code/23-next-phase-research-plan.md)
#
# De novo annotation of each purged ecotype assembly on raven. One chain per ecotype, the two
# chains in parallel at THREADS each (default 24, half of raven's 48):
#
#   1. RepeatModeler2 (-LTRStruct) on the purged assembly          -> repeats/{eco}/{eco}-families.fa
#   2. RepeatMasker with that library, soft-masked (-xsmall)       -> repeats/{eco}/{eco}.purged.fa.masked
#   3. HISAT2 index of the masked assembly; the ecotype's own pooled liver RNA-seq aligned to it
#      (lean -> lean, siscowet -> siscowet)                        -> rnaseq/{eco}.bam
#   4. BRAKER3 (RNA-seq BAM + OrthoDB Vertebrata + SaNama RefSeq proteins)
#                                                                  -> braker/{eco}/braker.{gtf,gff3,aa,codingseq}
#   5. BUSCO (protein mode, actinopterygii_odb12) on braker.aa     -> busco/{eco}/
#
# Every stage writes a .done marker and is skipped on re-run, so the script can be restarted after
# an interruption. A failing stage writes "<stage> FAILED" to logs/chain.{eco}.log.
# Note: GeneMark-ETP's masking-penalty step needs >~100 kb of reliable-gene sequence, so a smoke
# test on a contig with few genes can fail there ("input sequence size is too small"); this is a
# smoke-genome artifact, not a pipeline error. Inputs come from 23.17-denovo-setup.sh. Liftoff merge (TSEBRA / overlap rule),
# provenance tags and functional annotation of the BRAKER-unique genes follow in a later script.
#
# Run on raven (detached):
#   setsid nohup bash code/23.18-denovo-pipeline.sh > analyses/23-denovo-annotation/logs/pipeline.log 2>&1 &
#   (optional: bash code/23.18-denovo-pipeline.sh lean   -- one ecotype only)
# Smoke test on a genome subset: P4=<scratch dir> GENOME_DIR=<dir with {eco}.purged.fa> THREADS=8
#
# Author: Generated for project-lake-trout   Date: 2026-10-16
set -euo pipefail

BASE=/home/shared/16TB_HDD_01/sr320/github/project-lake-trout
P4=${P4:-$BASE/analyses/23-denovo-annotation}          # override for a smoke test
SYN=${GENOME_DIR:-$BASE/output/21.1-gene-anchored-synteny}
C=$BASE/analyses/23-denovo-annotation/containers
EVID=$BASE/analyses/23-denovo-annotation/evidence
TETOOLS=$C/tetools_2.00.sif
BRAKER=$C/braker3_v3.1.1.sif
BUSCO=$C/busco_6.1.0.sif
HISAT2=/home/shared/hisat2-2.2.1
SAMTOOLS=/usr/bin/samtools
THREADS=${THREADS:-24}
RM_PA=$(( THREADS / 4 ))          # RepeatMasker -pa: rmblast uses 4 threads per job
ECOS=${1:-"lean siscowet"}

export SINGULARITY_TMPDIR=$P4/tmp APPTAINER_TMPDIR=$P4/tmp TMPDIR=$P4/tmp
SX="singularity exec --cleanenv -B $BASE"
log() { echo "[$(date '+%F %T')] [$1] $2"; }
done_() { [ -e "$1.done" ]; }
mark() { touch "$1.done"; }

chain() {
  local eco=$1
  local fa=$SYN/$eco.purged.fa
  local R=$P4/repeats/$eco  RN=$P4/rnaseq  B=$P4/braker/$eco  Q=$P4/busco/$eco  L=$P4/logs
  mkdir -p "$R" "$RN" "$B" "$Q" "$L"

  # 1. RepeatModeler2 ---------------------------------------------------------------------------
  if ! done_ "$R/repeatmodeler"; then
    log $eco "RepeatModeler2 start"
    ( cd "$R" && $SX --pwd "$R" "$TETOOLS" BuildDatabase -name "$eco" "$fa" > "$L/$eco.builddb.log" 2>&1 \
      && $SX --pwd "$R" "$TETOOLS" RepeatModeler -database "$eco" -threads "$THREADS" -LTRStruct \
         > "$L/$eco.repeatmodeler.log" 2>&1 ) \
      || { log $eco "RepeatModeler FAILED (exit $?)"; return 1; }
    [ -s "$R/$eco-families.fa" ] || { log $eco "RepeatModeler produced no library"; return 1; }
    mark "$R/repeatmodeler"
  fi
  log $eco "repeat families: $(grep -c '^>' "$R/$eco-families.fa")"

  # 2. RepeatMasker -------------------------------------------------------------------------------
  local masked=$R/$eco.purged.fa.masked
  if ! done_ "$R/repeatmasker"; then
    log $eco "RepeatMasker start"
    $SX --pwd "$R" "$TETOOLS" RepeatMasker -pa "$RM_PA" -lib "$R/$eco-families.fa" -xsmall -gff \
      -dir "$R" "$fa" > "$L/$eco.repeatmasker.log" 2>&1 \
      || { log $eco "RepeatMasker FAILED (exit $?)"; return 1; }
    [ -s "$masked" ] || { log $eco "RepeatMasker produced no masked FASTA"; return 1; }
    mark "$R/repeatmasker"
  fi
  log $eco "masked: $(awk '!/^>/{n+=length($0); s+=gsub(/[acgtn]/,"")} END{printf "%.1f%% soft-masked", 100*s/n}' "$masked")"

  # 3. RNA-seq alignment ----------------------------------------------------------------------------
  local idx=$RN/$eco.masked.hisat2
  if ! done_ "$RN/$eco.index"; then
    log $eco "HISAT2 index"
    "$HISAT2/hisat2-build" -p "$THREADS" "$masked" "$idx" > "$L/$eco.hisat2-build.log" 2>&1
    mark "$RN/$eco.index"
  fi
  local bam=$RN/$eco.bam
  if ! done_ "$RN/$eco.align"; then
    log $eco "HISAT2 align own-ecotype RNA-seq"
    "$HISAT2/hisat2" -p "$THREADS" -x "$idx" -U "$EVID/rnaseq/$eco.fastq.gz" \
      --summary-file "$L/$eco.hisat2.summary.txt" 2> "$L/$eco.hisat2.log" \
      | "$SAMTOOLS" sort -@ 8 -m 2G -o "$bam" -
    "$SAMTOOLS" index "$bam"
    mark "$RN/$eco.align"
  fi

  # 4. BRAKER3 ---------------------------------------------------------------------------------------
  if ! done_ "$B/braker"; then
    log $eco "BRAKER3 start"
    if [ ! -d "$B/augustus_config" ]; then
      local cfg
      cfg=$($SX "$BRAKER" bash -c 'echo ${AUGUSTUS_CONFIG_PATH:-/opt/Augustus/config}')
      $SX "$BRAKER" cp -r "$cfg" "$B/augustus_config"
    fi
    $SX --pwd "$B" "$BRAKER" braker.pl --genome="$masked" --bam="$bam" \
      --prot_seq="$EVID/proteins/proteins.fa" --species="snamaycush_${eco}_p4" \
      --threads="$THREADS" --workingdir="$B" --AUGUSTUS_CONFIG_PATH="$B/augustus_config" --gff3 \
      > "$L/$eco.braker.log" 2>&1 \
      || { log $eco "BRAKER3 FAILED (exit $?); see $L/$eco.braker.log and $B/errors/"; return 1; }
    [ -s "$B/braker.aa" ] || { log $eco "BRAKER3 produced no proteins"; return 1; }
    mark "$B/braker"
  fi
  log $eco "BRAKER3 genes: $(awk '$3=="gene"' "$B/braker.gtf" | wc -l), proteins: $(grep -c '^>' "$B/braker.aa")"

  # 5. BUSCO -------------------------------------------------------------------------------------------
  if ! done_ "$Q/busco"; then
    log $eco "BUSCO protein mode"
    $SX --pwd "$Q" "$BUSCO" busco -i "$B/braker.aa" -m proteins -l actinopterygii_odb12 \
      --offline --download_path "$BASE/analyses/23-denovo-annotation/busco_downloads" -c "$THREADS" -o "$eco" --out_path "$Q" -f \
      > "$L/$eco.busco.log" 2>&1 || { log $eco "BUSCO FAILED (exit $?)"; return 1; }
    mark "$Q/busco"
  fi
  log $eco "BUSCO: $(grep -h 'C:' "$Q"/$eco/short_summary*.txt | head -1 | tr -s ' \t')"
  log $eco "chain done"
}

mkdir -p "$P4/logs" "$P4/tmp"
for f in "$TETOOLS" "$BRAKER" "$BUSCO" "$EVID/proteins/proteins.fa"; do
  [ -s "$f" ] || { echo "missing $f -- run code/23.17-denovo-setup.sh first"; exit 1; }
done
pids=()
for eco in $ECOS; do
  [ -s "$EVID/rnaseq/$eco.fastq.gz" ] || { echo "missing $eco RNA-seq"; exit 1; }
  chain "$eco" > "$P4/logs/chain.$eco.log" 2>&1 &
  pids+=($!)
done
status=0
for p in "${pids[@]}"; do wait "$p" || status=1; done
cat "$P4"/logs/chain.*.log
exit $status
