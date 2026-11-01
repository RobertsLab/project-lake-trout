#!/bin/bash
#SBATCH --job-name=p4-braker
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32
#SBATCH --mem=120G --time=4-00:00:00
#SBATCH --array=0-1
#SBATCH --output=%x-%A_%a.out --error=%x-%A_%a.err
#
# Phase 4 step 4 -- BRAKER per ecotype on its OWN masked assembly.
# Array task 0 = lean, 1 = siscowet.   ARM=protein (default) or ARM=rnaseq.
source "$(dirname "$0")/config.sh"

ECO="${ECOTYPES[$SLURM_ARRAY_TASK_ID]}"
MASKED="${REP}/${ECO}/${ECO}.purged.fa.masked"
OUT="${BRK}/${ARM}/${ECO}"
mkdir -p "$OUT"

[[ -s ~/.gm_key ]] || { echo "ERROR: GeneMark key missing/empty at ~/.gm_key"; exit 1; }

EVIDENCE=(--prot_seq="${IN}/proteins.combined.faa")
if [[ "$ARM" == "rnaseq" ]]; then
  EVIDENCE+=(--bam="${RNA}/${ECO}.merged.bam")
fi

apptainer exec \
  --bind "${SCRATCH}:${SCRATCH}" \
  --bind "${AUGUSTUS_CONFIG_PATH}:/augustus_config" \
  --env AUGUSTUS_CONFIG_PATH=/augustus_config \
  "${BRAKER_SIF}" braker.pl \
    --genome="${MASKED}" \
    "${EVIDENCE[@]}" \
    --species="snam_${ECO}_${ARM}" \
    --workingdir="${OUT}" \
    --threads=32 \
    --gff3

echo "DONE: ${OUT}/braker.gff3"
