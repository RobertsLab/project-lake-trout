#!/bin/bash
#SBATCH --job-name=p4-repmask
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32
#SBATCH --mem=64G --time=2-00:00:00
#SBATCH --array=0-1
#SBATCH --output=%x-%A_%a.out --error=%x-%A_%a.err
#
# Phase 4 step 2 -- soft-mask both assemblies against the SAME combined library.
# Array task 0 = lean, 1 = siscowet.
source "$(dirname "$0")/config.sh"
source "${CONDA%/bin/conda}/etc/profile.d/conda.sh"
conda activate "${PROGRAMS}/envs/repeatmasker-4.1.6-env"

ECO="${ECOTYPES[$SLURM_ARRAY_TASK_ID]}"
OUT="${REP}/${ECO}"
mkdir -p "$OUT"

# de novo library (step 1) + curated Salmoniformes from Dfam.
famdb.py -i "$(dirname "$(which RepeatMasker)")/Libraries/famdb" families \
  --format fasta_name --ancestors --descendants 'Salmoniformes' \
  > "${REP}/salmoniformes.fa" || echo "WARN: Dfam export failed; proceeding de-novo-only"
cat "${REP}/denovo_library.fa" "${REP}/salmoniformes.fa" > "${REP}/combined_library.fa"

RepeatMasker -pa 8 -xsmall -gff -lib "${REP}/combined_library.fa" \
  -dir "$OUT" "${IN}/${ECO}.purged.fa"
# -xsmall = SOFT mask (lowercase). BRAKER requires soft, not hard (N) masking.

echo "DONE: ${OUT}/${ECO}.purged.fa.masked"
