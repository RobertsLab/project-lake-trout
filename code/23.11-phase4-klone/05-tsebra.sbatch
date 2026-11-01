#!/bin/bash
#SBATCH --job-name=p4-tsebra
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=32G --time=02:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Phase 4 step 5 -- merge Liftoff (Track A) with BRAKER (Track B), tag provenance,
# and isolate the BRAKER-unique models.
#
# Liftoff models are the trusted backbone: same-species transfer carrying reference
# gene IDs, and Phases 2-3 are already keyed to those IDs. BRAKER contributes only
# what Liftoff structurally cannot -- genes absent from SaNama_1.0.
source "$(dirname "$0")/config.sh"
source "${CONDA%/bin/conda}/etc/profile.d/conda.sh"
conda activate "${PROGRAMS}/envs/tsebra-1.1.2-env"

for ECO in "${ECOTYPES[@]}"; do
  OUT="${MRG}/${ECO}"; mkdir -p "$OUT"
  LIFTOFF="${IN}/${ECO}.purged.liftoff.gff3_polished"
  BRAKER="${BRK}/${ARM}/${ECO}/braker.gff3"

  tsebra.py -g "${BRAKER}" -k "${LIFTOFF}" \
    -e "${BRK}/${ARM}/${ECO}/hintsfile.gff" \
    -c "${PROGRAMS}/envs/tsebra-1.1.2-env/share/tsebra/config/keep_ab_initio.cfg" \
    -o "${OUT}/merged.gff3"
  # -k keeps ALL Liftoff models unconditionally; BRAKER models join only on evidence.

  # BRAKER-unique = merged models with no positional overlap against Liftoff.
  conda activate "${PROGRAMS}/envs/bedtools-2.31-env"
  awk '$3=="gene"' "${OUT}/merged.gff3" | sort -k1,1 -k4,4n > "${OUT}/merged.genes.gff3"
  awk '$3=="gene"' "${LIFTOFF}"          | sort -k1,1 -k4,4n > "${OUT}/liftoff.genes.gff3"
  bedtools intersect -v -a "${OUT}/merged.genes.gff3" -b "${OUT}/liftoff.genes.gff3" \
    > "${OUT}/braker_unique.gff3"
  conda activate "${PROGRAMS}/envs/tsebra-1.1.2-env"

  echo "${ECO}: $(wc -l < "${OUT}/braker_unique.gff3") BRAKER-unique gene models"
done
