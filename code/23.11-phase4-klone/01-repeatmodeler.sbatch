#!/bin/bash
#SBATCH --job-name=p4-repmod
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32
#SBATCH --mem=120G --time=7-00:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Phase 4 step 1 -- de novo repeat library.
#
# ONE library for BOTH ecotypes: lean and siscowet are the same species, so their
# repeat content is effectively identical. Building twice doubles the most expensive
# step in the phase for no gain, and masking both assemblies against an IDENTICAL
# library also removes a confound where differential masking manufactures apparent
# gene-model differences.
#
# THIS JOB WILL LIKELY EXCEED ONE WALLTIME WINDOW. It is resubmit-safe: on restart it
# finds the newest RM_* round directory and resumes with -recoverDir. Just resubmit
# the same script until the library appears. On the preemptible `ckpt` partition this
# recovery logic is what makes the job survivable at all.
source "$(dirname "$0")/config.sh"
source "${CONDA%/bin/conda}/etc/profile.d/conda.sh"
conda activate "${PROGRAMS}/envs/repeatmodeler-2.0.5-env"

cd "$REP"
SEED=lean   # library built from lean, applied to both

if [[ ! -f "${SEED}.nsq" && ! -f "${SEED}.nin" ]]; then
  BuildDatabase -name "$SEED" "${IN}/${SEED}.purged.fa"
fi

# Resume if a previous round directory exists, otherwise start fresh.
RECOVER=""
LATEST=$(ls -1dt RM_* 2>/dev/null | head -1 || true)
if [[ -n "$LATEST" ]]; then
  echo "Resuming from ${LATEST}"
  RECOVER="-recoverDir ${LATEST}"
else
  echo "No prior round directory -- starting fresh."
fi

RepeatModeler -database "$SEED" -threads 32 -LTRStruct $RECOVER

cp "${SEED}-families.fa" "${REP}/denovo_library.fa"
echo "DONE: ${REP}/denovo_library.fa"
