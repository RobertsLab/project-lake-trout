#!/usr/bin/env bash
# Shared configuration for Phase 4 on Klone.
# Sourced by every 23.11-phase4-klone/*.sbatch script.
#
# VERIFY BEFORE FIRST SUBMISSION:
#   hyakalloc                    # account, partition, core/RAM ceilings
#   ls -l ~/.gm_key              # GeneMark licence must not be expired
set -euo pipefail

# --- Slurm -------------------------------------------------------------------
# `compute` is the allocated partition; `ckpt` is free but preemptible and REQUIRES
# the -recoverDir logic in 01-repeatmodeler.sbatch. Override per-job with --partition.
export ACCOUNT="${ACCOUNT:-srlab}"
export PARTITION="${PARTITION:-compute}"

# --- Paths -------------------------------------------------------------------
export SCRATCH="/gscratch/scrubbed/${USER}/lake-trout-phase4"   # 21-day auto-purge
export KEEP="/gscratch/srlab/${USER}/lake-trout-phase4/final"   # survives the purge
export PROGRAMS="/srlab/programs/miniforge3-24.7.1-0"
export CONDA="${PROGRAMS}/bin/conda"

export IN="${SCRATCH}/00_inputs"
export REP="${SCRATCH}/01_repeats"
export RNA="${SCRATCH}/02_rnaseq"
export BRK="${SCRATCH}/03_braker"
export MRG="${SCRATCH}/04_merge"
export SPC="${SCRATCH}/05_specific"
export FUN="${SCRATCH}/06_function"

# --- Source data on the local box (for 00-stage rsync) -----------------------
export LOCAL_HOST="${LOCAL_HOST:-sr320@raven}"   # adjust to the real hostname
export LOCAL_SYN="/home/shared/16TB_HDD_01/sr320/github/project-lake-trout/output/21.1-gene-anchored-synteny"

# --- Ecotypes ----------------------------------------------------------------
# Purged assemblies, NOT the raw hifiasm primaries -- see the plan doc for why.
export ECOTYPES=(lean siscowet)

# --- Evidence arm ------------------------------------------------------------
# protein : OrthoDB Actinopterygii + SaNama RefSeq proteins only.  RECOMMENDED --
#           PRJNA316738 is 36 bp single-end liver-only and supplies almost no
#           usable intron hints (see plan doc).
# rnaseq  : adds the PRJNA316738 BAMs.  Run as a comparison arm, not the backbone.
export ARM="${ARM:-protein}"

export BRAKER_SIF="${IN}/braker3.sif"
export AUGUSTUS_CONFIG_PATH="${IN}/augustus_config"   # MUST be writable; copied in 00-stage

mkdir -p "$IN" "$REP" "$RNA" "$BRK" "$MRG" "$SPC" "$FUN" "$KEEP"
