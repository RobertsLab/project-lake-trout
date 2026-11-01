#!/bin/bash
#SBATCH --job-name=p4-function
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=32
#SBATCH --mem=64G --time=1-00:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Phase 4 step 7 -- functional annotation of NEW models only, plus BUSCO QC.
# Lifted genes already have functions via analyses/18-annotation/gene_function_table.tsv;
# only BRAKER-unique models need a fresh pass (same approach as 18-build-gene-function-table.py).
source "$(dirname "$0")/config.sh"
source "${CONDA%/bin/conda}/etc/profile.d/conda.sh"

for ECO in "${ECOTYPES[@]}"; do
  OUT="${FUN}/${ECO}"; mkdir -p "$OUT"
  FAA="${SPC}/${ECO}/braker_unique.faa"

  conda activate "${PROGRAMS}/envs/eggnog-mapper-2.1.12-env"
  emapper.py -i "$FAA" -o "${OUT}/eggnog" --cpu 32 --itype proteins \
    --data_dir "${IN}/eggnog_db"

  conda activate "${PROGRAMS}/envs/interproscan-5.66-env"
  sed 's/\*//g' "$FAA" > "${OUT}/clean.faa"   # InterProScan rejects '*' stop characters
  interproscan.sh -i "${OUT}/clean.faa" -f TSV -goterms -cpu 32 -o "${OUT}/interpro.tsv"

  # BUSCO on the FULL merged proteome -- the number that tells you whether the
  # annotation is complete. Run on merged, not on the ecotype-specific subset.
  conda activate "${PROGRAMS}/envs/busco-5.7-env"
  agat_sp_extract_sequences.pl -g "${MRG}/${ECO}/merged.gff3" \
    -f "${IN}/${ECO}.purged.fa" -p -o "${OUT}/merged.faa"
  busco -i "${OUT}/merged.faa" -l actinopterygii_odb10 -m protein \
    -c 32 -o "busco_${ECO}" --out_path "$OUT"
done

# Persist everything worth keeping OUT of scrubbed before the 21-day purge.
for ECO in "${ECOTYPES[@]}"; do
  D="${KEEP}/${ECO}"; mkdir -p "$D"
  cp "${MRG}/${ECO}/merged.gff3"        "$D/" 2>/dev/null || true
  cp "${MRG}/${ECO}/braker_unique.gff3" "$D/" 2>/dev/null || true
  cp "${SPC}/${ECO}/ecotype_specific.ids" "$D/" 2>/dev/null || true
  cp "${FUN}/${ECO}/eggnog.emapper.annotations" "$D/" 2>/dev/null || true
  cp "${FUN}/${ECO}/interpro.tsv"       "$D/" 2>/dev/null || true
  cp "${FUN}/${ECO}/busco_${ECO}/short_summary"*.txt "$D/" 2>/dev/null || true
done
echo "Products copied to ${KEEP} -- safe from the scrubbed purge."
