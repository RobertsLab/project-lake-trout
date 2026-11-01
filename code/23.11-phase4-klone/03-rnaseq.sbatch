#!/bin/bash
#SBATCH --job-name=p4-rnaseq
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G --time=08:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Phase 4 step 3 -- PRJNA316738 fetch + spliced alignment.  SECONDARY ARM ONLY.
#
# CAVEAT, repeated here because it determines how much to trust the output:
# these are 72 x 36 bp SINGLE-end liver-only libraries (Illumina GA IIx, 2016).
# At 36 bp almost no read spans an intron with enough flank on both sides to yield
# a confident junction, so the intron hints BRAKER actually wants will be sparse.
# HISAT2 (not STAR) is used: it handles short reads better and needs far less RAM
# on a 2.6 Gb genome.
#
# Ecotype assignment comes from SraRunTable 'Ecotype': lean libraries -> lean
# assembly, siscowet -> siscowet. Never cross them; that is the whole point of
# native annotation.
source "$(dirname "$0")/config.sh"
source "${CONDA%/bin/conda}/etc/profile.d/conda.sh"
conda activate "${PROGRAMS}/envs/hisat2-2.2.1-env"

cd "$RNA"

# Run accessions per ecotype, derived from data/SraRunTable.csv (field 'Ecotype';
# note the upstream typo 'Sicowet' in the non-parasitized siscowet rows).
# Regenerate with:
#   python3 -I -c "import csv;[print(r['Run'],'lean' if 'Lean' in r['Ecotype'] else 'siscowet') for r in csv.DictReader(open('data/SraRunTable.csv'))]"
[[ -f runs.tsv ]] || { echo "ERROR: write runs.tsv (run<TAB>ecotype) first"; exit 1; }

while read -r RUN ECO; do
  [[ -f "${RUN}.fastq.gz" ]] || fasterq-dump --threads 16 --outdir . "$RUN"
  gzip -f "${RUN}.fastq" 2>/dev/null || true
done < runs.tsv

for ECO in "${ECOTYPES[@]}"; do
  MASKED="${REP}/${ECO}/${ECO}.purged.fa.masked"
  hisat2-build -p 16 "$MASKED" "${RNA}/${ECO}_idx"
  mapfile -t RUNS < <(awk -v e="$ECO" '$2==e{print $1}' runs.tsv)
  for RUN in "${RUNS[@]}"; do
    hisat2 -p 16 -x "${RNA}/${ECO}_idx" -U "${RUN}.fastq.gz" \
      | samtools sort -@ 4 -o "${RNA}/${ECO}.${RUN}.bam" -
  done
  samtools merge -@ 16 -f "${RNA}/${ECO}.merged.bam" "${RNA}/${ECO}".SRR*.bam
  samtools index "${RNA}/${ECO}.merged.bam"
  # Report the alignment + junction yield -- if junctions are as sparse as expected,
  # that is the quantitative basis for preferring the protein arm.
  echo "== ${ECO} ==" >> "${RNA}/alignment_summary.txt"
  samtools flagstat "${RNA}/${ECO}.merged.bam" >> "${RNA}/alignment_summary.txt"
done
