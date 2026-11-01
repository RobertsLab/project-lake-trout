#!/bin/bash
#SBATCH --job-name=p4-specific
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=64G --time=06:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Phase 4 step 6 -- reciprocal genomic-absence screen.  THIS IS A CORRECTNESS STEP.
#
# A BRAKER model present in lean and absent from siscowet's ANNOTATION is not
# evidence of an ecotype-specific gene -- it is equally consistent with prediction
# failure on the other assembly. The only defensible claim is GENOMIC absence.
#
# So: take each ecotype's BRAKER-unique proteins, map them back to the OTHER
# ecotype's genome with miniprot (splice-aware protein-to-genome), and keep only
# candidates with no credible locus there. Expect this to discard a large fraction
# of naive calls -- that is the screen working.
source "$(dirname "$0")/config.sh"
source "${CONDA%/bin/conda}/etc/profile.d/conda.sh"
conda activate "${PROGRAMS}/envs/miniprot-0.13-env"

IDENT_MIN=0.60   # protein identity above which we call the locus "present"
COV_MIN=0.60     # fraction of the query protein aligned

for ECO in "${ECOTYPES[@]}"; do
  OTHER=$([[ "$ECO" == "lean" ]] && echo siscowet || echo lean)
  OUT="${SPC}/${ECO}"; mkdir -p "$OUT"

  agat_sp_extract_sequences.pl -g "${MRG}/${ECO}/braker_unique.gff3" \
    -f "${IN}/${ECO}.purged.fa" -p -o "${OUT}/braker_unique.faa"

  miniprot -t 16 --gff "${IN}/${OTHER}.purged.fa" "${OUT}/braker_unique.faa" \
    > "${OUT}/vs_${OTHER}.gff"

  # Any query with a hit clearing both thresholds is PRESENT in the other genome
  # and therefore not ecotype-specific.
  awk -v i="$IDENT_MIN" -v c="$COV_MIN" '
    /^##PAF/ { ident=$13/$12; cov=($4-$3)/$2; if (ident>=i && cov>=c) print $1 }
  ' "${OUT}/vs_${OTHER}.gff" | sort -u > "${OUT}/present_in_other.ids"

  grep '^>' "${OUT}/braker_unique.faa" | sed 's/^>//;s/ .*//' | sort -u \
    > "${OUT}/all_candidates.ids"
  comm -23 "${OUT}/all_candidates.ids" "${OUT}/present_in_other.ids" \
    > "${OUT}/ecotype_specific.ids"

  echo "${ECO}: $(wc -l < "${OUT}/all_candidates.ids") candidates -> $(wc -l < "${OUT}/ecotype_specific.ids") survive genomic-absence screen"
done

# Cross-reference against Phase 2's independent reciprocal-PAV calls. Agreement
# between two independent methods is the >=2-lines-of-evidence standard Phase 5 applies.
echo "NEXT: intersect ecotype_specific.ids with analyses/23-reciprocal-pav/ calls."
