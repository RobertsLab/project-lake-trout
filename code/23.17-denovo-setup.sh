#!/bin/bash
# 23.17-denovo-setup.sh  (Phase 4, step 0, of code/23-next-phase-research-plan.md)
#
# Stages the tools and evidence for de novo annotation (BRAKER3) of the two purged ecotype
# assemblies on raven. Every step is idempotent (skipped when its output exists).
#
#   containers/  braker3_v3.1.1.sif   docker://teambraker/braker3:v3.1.1   (BRAKER3, GeneMark-ETP,
#                                                                         AUGUSTUS, ProtHint, TSEBRA)
#                tetools_2.00.sif     docker://dfam/tetools:2.00           (RepeatModeler2, RepeatMasker)
#                busco_6.1.0.sif      docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1
#   evidence/rnaseq/   PRJNA316738 liver RNA-seq (72 runs, 36 bp single-end GAIIx), pooled per ecotype:
#                      lean.fastq.gz (36 runs), siscowet.fastq.gz (36 runs)
#   evidence/proteins/ Vertebrata.fa (OrthoDB v12 partition, BRAKER team server)
#                      SaNama_1.0.refseq.protein.faa (NCBI Datasets, GCF_016432855.1)
#                      proteins.fa = both, the BRAKER3 --prot_seq input
#   busco_downloads/   actinopterygii_odb12 lineage (BUSCO QC)
#
# Run on raven:  bash code/23.17-denovo-setup.sh
# Author: Generated for project-lake-trout   Date: 2026-10-16
set -euo pipefail

BASE=/home/shared/16TB_HDD_01/sr320/github/project-lake-trout
P4=$BASE/analyses/23-denovo-annotation
SRA=/home/shared/sratoolkit.3.0.2-ubuntu64/bin
DATASETS=/home/shared/datasets
THREADS=16

mkdir -p "$P4"/{containers,evidence/rnaseq/sra,evidence/proteins,busco_downloads,logs,tmp}
export SINGULARITY_CACHEDIR=$P4/tmp/singularity APPTAINER_CACHEDIR=$P4/tmp/singularity
export SINGULARITY_TMPDIR=$P4/tmp APPTAINER_TMPDIR=$P4/tmp TMPDIR=$P4/tmp
log() { echo "[$(date '+%F %T')] $*"; }

# ---- containers -------------------------------------------------------------------------------
pull() {  # $1 sif  $2 uri
  if [ ! -s "$P4/containers/$1" ]; then
    log "pull $2"
    singularity pull --dir "$P4/containers" "$1" "$2"
  fi
}
pull braker3_v3.1.1.sif docker://teambraker/braker3:v3.1.1
pull tetools_2.00.sif   docker://dfam/tetools:2.00
pull busco_6.1.0.sif    docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_1

# ---- RNA-seq (PRJNA316738) ----------------------------------------------------------------------
RUNS=$P4/evidence/rnaseq/runs.tsv
if [ ! -s "$RUNS" ]; then
  python3 - "$BASE/data/SraRunTable.csv" > "$RUNS" <<'EOF'
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    eco = "lean" if r["Ecotype"].lower().startswith("lean") else "siscowet"
    print(f'{r["Run"]}\t{eco}\t{r["Ecotype"]}')
EOF
fi
for eco in lean siscowet; do
  out=$P4/evidence/rnaseq/$eco.fastq.gz
  [ -s "$out" ] && continue
  for run in $(awk -v e=$eco '$2==e{print $1}' "$RUNS"); do
    fq=$P4/evidence/rnaseq/sra/$run.fastq
    if [ ! -s "$fq" ] && [ ! -s "$fq.gz" ]; then
      log "fetch $run ($eco)"
      "$SRA/prefetch" -O "$P4/evidence/rnaseq/sra" "$run" > /dev/null
      "$SRA/fasterq-dump" -e $THREADS -t "$P4/tmp" -O "$P4/evidence/rnaseq/sra" \
        "$P4/evidence/rnaseq/sra/$run/$run.sra" > /dev/null
      gzip "$fq"
      rm -rf "$P4/evidence/rnaseq/sra/$run"
    fi
  done
  log "pool $eco RNA-seq"
  cat $(awk -v e=$eco -v d="$P4/evidence/rnaseq/sra" '$2==e{print d"/"$1".fastq.gz"}' "$RUNS") > "$out.tmp"
  mv "$out.tmp" "$out"
done

# ---- proteins -----------------------------------------------------------------------------------
PR=$P4/evidence/proteins
if [ ! -s "$PR/Vertebrata.fa" ]; then
  log "OrthoDB v12 Vertebrata"
  curl -sSfL -o "$PR/Vertebrata.fa.gz" https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Vertebrata.fa.gz
  gunzip "$PR/Vertebrata.fa.gz"
fi
if [ ! -s "$PR/SaNama_1.0.refseq.protein.faa" ]; then
  log "SaNama_1.0 RefSeq proteins"
  (cd "$P4/tmp" && "$DATASETS" download genome accession GCF_016432855.1 --include protein \
     --filename sanama_protein.zip > /dev/null && unzip -o -q sanama_protein.zip -d sanama_protein)
  cp "$P4"/tmp/sanama_protein/ncbi_dataset/data/GCF_016432855.1/protein.faa "$PR/SaNama_1.0.refseq.protein.faa"
fi
if [ ! -s "$PR/proteins.fa" ]; then
  cat "$PR/Vertebrata.fa" "$PR/SaNama_1.0.refseq.protein.faa" > "$PR/proteins.fa"
fi

# ---- BUSCO lineage ------------------------------------------------------------------------------
if [ ! -d "$P4/busco_downloads/lineages/actinopterygii_odb12" ]; then
  log "BUSCO actinopterygii_odb12"
  singularity exec -B "$P4" "$P4/containers/busco_6.1.0.sif" \
    busco --download actinopterygii_odb12 --download_path "$P4/busco_downloads" > "$P4/logs/busco_download.log" 2>&1
fi

# ---- summary ------------------------------------------------------------------------------------
{
  echo "containers:"; ls -la "$P4/containers" | awk 'NR>1{print "  "$5"\t"$9}'
  for eco in lean siscowet; do
    n=$(zcat "$P4/evidence/rnaseq/$eco.fastq.gz" | awk 'END{print NR/4}')
    echo "rnaseq $eco: $n reads"
  done
  echo "proteins: $(grep -c '^>' "$PR/proteins.fa") sequences (OrthoDB Vertebrata $(grep -c '^>' "$PR/Vertebrata.fa"), SaNama RefSeq $(grep -c '^>' "$PR/SaNama_1.0.refseq.protein.faa"))"
  echo "busco lineage: $(ls -d "$P4"/busco_downloads/lineages/* 2>/dev/null | xargs -n1 basename)"
} | tee "$P4/setup_summary.txt"
log "setup done"
