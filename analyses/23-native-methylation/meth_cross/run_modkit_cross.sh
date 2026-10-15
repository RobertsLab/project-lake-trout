#!/bin/bash
# 23.12 frame control: modkit on the Phase-2 (23.4) cross-ecotype alignments, so each ecotype frame
# has both ecotypes' reads on identical coordinates.
set -euo pipefail
cd /home/shared/16TB_HDD_01/sr320/github/project-lake-trout
MODKIT=code/04-pacbio/tools/modkit
SYN=output/21.1-gene-anchored-synteny
OUT=analyses/23-native-methylation/meth_cross
for s in bc2041 bc2068 bc2069 bc2070; do
  [ -s $OUT/$s.on_siscowet.CpG.bedmethyl ] || $MODKIT pileup analyses/23-reciprocal-pav/align/$s.to_siscowet.sorted.bam $OUT/$s.on_siscowet.CpG.bedmethyl --cpg --combine-strands --ref $SYN/siscowet.purged.fa --threads 38 > $OUT/$s.modkit.log 2>&1
done
for s in bc2071 bc2072 bc2073 bc2096; do
  [ -s $OUT/$s.on_lean.CpG.bedmethyl ] || $MODKIT pileup analyses/23-reciprocal-pav/align/$s.to_lean.sorted.bam $OUT/$s.on_lean.CpG.bedmethyl --cpg --combine-strands --ref $SYN/lean.purged.fa --threads 38 > $OUT/$s.modkit.log 2>&1
done
echo ALL_DONE
