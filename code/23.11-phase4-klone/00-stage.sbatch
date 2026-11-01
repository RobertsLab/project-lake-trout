#!/bin/bash
#SBATCH --job-name=p4-stage
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=16G --time=02:00:00
#SBATCH --output=%x-%j.out --error=%x-%j.err
#
# Phase 4 step 0 -- stage inputs onto Klone scratch and freeze them with md5s.
source "$(dirname "$0")/config.sh"

# Purged assemblies + existing Liftoff models from the local box.
for eco in "${ECOTYPES[@]}"; do
  rsync -av --progress \
    "${LOCAL_HOST}:${LOCAL_SYN}/${eco}.purged.fa" \
    "${LOCAL_HOST}:${LOCAL_SYN}/${eco}.purged.fa.fai" \
    "${LOCAL_HOST}:${LOCAL_SYN}/${eco}.purged.liftoff.gff3_polished" \
    "${IN}/"
done

# Reference proteins (same species -- strongest single protein evidence set).
wget -P "$IN" https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/432/855/GCF_016432855.1_SaNama_1.0/GCF_016432855.1_SaNama_1.0_protein.faa.gz
gunzip -f "${IN}/GCF_016432855.1_SaNama_1.0_protein.faa.gz"

# OrthoDB v11 Actinopterygii partition.
wget -P "$IN" https://bioinf.uni-greifswald.de/bioinf/partitioned_odb11/Actinopterygii.fa.gz
gunzip -f "${IN}/Actinopterygii.fa.gz"
cat "${IN}/Actinopterygii.fa" "${IN}/GCF_016432855.1_SaNama_1.0_protein.faa" > "${IN}/proteins.combined.faa"

# BRAKER3 container + a WRITABLE AUGUSTUS config (the usual first-run failure).
apptainer build "${BRAKER_SIF}" docker://teambraker/braker3:latest
apptainer exec "${BRAKER_SIF}" cp -r /opt/Augustus/config "${AUGUSTUS_CONFIG_PATH}"
chmod -R u+w "${AUGUSTUS_CONFIG_PATH}"

# Freeze. Compare these against the local .md5 files in 21.1-gene-anchored-synteny/.
cd "$IN" && md5sum ./*.fa ./*.gff3_polished ./*.faa > inputs.md5
echo "STAGED. Verify inputs.md5 matches the local box before proceeding."
