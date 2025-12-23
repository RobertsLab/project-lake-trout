# PacBio Analysis Outputs

This directory contains analysis outputs and reports from PacBio Revio sequencing data.

## Contents

### `sequencing_data_summary.md`

Comprehensive summary report of the PacBio Revio sequencing effort for Lake Trout genomics project. This report includes:

- Overview of sequencing data from https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/
- File inventory and organization
- Sequencing statistics (read counts, lengths, quality metrics)
- Sample information (Lean vs Siscowet subspecies)
- Technology overview (PacBio Revio HiFi sequencing)
- Recommended analysis workflows
- References and related resources

**To regenerate or update this report:**
```bash
cd ../../code/04-pacbio
python create_sequencing_summary.py --data-dir /path/to/data --analyze-bams
```

### `alignments/`

Directory for storing aligned BAM files and related outputs from pbmm2 alignment workflow.

## Related Scripts

Analysis scripts are located in `../../code/04-pacbio/`:
- `create_sequencing_summary.py` - Generate sequencing data summary reports
- `align_hifi_pbmm2.py` - Batch align HiFi reads using pbmm2
- See `../../code/04-pacbio/README.md` for detailed documentation

## Data Source

Primary sequencing data is hosted at:
https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/

## Project Context

This analysis is part of the Lake Trout (_Salvelinus namaycush_) genomics project comparing lean and siscowet subspecies. Related data and analyses:
- NCBI BioProject: [PRJNA674328](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA674328)
- Reference genome: GCF_016432855.1 (SaNama_1.0)
- RNAseq differential expression analysis: `../` (parent directory)
