# Quick Start Guide: PacBio Data Summary

This guide helps you quickly generate a summary report of PacBio Revio sequencing data.

## Prerequisites

- Python 3.11 or higher
- Optional: `samtools` for detailed BAM analysis
- Optional: `pysam` Python package for BAM analysis

## Quick Commands

### 1. View the Pre-Generated Summary

A template summary report has already been generated:

```bash
# View the report
cat analyses/04-pacbio/sequencing_data_summary.md

# Or open in your markdown viewer
```

### 2. Generate a New Summary Report

If you have access to the data or want to customize the report:

#### Option A: Generate Template Only (No Data Required)

```bash
cd code/04-pacbio
python create_sequencing_summary.py --generate-template
```

#### Option B: Analyze Local Data Directory

If you have downloaded the data locally:

```bash
cd code/04-pacbio
python create_sequencing_summary.py \
  --data-dir /path/to/pacbio/data \
  --output ../../analyses/04-pacbio/sequencing_data_summary.md
```

#### Option C: Full Analysis with BAM Statistics

For detailed read statistics (slower, requires samtools or pysam):

```bash
cd code/04-pacbio
python create_sequencing_summary.py \
  --data-dir /path/to/pacbio/data \
  --analyze-bams \
  --max-bams 5 \
  --output ../../analyses/04-pacbio/sequencing_data_summary.md
```

## Data Location

The PacBio Revio sequencing data is hosted at:
```
https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/
```

**Note**: This URL may require VPN or specific network access to reach.

## Downloading Data

If you have access to the data repository, you can download files using:

```bash
# Create data directory
mkdir -p data/pacbio-reads

# Download using wget (example)
wget -r -np -nH --cut-dirs=3 \
  https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/ \
  -P data/pacbio-reads/

# Or use curl to list files first
curl -s https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/ | \
  grep -oP 'href="\K[^"]+' | \
  while read file; do
    wget https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/$file \
      -P data/pacbio-reads/
  done
```

## Understanding the Report

The generated report includes:

### 1. **Overview Section**
- Project description
- Data source information
- Sequencing platform details

### 2. **File Inventory**
- List of all data files found
- Grouped by file type (BAM, XML, index files, etc.)

### 3. **Sequencing Statistics** (if BAMs analyzed)
- Total number of reads
- Mean/median read lengths
- N50 values
- Quality scores
- File sizes

### 4. **Sample Information**
- Lean lake trout subspecies details
- Siscowet lake trout subspecies details
- Biological context

### 5. **Analysis Workflows**
- Recommended commands for alignment
- Variant calling pipelines
- Methylation analysis approaches

### 6. **References**
- Links to related data (NCBI BioProject)
- Tool documentation
- Reference genome information

## Common Issues

### Cannot Access Data URL

If you cannot access the owl.fish.washington.edu URL:
- You may need VPN access to University of Washington network
- Contact the RobertsLab for access instructions
- Use `--generate-template` to create a report without data access

### Missing Dependencies

If you get import errors:
```bash
# Install pysam for Python-based BAM analysis
cd code/04-pacbio
uv add pysam

# Or install samtools system-wide
# Ubuntu/Debian:
sudo apt-get install samtools

# macOS:
brew install samtools
```

### Script Not Executable

If you get permission errors:
```bash
chmod +x code/04-pacbio/create_sequencing_summary.py
```

## Next Steps

After generating the summary:

1. **Review the report** to understand the data structure
2. **Download specific samples** you want to analyze
3. **Run alignment workflows** using `align_hifi_pbmm2.py`
4. **Perform downstream analysis** (variants, methylation, etc.)

For more details, see:
- `code/04-pacbio/README.md` - Detailed script documentation
- `analyses/04-pacbio/README.md` - Analysis outputs overview
- `code/05-pacbio-align.Rmd` - Example alignment workflow

## Support

For questions or issues:
- Open an issue in the project-lake-trout repository
- Contact the RobertsLab team
- Check existing documentation in the `notes/` directory
