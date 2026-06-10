# Lake Trout Genomics: Comparative Analysis of Lean and Siscowet Ecotypes

## Overview

Comprehensive genomic analysis of _Salvelinus namaycush_ (lake trout) comparing two distinct ecotypes: **lean** and **siscowet**. This repository contains multiple integrated analyses including:

- **RNAseq differential expression analysis** using parasitized/non-parasitized liver tissue
- **PacBio HiFi DNA methylation profiling** and differential methylation analysis
- **Presence-Absence Variation (PAV) analysis** to identify structural genomic variations
- **Interactive genome browser** for visualizing genomic features

### Reference Genome
- **Assembly**: GCF_016432855.1 (SaNama_1.0)
- **Species**: _Salvelinus namaycush_ (Lake Trout)
- **Source**: [NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_016432855.1/)

---

## Key Analyses

### 1. RNAseq Differential Expression

Analysis of liver RNAseq data from parasitized and non-parasitized samples ([NCBI BioProject PRJNA316738](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA316738)) to identify:
- Differentially expressed genes (DEGs) between subspecies
- Differentially expressed transcripts (DETs) and alternative isoforms
- Expression differences related to parasite status

**Key Results:**
- 202 differentially expressed transcripts (p < 0.05)
- Analysis performed using [Ballgown](https://github.com/alyssafrazee/ballgown)
- See [`analyses/README.md`](analyses/README.md) for detailed results

**Analysis Files:**
- [`code/01-ballgown-analysis.Rmd`](code/01-ballgown-analysis.Rmd) - Primary differential expression analysis
- [`code/02-gene-explore.qmd`](code/02-gene-explore.qmd) - Gene-specific exploration

### 2. PacBio HiFi DNA Methylation Analysis

Whole-genome DNA methylation profiling using PacBio HiFi sequencing with 5mC modification calling:
- Sample-level methylation profiles for both ecotypes
- Differential methylation analysis between lean and siscowet
- Identification of Differentially Methylated Regions (DMRs)

**Key Results:**
- 540,040 CpG sites tested
- 4,440 significant differentially methylated cytosines (DMCs, p < 0.05)
- 302 Differentially Methylated Regions (DMRs)
  - 20 hypermethylated in siscowet
  - 282 hypomethylated in siscowet

**Analysis Files:**
- [`code/04-pacbio/`](code/04-pacbio/) - PacBio workflow (alignment, QC, methylation calling)
- [`code/10-mCG-call.Rmd`](code/10-mCG-call.Rmd) - Methylation calling
- [`code/13.3-hifiasm-differential-methylation-plan.md`](code/13.3-hifiasm-differential-methylation-plan.md) - Plan for extending differential methylation analysis to ecotype-specific hifiasm assemblies
- [`code/14-diff-meth.Rmd`](code/14-diff-meth.Rmd) - Differential methylation analysis
- [`code/14-diff-meth.py`](code/14-diff-meth.py) - Python implementation for DMR identification

### 3. Presence-Absence Variation (PAV)

Genome-wide structural variation analysis identifying insertions and deletions between ecotypes:
- Coverage-based detection of absent regions (deletions)
- CIGAR-based detection of novel insertions
- Ecotype-specific and shared structural variants

**Key Results:**
- **Lean-specific**: 996,228 variants (770,891 insertions + 225,337 deletions)
- **Siscowet-specific**: 1,332,705 variants (1,086,799 insertions + 245,906 deletions)
- **Shared**: 878,372 variants common to both ecotypes

**Analysis Files:**
- [`code/11-pav.Rmd`](code/11-pav.Rmd) - PAV identification analysis
- [`code/12-pav.py`](code/12-pav.py) - Python implementation
- [`code/15-diff-pav.py`](code/15-diff-pav.py) - Differential PAV analysis

### 4. Interactive Genome Browser

Web-based genome browser for exploring PAV and methylation data across the genome.

**Features:**
- Visualize ecotype-specific insertions and deletions
- View differential methylation tracks
- Gene annotations with interactive navigation
- Mobile-responsive design

**Live Demo:** [https://sr320.github.io/project-lake-trout/genome-browser/](https://sr320.github.io/project-lake-trout/genome-browser/)

**Documentation:** [`genome-browser/README.md`](genome-browser/README.md)

---

## Repository Structure

```
project-lake-trout/
├── code/                    # Analysis scripts and notebooks
│   ├── 01-ballgown-analysis.Rmd       # RNAseq differential expression
│   ├── 02-gene-explore.qmd            # Gene exploration
│   ├── 04-pacbio/                     # PacBio HiFi analysis workflow
│   ├── 05-pacbio-align.Rmd            # PacBio alignment
│   ├── 07-pacbio-QC.Rmd               # PacBio quality control
│   ├── 10-mCG-call.Rmd                # Methylation calling
│   ├── 11-pav.Rmd                     # PAV analysis
│   ├── 13.3-hifiasm-differential-methylation-plan.md # Plan for hifiasm-based differential methylation
│   ├── 14-diff-meth.Rmd/py            # Differential methylation
│   └── 15-diff-pav.py                 # Differential PAV
├── data/                    # Raw data and metadata
│   ├── SraRunTable.csv                # RNAseq sample information
│   ├── ballgown-metadata.csv          # Ballgown metadata
│   └── *.bed                          # Gene annotations
├── analyses/                # Analysis outputs and results
│   ├── DEG-*.csv                      # Differentially expressed genes
│   ├── DET-*.csv                      # Differentially expressed transcripts
│   ├── 04-pacbio/                     # PacBio analysis outputs
│   ├── 14-diff-meth/                  # Methylation results
│   └── 15-diff-pav/                   # PAV results
├── genome-browser/          # Interactive genome browser
│   ├── index.html                     # Browser interface
│   ├── prepare_data.py                # Data preparation script
│   └── data/                          # Browser data files
└── figures/                 # Generated figures and plots
```

See README files in each subdirectory for detailed information about specific analyses.

---

## Sample Information

### RNAseq Samples
- **Lean Nonparasitized**: NPLL32, NPLL34, NPLL44, NPLL46, NPLL56, NPLL61
- **Lean Parasitized**: PLL20, PLL31, PLL43, PLL55, PLL59, PLL62
- **Siscowet Nonparasitized**: NPSL15, NPSL24, NPSL29, NPSL36, NPSL50, NPSL58
- **Siscowet Parasitized**: PSL13, PSL16, PSL35, PSL49, PSL53, PSL63

See [`data/SraRunTable.csv`](data/SraRunTable.csv) for complete RNAseq sample metadata.

### PacBio HiFi Samples (for methylation and PAV analysis)
- **Lean**: bc2041, bc2068, bc2069, bc2070
- **Siscowet**: bc2071, bc2072, bc2073, bc2096

---

## Pre-Analysis Data Processing

The following data processing steps were performed prior to analyses in this repository:

- [SRA RNAseq data retrieval, trimming, and QC](https://robertslab.github.io/sams-notebook/2022/07/06/SRA-Data-S.namaycush-SRA-BioProject-PRJNA674328-Download-and-QC.html)

- [RNAseq alignment and splice site identification using Hisat2 and StringTie](https://robertslab.github.io/sams-notebook/2022/08/10/Splice-Site-Identification-S.namaycush-Liver-Parasitized-and-Non-Parasitized-SRA-RNAseq-Using-Hisat2-Stingtie-with-Genome-GCF_016432855.1.html)

- [Convert NCBI genome GFF to BED](https://robertslab.github.io/sams-notebook/2022/08/18/Data-Wrangling-Convert-S.namaycush-NCBI-GFF-to-genes-only-BED-file-for-Use-in-Ballgown-Analysis.html)

---

## Technologies Used

- **R/RStudio**: Statistical analysis and visualization
  - Ballgown: Differential expression analysis
  - tidyverse: Data manipulation
- **Python**: Data processing and analysis pipelines
  - pysam, pandas, numpy: Data manipulation
  - modbampy: Modified base parsing
- **PacBio Tools**: HiFi sequencing analysis
  - pbmm2: Read alignment
  - pb-CpG-tools: Methylation calling
- **IGV.js**: Interactive genome visualization
- **Quarto/RMarkdown**: Reproducible analysis notebooks

---

## Citation

If you use data or methods from this repository, please cite:

- Lake Trout RNAseq data: [NCBI BioProject PRJNA316738](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA316738)
- Reference genome: GCF_016432855.1 (SaNama_1.0)

---

## Contact

**Roberts Lab**  
School of Aquatic and Fishery Sciences  
University of Washington

For questions or issues, please open a GitHub issue in this repository.