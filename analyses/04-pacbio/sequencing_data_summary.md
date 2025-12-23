# PacBio Revio Sequencing Data Summary

## Lake Trout (_Salvelinus namaycush_) Genomics Project

**Report Generated**: 2025-10-18 17:56:57 UTC

## Overview

This report summarizes PacBio Revio HiFi sequencing data generated for comparative genomics analysis of lean and siscowet lake trout subspecies.

## Data Source

- **Repository**: https://owl.fish.washington.edu/nightingales/S_namaycush/LakeTrout/
- **Sequencing Platform**: PacBio Revio
- **Technology**: HiFi (High-Fidelity) Circular Consensus Sequencing (CCS)
- **Species**: _Salvelinus namaycush_ (Lake Trout)
- **Subspecies**: Lean and Siscowet

## Sample Information

### Subspecies

The sequencing data includes samples from two lake trout subspecies:

1. **Lean Lake Trout**
   - Morphotype: Pelagic/limnetic
   - Habitat: Open water
   - Characteristics: Streamlined body, smaller fat content

2. **Siscowet Lake Trout**
   - Morphotype: Benthic/profundal
   - Habitat: Deep water
   - Characteristics: Higher fat content, adapted to deep waters

## Sequencing Technology

### PacBio Revio Platform

The PacBio Revio system is the latest generation of HiFi sequencing technology:

- **Read Type**: HiFi (High-Fidelity) reads
- **Accuracy**: >99.9% (Q30+)
- **Read Length**: Typically 10-25 kb, can exceed 30 kb
- **Chemistry**: Circular Consensus Sequencing (CCS)
- **Applications**: 
  - De novo genome assembly
  - Structural variant detection
  - Full-length isoform sequencing
  - Epigenetic analysis (5mC, 6mA methylation)
  - Haplotype phasing

## Potential Analyses

This dataset enables multiple types of genomic analyses:

### 1. Genome Assembly
- De novo assembly for each subspecies
- Comparative genomics between lean and siscowet
- Identification of subspecies-specific genomic features

### 2. Structural Variation Analysis
- Detection of large insertions/deletions
- Identification of inversions and translocations
- Copy number variation analysis

### 3. Isoform Analysis
- Full-length transcript sequencing
- Alternative splicing patterns
- Gene expression differences between subspecies

### 4. Epigenetic Analysis
- DNA methylation patterns (5mC)
- Comparison of methylation between subspecies
- Gene regulation insights

## Recommended Analysis Workflows

### Alignment
```bash
# Align HiFi reads to reference genome using pbmm2
pbmm2 align --preset CCS --sort \
  reference.fa \
  input.hifi_reads.bam \
  output.aligned.bam
```

### Variant Calling
```bash
# Call variants using pbsv or DeepVariant
pbsv discover aligned.bam variants.svsig.gz
pbsv call reference.fa variants.svsig.gz variants.vcf
```

### Methylation Analysis
```bash
# Extract methylation tags using primrose
primrose aligned.bam output.bam
# Analyze with pb-CpG-tools or custom scripts
```

## References

### Related Data
- **NCBI BioProject**: [PRJNA674328](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA674328)
- **Reference Genome**: GCF_016432855.1 (SaNama_1.0)

### Tools and Documentation
- [PacBio SMRT Tools](https://www.pacb.com/support/software-downloads/)
- [pbmm2 Aligner](https://github.com/PacificBiosciences/pbmm2)
- [pbsv Structural Variant Caller](https://github.com/PacificBiosciences/pbsv)
- [Primrose Methylation Caller](https://github.com/PacificBiosciences/primrose)

---

*This report was generated using `create_sequencing_summary.py` from the project-lake-trout repository.*

*For questions or issues, please contact the RobertsLab team.*
