`project-lake-trout/data`

---

- `20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed`: BED file generated from NCBI genome GFF. Documented in this notebook: [Convert NCBI genome GFF to BED](https://robertslab.github.io/sams-notebook/2022/08/18/Data-Wrangling-Convert-S.namaycush-NCBI-GFF-to-genes-only-BED-file-for-Use-in-Ballgown-Analysis.html)

- `ballgown-metadata.csv`: Metadata for individual [_Salvelinus namaycush_](lake trout) samples for use in Ballgown with the following columns:

  - `library`: SRA library name. Naming provides metadata info about each sample:

    - `NP` or `P`: non-parasitized or parasitized

    - `L` or `S`: lean or siscowet

    - `L`: liver

    - `NN`: Sample number
  
  - `subspecies`: Lean or Siscowet
  
  - `infection`: nonparisitized or parasitized
  
  - `strain`: Number classification to reflect subspecies and infection status:

    - `1`:lean parasitized

    - `2`: lean nonparasitized

    - `3`: siscowet parasitized

    - `4`: siscowet nonparasitized

- `fpkm-all_transcripts.csv`: Table of FPKM values for all transcripts.

- `SraRunTable.csv`: Full metadata table from NCBI SRA.

- `whole_tx_table.csv`: Ballgown table of all transcript data for all samples.