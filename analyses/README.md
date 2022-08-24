`project-lake-trout/analyses`

---

- `DEG-subspecies-pval_filtered_p-0.05.csv`: Differentially expressed genes (DEG) identified between subspecies (lean vs. siscowet), filtered for DEG having _only_ p-values <= 0.05. Contains the following columns:

  - `t_id`: Transcript ID (integer) assigned by `ballgown`.

  - `chr`: Chromosome ID from NCBI [_Salvelinus namaycush_ (lake trout)](https://en.wikipedia.org/wiki/Lake_trout) genome `GCF_016432855.1_SaNama_1.0_genomic.fna`.

  - `strand`: Originating DNA strand (+ or -).

  - `start`: One-based starting location of transcript on chromosome.

  - `end`: One-based ending location of transcript on chromosome.

  - `t_name`: Annotated transcript name taken from genome file.

  - `num_exons`: Number of exons comprising transcript.

  - `length`: Length of transcript (bp).

  - `gene_id`: Corresponding gene ID, if applicable.

  - `gene_name`: Corresponding gene name, if applicable.

  - `cov.*`: Sequencing read coverage for transcript for each sample.

  - `FPKM.*`: Relative expression levels for each transcript in each sample.

- `DEG_subspecies-pval_qval_filtered_p-q-0.05.csv`: Differentially expressed genes (DEG) identified between subspecies (lean vs. siscowet), filtered for DEG having p-values <= 0.05 _and_ q-values <= 0.05. NOTE: There are none!

- `DEG-subspecies-siscowet-pval_filtered_p-0.05.bed`: BED file of differentially expressed transcripts (DET) up-regulated in siscowet only, filtered for DET having _only_ p-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- `DET-subspecies-lean-pval_filtered_p-0.05.bed`: BED file of differentially expressed transcripts (DET) up-regulated in lean only, filtered for DET having _only_ p-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- `DET-subspecies-pval_filtered_p-0.05.bed`: BED file of differentially expressed transcripts (DET) up-regulated in lean only, filtered for DET having _only_ p-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.

- `DET-subspecies-pval_filtered_p-0.05.csv`: Differentially expressed transcripts (DET) identified between subspecies (lean vs. siscowet), filtered for DET having _only_ p-values <= 0.05. Contains the following columns:

  - `t_id`: Transcript ID (integer) assigned by `ballgown`.

  - `chr`: Chromosome ID from NCBI [_Salvelinus namaycush_ (lake trout)](https://en.wikipedia.org/wiki/Lake_trout) genome `GCF_016432855.1_SaNama_1.0_genomic.fna`.

  - `strand`: Originating DNA strand (+ or -).

  - `start`: One-based starting location of transcript on chromosome.

  - `end`: One-based ending location of transcript on chromosome.

  - `t_name`: Annotated transcript name taken from genome file.

  - `num_exons`: Number of exons comprising transcript.

  - `length`: Length of transcript (bp).

  - `gene_id`: Corresponding gene ID, if applicable.

  - `gene_name`: Corresponding gene name, if applicable.

  - `cov.*`: Sequencing read coverage for transcript for each sample.

  - `FPKM.*`: Relative expression levels for each transcript in each sample.

- `DET-subspecies-pval_qval_filtered_p-q-0.05.csv`: Differentially expressed transcripts (DET) identified between subspecies (lean vs. siscowet), filtered for DEG having p-values <= 0.05 _and_ q-values <= 0.05. NOTE: There are none!

- `DET-subspecies-siscowet-pval_filtered_p-0.05.bed`: BED file of differentially expressed transcripts (DET) up-regulated in siscowet only, filtered for DET having _only_ p-values <= 0.05. Includes optional
columns 4 (`name`), 5 (`score`), and 6 (`strand`). The `name` column is the transcript ID assigned by `ballgown`. The `score` column has been assigned an arbitrary value of 0.