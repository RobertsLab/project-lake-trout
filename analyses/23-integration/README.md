# Phase 2 integration — two-genome ecotype-divergence evidence

Gene-level integration of Phase-2 (native reciprocal PAV, ecotype-only gene presence, CNV) with Phase-1 SVs and the prior reference DMR/PAV tables. See `code/23.6-integrate-phase2.py`.

> **Superseded for interpretation by `refined/`** (code/23.7-refine-evidence-model.py): the ecotype_only / ref_pav / sv lines below are not independent, ecotype_only is Liftoff-sensitivity-inflated, and the GO tests are tandem-cluster driven. This file is kept as the Step-3 record.

## Evidence lines (genes per line)

- **reciprocal_pav**: 115
- **ecotype_only**: 14,647
- **cnv**: 2,359
- **sv**: 13,519
- **dmr**: 181
- **ref_pav**: 1,263

## Two-genome candidates (>=2 lines incl. a Phase-2 line): 2,895

A candidate is no longer a single-reference association: at least one of reciprocal_pav / ecotype_only / cnv (native, bidirectional) plus >=1 other line.

Top 20 by n_lines:

| gene | symbol | product | n_lines | lines |
|---|---|---|---|---|
| gene-LOC120017348 | LOC120017348 | A disintegrin and metalloproteinase with | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120017775 | LOC120017775 | ABC transporter G family member 23-like | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120017932 | LOC120017932 | histone H1 | 3 | cnv,sv,dmr |
| gene-LOC120017972 | LOC120017972 | histone H3 | 3 | cnv,sv,dmr |
| gene-LOC120017974 | LOC120017974 | galectin-related protein-like | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120018007 | LOC120018007 | histone H2A | 3 | cnv,sv,dmr |
| gene-LOC120018068 | LOC120018068 | histone H4 | 3 | cnv,sv,dmr |
| gene-LOC120018277 | LOC120018277 | uncharacterized LOC120018277 | 3 | cnv,sv,ref_pav |
| gene-LOC120020067 | LOC120020067 | corticotropin-releasing factor receptor  | 3 | cnv,sv,ref_pav |
| gene-LOC120022255 | LOC120022255 | NACHT, LRR and PYD domains-containing pr | 3 | cnv,sv,ref_pav |
| gene-LOC120022969 | LOC120022969 | vam6/Vps39-like protein | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120022978 | LOC120022978 |  | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120023423 | LOC120023423 | zinc finger protein 239-like | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120023487 | LOC120023487 | slit homolog 1 protein-like | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120024764 | LOC120024764 | ADAMTS-like protein 1 | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120024769 | LOC120024769 | pro-interleukin-16-like | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120024894 | LOC120024894 | cysteine--tRNA ligase, cytoplasmic-like | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120024897 | LOC120024897 | cGMP-inhibited 3',5'-cyclic phosphodiest | 3 | reciprocal_pav,ecotype_only,sv |
| gene-LOC120024918 | LOC120024918 | transmembrane channel-like protein 3 | 3 | ecotype_only,sv,ref_pav |
| gene-LOC120024927 | LOC120024927 | stAR-related lipid transfer protein 5-li | 3 | ecotype_only,sv,ref_pav |

## Phenotype survival

`phenotype_survival.tsv` lists lipid-metabolism and calcium-transport GO terms that remain enriched (FDR<0.25) in the two-genome sets. Total flagged: 13.
