# 18-annotation — Step 0: gene function table

Functional annotation backbone for the SaNama_1.0 (GCF_016432855.1) gene set, used to
attach function to differential-methylation and PAV gene hits
(see `code/18-diff-annotation-phenotype-plan.md`).

## Deliverable

`gene_function_table.tsv` — one row per gene/pseudogene, keyed on `gene_id`
(the same `gene-XXX` IDs as `data/20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed`).

| column | meaning |
|---|---|
| `gene_id` | join key (GFF `ID`, e.g. `gene-rbm24b`) |
| `symbol` | gene symbol (or `LOC######` if unnamed) |
| `geneid` | NCBI GeneID |
| `biotype` | protein_coding, lncRNA, tRNA, pseudogene, … |
| `chrom,start,end,strand` | genomic location |
| `is_named` | 1 if a real symbol, 0 if `LOC######` |
| `product` | RefSeq transcript product description |
| `n_go,go_bp,go_mf,go_cc` | GO term counts (total / BP / MF / CC) |
| `go_terms` | pipe-delimited GO IDs |

## Provenance

Built by `code/18-build-gene-function-table.py` from NCBI RefSeq annotation of
GCF_016432855.1 (downloaded to `raw/`):

- `*_genomic.gff.gz` — gene models + transcript `product=` names (join backbone)
- `*_gene_ontology.gaf.gz` — GO annotations by GeneID
- `*_feature_table.txt.gz` — retained for reference; not required by the build

Source: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/016/432/855/GCF_016432855.1_SaNama_1.0/

## QC (build run)

- 49,668 gene + pseudogene records; **all 46,359 BED gene IDs covered (0 missing)**.
- The 3,309 extra rows are pseudogenes (biotype=`pseudogene`), absent from the genes-only BED.
- 9,738 genes with a real symbol; 46,231 with a product description; 34,367 with ≥1 GO term.

## Steps 1-2 — feature-to-gene assignments

Built by `code/18.1-assign-features-to-genes.py` (pure stdlib; exon coords from the GFF).
Windows: **promoter = TSS ±2 kb** (strand-aware), **flanking = gene body ±5 kb**. Primary
`location_class` priority: promoter > exon > intron > upstream≤5kb > downstream≤5kb. Every
row carries symbol / product / biotype / GO count via the Step 0 join.

| file | what | rows |
|---|---|---|
| `dmr_gene_assignments.tsv` | DMR↔gene pairs | 293 pairs; **149/302 DMRs** within 5 kb of a gene (26 gene-body overlaps, matches the diff-meth report's 25) |
| `dmc_gene_assignments.tsv` | DMC↔gene pairs | 3,953 pairs; 2,186/4,440 DMCs near a gene |
| `pav_gene_assignments.tsv` | **stringent siscowet-specific deletion**↔gene pairs | 1,756 pairs; 1,543/3,465 deletions within 5 kb of a gene; **54 genes with an exon-overlapping deletion** (candidate copy-number/LOF in siscowet) |
| `pav_gene_burden.tsv` | per-gene lenient PAV burden (count + overlap bp) | 22,450 genes with ≥1 lenient PAV |

Distribution of DMR placements: 88 promoter, 137 upstream≤5kb, 52 downstream≤5kb, 11 intron,
5 exon. Stringent-PAV placements: 1,073 intron, 215 promoter, 187 up / 249 down, 32 exon.

**Reference-bias reminder (carried in `pav_gene_burden.tsv`):** the reference is lean, so the
siscowet-specific and deletion columns are inflated by divergence — see the caveats in
`code/18-diff-annotation-phenotype-plan.md`. Interpret PAV burden comparatively, not absolutely.

## Step 3 — integrated ranked candidates

Built by `code/18.2-integrate-candidates.py` → `integrated_candidate_genes.tsv`.
One row per gene with ≥1 DMR / DMC / stringent-PAV assignment within 5 kb (**2,036 genes**),
joined to lenient PAV burden and to liver RNAseq expression. `rank_score` rewards
convergence, promoter/exon placement, expression, and methylation↔expression concordance;
the ranked view deprioritizes repetitive ncRNA (`caution=repetitive_ncRNA`) so protein-coding
candidates surface first.

- **4 convergent genes** (a DMR *and* a stringent siscowet-specific deletion): `znf883-like`
  (LOC120032414 — exon DMR + exonic deletion, the top candidate), `XlCGF57.1-like`
  (LOC120040411), `septin-9-like` (LOC120043843), and one uncharacterized locus.
- Among protein-coding promoter-DMR genes: epoxide hydrolase 1-like, OCIA-domain 1-like,
  znf665-like, `rbm24b`, small G-protein signaling modulator 3-like.
- Exonic stringent siscowet deletions hit several **lipid-metabolism** genes — angiopoietin-
  related protein 5-like (`angptl5`), 2-acylglycerol O-acyltransferase 2-A-like (`mogat2`) —
  plus immune/adhesion genes (alpha-2-macroglobulin, CEACAM, DMBT1). These are the axes that
  distinguish the lean (lean-bodied) and siscowet (high-lipid) ecotypes; flagged for Step 5.

**Expression linkage is weak and that is expected:** `whole_tx_table.csv` is liver RNAseq from
a *separate parasite study* (12 lean + 12 siscowet, different individuals than the PacBio
fish), so most candidates are not liver-expressed (`log2fc=NA`) and only a handful yield a
concordance call. Treat expression as orthogonal support only, never confirmation.

## Step 4 — GO enrichment (over-representation)

Built by `code/18.3-go-enrichment.py` → `go_enrichment_{dmr,pav,union}.tsv`. Self-contained
hypergeometric ORA on the **local** GO annotations (so LOC genes are retained), with full
DAG ancestor propagation from `go-basic.obo`, BH-FDR within each set, and a `phenotype` flag
for ecotype-relevant terms (lipid / buoyancy / growth / ion transport / hypoxia / immune).
Background = 34,268 GO-annotated genes. Term-size filter 5–2000, ≥3 study genes.

**Read these results conservatively — three structural confounders dominate:**

1. **DMR set is a tandem-cluster artifact.** Only 30 of 181 DMR-genes carry GO; the 29
   "FDR<0.1" terms come almost entirely from **one histone gene cluster** (6 co-located
   histones → chromatin/nucleosome terms) and **three adjacent `znf883` paralogs** →
   transcription-repressor terms. This is ~2 loci, not broad functional convergence.
2. **PAV ORA has a gene-length bias.** Long genes accumulate deletions by chance; calcium/ion
   channels are among the longest genes in any genome, so the calcium-transport signal is
   partly length-driven, not necessarily selection. (A length-matched permutation null is the
   proper fix — deferred.)
3. **Reference bias.** PAV study set is siscowet-specific deletions on a *lean* reference, so
   it is divergence-inflated (see plan caveats).

**What survives as the most defensible signal (PAV set):**

| term | FDR | fold | genes (k) |
|---|---|---|---|
| calcium ion transport | **3.0e-3** | 3.0 | 19 |
| (calcium ion transmembrane transport / channel activity / channel complex) | <0.1 | — | — |
| P-type ion transporter activity | 0.10 | 4.1 | 5 |
| developmental / cell growth | 0.15 | 4.1 | 4 |
| lipid binding / phospholipid binding | 0.31 / 0.35 (ns) | 1.5 | 24 / 15 |

The top non-flagged PAV terms are **neuron projection development / axon guidance** and
**calcium channel** terms. **Lipid and growth terms are suggestive, not significant**, and
broad `lipid metabolic process` / `fatty acid metabolic process` are *depleted* (fold<1) —
so the lipid-phenotype thread from Step 3 is at the gene level, not a genome-wide GO signal.
The union set tracks the PAV set (PAV genes are 945/972 of it).

## Step 5 — phenotype synthesis

`code/18-diff-annotation-phenotype.Rmd` → render-ready HTML report (base R; needs pandoc, e.g.
RStudio) + `figures/18-top-candidates.png`. Ranked candidates + GO results turned into per-axis
**hypothesized** phenotype links (lipid/energy storage, body-shape/growth/muscle,
calcium/sensory-neural, immune, and a flagged chromatin tandem-cluster artifact), anchored to
the 17-landmark geometric-morphometric data in `data/measurements.xlsx`, with the full
limitation set (reference bias, no FDR-significant DMC, expression mismatch, gene-length bias,
associative-not-causal). Verified by knitting to markdown — all chunks execute.

## Notes / next steps

- GO terms here are RefSeq/InterPro-derived (mostly IEA). Step 4 enrichment can either use
  these directly or map symbols to zebrafish orthologs for a denser GO/KEGG background.
- ~36.6k genes remain `LOC######` (no symbol) but most carry an informative `product`
  (e.g. "septin-9-like"). Per the plan, only genes lacking both symbol and product need the
  optional eggNOG/DIAMOND ortholog pass; that is not required for Steps 1–3.
