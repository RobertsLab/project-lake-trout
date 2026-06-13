# 18 Differential Methylation & PAV — Gene Annotation and Phenotype Plan

## Goal

For the lean vs siscowet comparison on the **GCF_016432855.1 (SaNama_1.0)** reference,
determine **which genes** the differential methylation (DMRs) and differential PAVs fall
in/near, attach **functional annotation** to those genes, and reason about **plausible
phenotypic consequences** for the two ecotypes. The output is an interpretable, ranked
gene table plus a short biological narrative — not new variant calling.

## Inputs already in the repo

| Layer | File(s) | Notes |
|---|---|---|
| DMRs | `analyses/14-diff-meth/dmrs.bed`, `dmrs_annotated.csv`, `dmrs_genes_intersect.bed` | 302 DMRs; 25 already intersect genes. Direction in `dmrs.csv` (`hyper`/`hypo_siscowet`). |
| DMCs | `analyses/14-diff-meth/significant_dmcs.bed` | 4,440 single-CpG hits — wider net than DMRs. |
| Diff PAV (lenient) | `analyses/15-diff-pav/{lean,siscowet}_specific.{insertions,deletions}.bed` | Large, noisy. |
| Diff PAV (stringent) | `analyses/15-diff-pav/stringent.siscowet_specific.deletions.bed` | **3,465 deletions >100 bp, all-4-vs-none — the high-confidence set to lead with.** |
| Gene models | `data/20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed` | 46,359 genes, coords + name only (9,738 real symbols, rest `LOC######`). |
| RNAseq | `data/whole_tx_table.csv`, `fpkm-all_transcripts.csv` | Per-sample FPKM with `gene_name` — orthogonal expression support. |
| Phenotype | `data/measurements.xlsx` | Direct ecotype phenotype measures to anchor interpretation. |

## Key gap to close first

The repo has **no functional annotation** — only gene coordinates and symbols. Inferring
phenotype requires mapping each implicated gene to a product/description, GO/KEGG terms,
and (ideally) a vertebrate ortholog. This is the one new data dependency. **Step 0 below
acquires it; everything downstream reuses bedtools intersects the repo already does.**

---

## Plan of work

### Step 0 — Build a functional annotation table for SaNama_1.0 genes  ✅ DONE

> Implemented in `code/18-build-gene-function-table.py` →
> `analyses/18-annotation/gene_function_table.tsv` (see `analyses/18-annotation/README.md`).
> All 46,359 BED genes covered; 9,738 symbols, 46,231 products, 34,367 with GO. GO came
> straight from NCBI's `gene_ontology.gaf.gz`, so the eggNOG/DIAMOND ortholog pass below was
> not needed for first-pass coverage (revisit only for symbol-less + product-less LOC genes).


- Pull the NCBI annotation for GCF_016432855.1 (the `*_genomic.gff.gz` and the
  feature/`gene_ontology` or protein tables from the RefSeq `Annotation_comparison` /
  `*_feature_table.txt.gz`). Goal columns: `gene_symbol`, `LOC_id`, `product/description`,
  `GO_terms`, and a vertebrate ortholog (zebrafish/human) where available.
- For `LOC######` genes with no symbol, assign orthology by the protein FASTA against a
  reference proteome (zebrafish SwissProt) — eggNOG-mapper or a DIAMOND/BLAST + best-hit
  pass is sufficient. Many symbols are already zebrafish-style, so coverage should be good.
- Deliverable: `analyses/18-annotation/gene_function_table.tsv` keyed on the same
  `gene-XXX` IDs used in `data/...genes.bed`. This is the join key for all later steps.

### Step 1 — Assign DMRs (and DMCs) to genes with positional context  ✅ DONE

> `code/18.1-assign-features-to-genes.py` → `analyses/18-annotation/dmr_gene_assignments.tsv`
> + `dmc_gene_assignments.tsv`. 149/302 DMRs within 5 kb of a gene (26 gene-body, matches the
> report's 25); 88 DMRs land in promoters. Promoter=TSS±2kb, flanking=gene±5kb.


- Reuse the existing `bedtools intersect` pattern. Beyond exact overlap, add **proximity
  classes** because regulatory methylation often sits in promoters: intersect DMRs against
  gene bodies, and against promoter windows (TSS ±1–2 kb, strand-aware) and gene ±5 kb.
- Classify each DMR–gene pair as `promoter / gene_body / intergenic-near / distal`, and
  carry DMR `direction` (hyper/hypo in siscowet) and `mean_diff`.
- Run the same against `significant_dmcs.bed` to catch single-CpG signals near genes the
  302-DMR set misses.
- Deliverable: `analyses/18-annotation/dmr_gene_assignments.tsv`.

### Step 2 — Assign differential PAVs to genes, by overlap type  ✅ DONE

> `code/18.1-assign-features-to-genes.py` → `analyses/18-annotation/pav_gene_assignments.tsv`
> (stringent siscowet-specific deletions: 1,543/3,465 near a gene, **54 with exon overlap**)
> + `pav_gene_burden.tsv` (lenient lean/siscowet del+ins per-gene burden, 22,450 genes;
> reference-bias caveat carried in-file).


- **Lead with `stringent.siscowet_specific.deletions.bed`** (cleanest signal). Intersect
  with gene bodies and exons; classify as `exon-overlapping` (potential loss-of-function /
  copy change) vs `intronic` vs `flanking`.
- Repeat for the lenient lean- and siscowet-specific sets as a lower-confidence tier, and
  report a per-gene PAV burden (count + total bp of ecotype-specific deletion).
- Insertions: the stringent insertion set is empty, so treat the lenient insertion sets as
  exploratory only.
- Deliverable: `analyses/18-annotation/pav_gene_assignments.tsv` with a confidence tier.

### Step 3 — Integrate the two layers and add expression evidence  ✅ DONE

> `code/18.2-integrate-candidates.py` → `analyses/18-annotation/integrated_candidate_genes.tsv`
> (2,036 candidate genes; **4 convergent** DMR+stringent-PAV genes). Liver RNAseq joined on
> `gene_id`; expression linkage weak (different individuals/tissue) so used as support only.
> Ranked view deprioritizes repetitive ncRNA via a `caution` flag. Lipid-metabolism genes
> (`angptl5`, `mogat2`, epoxide hydrolase 1) flagged for Step 5.


- Join DMR-genes and PAV-genes on the Step 0 table to produce one master gene list with
  flags: `has_DMR`, `DMR_direction`, `in_promoter`, `has_stringent_PAV`, `PAV_exonic`.
- **Highlight convergent genes** hit by *both* methylation and PAV — strongest candidates.
- Cross-reference each candidate against the RNAseq FPKM table (`whole_tx_table.csv`) by
  `gene_name`: is the gene expressed in liver, and does lean-vs-siscowet expression move in
  the direction methylation/PAV would predict? Note the RNAseq is liver/parasite-focused,
  so treat as supporting, not confirmatory.
- Deliverable: `analyses/18-annotation/integrated_candidate_genes.tsv` (ranked).

### Step 4 — Functional enrichment  ✅ DONE (GO; KEGG deferred)

> `code/18.3-go-enrichment.py` → `analyses/18-annotation/go_enrichment_{dmr,pav,union}.tsv`.
> Self-contained hypergeometric ORA on local GO with DAG propagation (go-basic.obo), BH-FDR.
> **Robust signal: calcium ion transport/channels (FDR 3e-3) in the PAV set.** DMR enrichment
> is a histone + znf883 tandem-cluster artifact; lipid/growth terms only suggestive (FDR>0.1).
> Confounders flagged in README: gene-length bias (PAV ORA), reference bias, IEA annotations,
> gene clustering. KEGG deferred (needs KO/ortholog mapping; GO is the substantive result).


- Run GO/KEGG over-representation on the DMR-gene set, the PAV-gene set, and the union,
  each against the full annotated gene set as background (gprofiler2 / clusterProfiler with
  the zebrafish ortholog mapping from Step 0).
- Report enriched terms with FDR; flag any tied to **lipid/fatty-acid metabolism,
  buoyancy/swim-bladder, lipid storage, growth, depth/pressure or hypoxia response** —
  the axes that distinguish lean (shallow, lean-bodied) from siscowet (deep-water, high
  lipid content) ecotypes.
- Deliverable: enrichment tables + a couple of summary plots in `analyses/18-annotation/`.

### Step 5 — Phenotype interpretation  ✅ DONE

> `code/18-diff-annotation-phenotype.Rmd` (base-R, render-ready; renders where pandoc is
> available, e.g. RStudio) + figure `figures/18-top-candidates.png`. Per-axis evidence blocks
> (lipid/energy, body-shape/growth/muscle, calcium/sensory-neural, immune, chromatin-artifact),
> anchored to the morphometric `data/measurements.xlsx`, with the full causality/reference
> caveat set. Verified by knitting to markdown (all chunks execute, tables + figure populate).


- For the top convergent/enriched candidates, write a short evidence block each: gene,
  function, methylation/PAV/expression evidence, and the **specific lake-trout phenotype**
  it plausibly affects (lipid deposition, buoyancy, depth tolerance, growth rate, etc.).
- Anchor claims to `data/measurements.xlsx` where a measured ecotype trait lines up.
- Be explicit about causality limits: associations on a single reference, no functional
  validation; promoter-hypermethylation→silencing is a hypothesis, not a result.
- Deliverable: `code/18-diff-annotation-phenotype.Rmd` (or `.qmd`) rendering the ranked
  table + narrative, mirroring the `14-diff-meth.Rmd` style, and a top-candidates figure
  for the genome browser.

---

## Caveats / scope guards

- **The reference is a LEAN-background genome.** GCF_016432855.1 (SaNama_1.0) was assembled
  from a female gynogenetic **doubled haploid** derived from the **Seneca Lake (NY) brood
  stock** (Smith et al. 2022, Mol Ecol Resour; Iron River NFH). The paper gives no explicit
  ecotype, but the **Seneca Lake strain is a lean morphotype** (shallow <80 m, late-fall
  spawner) — so the reference is effectively lean, and from an Atlantic-drainage population
  distinct from Great Lakes fish, meaning even the lean samples are only an approximate match.
- **Directional PAV bias (lead with this caveat).** Siscowet diverges more from a lean
  reference, so siscowet reads map less completely and siscowet shows systematically *more*
  apparent variants — an artifact, not biology (consistent with the 1.33M vs 0.99M and the
  3,465 vs 0 stringent-deletion asymmetry already observed). Therefore:
    - Siscowet-specific vs lean-specific **counts are not magnitude-comparable**.
    - A "siscowet-specific deletion" is enriched for **lean-present/siscowet-absent** sequence
      and partly reflects the lean reference rather than true siscowet loss.
    - Lean-specific **deletions are suppressed** (can't delete what the reference lacks);
      **lean-specific *insertions*** (lean sequence absent from the reference) are the most
      trustworthy lean signal. Lead PAV interpretation with siscowet-specific deletions and
      lean-specific insertions — the two least-biased classes.
- **Methylation power bias.** Divergence-driven mapping loss costs siscowet covered CpGs in
  the most divergent regions, biasing DMR detection toward reference-concordant (lean-like)
  loci — conservative against exactly the regions where ecotypes differ most. Flag high-
  divergence regions as low-power.
- **Invisible siscowet gene content.** Genes absent/highly divergent in the lean reference are
  not in the gene BED, so a reference-only intersect cannot see siscowet-specific gene
  content. The per-ecotype hifiasm track (`code/13.3-...plan.md`) is the needed complement;
  keep the two analyses in distinct output folders.
- **Doubled-haploid reference**: a single fully homozygous haplotype stands in for the
  species — no allelic variation, not a pangenome. Affects the annotation backbone, not the
  per-sample methylation values (called on the bc20xx samples).
- **PAV direction semantics**: confirm the calling convention in `code/15.5-diff-pav.py`
  before asserting which ecotype "lost" a gene, and interpret it against the lean reference.
- **FDR reality**: 0 DMCs survive q<0.1 in the current report. Lead interpretation with the
  DMR-level and stringent-PAV sets; treat single-CpG and lenient-PAV hits as hypothesis-
  generating, and say so.
- **Annotation coverage**: phenotype reasoning is only as good as Step 0; `LOC` genes with
  no ortholog stay as "unknown function" rather than being forced into a story.

## Suggested implementation order

1. Step 0 (annotation table) — unblocks everything.
2. Steps 1 & 2 (intersects) — fast, reuse existing bedtools patterns.
3. Step 3 (integrate + expression).
4. Steps 4 & 5 (enrichment + write-up).
