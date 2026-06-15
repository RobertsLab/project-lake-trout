# 20 Annotation of Ecotype-Specific Genomes — Plan

## Goal

Produce structural (gene models) and functional (product / GO / ortholog) annotation for
the two ecotype-specific PacBio HiFi assemblies so that downstream ecotype comparisons
(differential methylation, PAV, expression) can be reported on each ecotype's own
coordinate system instead of only on the NCBI reference.

## Inputs

| Item | Path | Notes |
|---|---|---|
| Lean assembly | `data/genome/pb-hifiasm-lean-assembly.fa` | hifiasm **primary contigs**, 33,035 contigs, ~3.8 GB on disk |
| Siscowet assembly | `data/genome/pb-hifiasm-siscowet-assembly.fa` | hifiasm **primary contigs**, 24,184 contigs, ~4.1 GB on disk |
| Reference assembly + annotation | NCBI `GCF_016432855.1` (SaNama_1.0) | **Same species** — RefSeq GFF, protein FASTA, `gene_ontology.gaf.gz` |
| Reference gene BED | `data/20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed` | 46,359 genes already in repo |
| Functional table (reference) | `analyses/18-annotation/gene_function_table.tsv` | symbols/products/GO keyed to reference genes — built by `code/18-build-gene-function-table.py` |
| RNAseq | NCBI BioProject `PRJNA316738` (liver, lean/siscowet, ±parasite) | FASTQs not yet in repo; used in Ballgown (`code/01-ballgown-analysis.Rmd`) |

## Key facts that drive the strategy

1. **Same species as the reference.** The reference (GCF_016432855.1) is an annotated
   *S. namaycush* genome. That makes **annotation lift-over (Liftoff)** the highest-value,
   lowest-cost first pass — it transfers the existing 46k gene models (and their symbols,
   already in `gene_function_table.tsv`) directly onto each ecotype assembly with no gene
   prediction needed.
2. **These are contig-level, not chromosome-level.** Both are raw hifiasm primary contigs
   with no scaffolding. Lift-over works on contigs, but assembly QC and possible haplotype
   de-duplication should happen first (the on-disk sizes are well above the ~2.4–2.6 Gb
   expected for a salmonid haploid genome, suggesting retained haplotigs / duplicate content).
3. **De novo prediction adds the ecotype-specific genes Liftoff cannot.** Liftoff can only
   transfer genes that exist in the reference. Genes present in an ecotype but absent from
   SaNama_1.0 (relevant to the PAV story in `code/15-diff-pav.py`) need evidence-based
   prediction (BRAKER3 with RNAseq + protein evidence).
4. **Functional annotation is mostly already solved.** For lifted genes, reuse the existing
   `gene_function_table.tsv` join. Only genuinely *new* de novo gene models need a fresh
   functional pass (eggNOG-mapper / DIAMOND vs SwissProt + InterProScan), reusing the same
   approach as `code/18-build-gene-function-table.py`.

## Recommended approach

A **two-track hybrid**, run identically for each ecotype:

> **Track A (lift-over):** Liftoff transfer of the reference annotation → fast, gives
> directly comparable gene IDs and inherits the existing functional table.
> **Track B (de novo):** BRAKER3 (RNAseq + vertebrate protein evidence) → captures
> ecotype-specific / reference-absent genes. Merged with Track A; only Track-B-unique
> models get new functional annotation.

If compute or time is tight, **Track A alone is a defensible first deliverable** and can
be shipped before Track B is added.

---

## Plan of work

### Step 0 — Freeze inputs and assembly QC

- Treat the two FASTAs in `data/genome/` as frozen references; record md5s.
- Run **assembly stats** (contig count, total length, N50, longest contig) and **BUSCO**
  (`actinopterygii_odb10`, genome mode) on each assembly. Establishes completeness baseline
  and is the denominator for "how complete is the annotation."
- Inspect for retained haplotigs: if BUSCO **duplicated** is high or total length ≫ 2.6 Gb,
  run **purge_dups** and annotate the purged primary set. Keep the decision explicit and
  documented before annotating.
- Output: `analyses/20-ecotype-annotation/qc/{lean,siscowet}/`.

### Step 1 — Repeat identification and soft-masking

- Build a per-assembly repeat library with **RepeatModeler2**, then **RepeatMasker**
  (soft-mask: lowercase, do not hard-mask). Salmonids are highly repetitive, so this
  materially affects gene-prediction quality.
- Optionally combine the de novo library with a teleost/Dfam library.
- Output: soft-masked FASTA per ecotype + repeat summary tables under
  `analyses/20-ecotype-annotation/repeats/{lean,siscowet}/`.
- *Note:* Liftoff (Step 2) does not require masking, but BRAKER3 (Step 3) does — run masking
  before Track B.

### Step 2 — Track A: lift the reference annotation (Liftoff)

- Run **Liftoff** using `GCF_016432855.1` genome + its RefSeq GFF as the source, mapping onto
  each ecotype assembly. Enable `-polish` and `-copies` (capture duplicated/expanded genes).
- Produce per-ecotype GFF3 + a mapping of `ecotype gene ID → reference gene ID`.
- **Reuse the existing functional table:** join lifted genes back to
  `analyses/18-annotation/gene_function_table.tsv` on reference gene ID so lifted models
  inherit symbol / product / GO with no new computation.
- QC: report % of reference genes successfully lifted, single vs multi-copy, and unmapped
  genes (candidates for ecotype-specific loss / PAV).
- Output: `analyses/20-ecotype-annotation/liftoff/{lean,siscowet}/`.

### Step 3 — Track B: de novo evidence-based prediction (BRAKER3)

- Acquire RNAseq evidence: download `PRJNA316738` FASTQs (lean+siscowet liver), trim, and
  **HISAT2/STAR-align to the matching soft-masked assembly**. (Lean RNAseq → lean assembly,
  siscowet RNAseq → siscowet assembly.)
- Add protein evidence: a vertebrate/Actinopterygii protein set (OrthoDB partition) plus the
  SaNama_1.0 RefSeq proteins.
- Run **BRAKER3** (RNAseq BAM + protein DB) on each soft-masked assembly → de novo GFF3.
- Output: `analyses/20-ecotype-annotation/braker/{lean,siscowet}/`.

### Step 4 — Merge tracks and finalize gene set

- Merge Liftoff (Track A) and BRAKER3 (Track B) with **TSEBRA** (or a documented overlap
  rule): prefer the lifted model where they agree, add Track-B models that fall in regions
  with no lifted gene (the ecotype-specific set).
- Tag every final gene with provenance (`liftoff` / `braker` / `merged`).
- QC the final set with **BUSCO in protein mode** and gene-count / mono-exonic-fraction
  sanity checks; compare counts to the ~46k reference and to each other.
- Output: final per-ecotype GFF3 + protein/CDS FASTA in
  `analyses/20-ecotype-annotation/final/{lean,siscowet}/`.

### Step 5 — Functional annotation of new (Track-B-unique) genes

- For genes **not** inherited via Liftoff, run the same functional approach as
  `code/18-build-gene-function-table.py`: **eggNOG-mapper** (or DIAMOND vs SwissProt
  best-hit) for symbol/ortholog/GO, plus **InterProScan** for domains.
- Append to a per-ecotype `gene_function_table.tsv` with the same columns/keying as the
  reference table so all downstream joins are uniform across reference and ecotype genomes.
- Output: `analyses/20-ecotype-annotation/final/{lean,siscowet}/gene_function_table.tsv`.

### Step 6 — Deliverables, browser, and integration

- Per ecotype: soft-masked FASTA, final GFF3, protein/CDS FASTA, functional TSV, BUSCO +
  gene-count QC report, and a genes-only BED matching the
  `data/...genes.bed` convention.
- Add the new GFF3/BED as tracks to the existing genome browsers
  (`genome-browser/`, `jbrowse/`).
- This annotation is the prerequisite for putting the
  [`code/13.3-hifiasm-differential-methylation-plan.md`](13.3-hifiasm-differential-methylation-plan.md)
  DMR results onto ecotype coordinates and for resolving the PAV-implicated, reference-absent
  genes from `code/15-diff-pav.py`.

---

## Tooling summary

| Step | Primary tool | Conda env note |
|---|---|---|
| QC | `assembly-stats`/`seqkit`, BUSCO, purge_dups | follows `/srlab/programs/...` env pattern used in `13.1` |
| Repeats | RepeatModeler2 + RepeatMasker | de novo library per assembly |
| Track A | Liftoff (+ minimap2) | source = GCF_016432855.1 + RefSeq GFF |
| Track B | HISAT2/STAR + BRAKER3 | RNAseq from PRJNA316738 + OrthoDB proteins |
| Merge | TSEBRA / agat | provenance-tagged GFF3 |
| Function | eggNOG-mapper / DIAMOND + InterProScan | reuse `code/18-build-gene-function-table.py` logic |

## Open decisions before running

1. **Scope of first deliverable** — Liftoff-only (fast) vs full hybrid with BRAKER3.
2. **Haplotig handling** — purge_dups yes/no, pending Step 0 BUSCO duplication numbers.
3. **RNAseq reuse** — confirm `PRJNA316738` libraries are the right evidence set and which
   ecotype each maps to (see `data/ballgown-metadata.csv`).
4. **Compute target** — assumed Klone/Hyak (`/srlab/programs`, `/gscratch/scrubbed`), as in
   the `13.1` assembly notebooks; BRAKER3 + RepeatModeler are the heavy steps.
