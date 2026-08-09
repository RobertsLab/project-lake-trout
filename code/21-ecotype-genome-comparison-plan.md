# 21 Ecotype Genome-Level Comparison & Synteny — Plan

## Goal

Now that both ecotype assemblies (`lean`, `siscowet`) carry transferred gene models with
**reference-identical gene IDs** ([`code/20.2-liftoff-annotation.md`](20.2-liftoff-annotation.md)),
compare the two genomes at the whole-genome level. Concretely:

1. Place each ecotype's fragmented contigs onto the chromosome-scale reference so the two
   ecotypes can be compared in a shared chromosome coordinate frame.
2. Build **synteny** between lean and siscowet (and each vs. the reference): collinear blocks,
   rearrangements, inversions, and large structural variants.
3. Quantify **gene-level divergence**: shared vs. ecotype-specific genes, copy-number
   expansions/contractions, and where the PAV / differential-methylation signals fall in
   syntenic context.

## Inputs

| Item | Path / source | Notes |
|---|---|---|
| Lean assembly | `data/genome/pb-hifiasm-lean-assembly.fa` | **33,035 contigs**, 3.73 Gb, largest contig 5.3 Mb — contig-level, size-inflated |
| Siscowet assembly | `data/genome/pb-hifiasm-siscowet-assembly.fa` | **24,184 contigs**, 4.0 Gb, largest contig 4.8 Mb — contig-level, size-inflated |
| Chromosome-scale reference | `data/genome/ref-GCF_016432855.1/GCF_016432855.1_SaNama_1.0_genomic.fna` | **42 chromosomes** (`NC_*` >10 Mb), 2.35 Gb total — the anchor |
| Lean lifted genes | `output/20.1-liftoff-annotation/lean.liftoff.gff3_polished` / `.genes.bed` | 79,535 gene models, reference `gene-XXX` IDs |
| Siscowet lifted genes | `output/20.1-liftoff-annotation/siscowet.liftoff.gff3_polished` / `.genes.bed` | 82,691 gene models, reference `gene-XXX` IDs |
| Lifted proteins | `output/20.1-liftoff-annotation/{eco}.liftoff.proteins.faa` | for orthology cross-checks if needed |
| Unmapped reference features | `output/20.1-liftoff-annotation/{eco}.unmapped_features.txt` | lean 1,294 / siscowet 1,073 — loss / PAV candidates |
| Reference gene BED | `data/20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed` | 46,359 reference genes |
| PAV results | `code/15-diff-pav.py` outputs | to overlay on synteny |

## Key facts that drive the strategy

1. **The reference is chromosome-scale; the ecotype assemblies are not.** SaNama_1.0 has
   42 chromosome sequences. The ecotype hifiasm assemblies are tens of thousands of contigs
   with no scaffolding. **Direct lean-vs-siscowet whole-genome alignment on raw contigs would
   be uninterpretable** (no chromosome axis, contig-vs-contig hairball). The reference is the
   natural anchor: align/assign each ecotype's contigs to reference chromosomes first, then
   compare ecotypes *through* that shared frame.

2. **Orthology is already solved — for free.** Liftoff assigned every lifted gene the reference
   `gene-XXX` ID. A gene shared by lean and siscowet has the *same ID in both*. This eliminates
   the usual hardest synteny step (ortholog inference with OrthoFinder/Diamond). Gene-anchored
   synteny (MCScanX/GENESPACE-style) can use ID identity directly as the anchor pairs.

3. **Assembly size is inflated (~1.5× expected).** Salmonid haploid genome is ~2.4–2.6 Gb;
   both ecotype assemblies are 3.7–4.0 Gb → **retained haplotigs / duplicate content**. Raw
   assemblies will show spurious "duplications" in any alignment. **purge_dups before synteny**
   is strongly recommended (this is Step 0 of [`code/20-ecotype-genomes.md`](20-ecotype-genomes.md),
   not yet executed). At minimum, results must be interpreted knowing duplicate content is present.

4. **Two complementary resolutions.** Whole-genome **alignment** (nucmer/minimap2) catches
   sequence-level SVs and aligns intergenic/repeat regions; **gene-anchored synteny** catches
   collinearity and is robust to fragmentation because it follows gene order, not raw sequence.
   Run both; they cross-validate.

## Recommended approach

> **Anchor on the reference, compare ecotypes through it.** Layer two methods:
> **(A) whole-genome alignment** of each ecotype → reference (contig→chromosome assignment,
> dotplots, SV calls), and **(B) gene-anchored synteny** using the shared Liftoff gene IDs
> (collinear blocks lean↔reference↔siscowet, expansions, ecotype-specific genes).
> Then overlay PAV and DMR signals onto the syntenic map.

If time is tight, **(B) gene-anchored synteny is the highest-value first deliverable** — it is
cheap (operates on BED/GFF, not 4 Gb FASTAs), exploits the free orthology, and directly answers
"what is shared vs. ecotype-specific and where." (A) can follow.

---

## Plan of work

### Step 0 — (Recommended) purge haplotigs

- Run **purge_dups** per ecotype (uses the existing PacBio reads + self-alignment) to collapse
  retained haplotigs toward the expected ~2.5 Gb haploid size.
- Re-lift genes onto the purged assemblies *or* document that comparison uses raw assemblies with
  duplicate content acknowledged. **Decision point** — see open decisions.
- Output: `analyses/21-genome-comparison/purge/{lean,siscowet}/`.

### Step 1 — Contig → chromosome assignment (reference anchoring)

- Align each ecotype assembly to the reference with **minimap2 `asm5`** (or **nucmer**), then
  use **RagTag `scaffold`** (or a coverage/identity rule) to assign and order contigs onto the
  42 reference chromosomes.
- Produce a per-ecotype AGP / contig→chromosome map and chromosome-anchored coordinates.
- This gives both ecotypes a shared 42-chromosome axis for everything downstream.
- Output: `analyses/21-genome-comparison/anchor/{lean,siscowet}/`.

### Step 2 — Track A: whole-genome alignment & SVs

- **nucmer/`dnadiff`** (MUMmer; `nucmer` is available at `/usr/bin/nucmer`) lean↔reference,
  siscowet↔reference, and lean↔siscowet on the anchored sequences.
- Generate **dotplots** (mummerplot or **D-GENIES**) and `dnadiff` reports (alignment identity,
  coverage, inversions, translocations, indels).
- Call structural variants with **SyRI** (needs whole-genome alignment input) → syntenic regions,
  inversions, translocations, duplications, and SNPs/indels between ecotypes.
- Output: alignment deltas, dotplots, `*.syri.out` under `analyses/21-genome-comparison/wga/`.

### Step 3 — Track B: gene-anchored synteny

- Use the shared Liftoff gene IDs as anchors. Two equivalent routes:
  - **GENESPACE** (orthology + synteny + riparian plots) — feed it the three gene sets; or
  - **MCScanX / `jcvi` (MCscan)** — build anchor pairs directly from shared `gene-XXX` IDs
    (no BLAST needed for shared genes) and detect collinear blocks.
- Produce **riparian / synteny ribbon plots**: reference ↔ lean and reference ↔ siscowet,
  and lean ↔ siscowet, ordered by the 42 chromosomes.
- Identify: collinear blocks, breakpoints, lineage-specific inversions, and **gene
  expansions/contractions** (genes lifted in >1 copy via Liftoff `-copies`).
- Output: anchor tables, block tables, synteny plots under `analyses/21-genome-comparison/synteny/`.

### Step 4 — Shared vs. ecotype-specific gene set

- From the two `.genes.bed` files, compute the **gene-ID set intersection/difference**:
  shared genes, lean-only, siscowet-only (cross-reference each ecotype's
  `unmapped_features.txt` — a gene unmapped in one ecotype but present in the other is a strong
  presence/absence candidate).
- Build a **copy-number table** per gene per ecotype (from Liftoff `-copies` output) to flag
  CNV differences.
- Cross-check against [`code/15-diff-pav.py`](15-diff-pav.py): do PAV calls land on
  ecotype-specific / CNV-divergent genes here?
- Output: `analyses/21-genome-comparison/gene-sets/{shared,lean_only,siscowet_only,cnv}.tsv`.

### Step 5 — Overlay functional & methylation signal

- Map the **DMRs** from
  [`code/13.3-hifiasm-differential-methylation-plan.md`](13.3-hifiasm-differential-methylation-plan.md)
  onto the syntenic map: are differentially methylated regions in collinear blocks, near
  breakpoints, or in ecotype-specific regions?
- GO-enrich the ecotype-specific / CNV-divergent / DMR-overlapping gene sets reusing
  [`code/18.3-go-enrichment.py`](18.3-go-enrichment.py) and the inherited functional tables
  (`output/20.1-liftoff-annotation/{eco}.liftoff.gene_function_table.tsv`).
- Output: integrated table + enrichment results under `analyses/21-genome-comparison/integrate/`.

### Step 6 — Deliverables, browser, write-up

- Per comparison: dotplots, SyRI SV table, synteny ribbon plots, shared/specific gene tables,
  CNV table, and the integrated PAV/DMR overlay.
- Add SV / synteny-block tracks to the existing browsers (`genome-browser/`, `jbrowse/`).
- Lab-notebook post (`code/21.x-...md`) summarizing the comparison, mirroring the
  [`20.2`](20.2-liftoff-annotation.md) format.

---

## Tooling summary

| Step | Primary tool | Availability note |
|---|---|---|
| Purge | purge_dups | needs PacBio reads (already used in `13.1`) |
| Anchor | minimap2 (`asm5`) + RagTag | minimap2 in `liftoff_env`; RagTag via conda |
| WGA / SV | **nucmer/dnadiff** (MUMmer), D-GENIES, **SyRI** | `nucmer`/`mummer` at `/usr/bin`; SyRI via conda |
| Gene synteny | **GENESPACE** or MCScanX / `jcvi` | shared `gene-XXX` IDs = anchors, no ortholog step |
| Gene-set ops | `bedtools`, pandas | shared-ID set logic |
| Function / GO | reuse `code/18.3-go-enrichment.py` | inherited functional tables from `20.1` |
| Plots | D-GENIES dotplots, GENESPACE riparian, ggplot | — |

## Open decisions before running

1. **Purge haplotigs first?** Cleanest synteny needs it, but it means re-lifting genes onto
   purged assemblies (re-run [`20.1`](20.1-liftoff-annotation.Rmd)). Alternative: proceed on raw
   assemblies and flag duplicate content. **Recommendation: purge, then re-lift** — otherwise
   every "duplication" call is suspect.
2. **Anchoring method** — RagTag reference-scaffolding vs. lighter contig→chromosome label only.
   Scaffolding gives chromosome coordinates but introduces gap-filled joins; labeling is more
   conservative.
3. **Synteny engine** — GENESPACE (turnkey, nice plots, R-heavy) vs. MCScanX/`jcvi` (more manual,
   more control). Given the free shared IDs, either works; GENESPACE is the faster path to figures.
4. **SV caller** — SyRI requires near-chromosome-level, 1-to-1 input; confirm Step 1 anchoring is
   clean enough, else restrict SV calling to the best-anchored chromosomes.
5. **Compute target** — heavy steps are purge_dups and whole-genome alignment of 4 Gb assemblies;
   assume the same Klone/Hyak pattern as `13.1`.
```
