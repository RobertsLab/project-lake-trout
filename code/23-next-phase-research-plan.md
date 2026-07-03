# 23 — Next-Phase Research Plan: Two-Genome Ecotype Comparison

## Where we are

Both ecotype PacBio HiFi assemblies have been **purged of retained haplotigs** and now sit at
biologically sensible salmonid haploid sizes:

| Assembly | Contigs | Size | Lifted genes |
|---|---|---|---|
| `lean.purged.fa` | 12,700 | 2.58 Gb | 55,437 |
| `siscowet.purged.fa` | 8,159 | 2.84 Gb | 59,968 |

Done to date: Liftoff annotation on both purged assemblies with reference-identical gene IDs
([`20.1`](20.1-liftoff-annotation.Rmd)); gene-anchored MCScanX synteny — reference↔lean (1,897
blocks), reference↔siscowet (1,966), lean↔siscowet (4,113) — with dotplots and contig→chromosome
maps ([`21.1`](21.1-gene-anchored-synteny.Rmd)); a JBrowse 2 linear synteny view
([`22`](22-synteny-browser-integration-plan.md)).

**The central limitation this phase removes.** Every biological result so far — DMRs, PAVs,
candidate genes, GO enrichment ([`14`](14-diff-meth.py), [`15`](15-diff-pav.py),
[`18`](18-diff-annotation-phenotype-plan.md)) — was computed against the **lean-background
doubled-haploid reference** (GCF_016432855.1). That is the standing caveat in the project README:
siscowet-specific signal is divergence-inflated and not magnitude-comparable to lean. We now have
each ecotype's own genome, so the next phase re-derives the key comparisons on **native
coordinates** and cross-checks them **bidirectionally**, converting single-reference associations
into two-genome evidence.

## What has NOT been run (the gaps this plan fills)

1. **Whole-genome alignment & structural-variant calling** between the ecotypes — Track A of
   [`21`](21-ecotype-genome-comparison-plan.md) (nucmer/SyRI). Only gene-anchored synteny (Track B)
   was done. No inversions, translocations, or large-indel calls exist yet.
2. **Reciprocal PAV** — PAV was called only against the reference. Calling each ecotype against the
   *other's* genome removes the lean-reference bias directly.
3. **Ecotype-native differential methylation** — plan [`13.3`](13.3-hifiasm-differential-methylation-plan.md)
   was written but never executed. Methylation is still reference-anchored.
4. **De novo annotation (BRAKER3)** — Track B of [`20`](20-ecotype-genomes.md). Liftoff only
   transfers genes that exist in the reference; genuinely ecotype-specific genes (the core of the
   PAV story) are still invisible.

---

## Prioritized plan of work

Ordered by scientific value per unit effort. Phases 1–2 are the payoff of building two genomes and
should go first; 3 deepens methylation; 4 is the larger compute investment; 5 integrates and writes up.

### Phase 1 — Whole-genome SV map between ecotypes  *(highest value, lowest cost)*

Complete Track A of plan 21 on the purged assemblies.

- **Anchor** each purged assembly to the 42 reference chromosomes: `minimap2 -x asm5` →
  RagTag `scaffold`, producing a shared chromosome axis. (Contig→chromosome seed already exists in
  `output/21.1-gene-anchored-synteny/contig_to_chromosome.tsv`.)
- **Align & call SVs**: `nucmer`/`dnadiff` (MUMmer, at `/usr/bin/nucmer`) for lean↔ref,
  siscowet↔ref, and lean↔siscowet; then **SyRI** on the anchored alignments → syntenic regions,
  inversions, translocations, duplications, and SNPs/indels.
- **Cross-validate** SyRI SVs against the existing PAV calls ([`15`](15-diff-pav.py)): reference-based
  PAVs that fall inside SyRI-confirmed ecotype SVs are the high-confidence set.
- Output: `analyses/23-genome-sv/{anchor,wga,syri}/`, dotplots, and a merged SV table.

### Phase 2 — Reciprocal PAV + shared/ecotype-specific gene set

Turn single-reference PAV into a bidirectional, gene-aware call set.

- **Reciprocal PAV**: run the [`15`](15-diff-pav.py) coverage/CIGAR approach both directions
  (lean reads → siscowet assembly, and vice-versa). A region absent in one direction *and* present
  in the other is a confident presence/absence event, no longer inflated by lean-reference distance.
- **Gene-set ops** (Step 4 of plan 21, still open): from the two `*.purged.liftoff.genes.bed` files
  compute shared / lean-only / siscowet-only gene-ID sets, cross-referencing each
  `*.purged.unmapped_features.txt`. Build a per-gene **copy-number table** from Liftoff `-copies`.
- **Integrate**: intersect reciprocal PAV and CNV-divergent genes with Phase-1 SVs and the existing
  DMR/candidate tables; re-run GO enrichment ([`18.3`](18.3-go-enrichment.py)) on the two-genome
  gene sets. Test whether the lipid-metabolism / calcium-transport signals survive on native coords.
- Output: `analyses/23-reciprocal-pav/`, `analyses/23-gene-sets/{shared,lean_only,siscowet_only,cnv}.tsv`.

### Phase 3 — Ecotype-native differential methylation

Execute plan [`13.3`](13.3-hifiasm-differential-methylation-plan.md) end to end.

- Re-align each ecotype's PacBio HiFi reads to its **own** purged assembly with
  [`04-pacbio/align_hifi_pbmm2.py`](04-pacbio/align_hifi_pbmm2.py); call 5mC with the
  `modkit pileup` workflow from [`10-mCG-call.Rmd`](10-mCG-call.Rmd).
- **Coordinate harmonization** (the key constraint in 13.3): use the Liftoff gene-ID bridge in
  `*.gene_function_table.tsv` to project per-CpG methylation from each native assembly onto shared
  reference/gene coordinates before the cross-ecotype test, replacing the current `chrom:pos` merge
  in [`14-diff-meth.py`](14-diff-meth.py).
- Re-run the DMR test on harmonized coordinates and compare against the reference-based DMRs — which
  of the 302 original DMRs are confirmed vs. reference artifacts.
- Output: `analyses/23-native-methylation/`.

### Phase 4 — De novo annotation of ecotype-specific genes (BRAKER3)

Track B of plan [`20`](20-ecotype-genomes.md) — the larger compute item.

- RepeatModeler2 + RepeatMasker soft-masking per assembly (salmonids are highly repetitive).
- Acquire `PRJNA316738` RNAseq, align lean→lean / siscowet→siscowet; add an Actinopterygii OrthoDB
  protein set + SaNama RefSeq proteins; run **BRAKER3** per ecotype.
- Merge with Liftoff via **TSEBRA**; keep only BRAKER-unique models (the ecotype-specific set),
  tag provenance, and functionally annotate the new models with eggNOG-mapper/InterProScan
  (same approach as [`18-build-gene-function-table.py`](18-build-gene-function-table.py)).
- QC with BUSCO (`actinopterygii_odb10`, protein mode).
- Output: `analyses/23-denovo-annotation/final/{lean,siscowet}/`.

### Phase 5 — Integration, browser, manuscript

- Add SV, reciprocal-PAV, and native-DMR tracks to the browsers (`genome-browser/`, `jbrowse/`),
  extending the synteny view from [`22`](22-synteny-browser-integration-plan.md).
- Assemble the two-genome ecotype-divergence narrative (lipid/buoyancy and depth-adaptation
  hypotheses), stating for each candidate whether it is now supported by ≥2 independent lines
  (native DMR, reciprocal PAV, SV, expression) rather than single-reference association.
- Draft the assembly + comparative-genomics manuscript; the assemblies also warrant NCBI submission.

---

## Recommended sequencing & decision points

- **Do Phase 1 first.** It is cheap (MUMmer + SyRI on ~2.5 Gb assemblies), and its SV map is the
  scaffold every other phase overlays onto.
- **Phases 1–3 are the scientific core** and directly retire the single-reference caveat. If effort
  must be capped, ship these three and defer Phase 4.
- **Phase 4 (BRAKER3) is the compute-heavy fork.** Only launch it once RNAseq is staged and a GPU/
  large-CPU node is available; it can proceed in parallel with Phases 2–3.
- **Open question — a chromosome-scale ecotype assembly?** The purged assemblies are still
  contig-level (thousands of contigs). If Hi-C or a linkage map is obtainable, scaffolding at least
  one ecotype to chromosomes would make it a standalone reference and sharpen every comparison here.
  Absent that, reference-anchoring (Phase 1) is the pragmatic substitute — decide before Phase 4.
