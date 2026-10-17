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

---

## Progress log

### Phase 1 — DONE (see [`23.1`](23.1-genome-sv-map.Rmd), [`23.2`](23.2-sv-pav-crossvalidate.py))
SV map + PAV cross-validation for all three comparisons. A repaired step is noted below.

### Phase 2 — DONE (2026-07)

**Step 0 — repaired the Phase-1 siscowet SV gap** ([`23.3`](23.3-sisco-sv-repair.py)). `sisco_vs_ref`
SyRI had aborted ("Unequal number of chromosomes"): RagTag had scaffolded a 15 kb siscowet contig
onto the reference **mitochondrion** `NC_036392.1`, giving 43 anchored scaffolds vs the reference's
42 nuclear chromosomes. Fixed by trimming the mito from `siscowet.anchored.fa` and dropping its
record from the reused WGA BAM, then re-running SyRI. Regenerated `merged_sv_table.tsv` /
`sv_summary.tsv` (now symmetric across all three comparisons) and re-ran [`23.2`](23.2-sv-pav-crossvalidate.py):
siscowet PAV is now SV-confirmed (insertions 17.6 %, deletions 18.8 %; cf. lean 18.3 % / 20.4 %).
The [`23.1`](23.1-genome-sv-map.Rmd) anchor chunk was patched (size filter on `NC_*_RagTag`
scaffolds) so the mito cannot re-enter.

**Step 1 — reciprocal PAV** ([`23.4`](23.4-reciprocal-pav.py), `analyses/23-reciprocal-pav/`). Each
ecotype's 4 HiFi samples aligned (pbmm2, CCS) to the *other* ecotype's purged assembly (8
alignments, ~100 % mapped). `11-pav`-style mosdepth + CIGAR calling, then 4-sample consensus on
native coordinates. Interior consensus-absent regions (present in the target ecotype's own genome,
absent from all 4 querying-ecotype read sets) give ecotype-specific presence without lean-reference
inflation:

| Direction | target-specific present | genes | querying-specific novel seq |
|---|---|---|---|
| lean reads → siscowet | 340 regions / 0.70 Mb (siscowet) | 71 | 115,171 lean clusters / 10.3 Mb |
| siscowet reads → lean | 343 regions / 0.87 Mb (lean) | 44 | 175,748 siscowet clusters / 17.2 Mb |

**Step 2 — gene sets + CNV** ([`23.5`](23.5-gene-set-ops.py), `analyses/23-gene-sets/`).
50,379 shared, 5,058 lean-only, 9,589 siscowet-only genes; 2,359 CNV-divergent (copy count dominated
by tRNA/repeat families). Ecotype-"only" calls partly corroborated by the other assembly's
`unmapped_features` (1,474/5,058 lean-only; 3,040/9,589 siscowet-only — capped by the unmapped list
sizes).

**Step 3 — integration + GO re-run** ([`23.6`](23.6-integrate-phase2.py), `analyses/23-integration/`).
Gene-level evidence matrix over 6 lines (reciprocal_pav 115, ecotype_only 14,647, cnv 2,359, sv
13,519, dmr 181, ref_pav 1,263). **2,895 two-genome candidates** carry ≥2 independent lines including
a Phase-2 line (109 with ≥3). GO ORA re-run (reusing [`18.3`](18.3-go-enrichment.py) machinery):
the **calcium-transport / depth-adaptation signal survives strongly on native two-genome
coordinates** — lean-only genes are enriched for voltage-gated calcium channel activity (FDR 9e-4)
and calcium ion transport (FDR 9e-3); reciprocal-PAV genes for calcium ion binding (FDR 0.09).
**The lipid-metabolism / buoyancy signal did NOT reach FDR<0.25 in any two-genome set**, i.e. it was
more reference-dependent than the calcium signal and should be treated cautiously pending Phase 3
(native methylation) and Phase 4 (de novo genes). See `analyses/23-integration/phenotype_survival.tsv`.

**Step 4 — evidence-model refinement (2026-09)** ([`23.7`](23.7-refine-evidence-model.py),
`analyses/23-integration/refined/`). Review of the Step-3 matrix found three problems: (i) the
`ecotype_only`, `ref_pav` and `sv` lines all measure distance from the lean reference and co-occur
heavily (1,565 genes share ecotype_only+ref_pav), so `n_lines` over-stated independence; (ii)
`ecotype_only` (14,647) is dominated by Liftoff sensitivity between two contig-level assemblies —
only ~30 % of each "only" set is corroborated by the other assembly's `unmapped_features`; (iii) the
GO signals were tandem-cluster driven (18 hemoglobin genes at one NC_052347.1 locus lead lean_only;
histone / tRNA / snRNA / protocadherin arrays dominate cnv). 23.7 collapses the three reference-
divergence lines into one, splits ecotype_only into corroborated (kept as its own line) and
Liftoff-only (demoted), and collapses same-family genes within 50 kb to one representative in study
and background before the GO test.

| line (refined) | genes |
|---|---|
| reciprocal_pav | 114 |
| ecotype_only_corroborated | 4,382 (lean 1,434 / siscowet 2,948; unique reference IDs after collapsing Liftoff `_N` copy suffixes — cf. 1,474 / 3,040 raw rows in 23.5) |
| cnv | 2,359 |
| dmr | 181 |
| ref_divergence (sv ∪ ref_pav ∪ Liftoff-only) | 16,068 |

Tiers: **A** (reciprocal_pav + ≥1) = 77 genes (66 flag-free; 48 siscowet-present, 29 lean-present);
**B** (corroborated-only or cnv + ≥1) = 4,011; **C** (dmr + ref_divergence only) = 15. The 23.6
figure of 2,895 "two-genome candidates" is superseded by tier A+B = 4,088, of which only tier A has
read-level native support.

Outcome for the hypotheses: the **calcium signal survives collapsing** (lean_only_corroborated:
voltage-gated calcium channel activity fold 4.5, FDR 0.028 collapsed vs 8e-4 uncollapsed; calcium
ion transport FDR 0.062; reciprocal_pav calcium ion binding FDR 0.10). **No lipid term reaches
FDR 0.1** in any set (best: long-chain fatty acid metabolic process, FDR 0.24, 4 clusters). Tier-A
GO is small-n but points at insulin response / muscle contraction (gpia, mitochondrial PEPCK, two
PDEs, grb10-like, ROR-alpha paralogs partitioned one per ecotype, calsequestrin-1, troponin T).
Phase-18 headline fate (`phase18_candidate_fate.tsv`): 3 of 4 convergent genes are tier C
(reference-only), LOC120040411 is tier B via cnv inside a 5-paralog cluster; `angptl5` is tier B
(corroborated siscowet-only, no read support), `mogat2` and epoxide hydrolase 1 have no native line.

**Step 5 — read-coverage filter on corroborated ecotype-only genes (2026-10-14)** (23.7 updated).
Checking `angptl5` showed lean reads cover it at genome depth in all four lean fish when aligned to
the siscowet assembly; the whole ~280 kb block (LOC120049129–LOC120050016, ref
NC_052312.1:40.30–40.54 Mb) is missing from the lean *assembly*, not lean fish. Two
assemblies agreeing on absence is therefore not enough. Each corroborated gene is now checked
against the 23.4 mosdepth windows of the other ecotype's reads on the assembly that carries it;
genes with ≥80 % of their length at ≥25 % of genome-mean depth in ≥3 of 4 fish are `read_covered`
(assembly gaps) and lose the line. The coverage distribution is bimodal; the cutoff sits in the gap.

Result: **3,891 of 4,382 (89 %) corroborated genes are read-covered**; 491 remain (lean 147 /
siscowet 344). Tier A unchanged (77 / 66 flag-free); tier B 4,011 → 2,509. All 18 hemoglobin genes
and all 9 genes behind the lean_only voltage-gated calcium channel enrichment are read-covered, so
**the calcium-channel signal was an assembly-gap effect** and is withdrawn. Remaining calcium
evidence: reciprocal_pav calcium ion binding (FDR 0.10, 8 genes) and individual tier-A genes.
Lipid unchanged (no term at FDR 0.1). `angptl5` drops to reference divergence only. Caveat:
coverage counts MAPQ-0 reads, so recent paralogs can look covered; BAM-level MAPQ filtering on the
Linux box would tighten this.
The hemoglobin cluster absent from the siscowet assembly is not flagged by reciprocal PAV, consistent
with an assembly gap. Landing page and README updated to this model.

### Phase 3 — pass 1 (2026-07, run on raven; committed there 2026-07-09 but not pushed until 2026-10)

Ecotype-native differential methylation, gene-anchored ([`23.8`](23.8-native-meth-pipeline.py),
[`23.9`](23.9-native-meth-dmg.py), `analyses/23-native-methylation/`).

**Step A — native calling** ([`23.8`](23.8-native-meth-pipeline.py)). Each ecotype's 4 HiFi samples
aligned (pbmm2 CCS, MM/ML preserved, reusing the Phase-2 purged indexes) to its OWN purged assembly;
5mC called with `modkit pileup --cpg --combine-strands`. Tooling note: modkit was not installed, and
the current **v0.6.4 rejects these reads** (no `MN` tag + dual 6mA/5mC mods); dropped to **v0.5.0**
(the version 10-mCG-call used) which processes them cleanly. Per-sample CpGs at cov≥10 ranged widely
(bc2069 0.68 M → bc2096 27 M), reflecting uneven per-sample depth.

**Step B — gene-anchored DMG** ([`23.9`](23.9-native-meth-dmg.py)). Coordinate harmonization by the
shared Liftoff gene ID (per approved plan): each CpG assigned to a gene on its native assembly,
coverage-weighted mean methylation per gene per sample (`bedtools map`, ≥3 CpGs/gene), then
Mann-Whitney lean vs siscowet per gene (14-diff-meth recipe: 15% diff, BH-FDR).

Result: **10,838 genes testable** (limited by native coverage in ≥3 samples/ecotype);
**63 nominal DMGs** (p<0.05, |Δ|≥15%; 37 hyper / 26 hypo in siscowet) but **0 at FDR<0.10** — with
only 4 vs 4 samples the Mann-Whitney floor is p≈0.029, so nothing survives multiple testing. Of the
181 reference-DMR genes, only **45 were testable natively** and **0 were confirmed** as native DMGs
(chance overlap ≈0.3), i.e. the reference-based methylation differences do **not** reproduce on
native two-genome coordinates.

*Superseded by pass 2 below — the "0 confirmed" result was a power artifact.*

**Interpretation at the time (with caveats).** This corroborates the Phase-2 theme that single-reference
signals are fragile — the reference DMRs show no native confirmation. But the result is strongly
**power- and coverage-limited**: n=4/group caps sensitivity, gene-level aggregation is coarser than
per-CpG DMRs, and only ~1/5 of shared genes (and 45/181 reference-DMR genes) had enough native
coverage to test. So this is best read as "reference DMRs are not corroborated natively at current
depth," not "definitively artifacts." A deeper or CpG-level native test (or more samples) is the way
to settle it. Outputs: `native_methylation_by_gene.tsv`, `dmg_results.tsv`, `dmg_significant.tsv`,
`reference_dmr_confirmation.tsv`, `native_dmg_report.txt`.

### Phase 3 — pass 2 (2026-10-14): count-based test, reference-DMR placement, frame control

Pass 1 kept only CpGs with ≥10 reads per fish (native depth is ~8–20×, so bc2069 kept 0.68 M CpGs)
and used Mann-Whitney on 4 vs 4 (smallest possible p = 0.029). Pass 2 reuses the 23.8 alignments
and bedMethyl files; all compute ran on raven (Klone was not needed — data and indexes were already
there, and no realignment was required).

**Step C — count-based features** ([`23.11`](23.11-native-meth-counts.py)). Every 5mC CpG call
(Nvalid ≥ 1) pooled into promoter (2 kb upstream of TSS) and gene-body features of the 36,668 genes
single-copy in both Liftoff annotations, paired by gene ID. Quasi-binomial test (closed-form
two-group GLM, Pearson dispersion, t on n−2 df), BH per feature class. 34,740 bodies / 34,442
promoters testable (vs 10,838 genes in pass 1). 1,031 promoters and 245 bodies at FDR < 0.10 and
|Δ| ≥ 10 points. Label permutation (all 34 other balanced 4/4 splits, including splits by
sequencing run) gives 0–1 hits, so the test is calibrated and run batch is not the driver.

**Step C′ — reference DMRs on native coordinates** (23.11 Part B). Each of the 302 reference DMRs
± 500 bp was placed on both purged assemblies with minimap2 `asm5`. **About half are MAPQ 0 — present
in ≥ 2 near-identical copies — in the native assemblies (lean 47 %, siscowet 54 %), against 3–5 %
of 1,000 random reference windows; only 9 % reach MAPQ 60 vs 82–86 % of random windows**
(`ref_dmr_placement_mapq.tsv`). The reference DMRs are concentrated in duplicated / repetitive
sequence where reference read placement is ambiguous. The 23 placed uniquely on both assemblies
mostly reproduce (sign agrees 19/22, Spearman ρ = 0.68).

**Step D — frame control** ([`23.12`](23.12-native-meth-frame-control.py)). In Step C ecotype is
confounded with assembly. modkit was run on the 23.4 cross alignments (MM/ML preserved;
`meth_cross/run_modkit_cross.sh`), so each assembly holds both ecotypes' reads on identical
coordinates, and every feature was re-tested inside the lean frame and inside the siscowet frame.

| | lean frame | siscowet frame |
|---|---|---|
| features tested | 69,635 | 69,720 |
| p < 0.05 | 5.2 % | 5.2 % |
| π0 estimate | 0.96 | 0.97 |
| FDR < 0.10 | 0 (min q 0.91) | 0 (min q 0.91) |
| reference DMRs testable / sign agrees / ρ | 48 / 47 / 0.71 | 38 / 34 / 0.61 |
| reference DMRs same sign, p < 0.05, \|Δ\| ≥ 10 | 32 | 19 |

**Result.** (1) **The Step-C feature hits are almost all coordinate-harmonization artifacts**: e.g.
the `tmeff2a` promoter reads ~65 % methylated in all eight fish on the lean assembly and ~10 % in all
eight on the siscowet assembly — the two Liftoff "promoters" are different sequence, not different
methylation. Correlation of Step-C differences with either frame is ρ ≈ 0.1. Only **63 features (57
promoters, 6 bodies) keep sign, p < 0.05 and |Δ| ≥ 10 in both frames** (`dmf_robust.tsv`); these are
nominal candidates, none in tier A/B, 19 carrying only ref_divergence. They include individual genes
of interest — a DGAT2-like promoter (+21 points in siscowet), hippocalcin-like 1, pro-opiomelanocortin
A, `sik2b`, carbonic anhydrase 4-like — but GO over-representation gives only generic terms
(FDR 0.09) and **no lipid or calcium term**. (2) **Genome-wide, promoter and gene-body methylation
does not differ detectably between ecotypes** (null p-value distribution in both frames). (3) **The
reference DMRs that can be placed uniquely are real**: they replicate inside both frames, including
the siscowet frame where no read touches the lean reference. Pass 1's "0 confirmed" was a power
artifact. The other ~half of reference DMRs lie in multi-copy sequence and cannot be adjudicated with
these assemblies.

Methodological lesson: pairing features across two assemblies by Liftoff gene ID is not a valid
basis for a quantitative cross-ecotype comparison; a shared frame (both ecotypes' reads on one
assembly), tested in both directions, is. Tables: `dmf_frame_{lean,siscowet}.tsv.gz`,
`dmf_robust.tsv`, `ref_dmr_frame_tests.tsv`, `ref_dmr_native_loci.tsv`, reports
`native_meth_counts_report.txt` and `native_meth_frame_report.txt`.

**Step E — CpG-level DMR scan inside each frame (2026-10-15)**
([`23.13`](23.13-native-meth-cpg-scan.py), `analyses/23-native-methylation/cpg_scan/`). The 14-diff-meth
design (DMCs chained within 250 bp, same direction) repeated inside each frame on CpGs with ≥ 5
reads in ≥ 7 of 8 fish (label-free filter; 10.05 M CpGs lean frame, 9.41 M siscowet frame), per-CpG
quasi-binomial test, and **empirical FDR from the 18 fully balanced relabelings** (2 lean + 2
siscowet per group) over a threshold grid (`dmr_scan_grid.tsv`):

| setting (DMC p, \|Δ\|, min DMCs) | lean frame DMRs (eFDR) | siscowet frame DMRs (eFDR) |
|---|---|---|
| 0.05, 15, 3 — the 14-diff-meth recipe | 2,533 (0.43) | 2,535 (0.39) |
| **0.01, 15, 3 — primary** | **268 (0.07)** | **304 (0.05)** |
| 0.001, 25, 3 | 12 (0.00) | 15 (0.00) |

At the 14-diff-meth thresholds ~40 % of DMRs are expected false, which helps explain why the 302
reference DMRs were a mix of real and spurious loci.

**Result.** Primary DMRs are short (median ~200 bp, 5 DMCs) with large effects (median |Δ| 30
points), and mostly **hypomethylated in siscowet** (206/268 and 248/304; at those DMRs lean fish
average 59–65 % and siscowet fish 27–31 %). About 43 % place uniquely on the other assembly; of
those, **101/111 (lean) and 119/131 (siscowet) replicate in the other frame** (same sign, p < 0.05,
|Δ| ≥ 10; Spearman ρ 0.97 / 0.93). 101 loci were called independently in both frames; ~137 unique
replicated loci, 113 genes. Robustness: lean and siscowet fish are completely separated at 265/268
and 297/304 DMRs; dropping any one fish never flips a DMR or takes it below 10 points; 94 % keep
their sign within each of the three sequencing runs. 63 % are intergenic, ~20 % promoter, ~14 %
body; 10 per frame overlap a placed reference DMR. GO on replicated-DMR genes: transcription-factor
/ cis-regulatory DNA binding (FDR 3e-4), but 12 of the 16 genes are dispersed C2H2 zinc-finger
genes (a repeat-rich family not merged by 50 kb collapsing); the others are en1a, a forkhead E1-like,
a HES-4-like and an ATF-like 3. No lipid term; calcium ion binding only at FDR 0.23 (3 genes). One
tier-A gene (LOC120041427) and 14 tier-B genes (cnv + ref_divergence, mostly uncharacterized) carry
a replicated DMR; 4 genes overlap the 23.12 robust features (AFG3-like 1, nudt19, a NET1-like,
LOC120034056).

**Phase 3 conclusion.** Ecotype methylation differences are **localized, large and reproducible**
on native coordinates — a few hundred short regions, mostly siscowet hypomethylation, often
intergenic — while promoter- and gene-wide methylation does not differ. The fixed, all-fish
consistency of these differences is compatible with cis-genetic control as well as with a stable
epigenetic ecotype difference; n = 4 wild-caught fish per ecotype cannot separate the two.

**Step F — native DMRs as evidence line and browser tracks (2026-10-15)**
([`23.14`](23.14-native-dmr-tracks.py); [`23.7`](23.7-refine-evidence-model.py) rerun). Each 23.13
DMR gets a status: replicated (220), single_frame (330, not placeable / testable on the other
assembly), not_replicated (22). Genes within 5 kb on the DMR's own assembly are listed with a
strand-aware location class (`cpg_scan/native_dmr_gene_assignments.tsv`). **`native_dmr` (71 genes
near a replicated DMR) replaces the reference `dmr` line** in 23.7; the reference set is kept as the
`ref_dmr` column. Tier B now also admits native_dmr; tier C (reference DMR + reference divergence)
is retired. Tiers: A 77 → 76 (a 28S rRNA lost its reference-DMR line), B 2,509 → 2,527 (+22 via
native_dmr, −4 that relied on the reference DMR), former tier C (15) → single-line. 30 tier-B genes carry a
native DMR; none in tier A. native_dmr GO is the same C2H2 zinc-finger repressor signal as 23.13;
no lipid or calcium term reaches significance.

Tracks (BED9, itemRgb; blue = lower in siscowet, red = higher; dark = replicated, light = single
frame): all 572 DMRs on the lean and siscowet purged assemblies in JBrowse; 243 DMRs placed uniquely
on SaNama_1.0 (152 of 220 replicated), merged into 160 reference loci (85 replicated) for IGV and
JBrowse. Reference DMR tracks were renamed "Reference DMRs".

**Step G — sequence differences near native DMRs (2026-10-15)**
([`23.15`](23.15-dmr-fixed-variants.py), `analyses/23-native-methylation/fixed_variants/`). Inside each
frame, all eight fish were genotyped (bcftools 1.14 `-X pacbio-ccs`, MAPQ ≥ 20, QUAL ≥ 20) over every
DMR ± 5 kb and over 10 matched control windows per DMR (same length, tested-CpG count within 30 %,
≥ 20 kb from any DMR). Callable = ≥ 4 reads in all 8 fish. Window-level tests (fraction of windows
hit, matched resampling), so one divergent block cannot dominate.

| lean frame / siscowet frame | DMR windows | matched controls | p |
|---|---|---|---|
| fixed difference (all 4 vs all 4) within 1 kb | 0 / 0 % | 0.1 / 0.2 % | — |
| fixed difference 1–5 kb | 1.5 / 1.6 % | 0.3 / 0.3 % | 0.008 / 0.0005 |
| SNV \|ΔAF\| ≥ 0.5 within 1 kb | 14.2 / 13.2 % | 7.2 / 7.6 % | 0.0001 / 0.0003 |
| same, replicated DMRs | 22.8 / 18.5 % | 6.3 / 6.1 % | 0.0001 / 0.0001 |
| ecotype depth ratio ≥ 1.4-fold | 26.9 / 33.2 % | 14.0 / 15.5 % | 0.0001 / 0.0001 |

**Result.** (1) **No DMR has a fixed sequence difference within 1 kb**; fixed differences are rare
everywhere in these eight fish (62–75 per frame in ~16 Mb callable, ~1 per 300–500 kb). The 1–5 kb
excess comes mostly from one interferon cluster (interferon a3-like / alpha-6-like; LDMR_162/163 =
SDMR_114/115), where a block of 12–18 fixed differences flanks replicated hypomethylated DMRs in
both frames. (2) **Differentiated SNVs are 2–3× enriched within 1 kb of DMRs**, most strongly at
replicated DMRs and consistently in both frames, so a substantial minority of DMRs (~15–20 %) sit in
genetically differentiated haplotypes. (3) **About a third of DMRs show an ecotype read-depth
difference** (vs ~15 % of controls), in both directions and not tied to methylation direction: copy-
number variation of collapsed paralogs (salmonid WGD) is a second, partly separate contributor. DMRs
with differentiated SNVs are not more depth-deviant (21–25 % vs 28–35 %), so (2) is not just (3).
QC: 8–12 % of covered variant sites are heterozygous-looking in all 8 fish (collapsed paralogs).

**Interpretation.** Most native DMRs have neither a nearby fixed difference nor a differentiated
SNV, so a simple "fixed cis variant" model does not explain them; but genetic differentiation and
copy-number differences are each enriched at DMRs. A differentiated SNV within 5 kb or a depth
difference marks 47–52 % of DMR windows vs 28–30 % of matched controls, an excess of ~20 points. Genetic
and epigenetic ecotype divergence are therefore entangled, and n = 4 wild fish per ecotype cannot
resolve causality. Allele-specific methylation in heterozygous fish (HiFi reads carry both) is the
direct next test.

**Step H — allele-specific methylation at differentiated DMRs (2026-10-16)**
([`23.16`](23.16-dmr-allele-specific-meth.py), `analyses/23-native-methylation/asm/`). In fish
heterozygous for a differentiated SNV (|ΔAF| ≥ 0.5, within 5 kb), HiFi reads spanning the SNV and the
DMR core were split by allele and 5mC (MM/ML, ML ≥ 128) pooled per allele at the core CpGs.
ASM = siscowet-allele % − lean-allele %, same fish. Cis predicts ASM ≈ DMR difference (slope 1);
trans / environment predicts ASM ≈ 0. Controls: the 23.15 matched windows with a differentiated SNV.

| | lean frame | siscowet frame |
|---|---|---|
| DMRs testable (phased heterozygotes) | 55 of 268 | 66 of 304 |
| \|ASM\| median, DMR vs control | 34 vs 6 points | 20 vs 5 points |
| slope ASM ~ DMR difference (95 % CI) | 0.79 (0.60–0.98) | 0.72 (0.53–0.91) |
| sign agreement | 47/55 | 51/66 |
| excluding ecotype depth-different windows | 0.84 (0.63–1.05), 38/43 | 0.78 (0.57–1.01), 37/48 |
| replicated DMRs only | 0.81 (0.55–1.05), 27/30 | 0.82 (0.60–1.05), 29/33 |
| slope in lean / siscowet heterozygotes | 0.82 / 0.65 | 0.90 / 0.55 |
| cis-like / intermediate / trans-like | 32 / 19 / 4 | 27 / 32 / 7 |

**Result.** At the DMRs that can be phased, **methylation follows the haplotype**: the two alleles
inside one fish differ by most of the ecotype difference (slope ~0.8), in lean and in siscowet
heterozygotes, robust to excluding copy-number-suspect windows and to restricting to replicated DMRs.
Cis-like DMRs include `aldh1a3`, talin-2-like, FMR1-like, transmembrane protein 271-like, mitochondrial
Rho GTPase 1-A, an ATF-like bZIP, several C2H2 zinc-finger genes and interferon a3-like. The slope is
lower in siscowet heterozygotes (0.55–0.67) than in lean (0.82–0.94), so part of the difference does
not follow the allele in siscowet fish (a trans or environmental component, or incomplete phasing).
Scope: this covers the ~20 % of DMRs with a differentiated SNV and a heterozygous fish; the rest cannot
be phased with these data, and the two frames share many loci.

**Phase 3 conclusion (updated).** Where it can be tested, ecotype methylation divergence at the native
DMRs is largely **cis-genetic** — carried by differentiated haplotypes rather than set by ecotype
environment — with a minority trans / environmental component, larger in siscowet.

**Next:** Phase 4 (BRAKER3) on raven; for methylation, a population sample (more fish, both
ecotypes) would test the cis haplotypes and the untestable ~80 % of DMRs.

### Phase 4 — de novo annotation (BRAKER3) — RUNNING (started 2026-10-16 on raven)

**Setup** ([`23.17`](23.17-denovo-setup.sh), `analyses/23-denovo-annotation/`): containers
`teambraker/braker3:v3.1.1` (BRAKER 3.0.8, GeneMark-ETP, AUGUSTUS, ProtHint, TSEBRA),
`dfam/tetools:2.00` (RepeatModeler 2.0.9, RepeatMasker 4.2.4) and BUSCO 6.1.0
(`actinopterygii_odb12`); RNA-seq PRJNA316738 (72 liver runs, 36 bp single-end GAIIx) pooled per
ecotype (lean 257 M reads, siscowet 237 M); proteins = OrthoDB v12 Vertebrata (19.39 M) + SaNama_1.0
RefSeq (58 k).

**Pipeline** ([`23.18`](23.18-denovo-pipeline.sh)), one resumable chain per ecotype, both in parallel at
24 threads (niced): RepeatModeler2 `-LTRStruct` → RepeatMasker soft-mask with the ecotype's own
library → HISAT2 own-ecotype RNA-seq → BRAKER3 (RNA + proteins) → BUSCO protein mode.

**Smoke test** (one ~5 Mb gene-dense contig per assembly, 12 threads): every stage ran; lean contig
44 % soft-masked, 149 BRAKER genes vs 106 Liftoff genes. Timing: BRAKER 83 min, of which DIAMOND
vs the 19 M proteins was 3 min and ProtHint / spaln ~40 min (scales with gene number). Extrapolated
full run: RepeatModeler ~3–4 days, BRAKER ~6–8 days per ecotype, ~1.5–2 weeks wall time.
Caveat: the only RNA-seq is 36 bp single-end liver, so splice support is weak and non-liver genes
rely on protein evidence.

**Merge script ready** ([`23.19`](23.19-merge-annotation.py), written while BRAKER runs; tested on the
lean smoke contig). BRAKER genes are classed by bedtools overlap with the polished Liftoff GFF3:
`shared` (same-strand CDS overlap; Liftoff model kept, CDS Jaccard recorded), `liftoff_locus`
(inside a Liftoff gene of any biotype / strand; reported only) or `new` (no Liftoff gene; added with
provenance=braker). New genes get hq / lq from protein length ≥ 50 aa, CDS repeat fraction < 0.5 and
evidence (≥ 1 intron matching an RNA-seq b2h or ProtHint intron hint, or ProtHint CDSpart coverage
≥ 0.5), and a reference-homology status from DIAMOND vs SaNama RefSeq proteins (ref_homolog ≥ 90 %
identity = likely Liftoff miss; ref_paralog; no_ref_hit = novel candidate). Outputs per ecotype in
`final/{eco}/`: merged GFF3 (gffread-valid), genes BED, new-gene GFF3 / protein / CDS FASTA, merged
proteins for BUSCO, per-gene class table, summary. Smoke contig (lean, 4.8 Mb): 149 BRAKER genes =
67 shared (median CDS Jaccard 0.96) / 22 liftoff_locus / 60 new (16 hq; 37 no_ref_hit, 3 of them
hq); 59 % of Liftoff protein-coding genes matched by a BRAKER CDS.

Smoke note: the siscowet smoke BRAKER run failed in GeneMark-ETP's masking-penalty step ("input
sequence size is too small", 11 reliable genes on the 5 Mb contig), which is a small-genome
artifact. 23.18 now logs stage failures to the chain log; this change is committed but **not**
deployed to raven while the pipeline that reads it is running.

**Remaining after BRAKER:** run 23.19 on both ecotypes, BUSCO on the merged proteins, functional
annotation of new hq genes (eggNOG / DIAMOND + InterProScan as in 18-build-gene-function-table.py),
browser tracks.
