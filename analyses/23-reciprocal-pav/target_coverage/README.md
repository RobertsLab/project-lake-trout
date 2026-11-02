# Target-ecotype read coverage of the reciprocal-PAV regions (23.23)

Built by `code/23.23-reciprocal-pav-target-coverage.py` on 2026-11-02 from the Phase 3 own-assembly
alignments (`analyses/23-native-methylation/align/{sample}.to_own.sorted.bam`, not committed).

## Question

23.4 calls a region of ecotype X's purged assembly "X-specific present" when all four fish of the
other ecotype have zero coverage there and the flanks are covered. It never checks whether X's own
four fish cover the region. Both assemblies pool all four fish of their ecotype, so a region could be
fixed in X, segregating within X, a retained haplotig, or an assembly artifact. Every tier-A gene in
23.7 requires a reciprocal-PAV call.

## Method

Per region and per target-ecotype fish: mean depth, fraction of bases at ≥1x / ≥3x, and the ratio of
region depth to the mean of the two 2 kb flanks on the same assembly (local single-copy baseline).
Per fish: `present` (frac ≥1x ≥ 0.80 and ratio ≥ 0.60), `absent` (frac ≥1x < 0.20), else `low`.
Per region: `fixed_present` (4/4 present), `polymorphic` (≥1 absent and ≥1 present), `low_all`,
`absent_all`, `mixed`. Null: the same regions shuffled three times to random positions on the same
contigs (bedtools shuffle), classified identically, to measure the chance rate at 2–9x per-fish depth.

## Result

| | siscowet-present (n = 340) | shuffled null (mean of 3) | lean-present (n = 343) | shuffled null (mean of 3) |
|---|---|---|---|---|
| fixed_present | **32 (9 %)** | 247 (73 %) | **21 (6 %)** | 193 (56 %) |
| polymorphic | 192 (57 %) | 43 (13 %) | 231 (67 %) | 74 (22 %) |
| mixed | 110 (32 %) | 49 (14 %) | 86 (25 %) | 75 (22 %) |
| low_all | 6 | 1 | 5 | 2 |
| absent_all | 0 | 0 | 0 | 0 |
| per-fish absent rate | 19–36 % | 3–7 % | 20–48 % | 3–14 % |
| pooled 4-fish depth ≥ 0.6 × flank | 155 (46 %) | | 178 (52 %) | |
| pooled 0.25–0.6 × flank | 145 (43 %) | | 142 (41 %) | |
| pooled < 0.25 × flank | 40 (12 %) | | 23 (7 %) | |

Per-fish median flank depth: siscowet 3.7–9.4x, lean 2.2–6.3x. At that depth a random region is
missed by one fish 3–14 % of the time, so `polymorphic` has a chance floor of 13–22 %; the observed
57–67 % is three to five times that, and the per-fish absent rates are four to six times the null.
Absent fish are spread across all four fish of each ecotype, not concentrated in the shallowest.

**Tier A (76 genes):** fixed_present 7, mixed 21, polymorphic 47, low_all 1. The seven fixed genes
are `gpia`, `LOC120029476` (KH-domain RNA-binding), `LOC120031768` (myotubularin-related 5-like),
`LOC120048669` (CEACAM5-like), `LOC120053938` (ROR-alpha A-like) in siscowet and `nek7`, `tmeff2a`
in lean. Median best pooled region/flank ratio: fixed 0.92, polymorphic 0.58, mixed 0.52.

## Interpretation

Most reciprocal-PAV regions are **not fixed in the ecotype whose assembly carries them**. About
two-thirds are segregating within that ecotype (at least one of its own fish lacks the sequence),
and roughly 40 % sit at half depth in the pooled reads (heterozygous in the assembled fish, or a
retained haplotig). The 23.4 call therefore mostly marks large indels that are polymorphic in both
ecotypes and happened to be absent from the four querying fish: with eight haplotypes sampled, an
allele at frequency 0.3 is absent from all four fish 6 % of the time, which across thousands of
segregating indels yields hundreds of apparent "ecotype-specific" regions.

Consequences for the evidence model:
- `reciprocal_pav` as used in 23.7 is a within-sample absence, not an ecotype-level one. Only the
  `fixed_present` subset (53 regions, 9 genes; 7 tier-A genes) has the property the line was meant to
  carry, and even that is fixed only among four fish.
- The 66 tier-A genes whose second line is `ref_divergence` have no independent second line (that
  line co-occurs with any PAV by construction), so tier A should be restated as: 7 genes with
  fixed-in-four-fish, bidirectional read evidence; the rest are within-ecotype polymorphisms.
- A population sample (tens of fish per ecotype) is the only way to turn these into allele-frequency
  differences; the current design cannot.

## Files

- `{eco}.region_coverage.tsv` — per region: per-fish depth / frac1x / frac3x / flank / ratio / status, pooled columns, class
- `{eco}.null{1,2,3}.region_coverage.tsv` — the same for shuffled regions
- `{eco}.gene_coverage.tsv` — genes overlapping regions (bedtools intersect with the native Liftoff genes), best class
- `tierA_target_coverage.tsv` — the 76 tier-A genes with their class
- `target_coverage_summary.txt` — the counts above; `run.log`, `run.pass1.log` — run logs
