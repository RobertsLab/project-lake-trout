# Lean vs published lean as a baseline (23.17)

Built by `code/23.17-lean-reference-baseline.py` from the three existing SyRI runs (23.1 / 23.3).
No new alignment.

**The reference fish.** SaNama_1.0 (GCF_016432855.1, BioSample SAMN16985724) is a Seneca-strain,
gynogenetic doubled-haploid lean female, so it carries one homozygous lean haplotype. Our assemblies are
hifiasm primary contigs from four pooled fish per ecotype, purged and RagTag-anchored to that reference.
Lean vs reference therefore measures within-lean variation (strain / population) plus assembly artifacts.
If the two are not the same population, it also includes population divergence, which makes it a generous
baseline.

## 1. Genome-wide scale (`genome_scale.tsv`)

| | lean vs ref | siscowet vs ref | lean vs siscowet | lean-ref / lean-sisco |
|---|---|---|---|---|
| SNPs | 1.14 M | 1.23 M | 1.24 M | 0.92 |
| deletions (count) | 232 k | 258 k | 224 k | 1.03 |
| insertions (count) | 250 k | 281 k | 225 k | 1.11 |
| highly diverged (Mb) | 231 | 175 | 197 | 1.17 |
| syntenic (Mb, ref side) | 1,679 | 1,625 | 1,519 | — |
| inversions (Mb) | 9.9 | 21.3 | 29.3 | 0.34 |
| translocations (Mb) | 62.6 | 107.9 | 120.8 | 0.52 |
| duplications, query (Mb) | 15.3 | 38.7 | 33.7 | 0.45 |
| not aligned, query (Mb) | 230 | 380 | 491 | 0.47 |

At the sequence level (SNPs, indels, diverged blocks), two lean genomes differ about as much as lean and
siscowet do. Only large rearrangements are 2–3× larger between ecotypes. Siscowet vs reference shows most
of that excess on its own, and the siscowet assembly is larger (2.84 vs 2.58 Gb), so assembly differences
can't be ruled out as the cause.

## 2. Which side of each lean–siscowet SV differs from the published lean (`polarized_sv.tsv`)

100,671 lean–siscowet events of 50 bp or more. 83,552 have both sides placed; the rest are one-sided
unaligned blocks. A side "differs" when at least 50% of it (padded 50 bp) lies in that assembly's own SVs
against the reference. These SVs cover 29% of the lean genome and 36% of the siscowet genome. The null
shuffles positions within chromosomes (3 seeds); it ignores SV clustering, so it is lenient.

| class | events | siscowet side | lean side | both | neither |
|---|---|---|---|---|---|
| all | 83,552 | 0.30 (null 0.26) | 0.25 (0.18) | 0.36 (0.11) | 0.10 (0.45) |
| indel | 60,866 | 0.31 | 0.28 | 0.34 | 0.07 |
| diverged | 10,913 | 0.22 | 0.21 | 0.45 | 0.12 |
| translocation | 3,558 | 0.22 | 0.13 | 0.51 | 0.15 |
| duplication | 6,854 | 0.33 | 0.06 | 0.32 | 0.29 |
| inversion | 335 | 0.39 | 0.24 | 0.37 | 0.01 |

Fractions of events; the bp-weighted fractions are in the TSV.

- **lean side** means siscowet matches the published lean and our lean assembly is the odd one out. A
  quarter of lean–siscowet SVs (99 Mb) are of this kind, so they are not lean-ecotype differences.
- For **indels and diverged blocks**, lean side and siscowet side are about equal. This is the pattern
  expected from polymorphism or assembly noise, not from an ecotype axis.
- **Duplications, translocations and inversions** skew toward the siscowet side, which fits either real
  siscowet structure or the larger siscowet assembly.
- **both** (all three genomes differ) is three times the null. These are hotspots, most likely repeats and
  regions that are hard to assemble.

## 3. Reciprocal-PAV genes (`reciprocal_gene_absent_fraction.tsv`, `gene_alignment*.tsv`)

All reciprocal-PAV and corroborated ecotype-only genes are Liftoff-lifted reference genes, so every one is
annotated in the published lean genome. That includes the 71 siscowet-present genes absent from all four of
our lean fish.

A reciprocal-PAV gene is usually **not a whole-gene absence.** It is a large gene (median 56–62 kb)
containing a short native presence/absence segment: the median absent fraction is 3% of the gene, and only
6 of 115 genes are at least 50% absent. Five of those six are extra Liftoff copies (`_N` suffix), and the
sixth is an uncharacterized LOC (LOC120056879, tier A). This fits the SyRI check: the reference interval of these genes is
almost fully aligned by both assemblies (median aligned fraction 0.99–1.0). Corroborated ecotype-only genes
behave as expected: the missing assembly does not align across them (median 0.0 for siscowet-only genes in
lean).

**What this means for tier A.** Reciprocal PAV is real, read-level and bidirectional, but at gene level it
mostly marks sub-gene segments (1 kb resolution). Whether these segments hit exons, introns or flanks has
not been tested, and that is the next check before any tier-A gene is described as present/absent.

## Files

| file | content |
|---|---|
| genome_scale.tsv | SyRI summary rows × comparison, with ratios to lean vs siscowet |
| polarized_sv.tsv / polarized_sv_null.tsv | counts and bp per SV class × polarity, observed / shuffled |
| mask_fractions.tsv | genome fraction in each assembly's SVs against the reference |
| gene_alignment.tsv / gene_alignment_summary.tsv | aligned fraction of each gene's reference interval, per assembly |
| reciprocal_gene_absent_fraction.tsv | per reciprocal-PAV gene, length and bp inside a native present-only region (computed locally from 23.4 BEDs) |
