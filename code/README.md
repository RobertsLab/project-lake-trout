`project-lake-trout/code`

---

Scripts are numbered in execution order. Plans and progress logs: `18-diff-annotation-phenotype-plan.md`,
`20-ecotype-genomes.md`, `21-ecotype-genome-comparison-plan.md`, `22-synteny-browser-integration-plan.md`,
`23-next-phase-research-plan.md` (Phases 1–5 with results).

| range | what | key outputs |
|---|---|---|
| `01`–`02` | RNA-seq differential expression (Ballgown) | `analyses/DEG-*`, `analyses/DET-*` |
| `03`–`10`, `16`–`17` | PacBio HiFi alignment, QC, 5mC calling (`04-pacbio/` holds the pbmm2 workflow) | `analyses/04-pacbio`, `08-*`, `09-*` |
| `11`–`12`, `15` | Reference-based PAV and differential PAV | `analyses/15-diff-pav` |
| `14` | Reference-based differential methylation (302 DMRs) | `analyses/14-diff-meth` |
| `18` | RefSeq function table, feature→gene assignment, candidates, GO | `analyses/18-annotation` (superseded by 23) |
| `13.1`–`13.3`, `20`–`22` | hifiasm assemblies, Liftoff, purge_dups + MCScanX synteny, browser tracks | `output/13.1-*`, `output/20.1-*`, `output/21.1-*`, `output/22-*` |
| `23.1`–`23.3` | Phase 1: SyRI structural-variant map, PAV cross-validation | `analyses/23-genome-sv` (not committed) |
| `23.4`–`23.7` | Phase 2: reciprocal PAV, gene sets, integration, refined evidence model | `analyses/23-reciprocal-pav`, `23-gene-sets`, `23-integration/refined` |
| `23.8`–`23.16` | Phase 3: native methylation, CpG scan, DMR tracks, fixed variants, allele-specific methylation | `analyses/23-native-methylation` |
| `23.17-lean-reference-baseline.py`, `23.18-reciprocal-pav-exon.py` | Lean vs published-lean baseline; exon placement of reciprocal-PAV segments | `analyses/23-lean-baseline` |
| `23.17-denovo-setup.sh` … `23.22` | Phase 4: RepeatModeler2 + BRAKER3, Liftoff merge, function, browser tracks, reciprocal absence screen | `analyses/23-denovo-annotation` |

Numbering note: `23.17` and `23.18` are each used twice (the lean-baseline Python scripts and the
Phase 4 shell scripts were written in parallel). The file names are unambiguous; references in the plan
log say which.
