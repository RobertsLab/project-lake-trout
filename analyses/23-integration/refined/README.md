# Refined two-genome evidence model (23.7)

Re-scoring of `../gene_evidence_matrix.tsv` with independent evidence lines, a corroboration tier for ecotype-only genes, and tandem-cluster-collapsed GO tests. Built by `code/23.7-refine-evidence-model.py`; rationale in the docstring and in `code/23-next-phase-research-plan.md`.

## Evidence lines (genes per line)

| line | genes | basis |
|---|---|---|
| reciprocal_pav | 114 | native, bidirectional read coverage (Phase 2) |
| ecotype_only_corroborated | 4,382 | lifted into one assembly only AND in the other's unmapped list (lean 1,434 / siscowet 2,948) |
| cnv | 2,359 | Liftoff copy-number divergence |
| dmr | 181 | reference-based DMR within 5 kb |
| ref_divergence | 16,068 | SyRI SV (12,392) ∪ reference PAV (1,263) ∪ Liftoff-only ecotype_only (3,449) |

SV line source: reused from analyses/23-integration/gene_evidence_matrix.tsv (SyRI outputs untracked).

## Candidate tiers

| tier | definition | genes | of which no caution flag |
|---|---|---|---|
| A | reciprocal_pav + ≥1 other line | 77 | 66 |
| B | corroborated ecotype_only or cnv + ≥1 other line, no reciprocal PAV | 4,011 | 2,613 |
| C | dmr + ref_divergence only (reference-only association) | 15 | 6 |

Compare 23.6: 2,895 "two-genome candidates" at ≥2 of 6 lines; refined tier A+B = 4,088, tier A = 77.

### Tier A (native-anchored), protein-coding, no caution flag

| gene | symbol | product | present in | lines |
|---|---|---|---|---|
| gene-LOC120022070 | LOC120022070 | growth factor receptor-bound protein 10-like | lean | reciprocal_pav,cnv,ref_divergence |
| gene-LOC120024897 | LOC120024897 | cGMP-inhibited 3',5'-cyclic phosphodiesterase B-li | siscowet | reciprocal_pav,ecotype_only_corroborated,ref_divergence |
| gene-LOC120034415 | LOC120034415 | inactive rhomboid protein 1-like | lean | reciprocal_pav,ecotype_only_corroborated,ref_divergence |
| gene-LOC120050135 | LOC120050135 | nucleosome assembly protein 1-like 1-A | siscowet | reciprocal_pav,ecotype_only_corroborated,ref_divergence |
| gene-gpia | gpia | glucose-6-phosphate isomerase a | siscowet | reciprocal_pav,cnv,ref_divergence |
| gene-ttll9 | ttll9 | tubulin tyrosine ligase-like family, member 9 | siscowet | reciprocal_pav,cnv,ref_divergence |
| gene-LOC120020198 | LOC120020198 | synapsin-2-like | lean | reciprocal_pav,ref_divergence |
| gene-LOC120021129 | LOC120021129 | CMP-N-acetylneuraminate-beta-galactosamide-alpha-2 | lean | reciprocal_pav,ref_divergence |
| gene-LOC120024458 | LOC120024458 | sodium/potassium-transporting ATPase subunit beta- | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120024920 | LOC120024920 | major histocompatibility complex class I-related g | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120025009 | LOC120025009 | protein kinase C-binding protein NELL1-like | lean | reciprocal_pav,ref_divergence |
| gene-LOC120025522 | LOC120025522 | ankyrin-1-like | lean | reciprocal_pav,ref_divergence |
| gene-LOC120025537 | LOC120025537 | astrotactin-2-like | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120025549 | LOC120025549 | leucine-rich repeat transmembrane neuronal protein | lean | reciprocal_pav,ref_divergence |
| gene-LOC120027905 | LOC120027905 | band 4.1-like protein 3 | lean | reciprocal_pav,ref_divergence |
| gene-LOC120028089 | LOC120028089 | phosphoenolpyruvate carboxykinase [GTP], mitochond | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120029476 | LOC120029476 | KH domain-containing, RNA-binding, signal transduc | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120030672 | LOC120030672 | calsequestrin-1-like | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120031305 | LOC120031305 | opioid-binding protein/cell adhesion molecule-like | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120031768 | LOC120031768 | myotubularin-related protein 5-like | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120032210 | LOC120032210 | kinesin-like protein KIF21A | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120032728 | LOC120032728 | E3 ubiquitin-protein ligase TRIM39-like | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120038650 | LOC120038650 | prickle-like protein 1 | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120040161 | LOC120040161 | sperm-associated antigen 1A-like | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120040412 | LOC120040412 | plexin-B2-like | lean | reciprocal_pav,ref_divergence |
| gene-LOC120042766 | LOC120042766 | tumor necrosis factor receptor superfamily member  | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120043213 | LOC120043213 | extracellular serine/threonine protein kinase FAM2 | lean | reciprocal_pav,ref_divergence |
| gene-LOC120048047 | LOC120048047 | alpha-1,3-mannosyl-glycoprotein 4-beta-N-acetylglu | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120048669 | LOC120048669 | carcinoembryonic antigen-related cell adhesion mol | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120049548 | LOC120049548 | glucosidase 2 subunit beta-like | lean | reciprocal_pav,ref_divergence |
| gene-LOC120053194 | LOC120053194 | NT-3 growth factor receptor-like | lean | reciprocal_pav,ref_divergence |
| gene-LOC120053598 | LOC120053598 | lisH domain-containing protein ARMC9 | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120053718 | LOC120053718 | cAMP-specific 3',5'-cyclic phosphodiesterase 4D-li | lean | reciprocal_pav,ref_divergence |
| gene-LOC120053938 | LOC120053938 | nuclear receptor ROR-alpha A-like | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120054371 | LOC120054371 | lipoma-preferred partner homolog | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120054763 | LOC120054763 | DNA-directed RNA polymerase, mitochondrial-like | lean | reciprocal_pav,ref_divergence |
| gene-LOC120055734 | LOC120055734 | LARGE xylosyl- and glucuronyltransferase 2-like | lean | reciprocal_pav,ref_divergence |
| gene-LOC120056781 | LOC120056781 | kin of IRRE-like protein 3 | lean | reciprocal_pav,ref_divergence |
| gene-LOC120056805 | LOC120056805 | troponin T, fast skeletal muscle isoforms-like | siscowet | reciprocal_pav,ref_divergence |
| gene-LOC120060868 | LOC120060868 | agrin-like | siscowet | reciprocal_pav,ref_divergence |

## GO after tandem-cluster collapsing

Genome-wide, 2,965 clusters of ≥2 same-family genes within 50 kb were collapsed to one representative each (background 34,266 genes → 31,285 clusters). Lipid / calcium terms at FDR<0.25 in any set, collapsed vs uncollapsed: see `phenotype_survival_refined.tsv`. Surviving after collapsing: 15 terms (calcium 14, lipid 1).

| set | group | term | k (clusters) | FDR collapsed | FDR uncollapsed |
|---|---|---|---|---|---|
| cnv | calcium | calcium sensitive guanylate cyclase activator activity |  | — | 1.578e-01 |
| lean_only_corroborated | calcium | calcium ion transmembrane transporter activity | 14 | 1.290e-01 | 2.621e-02 |
| lean_only_corroborated | calcium | calcium channel activity | 12 | 1.290e-01 | 2.606e-02 |
| lean_only_corroborated | calcium | voltage-gated calcium channel complex | 6 | 1.536e-01 | 5.291e-02 |
| lean_only_corroborated | calcium | calcium ion transmembrane transport | 11 | 1.648e-01 | 6.723e-02 |
| lean_only_corroborated | calcium | high voltage-gated calcium channel activity | 4 | 1.884e-01 | 1.631e-01 |
| lean_only_corroborated | calcium | calcium ion import | 5 | 2.207e-01 | 2.096e-01 |
| lean_only_corroborated | calcium | calcium ion transmembrane import into cytosol | 5 | 2.207e-01 | 2.096e-01 |
| lean_only_corroborated | calcium | regulation of voltage-gated calcium channel activity | 3 | 2.474e-01 | 2.456e-01 |
| lean_only_corroborated | calcium | calcium channel complex | 6 | 2.803e-01 | 1.631e-01 |
| lean_only_corroborated | calcium | voltage-gated calcium channel activity | 9 | 2.835e-02 | 8.363e-04 |
| lean_only_corroborated | calcium | calcium ion transport | 15 | 6.214e-02 | 9.246e-03 |
| lean_only_corroborated | lipid | long-chain fatty acid metabolic process | 4 | 2.430e-01 | 3.561e-01 |
| reciprocal_pav | calcium | calcium ion binding | 8 | 1.049e-01 | 9.357e-02 |
| siscowet_only_corroborated | calcium | regulation of calcium ion-dependent exocytosis | 9 | 1.324e-01 | 4.779e-01 |
| siscowet_only_corroborated | calcium | cellular response to calcium ion | 10 | 2.258e-01 | 5.129e-01 |
| tierA | calcium | calcium ion binding | 5 | 2.072e-01 | 2.232e-01 |
| tierAB | lipid | short-chain fatty acid catabolic process | 3 | 3.660e-01 | 2.278e-01 |

### Largest collapsed clusters per study set

| set | locus | family | study genes | cluster genes | example |
|---|---|---|---|---|---|
| cnv | NC_052308.1:33,025,163-33,087,353 | spliceosomal | 33 | 46 | U4 spliceosomal RNA |
| cnv | NC_052311.1:29,716,908-30,037,076 | protocadherin | 27 | 46 | protocadherin alpha-4-like |
| cnv | NC_052335.1:26,147,236-26,155,916 | small | 12 | 36 | small nucleolar RNA SNORD50 |
| cnv | NC_052308.1:72,299,678-72,441,127 | histone | 11 | 19 | histone H2B |
| cnv | NW_024058168.1:637-20,029 | small | 11 | 20 | small nucleolar RNA U3 |
| cnv | NC_052314.1:2,399,403-2,423,268 | small | 11 | 18 | small nucleolar RNA U3 |
| cnv | NC_052319.1:22,927,602-23,120,998 | gamma | 11 | 18 | gamma-crystallin M2-like |
| cnv | NW_024057805.1:9,773-83,599 | spliceosomal | 10 | 10 | U1 spliceosomal RNA |
| cnv | NW_024057745.1:22,793-87,701 | spliceosomal | 9 | 15 | U1 spliceosomal RNA |
| cnv | NW_024058207.1:1,536-62,831 | ribosomal | 9 | 26 | 5S ribosomal RNA |
| cnv | NC_052308.1:42,960,918-43,157,673 | loc | 8 | 10 | uncharacterized LOC120063518 |
| cnv | NW_024061618.1:20,152-68,241 | spliceosomal | 7 | 9 | U1 spliceosomal RNA |
| cnv | NC_052308.1:65,159,504-65,440,296 | keratin | 7 | 15 | keratin, type I cytoskeletal 18-A-like |
| cnv | NW_024058286.1:142,629-172,621 | spliceosomal | 7 | 8 | U1 spliceosomal RNA |
| cnv | NC_052332.1:142,268-144,191 | small | 6 | 6 | small nucleolar RNA SNORD36 |
| cnv | NW_024057745.1:1,055-63,459 | ribosomal | 6 | 26 | 5S ribosomal RNA |
| cnv | NC_052308.1:43,215,320-43,364,997 | loc | 6 | 10 | uncharacterized LOC120063534 |
| cnv | NC_052309.1:63,585,120-63,890,752 | interferon | 5 | 22 | interferon alpha-1-like |
| cnv | NC_052326.1:37,183,416-37,253,971 | green | 5 | 5 | green-sensitive opsin-like |
| cnv | NC_052322.1:2,765,268-2,772,096 | ribosomal | 5 | 5 | 5S ribosomal RNA |
| cnv | NW_024060776.1:345-20,138 | spliceosomal | 5 | 9 | U4 spliceosomal RNA |
| cnv | NW_024061248.1:94,514-110,035 | spliceosomal | 5 | 6 | U1 spliceosomal RNA |
| cnv | NC_052348.1:20,187,668-20,224,456 | zinc | 5 | 5 | zinc finger protein 629-like |
| cnv | NC_052331.1:33,568,406-33,571,360 | small | 5 | 5 | small nucleolar RNA SNORD2 |
| cnv | NW_024058273.1:4,572-36,436 | spliceosomal | 4 | 6 | U1 spliceosomal RNA |

## Fate of the Phase-18 headline candidates

| gene | symbol | product | Phase-18 role | refined tier | lines | native support |
|---|---|---|---|---|---|---|
| gene-LOC120032414 | LOC120032414 | zinc finger protein 883-like | convergent | C | dmr,ref_divergence | N |
| gene-LOC120040411 | LOC120040411 | gastrula zinc finger protein XlCGF57.1-like | convergent | B | cnv,dmr,ref_divergence | Y |
| gene-LOC120043843 | LOC120043843 | septin-9-like | convergent | C | dmr,ref_divergence | N |
| gene-LOC120039781 | LOC120039781 |  | convergent | C | dmr,ref_divergence | N |
| gene-LOC120050008 | LOC120050008 | angiopoietin-related protein 5-like | lipid_axis | B | ecotype_only_corroborated,ref_divergence | Y |
| gene-LOC120041635 | LOC120041635 | 2-acylglycerol O-acyltransferase 2-A-like | lipid_axis | none | ref_divergence | N |
| gene-LOC120029926 | LOC120029926 | epoxide hydrolase 1-like | lipid_axis | C | dmr,ref_divergence | N |

## Reading guide

- Only tier A carries read-level, bidirectional evidence on native coordinates. Tier B rests on two Liftoff annotations agreeing; tier C is the Phase-18 style single-reference association and is retained for continuity only.
- `ref_divergence` is one line no matter how many of SV / reference PAV / Liftoff-only presence a gene carries: they share a cause (distance from the lean reference).
- `caution` flags noncoding biotypes, repeat / tandem families, and clusters of ≥5 genes; rank on the flag-free subset first.
- Phase 3 (native methylation) will add the first native epigenetic line and can promote tier-B genes; Phase 4 (BRAKER3) addresses genes Liftoff cannot see at all.
