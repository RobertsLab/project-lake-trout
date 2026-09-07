#!/usr/bin/env python3
"""
23.7-refine-evidence-model.py  (Phase 2 follow-up; see progress log in
code/23-next-phase-research-plan.md)

WHY. The 23.6 evidence matrix counted six "independent" lines per gene, but three of them --
`ecotype_only`, `ref_pav`, `sv` -- all measure divergence from the lean-background reference and
co-occur heavily (1,565 genes carry ecotype_only+ref_pav; 543 carry ref_pav+sv), so `n_lines`
over-stated independence. Two further problems: `ecotype_only` (14,647 genes) is dominated by
Liftoff sensitivity differences between two contig-level assemblies (only ~30 % of each "only" set
is corroborated by the other assembly's unmapped list), and the GO tests were driven by tandem
clusters (18 hemoglobin genes at one locus on NC_052347.1 lead lean_only; histone and tRNA
clusters dominate cnv), the same artifact already flagged for the DMR set in Phase 18.

WHAT THIS DOES.
  1. Collapses ecotype_only / ref_pav / sv into ONE line, `ref_divergence`.
  2. Splits ecotype_only into two tiers: `corroborated` (gene lifted into only one assembly AND
     explicitly listed in the other assembly's Liftoff unmapped_features -- both assemblies agree)
     keeps its own line, `ecotype_only_corroborated`; the Liftoff-only remainder is demoted into
     `ref_divergence`.
  3. Re-scores every gene on FIVE lines that are independent in data type or direction:
        reciprocal_pav            native, bidirectional read evidence         (Phase 2)
        ecotype_only_corroborated two-assembly agreement on gene absence      (Phase 2)
        cnv                       Liftoff copy-number divergence              (Phase 2)
        dmr                       reference-based methylation                 (Phase 14/18)
        ref_divergence            SyRI SV | reference PAV | Liftoff-only ecotype_only
     and assigns a tier:
        A  reciprocal_pav + >=1 other line          (native-anchored two-genome candidate)
        B  corroborated ecotype_only or cnv + >=1   (two-assembly candidate, no read support yet)
        C  >=2 lines but only dmr + ref_divergence  (reference-only association, Phase-18 style)
  4. Re-runs GO over-representation (18.3 machinery) with TANDEM-CLUSTER COLLAPSING: genes of the
     same product family (first product token) within 50 kb of each other on the same chromosome
     are collapsed genome-wide to one representative carrying the union of their GO terms, for the
     study set AND the background. Uncollapsed results are kept alongside for comparison.
  5. Reports the fate of the Phase-18 headline candidates (4 convergent genes, lipid-axis genes)
     under the refined model, and a cluster-diagnostics table naming the collapsed clusters.

INPUTS (all in repo):
  analyses/23-integration/gene_evidence_matrix.tsv   (sv flag; SyRI raw outputs are not tracked)
  analyses/23-reciprocal-pav/{lean,siscowet}_specific.genes.tsv
  analyses/23-gene-sets/{lean_only,siscowet_only,cnv}.tsv
  analyses/18-annotation/{gene_function_table,dmr_gene_assignments,pav_gene_assignments,
                          integrated_candidate_genes}.tsv
  analyses/18-annotation/raw/go-basic.obo   (re-download: http://purl.obolibrary.org/obo/go/go-basic.obo)

OUTPUTS (analyses/23-integration/refined/):
  gene_evidence_refined.tsv          every gene with >=1 line, refined flags, tier, cluster
  two_genome_candidates_refined.tsv  tier A + B
  go_enrichment_{set}.collapsed.tsv / .uncollapsed.tsv
  cluster_diagnostics.tsv            tandem clusters (>=3 study genes) collapsed per study set
  phenotype_survival_refined.tsv     lipid / calcium terms, collapsed vs uncollapsed FDR
  phase18_candidate_fate.tsv         what happened to the Phase-18 headline genes
  README.md

Author: Generated for project-lake-trout
Date: 2026-09-06
"""

import csv
import importlib.util
import re
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
INT = BASE / "analyses" / "23-integration"
RECIP = BASE / "analyses" / "23-reciprocal-pav"
GENESETS = BASE / "analyses" / "23-gene-sets"
ANN = BASE / "analyses" / "18-annotation"
OUT = INT / "refined"

CLUSTER_GAP = 50_000          # bp between consecutive same-family genes to chain a tandem cluster
FAMILY_STOP = {"protein", "uncharacterized", "putative", "probable", "hypothetical", "the", "and"}
PHENO_LIPID = re.compile(
    r"lipid|fatty.?acid|triglyceride|triacylglycerol|sterol|cholesterol|lipo|"
    r"buoyan|swim.?bladder|gas bladder|adipos", re.I)
PHENO_CALCIUM = re.compile(r"calcium|\bca2\b|calcium.ion|voltage.gated calcium", re.I)
REPEAT_FAMILY = re.compile(r"histone|tRNA|zinc finger|NACHT|immunoglobulin|hemoglobin|"
                           r"uncharacterized|transposon|retrotranspos|gag-pol|reverse transcriptase", re.I)
NONCODING = {"tRNA", "rRNA", "snRNA", "snoRNA", "misc_RNA", "ncRNA", "lncRNA", "pseudogene"}

PHASE18_CONVERGENT = ["gene-LOC120032414", "gene-LOC120040411",
                      "gene-LOC120043843", "gene-LOC120039781"]
PHASE18_LIPID = re.compile(r"angiopoietin-related protein 5|acylglycerol O-acyltransferase|"
                           r"epoxide hydrolase 1", re.I)


def load_module(fname, modname):
    spec = importlib.util.spec_from_file_location(modname, BASE / "code" / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def canon(gid):
    """Liftoff -copies suffixes (gene-LOC1200_1) collapse onto the reference gene ID."""
    return re.sub(r"_\d+$", "", gid)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def family_key(product):
    p = product.lower().replace("-like", "")
    toks = [t for t in re.split(r"[^a-z]+", p) if len(t) >= 3 and t not in FAMILY_STOP]
    return toks[0] if toks else ""


# ------------------------------------------------------------------- tandem clusters
def build_clusters(fn):
    """Genome-wide: chain same-family genes within CLUSTER_GAP on a chromosome. Returns
    rep_of[gid] -> representative gid, members[rep] -> [gids]."""
    by_key = defaultdict(list)
    for gid, r in fn.items():
        k = family_key(r["product"])
        if not k or not r["chrom"]:
            continue
        by_key[(r["chrom"], k)].append((int(r["start"]), int(r["end"]), gid))
    rep_of, members = {}, {}
    for _, lst in by_key.items():
        lst.sort()
        cur, cur_end = [], -1
        for s, e, g in lst:
            if cur and s - cur_end > CLUSTER_GAP:
                rep = cur[0]
                members[rep] = list(cur)
                for m in cur:
                    rep_of[m] = rep
                cur = []
            cur.append(g)
            cur_end = max(cur_end, e)
        if cur:
            rep = cur[0]
            members[rep] = list(cur)
            for m in cur:
                rep_of[m] = rep
    for gid in fn:                      # singletons / no product
        rep_of.setdefault(gid, gid)
        members.setdefault(gid, [gid])
    return rep_of, members


# ------------------------------------------------------------------- main
def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- annotation ----
    fn = {}
    for r in read_tsv(ANN / "gene_function_table.tsv"):
        fn[r["gene_id"]] = r
    rep_of, members = build_clusters(fn)
    n_multi = sum(1 for m in members.values() if len(m) > 1)
    print(f"  tandem clusters (>=2 genes): {n_multi:,}; "
          f"genes in clusters: {sum(len(m) for m in members.values() if len(m) > 1):,}")

    # ---- evidence sets (canonical reference gene IDs) ----
    recip_dir = {}
    for side, fname in (("lean", "lean_specific.genes.tsv"),
                        ("siscowet", "siscowet_specific.genes.tsv")):
        for r in read_tsv(RECIP / fname):
            recip_dir[canon(r["gene_id"])] = side       # gene present only in <side>
    recip = set(recip_dir)

    eco_tier, eco_side, copies = {}, {}, {}
    for side, fname, ccol, ucol in (("lean", "lean_only.tsv", "copies_lean", "in_siscowet_unmapped"),
                                    ("siscowet", "siscowet_only.tsv", "copies_siscowet", "in_lean_unmapped")):
        for r in read_tsv(GENESETS / fname):
            g = canon(r["gene_id"])
            eco_side[g] = side
            eco_tier[g] = "corroborated" if r.get(ucol, "").strip() == "Y" else "liftoff_only"
    eco_corr = {g for g, t in eco_tier.items() if t == "corroborated"}
    eco_lift = {g for g, t in eco_tier.items() if t == "liftoff_only"}

    cnv = set()
    for r in read_tsv(GENESETS / "cnv.tsv"):
        g = canon(r["gene_id"])
        cnv.add(g)
        copies[g] = (r["copies_lean"], r["copies_siscowet"])

    dmr = {canon(r["gene_id"]) for r in read_tsv(ANN / "dmr_gene_assignments.tsv") if r["gene_id"]}
    ref_pav = {canon(r["gene_id"]) for r in read_tsv(ANN / "pav_gene_assignments.tsv") if r["gene_id"]}

    # sv: recompute from SyRI if the (untracked) outputs are present, else reuse the 23.6 matrix
    old = {canon(r["gene_id"]): r for r in read_tsv(INT / "gene_evidence_matrix.tsv")}
    syri_dir = BASE / "analyses" / "23-genome-sv" / "syri"
    if (syri_dir / "lean_vs_ref.syri.out").exists() and (syri_dir / "sisco_vs_ref.syri.out").exists():
        m236 = load_module("23.6-integrate-phase2.py", "integrate_236")
        sv = {canon(g) for g in m236.sv_affected_genes()}
        sv_src = "recomputed from SyRI"
    else:
        sv = {g for g, r in old.items() if r.get("sv") == "1"}
        sv_src = "reused from analyses/23-integration/gene_evidence_matrix.tsv (SyRI outputs untracked)"
    print(f"  sv line: {len(sv):,} genes ({sv_src})")

    ref_div = sv | ref_pav | eco_lift

    lines = {
        "reciprocal_pav": recip,
        "ecotype_only_corroborated": eco_corr,
        "cnv": cnv,
        "dmr": dmr,
        "ref_divergence": ref_div,
    }
    for k, v in lines.items():
        print(f"  {k:26s}: {len(v):,} genes")

    # Phase-18 rank scores for context
    p18 = {r["gene_id"]: r for r in read_tsv(ANN / "integrated_candidate_genes.tsv")}

    # ---- refined matrix ----
    all_genes = set().union(*lines.values())
    flag_cols = list(lines)
    rows = []
    for g in all_genes:
        flags = {c: int(g in lines[c]) for c in flag_cols}
        n = sum(flags.values())
        if flags["reciprocal_pav"] and n >= 2:
            tier = "A"
        elif (flags["ecotype_only_corroborated"] or flags["cnv"]) and n >= 2:
            tier = "B"
        elif n >= 2:
            tier = "C"
        else:
            tier = ""
        a = fn.get(g, {})
        rep = rep_of.get(g, g)
        csize = len(members.get(rep, [g]))
        product = a.get("product", "")
        biotype = a.get("biotype", "")
        caution = []
        if biotype in NONCODING:
            caution.append("noncoding")
        if REPEAT_FAMILY.search(product):
            caution.append("repeat_or_tandem_family")
        if csize >= 5:
            caution.append(f"tandem_cluster_n{csize}")
        cl, cs = copies.get(g, ("", ""))
        rows.append({
            "gene_id": g, "symbol": a.get("symbol", ""), "product": product, "biotype": biotype,
            "chrom": a.get("chrom", ""), "start": a.get("start", ""), "end": a.get("end", ""),
            "n_go": a.get("n_go", "0"),
            "tier": tier, "n_lines_refined": n,
            "lines_refined": ",".join(c for c in flag_cols if flags[c]),
            **flags,
            "reciprocal_present_in": recip_dir.get(g, ""),
            "ecotype_only_tier": eco_tier.get(g, ""),
            "ecotype_only_side": eco_side.get(g, ""),
            "copies_lean": cl, "copies_siscowet": cs,
            "sv": int(g in sv), "ref_pav": int(g in ref_pav),
            "n_lines_23_6": old.get(g, {}).get("n_lines", ""),
            "cluster_rep": rep, "cluster_size": csize,
            "caution": ";".join(caution),
            "phase18_convergent": int(g in PHASE18_CONVERGENT),
            "phase18_rank_score": p18.get(g, {}).get("rank_score", ""),
        })
    tier_order = {"A": 0, "B": 1, "C": 2, "": 3}
    rows.sort(key=lambda r: (tier_order[r["tier"]], -r["n_lines_refined"],
                             r["caution"] != "", r["gene_id"]))
    cols = list(rows[0].keys())

    def write_rows(path, rs):
        with open(path, "w") as o:
            o.write("\t".join(cols) + "\n")
            for r in rs:
                o.write("\t".join(str(r[c]) for c in cols) + "\n")

    write_rows(OUT / "gene_evidence_refined.tsv", rows)
    cand = [r for r in rows if r["tier"] in ("A", "B")]
    write_rows(OUT / "two_genome_candidates_refined.tsv", cand)
    tiers = defaultdict(int)
    for r in rows:
        tiers[r["tier"] or "single"] += 1
    tierA = [r for r in rows if r["tier"] == "A"]
    tierA_clean = [r for r in tierA if not r["caution"]]
    print(f"  tiers: {dict(tiers)}; tier-A without caution flags: {len(tierA_clean)}")

    # ---- GO: collapsed vs uncollapsed ----
    go = load_module("18.3-go-enrichment.py", "go_enrichment_183")
    print("  loading GO DAG ...")
    name, ns, parents, alt = go.parse_obo(go.OBO)
    _, anc_fn = go.build_ancestors(parents)
    gene_go = go.load_gene_go(name, alt, anc_fn)

    def term_index(gg):
        t2g = defaultdict(set)
        for gid, terms in gg.items():
            for t in terms:
                t2g[t].add(gid)
        return t2g

    universe_u = set(gene_go)
    t2g_u = term_index(gene_go)

    rep_go = defaultdict(set)
    for gid, terms in gene_go.items():
        rep_go[rep_of.get(gid, gid)] |= terms
    universe_c = set(rep_go)
    t2g_c = term_index(rep_go)
    print(f"  GO background: {len(universe_u):,} genes -> {len(universe_c):,} collapsed clusters")

    study_sets = {
        "reciprocal_pav": recip,
        "tierA": {r["gene_id"] for r in tierA},
        "tierAB": {r["gene_id"] for r in cand},
        "lean_only_corroborated": {g for g in eco_corr if eco_side[g] == "lean"},
        "siscowet_only_corroborated": {g for g in eco_corr if eco_side[g] == "siscowet"},
        "cnv": cnv,
    }
    pheno = {}      # (set, go_id) -> row
    diag = []
    for key, study in study_sets.items():
        # uncollapsed
        rows_u, n_u, _ = go.ora(study, t2g_u, universe_u, name, ns)
        fdr_u = {r["go_id"]: r["fdr"] for r in rows_u}
        # collapsed
        study_c = {rep_of.get(g, g) for g in study}
        rows_c, n_c, _ = go.ora(study_c, t2g_c, universe_c, name, ns)
        sig_u = sum(1 for r in rows_u if r["fdr"] < 0.1)
        sig_c = sum(1 for r in rows_c if r["fdr"] < 0.1)
        print(f"  [{key:27s}] genes {n_u:5d} -> clusters {n_c:5d} | "
              f"terms FDR<0.1: {sig_u:3d} -> {sig_c:3d}")
        for r in rows_c:
            grp = ("lipid" if PHENO_LIPID.search(r["name"]) else
                   "calcium" if PHENO_CALCIUM.search(r["name"]) else None)
            if grp and (r["fdr"] < 0.25 or fdr_u.get(r["go_id"], 1) < 0.25):
                pheno[(key, r["go_id"])] = {
                    "set": key, "group": grp, "go_id": r["go_id"], "name": r["name"],
                    "fold_collapsed": r["fold_enrichment"], "k_collapsed": r["study_k"],
                    "fdr_collapsed": f"{r['fdr']:.3e}",
                    "fdr_uncollapsed": f"{fdr_u[r['go_id']]:.3e}" if r["go_id"] in fdr_u else "",
                    "survives_collapsed_fdr0.25": "Y" if r["fdr"] < 0.25 else "N",
                    "genes": r["genes"]}
        for r in rows_u:            # terms that vanish entirely after collapsing
            grp = ("lipid" if PHENO_LIPID.search(r["name"]) else
                   "calcium" if PHENO_CALCIUM.search(r["name"]) else None)
            if grp and r["fdr"] < 0.25 and (key, r["go_id"]) not in pheno:
                pheno[(key, r["go_id"])] = {
                    "set": key, "group": grp, "go_id": r["go_id"], "name": r["name"],
                    "fold_collapsed": "", "k_collapsed": "", "fdr_collapsed": "",
                    "fdr_uncollapsed": f"{r['fdr']:.3e}",
                    "survives_collapsed_fdr0.25": "N (term dropped: <3 clusters)",
                    "genes": r["genes"]}
        go.write(rows_u, OUT / f"go_enrichment_{key}.uncollapsed.tsv")
        go.write(rows_c, OUT / f"go_enrichment_{key}.collapsed.tsv")

        # cluster diagnostics: clusters contributing >=3 study genes
        per_rep = defaultdict(list)
        for g in study:
            per_rep[rep_of.get(g, g)].append(g)
        for rep, gs in per_rep.items():
            if len(gs) >= 3:
                a = fn.get(rep, {})
                mem = members[rep]
                diag.append({
                    "set": key, "chrom": a.get("chrom", ""),
                    "start": min(int(fn[m]["start"]) for m in mem if m in fn),
                    "end": max(int(fn[m]["end"]) for m in mem if m in fn),
                    "family": family_key(a.get("product", "")),
                    "n_study_genes": len(gs), "n_cluster_genes": len(mem),
                    "example_product": a.get("product", "")[:60],
                    "genes": ",".join(sorted(g.replace("gene-", "") for g in gs)[:12]),
                })
    diag.sort(key=lambda d: (d["set"], -d["n_study_genes"]))
    with open(OUT / "cluster_diagnostics.tsv", "w") as o:
        dc = ["set", "chrom", "start", "end", "family", "n_study_genes", "n_cluster_genes",
              "example_product", "genes"]
        o.write("\t".join(dc) + "\n")
        for d in diag:
            o.write("\t".join(str(d[c]) for c in dc) + "\n")

    pc = ["set", "group", "go_id", "name", "fold_collapsed", "k_collapsed", "fdr_collapsed",
          "fdr_uncollapsed", "survives_collapsed_fdr0.25", "genes"]
    pheno_rows = sorted(pheno.values(), key=lambda x: (x["set"], x["group"],
                                                        x["fdr_collapsed"] or "9"))
    with open(OUT / "phenotype_survival_refined.tsv", "w") as o:
        o.write("\t".join(pc) + "\n")
        for r in pheno_rows:
            o.write("\t".join(str(r[c]) for c in pc) + "\n")
    surv = [r for r in pheno_rows if r["survives_collapsed_fdr0.25"] == "Y"]
    print(f"  phenotype terms FDR<0.25 after collapsing: {len(surv)} "
          f"(calcium {sum(1 for r in surv if r['group']=='calcium')}, "
          f"lipid {sum(1 for r in surv if r['group']=='lipid')})")

    # ---- Phase-18 headline fate ----
    fate_ids = list(PHASE18_CONVERGENT)
    for r in read_tsv(ANN / "pav_gene_assignments.tsv"):
        if r["gene_id"] and PHASE18_LIPID.search(r.get("product", "")) and r.get("exon_overlap") == "1":
            fate_ids.append(canon(r["gene_id"]))
    for r in read_tsv(ANN / "dmr_gene_assignments.tsv"):
        if r["gene_id"] and PHASE18_LIPID.search(r.get("product", "")):
            fate_ids.append(canon(r["gene_id"]))
    seen, fate = set(), []
    byid = {r["gene_id"]: r for r in rows}
    for g in fate_ids:
        if g in seen:
            continue
        seen.add(g)
        a = fn.get(g, {})
        r = byid.get(g)
        fate.append({
            "gene_id": g, "symbol": a.get("symbol", ""), "product": a.get("product", "")[:60],
            "phase18_role": "convergent" if g in PHASE18_CONVERGENT else "lipid_axis",
            "phase18_rank_score": p18.get(g, {}).get("rank_score", ""),
            "refined_tier": (r["tier"] or "none") if r else "none",
            "n_lines_refined": r["n_lines_refined"] if r else 0,
            "lines_refined": r["lines_refined"] if r else "",
            "native_support": "Y" if r and (r["reciprocal_pav"] or r["ecotype_only_corroborated"]
                                            or r["cnv"]) else "N",
        })
    fc = list(fate[0].keys())
    with open(OUT / "phase18_candidate_fate.tsv", "w") as o:
        o.write("\t".join(fc) + "\n")
        for r in fate:
            o.write("\t".join(str(r[c]) for c in fc) + "\n")

    # ---- README ----
    def n(set_):
        return f"{len(set_):,}"
    with open(OUT / "README.md", "w") as o:
        o.write("# Refined two-genome evidence model (23.7)\n\n")
        o.write("Re-scoring of `../gene_evidence_matrix.tsv` with independent evidence lines, a "
                "corroboration tier for ecotype-only genes, and tandem-cluster-collapsed GO tests. "
                "Built by `code/23.7-refine-evidence-model.py`; rationale in the docstring and in "
                "`code/23-next-phase-research-plan.md`.\n\n")
        o.write("## Evidence lines (genes per line)\n\n")
        o.write("| line | genes | basis |\n|---|---|---|\n")
        o.write(f"| reciprocal_pav | {n(recip)} | native, bidirectional read coverage (Phase 2) |\n")
        o.write(f"| ecotype_only_corroborated | {n(eco_corr)} | lifted into one assembly only AND "
                f"in the other's unmapped list (lean {sum(1 for g in eco_corr if eco_side[g]=='lean'):,} / "
                f"siscowet {sum(1 for g in eco_corr if eco_side[g]=='siscowet'):,}) |\n")
        o.write(f"| cnv | {n(cnv)} | Liftoff copy-number divergence |\n")
        o.write(f"| dmr | {n(dmr)} | reference-based DMR within 5 kb |\n")
        o.write(f"| ref_divergence | {n(ref_div)} | SyRI SV ({n(sv)}) ∪ reference PAV ({n(ref_pav)}) ∪ "
                f"Liftoff-only ecotype_only ({n(eco_lift)}) |\n\n")
        o.write(f"SV line source: {sv_src}.\n\n")
        o.write("## Candidate tiers\n\n")
        o.write("| tier | definition | genes | of which no caution flag |\n|---|---|---|---|\n")
        for t, d in (("A", "reciprocal_pav + ≥1 other line"),
                     ("B", "corroborated ecotype_only or cnv + ≥1 other line, no reciprocal PAV"),
                     ("C", "dmr + ref_divergence only (reference-only association)")):
            tr = [r for r in rows if r["tier"] == t]
            o.write(f"| {t} | {d} | {len(tr):,} | {sum(1 for r in tr if not r['caution']):,} |\n")
        o.write(f"\nCompare 23.6: 2,895 \"two-genome candidates\" at ≥2 of 6 lines; refined tier A+B = "
                f"{len(cand):,}, tier A = {len(tierA):,}.\n\n")
        o.write("### Tier A (native-anchored), protein-coding, no caution flag\n\n")
        o.write("| gene | symbol | product | present in | lines |\n|---|---|---|---|---|\n")
        for r in [x for x in tierA_clean if x["biotype"] == "protein_coding"][:40]:
            o.write(f"| {r['gene_id']} | {r['symbol']} | {r['product'][:50]} | "
                    f"{r['reciprocal_present_in']} | {r['lines_refined']} |\n")
        o.write("\n## GO after tandem-cluster collapsing\n\n")
        o.write(f"Genome-wide, {n_multi:,} clusters of ≥2 same-family genes within {CLUSTER_GAP//1000} kb "
                f"were collapsed to one representative each (background {len(universe_u):,} genes → "
                f"{len(universe_c):,} clusters). Lipid / calcium terms at FDR<0.25 in any set, "
                f"collapsed vs uncollapsed: see `phenotype_survival_refined.tsv`. "
                f"Surviving after collapsing: {len(surv)} terms "
                f"(calcium {sum(1 for r in surv if r['group']=='calcium')}, "
                f"lipid {sum(1 for r in surv if r['group']=='lipid')}).\n\n")
        o.write("| set | group | term | k (clusters) | FDR collapsed | FDR uncollapsed |\n|---|---|---|---|---|---|\n")
        for r in pheno_rows:
            o.write(f"| {r['set']} | {r['group']} | {r['name']} | {r['k_collapsed']} | "
                    f"{r['fdr_collapsed'] or '—'} | {r['fdr_uncollapsed'] or '—'} |\n")
        o.write("\n### Largest collapsed clusters per study set\n\n")
        o.write("| set | locus | family | study genes | cluster genes | example |\n|---|---|---|---|---|---|\n")
        for d in diag[:25]:
            o.write(f"| {d['set']} | {d['chrom']}:{d['start']:,}-{d['end']:,} | {d['family']} | "
                    f"{d['n_study_genes']} | {d['n_cluster_genes']} | {d['example_product']} |\n")
        o.write("\n## Fate of the Phase-18 headline candidates\n\n")
        o.write("| gene | symbol | product | Phase-18 role | refined tier | lines | native support |\n"
                "|---|---|---|---|---|---|---|\n")
        for r in fate:
            o.write(f"| {r['gene_id']} | {r['symbol']} | {r['product']} | {r['phase18_role']} | "
                    f"{r['refined_tier']} | {r['lines_refined'] or '—'} | {r['native_support']} |\n")
        o.write("\n## Reading guide\n\n")
        o.write("- Only tier A carries read-level, bidirectional evidence on native coordinates. "
                "Tier B rests on two Liftoff annotations agreeing; tier C is the Phase-18 style "
                "single-reference association and is retained for continuity only.\n")
        o.write("- `ref_divergence` is one line no matter how many of SV / reference PAV / "
                "Liftoff-only presence a gene carries: they share a cause (distance from the "
                "lean reference).\n")
        o.write("- `caution` flags noncoding biotypes, repeat / tandem families, and clusters of ≥5 "
                "genes; rank on the flag-free subset first.\n")
        o.write("- Phase 3 (native methylation) will add the first native epigenetic line and can "
                "promote tier-B genes; Phase 4 (BRAKER3) addresses genes Liftoff cannot see at all.\n")
    print(f"\nDone. Outputs -> {OUT}")


if __name__ == "__main__":
    main()
