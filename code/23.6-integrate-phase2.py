#!/usr/bin/env python3
"""
23.6-integrate-phase2.py  (Phase 2, Step 3 of code/23-next-phase-research-plan.md)

Integrate the two-genome evidence produced in Phase 2 and re-test the ecotype-divergence
hypotheses on native / bidirectional data instead of single-reference association.

Gene-level evidence collected per reference gene ID (gene-LOC######; the shared namespace across
Liftoff, PAV, DMR and candidate tables):

  reciprocal_pav  gene overlaps a native reciprocal-PAV specific-present region
                  (analyses/23-reciprocal-pav/{lean,siscowet}_specific.genes.tsv)   [Phase 2, native]
  ecotype_only    gene present in only one ecotype's Liftoff annotation
                  (analyses/23-gene-sets/{lean,siscowet}_only.tsv)                   [Phase 2]
  cnv             copy-number differs between ecotypes
                  (analyses/23-gene-sets/cnv.tsv)                                    [Phase 2]
  sv              gene overlaps a SyRI structural SV vs reference
                  (analyses/23-genome-sv/syri/{lean,sisco}_vs_ref.syri.out x ref gene bed) [Phase 1]
  dmr             gene assigned a DMR (analyses/18-annotation/dmr_gene_assignments.tsv)     [prior]
  ref_pav         gene assigned a reference PAV (analyses/18-annotation/pav_gene_assignments.tsv) [prior]

`n_lines` = number of independent evidence categories. The two-genome deliverable is the set of
genes supported by >= 2 lines, at least one of which is a NEW Phase-2 line (reciprocal_pav /
ecotype_only / cnv) -- i.e. no longer single-reference association.

We then re-run the GO over-representation analysis of code/18.3-go-enrichment.py (reusing its OBO /
DAG / hypergeometric machinery, imported directly) on the two-genome gene sets, and report whether
the lipid-metabolism / calcium-transport phenotype terms survive on native coordinates.

Outputs (analyses/23-integration/):
  gene_evidence_matrix.tsv     every gene with >=1 line, all flags + n_lines + annotation
  two_genome_candidates.tsv    genes with >=2 lines incl. a Phase-2 line, ranked
  go_enrichment_{reciprocal_pav,lean_only,siscowet_only,cnv}.tsv
  phenotype_survival.tsv       lipid/calcium terms and whether they remain enriched
  README.md

Author: Generated for project-lake-trout
Date: 2026-07-04
"""

import csv
import importlib.util
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
RECIP = BASE / "analyses" / "23-reciprocal-pav"
GENESETS = BASE / "analyses" / "23-gene-sets"
SV = BASE / "analyses" / "23-genome-sv" / "syri"
ANN = BASE / "analyses" / "18-annotation"
OUT = BASE / "analyses" / "23-integration"
REF_GENE_BED = BASE / "data" / "20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed"
BEDTOOLS = "/usr/bin/bedtools"

# Structural SV types (SyRI col 11) that should intersect genes for the `sv` line.
SV_STRUCTURAL = {"INV", "TRANS", "INVTR", "DUP", "INVDP", "CPG", "CPL", "HDR", "TDM"}

PHENO_LIPID = re.compile(
    r"lipid|fatty.?acid|triglyceride|triacylglycerol|sterol|cholesterol|lipo|"
    r"buoyan|swim.?bladder|gas bladder|adipos", re.I)
PHENO_CALCIUM = re.compile(r"calcium|\bca2\b|calcium.ion|voltage.gated calcium", re.I)


# ------------------------------------------------------------ import 18.3 as a module
def load_go_module():
    path = BASE / "code" / "18.3-go-enrichment.py"
    spec = importlib.util.spec_from_file_location("go_enrichment_183", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------ evidence readers
def col_ids(path, col="gene_id", sep="\t"):
    s = set()
    if not Path(path).exists():
        print(f"  WARNING missing {path}")
        return s
    with open(path) as f:
        r = csv.DictReader(f, delimiter=sep)
        for row in r:
            if row.get(col):
                s.add(row[col])
    return s


def bed_col4_ids(path):
    s = set()
    if not Path(path).exists():
        print(f"  WARNING missing {path}")
        return s
    with open(path) as f:
        for line in f:
            if line.startswith(("chrom\t", "#", "track")):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                s.add(p[3])
    return s


def sv_affected_genes():
    """Reference genes overlapping a structural SyRI SV in either reference comparison."""
    OUT.mkdir(parents=True, exist_ok=True)
    sv_bed = OUT / "_ref_structural_sv.bed"
    with open(sv_bed, "w") as o:
        for tag in ("lean_vs_ref", "sisco_vs_ref"):
            syri = SV / f"{tag}.syri.out"
            if not syri.exists():
                print(f"  WARNING missing {syri}")
                continue
            with open(syri) as f:
                for line in f:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 11 or p[10] not in SV_STRUCTURAL:
                        continue
                    try:
                        a, b = int(float(p[1])), int(float(p[2]))
                    except ValueError:
                        continue
                    if b < a:
                        a, b = b, a
                    if b == a:
                        b = a + 1
                    o.write(f"{p[0]}\t{a}\t{b}\t{p[10]}\n")
    sh(f"sort -k1,1 -k2,2n '{sv_bed}' -o '{sv_bed}'")
    hit = OUT / "_sv_genes.bed"
    sh(f"'{BEDTOOLS}' intersect -a '{REF_GENE_BED}' -b '{sv_bed}' -u > '{hit}'")
    genes = bed_col4_ids(hit)
    sv_bed.unlink(missing_ok=True)
    hit.unlink(missing_ok=True)
    return genes


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True)


# ------------------------------------------------------------ main
def main():
    OUT.mkdir(parents=True, exist_ok=True)

    # ---- collect evidence sets (all in gene-LOC###### namespace) ----
    recip = bed_col4_ids(RECIP / "lean_specific.genes.tsv") | \
        bed_col4_ids(RECIP / "siscowet_specific.genes.tsv")
    lean_only = col_ids(GENESETS / "lean_only.tsv")
    sisco_only = col_ids(GENESETS / "siscowet_only.tsv")
    ecotype_only = lean_only | sisco_only
    cnv = col_ids(GENESETS / "cnv.tsv")
    dmr = col_ids(ANN / "dmr_gene_assignments.tsv")
    ref_pav = col_ids(ANN / "pav_gene_assignments.tsv")
    sv_genes = sv_affected_genes()

    lines = {
        "reciprocal_pav": recip,
        "ecotype_only": ecotype_only,
        "cnv": cnv,
        "sv": sv_genes,
        "dmr": dmr,
        "ref_pav": ref_pav,
    }
    phase2_lines = {"reciprocal_pav", "ecotype_only", "cnv"}
    for k, v in lines.items():
        print(f"  {k:15s}: {len(v):,} genes")

    all_genes = set().union(*lines.values())

    # ---- gene annotation from the function table ----
    fn = {}
    with open(ANN / "gene_function_table.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fn[r["gene_id"]] = (r.get("symbol", ""), r.get("product", ""),
                                r.get("biotype", ""), r.get("n_go", "0"))

    # ---- evidence matrix ----
    flag_cols = ["reciprocal_pav", "ecotype_only", "cnv", "sv", "dmr", "ref_pav"]
    rows = []
    for g in all_genes:
        flags = {c: (1 if g in lines[c] else 0) for c in flag_cols}
        n_lines = sum(flags.values())
        has_p2 = any(flags[c] for c in phase2_lines)
        s, p, bt, ng = fn.get(g, ("", "", "", "0"))
        rows.append({"gene_id": g, "symbol": s, "product": p, "biotype": bt,
                     "n_go": ng, "n_lines": n_lines, "has_phase2_line": int(has_p2),
                     **flags})
    rows.sort(key=lambda r: (-r["n_lines"], -r["has_phase2_line"], r["gene_id"]))

    mat_cols = (["gene_id", "symbol", "product", "biotype", "n_go", "n_lines",
                 "has_phase2_line"] + flag_cols)
    with open(OUT / "gene_evidence_matrix.tsv", "w") as o:
        o.write("\t".join(mat_cols) + "\n")
        for r in rows:
            o.write("\t".join(str(r[c]) for c in mat_cols) + "\n")

    cand = [r for r in rows if r["n_lines"] >= 2 and r["has_phase2_line"]]
    with open(OUT / "two_genome_candidates.tsv", "w") as o:
        o.write("\t".join(mat_cols) + "\n")
        for r in cand:
            o.write("\t".join(str(r[c]) for c in mat_cols) + "\n")
    print(f"\n  two-genome candidates (>=2 lines incl. Phase-2): {len(cand):,}")

    # ---- GO re-run on the two-genome sets (reuse 18.3 machinery) ----
    go = load_go_module()
    print("\n  loading GO DAG (18.3 machinery) ...")
    name, ns, parents, alt = go.parse_obo(go.OBO)
    _, anc_fn = go.build_ancestors(parents)
    gene_go = go.load_gene_go(name, alt, anc_fn)
    universe = set(gene_go)
    term2genes = defaultdict(set)
    for gid, terms in gene_go.items():
        for t in terms:
            term2genes[t].add(gid)
    print(f"  GO background: {len(universe):,} genes")

    study_sets = {
        "reciprocal_pav": recip,
        "lean_only": lean_only,
        "siscowet_only": sisco_only,
        "cnv": cnv,
    }
    pheno_rows = []
    for key, study in study_sets.items():
        ora_rows, n, N = go.ora(study, term2genes, universe, name, ns)
        # Read numeric fdr BEFORE go.write() -- write() formats fdr/p_value in place to strings.
        sig = sum(1 for r in ora_rows if r["fdr"] < 0.1)
        for r in ora_rows:
            nm = r["name"]
            grp = ("lipid" if PHENO_LIPID.search(nm) else
                   "calcium" if PHENO_CALCIUM.search(nm) else None)
            if grp and r["fdr"] < 0.25:
                pheno_rows.append({
                    "set": key, "group": grp, "go_id": r["go_id"], "name": nm,
                    "fold_enrichment": r["fold_enrichment"], "study_k": r["study_k"],
                    "fdr": f"{r['fdr']:.3e}", "genes": r["genes"],
                })
        go.write(ora_rows, OUT / f"go_enrichment_{key}.tsv")
        print(f"  [{key}] study={n} (of {len(study)}), bg={N}, "
              f"terms={len(ora_rows)}, FDR<0.1={sig}")

    with open(OUT / "phenotype_survival.tsv", "w") as o:
        cols = ["set", "group", "go_id", "name", "fold_enrichment", "study_k", "fdr", "genes"]
        o.write("\t".join(cols) + "\n")
        for r in sorted(pheno_rows, key=lambda x: (x["set"], x["group"])):
            o.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"  phenotype-flagged enriched terms (FDR<0.25): {len(pheno_rows)}")

    # ---- README ----
    top = cand[:20]
    with open(OUT / "README.md", "w") as o:
        o.write("# Phase 2 integration — two-genome ecotype-divergence evidence\n\n")
        o.write("Gene-level integration of Phase-2 (native reciprocal PAV, ecotype-only gene "
                "presence, CNV) with Phase-1 SVs and the prior reference DMR/PAV tables. "
                "See `code/23.6-integrate-phase2.py`.\n\n")
        o.write("> **Superseded for interpretation by `refined/`** (code/23.7-refine-evidence-model.py): the ecotype_only / ref_pav / sv lines below are not independent, ecotype_only is Liftoff-sensitivity-inflated, and the GO tests are tandem-cluster driven. This file is kept as the Step-3 record.\n\n")
        o.write("## Evidence lines (genes per line)\n\n")
        for k in flag_cols:
            o.write(f"- **{k}**: {len(lines[k]):,}\n")
        o.write(f"\n## Two-genome candidates (>=2 lines incl. a Phase-2 line): "
                f"{len(cand):,}\n\n")
        o.write("A candidate is no longer a single-reference association: at least one of "
                "reciprocal_pav / ecotype_only / cnv (native, bidirectional) plus >=1 other line.\n\n")
        o.write("Top 20 by n_lines:\n\n")
        o.write("| gene | symbol | product | n_lines | lines |\n|---|---|---|---|---|\n")
        for r in top:
            act = ",".join(c for c in flag_cols if r[c])
            o.write(f"| {r['gene_id']} | {r['symbol']} | {r['product'][:40]} | "
                    f"{r['n_lines']} | {act} |\n")
        o.write("\n## Phenotype survival\n\n")
        o.write("`phenotype_survival.tsv` lists lipid-metabolism and calcium-transport GO terms "
                "that remain enriched (FDR<0.25) in the two-genome sets. "
                f"Total flagged: {len(pheno_rows)}.\n")

    print(f"\nDone. Outputs -> {OUT}")


if __name__ == "__main__":
    main()
