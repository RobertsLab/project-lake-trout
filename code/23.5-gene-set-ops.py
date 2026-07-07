#!/usr/bin/env python3
"""
23.5-gene-set-ops.py  (Phase 2, Step 2 of code/23-next-phase-research-plan.md)

Turn the two native-coordinate Liftoff annotations into gene-set and copy-number evidence:

  * shared / lean-only / siscowet-only gene-ID sets, from the two *.purged.liftoff.genes.bed files.
    Presence = >= 1 lifted copy of that reference gene in that ecotype's genome. Because both
    annotations were lifted from the SAME reference (identical gene-ID namespace), the IDs are
    directly comparable. An ecotype-"only" gene should also appear in the OTHER ecotype's
    *.purged.unmapped_features.txt (it failed to lift there) -- we cross-reference and flag that.
  * a per-gene COPY-NUMBER table from Liftoff's copy annotations (each mapped copy is a separate
    `gene` feature in the gff; total copies = number of gene features sharing a Name=). Genes whose
    copy count differs between ecotypes are the CNV-divergent set.

All output genes are annotated with symbol / product / GO count from the Step-0 gene function table
(analyses/18-annotation/gene_function_table.tsv) so downstream GO/phenotype work (23.6) can join.

Outputs (analyses/23-gene-sets/):
  shared.tsv           gene_id, symbol, product, copies_lean, copies_siscowet, n_go
  lean_only.tsv        gene_id, symbol, product, copies_lean, in_siscowet_unmapped, n_go
  siscowet_only.tsv    gene_id, symbol, product, copies_siscowet, in_lean_unmapped, n_go
  cnv.tsv              gene_id, symbol, product, copies_lean, copies_siscowet, delta, n_go  (CNV-divergent)
  gene_set_summary.csv counts per category

Author: Generated for project-lake-trout
Date: 2026-07-04
"""

import csv
import re
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
ANN = BASE / "analyses" / "18-annotation"
OUT = BASE / "analyses" / "23-gene-sets"

ECO = {
    "lean":     {"bed": SYN / "lean.purged.liftoff.genes.bed",
                 "gff": SYN / "lean.purged.liftoff.gff3",
                 "unmapped": SYN / "lean.purged.unmapped_features.txt"},
    "siscowet": {"bed": SYN / "siscowet.purged.liftoff.genes.bed",
                 "gff": SYN / "siscowet.purged.liftoff.gff3",
                 "unmapped": SYN / "siscowet.purged.unmapped_features.txt"},
}

NAME_RE = re.compile(r"(?:^|;)Name=([^;]+)")


def bed_gene_ids(path):
    ids = set()
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                ids.add(p[3])
    return ids


def gff_copy_numbers(path):
    """gene_id ('gene-'+Name) -> total copy count (number of `gene` features with that Name)."""
    copies = defaultdict(int)
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            m = NAME_RE.search(p[8])
            if not m:
                continue
            copies["gene-" + m.group(1)] += 1
    return copies


def load_unmapped(path):
    return {l.strip() for l in open(path) if l.strip()}


def load_function_table():
    fn = {}
    with open(ANN / "gene_function_table.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fn[r["gene_id"]] = (r.get("symbol", ""), r.get("product", ""), r.get("n_go", "0"))
    return fn


def annot(gid, fn):
    return fn.get(gid, ("", "", "0"))


def main():
    OUT.mkdir(parents=True, exist_ok=True)

    lean_ids = bed_gene_ids(ECO["lean"]["bed"])
    sisco_ids = bed_gene_ids(ECO["siscowet"]["bed"])
    print(f"lean genes    : {len(lean_ids):,}")
    print(f"siscowet genes: {len(sisco_ids):,}")

    cn_lean = gff_copy_numbers(ECO["lean"]["gff"])
    cn_sisco = gff_copy_numbers(ECO["siscowet"]["gff"])

    lean_unmapped = load_unmapped(ECO["lean"]["unmapped"])
    sisco_unmapped = load_unmapped(ECO["siscowet"]["unmapped"])
    fn = load_function_table()

    shared = lean_ids & sisco_ids
    lean_only = lean_ids - sisco_ids
    sisco_only = sisco_ids - lean_ids
    print(f"shared={len(shared):,}  lean_only={len(lean_only):,}  "
          f"siscowet_only={len(sisco_only):,}")

    # shared.tsv
    with open(OUT / "shared.tsv", "w") as o:
        o.write("gene_id\tsymbol\tproduct\tcopies_lean\tcopies_siscowet\tn_go\n")
        for g in sorted(shared):
            s, p, ng = annot(g, fn)
            o.write(f"{g}\t{s}\t{p}\t{cn_lean.get(g,1)}\t{cn_sisco.get(g,1)}\t{ng}\n")

    # lean_only.tsv  (+ confirm it failed to lift into siscowet)
    conf_lo = 0
    with open(OUT / "lean_only.tsv", "w") as o:
        o.write("gene_id\tsymbol\tproduct\tcopies_lean\tin_siscowet_unmapped\tn_go\n")
        for g in sorted(lean_only):
            s, p, ng = annot(g, fn)
            unm = "Y" if g in sisco_unmapped else ""
            conf_lo += unm == "Y"
            o.write(f"{g}\t{s}\t{p}\t{cn_lean.get(g,1)}\t{unm}\t{ng}\n")

    # siscowet_only.tsv
    conf_so = 0
    with open(OUT / "siscowet_only.tsv", "w") as o:
        o.write("gene_id\tsymbol\tproduct\tcopies_siscowet\tin_lean_unmapped\tn_go\n")
        for g in sorted(sisco_only):
            s, p, ng = annot(g, fn)
            unm = "Y" if g in lean_unmapped else ""
            conf_so += unm == "Y"
            o.write(f"{g}\t{s}\t{p}\t{cn_sisco.get(g,1)}\t{unm}\t{ng}\n")

    # cnv.tsv  (shared genes whose copy count differs)
    cnv = []
    for g in shared:
        cl, cs = cn_lean.get(g, 1), cn_sisco.get(g, 1)
        if cl != cs:
            cnv.append((g, cl, cs, cs - cl))
    cnv.sort(key=lambda x: -abs(x[3]))
    with open(OUT / "cnv.tsv", "w") as o:
        o.write("gene_id\tsymbol\tproduct\tcopies_lean\tcopies_siscowet\tdelta_siscowet_minus_lean\tn_go\n")
        for g, cl, cs, d in cnv:
            s, p, ng = annot(g, fn)
            o.write(f"{g}\t{s}\t{p}\t{cl}\t{cs}\t{d}\t{ng}\n")

    # summary
    with open(OUT / "gene_set_summary.csv", "w") as o:
        o.write("category,n_genes,note\n")
        o.write(f"lean_total,{len(lean_ids)},genes lifted into lean\n")
        o.write(f"siscowet_total,{len(sisco_ids)},genes lifted into siscowet\n")
        o.write(f"shared,{len(shared)},lifted into both\n")
        o.write(f"lean_only,{len(lean_only)},in lean not siscowet "
                f"({conf_lo} also in siscowet unmapped_features)\n")
        o.write(f"siscowet_only,{len(sisco_only)},in siscowet not lean "
                f"({conf_so} also in lean unmapped_features)\n")
        o.write(f"cnv_divergent,{len(cnv)},shared genes with differing copy number\n")

    print(f"\nlean_only confirmed in siscowet unmapped_features : {conf_lo}/{len(lean_only)}")
    print(f"siscowet_only confirmed in lean unmapped_features : {conf_so}/{len(sisco_only)}")
    print(f"CNV-divergent genes                               : {len(cnv):,}")
    print(f"\noutputs -> {OUT}")


if __name__ == "__main__":
    main()
