#!/usr/bin/env python3
"""
18.1-assign-features-to-genes.py  (Steps 1-2 of code/18-diff-annotation-phenotype-plan.md)

Step 1 - assign DMRs and DMCs to genes with positional context (promoter / exon /
         intron / flanking), strand-aware.
Step 2 - assign differential PAVs to genes:
           * stringent siscowet-specific deletions -> per-overlap classification
           * lenient lean/siscowet specific del+ins -> per-gene burden (count + bp)

Pure stdlib (no bedtools). Joins straight onto Step 0's gene_function_table.tsv so every
output row already carries symbol / product / biotype / GO count.

Coordinates: everything treated numerically with half-open overlap (qs < end & qe > start).
Gene-model coords come from the GFF (1-based); BED features are 0-based. The <=1 bp slack
this introduces is immaterial at the promoter/gene scale used here.

Windows:  PROMOTER = TSS +/- 2 kb (strand-aware);  NEAR = gene body +/- 5 kb (flanking).
Primary location_class priority:  promoter > exon > intron > upstream<=5kb > downstream<=5kb.
"""

import csv
import gzip
import re
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
ANN = BASE / "analyses" / "18-annotation"
DM = BASE / "analyses" / "14-diff-meth"
PV = BASE / "analyses" / "15-diff-pav"
GFF = ANN / "raw" / "GCF_016432855.1_SaNama_1.0_genomic.gff.gz"
GENE_TBL = ANN / "gene_function_table.tsv"

PROMOTER = 2000
NEAR = 5000

# ---------------------------------------------------------------- load genes
def load_genes():
    genes = {}                       # gene_id -> dict
    by_chrom = defaultdict(list)     # chrom -> list of gene dicts (sorted by start)
    with open(GENE_TBL) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                s, e = int(row["start"]), int(row["end"])
            except ValueError:
                continue
            g = {
                "gene_id": row["gene_id"], "symbol": row["symbol"],
                "product": row["product"], "biotype": row["biotype"],
                "n_go": row["n_go"], "chrom": row["chrom"],
                "start": s, "end": e, "strand": row["strand"],
            }
            genes[g["gene_id"]] = g
            by_chrom[g["chrom"]].append(g)
    for c in by_chrom:
        by_chrom[c].sort(key=lambda g: g["start"])
    # parallel arrays of starts for bisect, plus running max-end for back-scan bound
    starts = {c: [g["start"] for g in lst] for c, lst in by_chrom.items()}
    return genes, by_chrom, starts


# --------------------------------------------------- exon index from the GFF
def load_exons():
    exons = defaultdict(list)        # gene_id -> list of (start,end)
    with gzip.open(GFF, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "exon":
                continue
            m = re.search(r"(?:^|;)gene=([^;]+)", c[8])
            if not m:
                continue
            gid = "gene-" + m.group(1)
            exons[gid].append((int(c[3]), int(c[4])))
    return exons


# --------------------------------------------------------------- classify
def candidate_genes(chrom, qs, qe, by_chrom, starts):
    """Genes whose body+/-NEAR could touch [qs,qe). Back-scan handles long introns."""
    lst = by_chrom.get(chrom)
    if not lst:
        return []
    st = starts[chrom]
    hi = bisect_right(st, qe + NEAR)           # genes starting after this can't reach
    out = []
    # scan backward; stop once gene starts are far enough left that even the
    # longest plausible gene can't reach qs-NEAR. Genes are large, so bound generously.
    limit_left = qs - NEAR - 3_000_000          # > longest lake trout gene
    for i in range(hi - 1, -1, -1):
        g = lst[i]
        if g["start"] < limit_left:
            break
        if g["end"] + NEAR >= qs:               # body+/-NEAR overlaps query
            out.append(g)
    return out


def classify(g, qs, qe, exons):
    body = qs < g["end"] and qe > g["start"]
    tss = g["start"] if g["strand"] == "+" else g["end"]
    p0, p1 = tss - PROMOTER, tss + PROMOTER
    in_prom = qs < p1 and qe > p0
    exon_ov = False
    if body:
        for (es, ee) in exons.get(g["gene_id"], ()):
            if qs < ee and qe > es:
                exon_ov = True
                break
    if body:
        dist = 0
    elif qe <= g["start"]:
        dist = g["start"] - qe
    else:
        dist = qs - g["end"]
    # primary class
    if in_prom:
        cls = "promoter"
    elif exon_ov:
        cls = "exon"
    elif body:
        cls = "intron"
    else:
        upstream = (g["strand"] == "+" and qe <= g["start"]) or \
                   (g["strand"] == "-" and qs >= g["end"])
        cls = "upstream<=5kb" if upstream else "downstream<=5kb"
    return cls, int(in_prom), int(exon_ov), int(body), dist


# --------------------------------------------------------------- Step 1: DMR
def assign_dmrs(genes, by_chrom, starts, exons):
    out = ANN / "dmr_gene_assignments.tsv"
    hits = 0
    dmrs_with_gene = set()
    with open(DM / "dmrs.csv") as fh, open(out, "w") as o:
        o.write("\t".join([
            "dmr_id", "chrom", "dmr_start", "dmr_end", "direction", "meth_diff",
            "n_cpgs", "gene_id", "symbol", "product", "biotype", "location_class",
            "in_promoter", "exon_overlap", "gene_body_overlap", "distance_bp", "n_go",
        ]) + "\n")
        for row in csv.DictReader(fh):
            chrom = row["chrom"]
            qs, qe = int(row["start"]), int(row["end"])
            for g in candidate_genes(chrom, qs, qe, by_chrom, starts):
                cls, ip, ex, bd, dist = classify(g, qs, qe, exons)
                hits += 1
                dmrs_with_gene.add(row["dmr_id"])
                o.write("\t".join(str(x) for x in [
                    row["dmr_id"], chrom, qs, qe, row["direction"],
                    f"{float(row['mean_diff']):.2f}", row["n_cpgs"],
                    g["gene_id"], g["symbol"], g["product"], g["biotype"],
                    cls, ip, ex, bd, dist, g["n_go"],
                ]) + "\n")
    return out, hits, len(dmrs_with_gene)


# --------------------------------------------------------------- Step 1: DMC
def assign_dmcs(genes, by_chrom, starts, exons):
    out = ANN / "dmc_gene_assignments.tsv"
    hits = 0
    dmcs_with_gene = set()
    with open(DM / "significant_dmcs.csv") as fh, open(out, "w") as o:
        o.write("\t".join([
            "site_id", "chrom", "pos", "meth_diff", "direction", "gene_id", "symbol",
            "product", "biotype", "location_class", "in_promoter", "exon_overlap",
            "distance_bp",
        ]) + "\n")
        for row in csv.DictReader(fh):
            chrom = row["chrom"]
            pos = int(row["pos"])
            qs, qe = pos, pos + 1
            for g in candidate_genes(chrom, qs, qe, by_chrom, starts):
                cls, ip, ex, bd, dist = classify(g, qs, qe, exons)
                hits += 1
                dmcs_with_gene.add(row["site_id"])
                o.write("\t".join(str(x) for x in [
                    row["site_id"], chrom, pos, f"{float(row['meth_diff']):.2f}",
                    row["direction"], g["gene_id"], g["symbol"], g["product"],
                    g["biotype"], cls, ip, ex, dist,
                ]) + "\n")
    return out, hits, len(dmcs_with_gene)


# ------------------------------------------------- Step 2: stringent PAV
def assign_stringent_pav(genes, by_chrom, starts, exons):
    out = ANN / "pav_gene_assignments.tsv"
    src = PV / "stringent.siscowet_specific.deletions.bed"
    hits = 0
    pav_with_gene = 0
    total = 0
    with open(src) as fh, open(out, "w") as o:
        o.write("\t".join([
            "pav_id", "chrom", "start", "end", "size_bp", "ecotype_specific",
            "variant_type", "gene_id", "symbol", "product", "biotype",
            "location_class", "exon_overlap", "gene_body_overlap", "distance_bp",
            "n_go", "confidence",
        ]) + "\n")
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 3:
                continue
            total += 1
            chrom, qs, qe = c[0], int(c[1]), int(c[2])
            pav_id = f"{chrom}:{qs}-{qe}"
            found = False
            for g in candidate_genes(chrom, qs, qe, by_chrom, starts):
                cls, ip, ex, bd, dist = classify(g, qs, qe, exons)
                hits += 1
                found = True
                o.write("\t".join(str(x) for x in [
                    pav_id, chrom, qs, qe, qe - qs, "siscowet", "deletion",
                    g["gene_id"], g["symbol"], g["product"], g["biotype"],
                    cls, ex, bd, dist, g["n_go"], "stringent",
                ]) + "\n")
            if found:
                pav_with_gene += 1
    return out, hits, pav_with_gene, total


# ------------------------------------------- Step 2: lenient burden (sweep)
def overlap_join(bed_path, by_chrom):
    """Sweep-line: per gene_id accumulate (count, overlap_bp clipped to gene body)."""
    # bucket PAV intervals by chrom, sorted by start
    pav = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            c = line.split("\t")
            if len(c) < 3:
                continue
            try:
                pav[c[0]].append((int(c[1]), int(c[2])))
            except ValueError:
                continue
    cnt = defaultdict(int)
    bp = defaultdict(int)
    for chrom, plist in pav.items():
        glist = by_chrom.get(chrom)
        if not glist:
            continue
        plist.sort()
        # active genes sorted by end; two-pointer
        gi = 0
        ng = len(glist)
        import heapq
        active = []  # heap of (end, idx)
        for (ps, pe) in plist:
            while gi < ng and glist[gi]["start"] < pe:
                g = glist[gi]
                heapq.heappush(active, (g["end"], gi))
                gi += 1
            while active and active[0][0] <= ps:
                heapq.heappop(active)
            for (gend, idx) in active:
                g = glist[idx]
                ov = min(pe, g["end"]) - max(ps, g["start"])
                if ov > 0:
                    cnt[g["gene_id"]] += 1
                    bp[g["gene_id"]] += ov
    return cnt, bp


def lenient_burden(genes, by_chrom):
    sets = {
        "lean_del": PV / "lean_specific.deletions.bed",
        "sisco_del": PV / "siscowet_specific.deletions.bed",
        "lean_ins": PV / "lean_specific.insertions.bed",
        "sisco_ins": PV / "siscowet_specific.insertions.bed",
    }
    res = {}
    for key, path in sets.items():
        print(f"  burden: {key} ...", file=sys.stderr)
        res[key] = overlap_join(path, by_chrom)
    out = ANN / "pav_gene_burden.tsv"
    cols = ["gene_id", "symbol", "product", "biotype",
            "lean_del_n", "lean_del_bp", "sisco_del_n", "sisco_del_bp",
            "lean_ins_n", "sisco_ins_n", "ecotype_skew_note"]
    n_rows = 0
    with open(out, "w") as o:
        o.write("\t".join(cols) + "\n")
        for gid, g in genes.items():
            ld_n = res["lean_del"][0].get(gid, 0)
            ld_bp = res["lean_del"][1].get(gid, 0)
            sd_n = res["sisco_del"][0].get(gid, 0)
            sd_bp = res["sisco_del"][1].get(gid, 0)
            li_n = res["lean_ins"][0].get(gid, 0)
            si_n = res["sisco_ins"][0].get(gid, 0)
            if (ld_n + sd_n + li_n + si_n) == 0:
                continue
            n_rows += 1
            o.write("\t".join(str(x) for x in [
                gid, g["symbol"], g["product"], g["biotype"],
                ld_n, ld_bp, sd_n, sd_bp, li_n, si_n,
                "lean-reference: siscowet-specific & deletion counts inflated by divergence",
            ]) + "\n")
    return out, n_rows


# --------------------------------------------------------------- main
def main():
    print("loading genes + exons ...", file=sys.stderr)
    genes, by_chrom, starts = load_genes()
    exons = load_exons()
    print(f"  {len(genes)} genes, exons for {len(exons)} genes", file=sys.stderr)

    print("Step 1: DMRs -> genes ...", file=sys.stderr)
    f, h, n = assign_dmrs(genes, by_chrom, starts, exons)
    print(f"  {h} DMR-gene pairs; {n}/302 DMRs have >=1 gene within 5kb -> {f.name}",
          file=sys.stderr)

    print("Step 1: DMCs -> genes ...", file=sys.stderr)
    f, h, n = assign_dmcs(genes, by_chrom, starts, exons)
    print(f"  {h} DMC-gene pairs; {n} DMCs have >=1 gene within 5kb -> {f.name}",
          file=sys.stderr)

    print("Step 2: stringent siscowet-specific deletions -> genes ...", file=sys.stderr)
    f, h, pg, tot = assign_stringent_pav(genes, by_chrom, starts, exons)
    print(f"  {h} PAV-gene pairs; {pg}/{tot} stringent deletions hit a gene+/-5kb -> {f.name}",
          file=sys.stderr)

    print("Step 2: lenient PAV per-gene burden (sweep) ...", file=sys.stderr)
    f, n = lenient_burden(genes, by_chrom)
    print(f"  {n} genes with >=1 lenient PAV overlap -> {f.name}", file=sys.stderr)
    print("done.", file=sys.stderr)


if __name__ == "__main__":
    main()
