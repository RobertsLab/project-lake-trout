#!/usr/bin/env python3
"""
18-build-gene-function-table.py  (Step 0 of code/18-diff-annotation-phenotype-plan.md)

Build a functional annotation table for the SaNama_1.0 (GCF_016432855.1) gene set,
keyed on the same `gene-XXX` IDs used in
  data/20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed
so DMR/PAV gene assignments (Steps 1-2) can join straight onto function.

Inputs (downloaded from NCBI RefSeq into analyses/18-annotation/raw/):
  - *_genomic.gff.gz          gene models + transcript `product=` names  (join backbone)
  - *_gene_ontology.gaf.gz    GO annotations by NCBI GeneID
  - *_feature_table.txt.gz    (not required; GFF + GAF are sufficient)

Output:
  - analyses/18-annotation/gene_function_table.tsv

Columns: gene_id, symbol, geneid, biotype, chrom, start, end, strand,
         is_named (symbol is a real gene symbol, not a LOC id),
         product, n_go, go_bp, go_mf, go_cc, go_terms
"""

import gzip
import re
import sys
from pathlib import Path
from collections import defaultdict

BASE = Path(__file__).resolve().parent.parent
RAW = BASE / "analyses" / "18-annotation" / "raw"
OUT = BASE / "analyses" / "18-annotation" / "gene_function_table.tsv"
GENES_BED = BASE / "data" / "20220818-snam-GCF_016432855.1_SaNama_1.0_genes.bed"

GFF = RAW / "GCF_016432855.1_SaNama_1.0_genomic.gff.gz"
GAF = RAW / "GCF_016432855.1_SaNama_1.0_gene_ontology.gaf.gz"

# transcript-level feature types that carry a `product=` worth surfacing
TRANSCRIPT_TYPES = {
    "mRNA", "lnc_RNA", "transcript", "tRNA", "rRNA", "snRNA", "snoRNA",
    "guide_RNA", "antisense_RNA", "primary_transcript", "V_gene_segment",
    "C_gene_segment", "ncRNA",
}


def attrs(field):
    """Parse a GFF column-9 attribute string into a dict (URL-decoded values)."""
    d = {}
    for kv in field.rstrip(";").split(";"):
        if "=" not in kv:
            continue
        k, v = kv.split("=", 1)
        # GFF3 percent-encoding: %2C -> comma, %25 -> %, %3B -> ;, etc.
        v = re.sub(r"%([0-9A-Fa-f]{2})", lambda m: chr(int(m.group(1), 16)), v)
        d[k] = v
    return d


def parse_gff(path):
    """Return (genes dict keyed by gene_id, gene_id->product, rna_id->gene_id)."""
    genes = {}
    product = {}            # gene_id -> representative product string
    product_seen = defaultdict(set)
    rna2gene = {}           # rna ID -> gene_id  (to attribute CDS products if needed)
    opener = gzip.open(path, "rt")
    with opener as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            ftype = c[2]
            if ftype == "gene" or ftype == "pseudogene":
                a = attrs(c[8])
                gid = a.get("ID")
                if not gid:
                    continue
                geneid = ""
                for x in a.get("Dbxref", "").split(","):
                    if x.startswith("GeneID:"):
                        geneid = x.split(":", 1)[1]
                        break
                symbol = a.get("Name", a.get("gene", gid))
                genes[gid] = {
                    "gene_id": gid,
                    "symbol": symbol,
                    "geneid": geneid,
                    "biotype": a.get("gene_biotype", "pseudogene" if ftype == "pseudogene" else ""),
                    "chrom": c[0],
                    "start": c[3],
                    "end": c[4],
                    "strand": c[6],
                }
            elif ftype in TRANSCRIPT_TYPES:
                a = attrs(c[8])
                parent = a.get("Parent", "")
                rid = a.get("ID", "")
                if parent.startswith("gene-"):
                    if rid:
                        rna2gene[rid] = parent
                    prod = a.get("product")
                    if prod:
                        product_seen[parent].add(prod)
    # collapse multiple isoform products: prefer the shortest cleanest, else join
    for gid, prods in product_seen.items():
        prods = sorted(prods)
        if len(prods) == 1:
            product[gid] = prods[0]
        else:
            # drop "%, transcript variant X" style suffixes for the representative
            base = sorted({re.sub(r",?\s*transcript variant.*$", "", p) for p in prods})
            product[gid] = base[0] if len(base) == 1 else "; ".join(prods[:3])
    return genes, product, rna2gene


def parse_gaf(path):
    """Return geneid(str) -> {'P':set,'F':set,'C':set} of GO ids."""
    go = defaultdict(lambda: {"P": set(), "F": set(), "C": set()})
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            geneid = c[1]          # NCBIGene id
            goid = c[4]
            aspect = c[8]          # P / F / C
            if aspect in ("P", "F", "C"):
                go[geneid][aspect].add(goid)
    return go


def main():
    for p in (GFF, GAF, GENES_BED):
        if not p.exists():
            sys.exit(f"missing input: {p}")

    print("parsing GFF ...", file=sys.stderr)
    genes, product, _ = parse_gff(GFF)
    print(f"  {len(genes)} gene/pseudogene records", file=sys.stderr)

    print("parsing GAF ...", file=sys.stderr)
    go = parse_gaf(GAF)
    print(f"  {len(go)} genes with GO", file=sys.stderr)

    loc_re = re.compile(r"^LOC\d+$")
    cols = ["gene_id", "symbol", "geneid", "biotype", "chrom", "start", "end",
            "strand", "is_named", "product", "n_go", "go_bp", "go_mf", "go_cc",
            "go_terms"]

    n_named = n_prod = n_go = 0
    with open(OUT, "w") as out:
        out.write("\t".join(cols) + "\n")
        for gid in sorted(genes):
            g = genes[gid]
            sym = g["symbol"]
            is_named = "0" if loc_re.match(sym) else "1"
            if is_named == "1":
                n_named += 1
            prod = product.get(gid, "")
            if prod:
                n_prod += 1
            gg = go.get(g["geneid"], {"P": set(), "F": set(), "C": set()})
            allgo = sorted(gg["P"] | gg["F"] | gg["C"])
            if allgo:
                n_go += 1
            row = [
                gid, sym, g["geneid"], g["biotype"], g["chrom"], g["start"],
                g["end"], g["strand"], is_named, prod,
                str(len(allgo)), str(len(gg["P"])), str(len(gg["F"])),
                str(len(gg["C"])), "|".join(allgo),
            ]
            out.write("\t".join(row) + "\n")

    # ---- QC against the gene BED join key ----
    bed_ids = set()
    with open(GENES_BED) as fh:
        for line in fh:
            bed_ids.add(line.split("\t")[3])
    tbl_ids = set(genes)
    missing = bed_ids - tbl_ids
    extra = tbl_ids - bed_ids

    print("\n=== QC ===", file=sys.stderr)
    print(f"genes in table              : {len(genes)}", file=sys.stderr)
    print(f"  with real symbol          : {n_named}", file=sys.stderr)
    print(f"  with product name         : {n_prod}", file=sys.stderr)
    print(f"  with >=1 GO term          : {n_go}", file=sys.stderr)
    print(f"gene IDs in BED             : {len(bed_ids)}", file=sys.stderr)
    print(f"  BED IDs missing from table: {len(missing)}", file=sys.stderr)
    print(f"  table IDs not in BED      : {len(extra)}", file=sys.stderr)
    if missing:
        print("  e.g. missing:", list(sorted(missing))[:5], file=sys.stderr)
    if extra:
        print("  e.g. extra:", list(sorted(extra))[:5], file=sys.stderr)
    print(f"\nwrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
