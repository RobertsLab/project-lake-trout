#!/usr/bin/env python3
"""
22 - Build reference-anchored synteny + functional-annotation browser tracks (Phase 1).

Projects the gene-anchored ecotype synteny (code/21.1) and the liftoff functional
annotation (code/20.1) onto the reference SaNama_1.0 coordinate system so they can be
served as ordinary tracks in the existing IGV.js and JBrowse 2 browsers.

Outputs (written to output/22-synteny-browser-tracks/, then staged into both browsers):
  - genes_annotated.gff3        reference genes + symbol/product/biotype/GO + ecotype presence
  - synteny_lean_blocks.bed     reference footprint of each reference<->lean synteny block
  - synteny_siscowet_blocks.bed reference footprint of each reference<->siscowet synteny block

Coordinate note: MCScanX sanitized chromosome names (rfNC0523351) are un-sanitized back to
reference accessions (NC_052335.1) so features line up with the browser sequence.

Usage:  python code/22-synteny-browser-tracks.py
"""

import csv
import hashlib
import os
import re
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYN_DIR = os.path.join(ROOT, "output", "21.1-gene-anchored-synteny")
ANNOT_DIR = os.path.join(ROOT, "output", "20.1-liftoff-annotation")
GENES_BED = os.path.join(ROOT, "genome-browser", "data", "annotations", "genes.bed")
OUT_DIR = os.path.join(ROOT, "output", "22-synteny-browser-tracks")

# Where the finished tracks get copied so both browsers can serve them.
STAGE_DIRS = [
    os.path.join(ROOT, "genome-browser", "data", "synteny"),
    os.path.join(ROOT, "jbrowse", "data", "lake-trout"),
]


def strip_copy_suffix(gene_id):
    """gene-LOC120041454_1 -> gene-LOC120041454 (Liftoff -copies suffix)."""
    return re.sub(r"_\d+$", "", gene_id)


def unsanitize_ref_chr(name):
    """rfNC0523351 -> NC_052335.1 ; rfNW0240584641 -> NW_024058464.1.

    MCScanX needs short alnum chr ids, so 21.1 stripped '_' and '.' and prefixed 'rf'.
    The trailing digit is the accession version; everything before it is the number.
    """
    m = re.match(r"^rf([A-Z]{2})(\d+)(\d)$", name)
    if not m:
        return None
    return f"{m.group(1)}_{m.group(2)}.{m.group(3)}"


def gff3_escape(value):
    """URL-encode GFF3 attribute reserved characters."""
    out = []
    for ch in value:
        if ch in ";=&,\t\n\r%":
            out.append("%{:02X}".format(ord(ch)))
        else:
            out.append(ch)
    return "".join(out)


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# A1. Functionally-annotated reference gene track (GFF3)
# ---------------------------------------------------------------------------

def load_function_table(path):
    """gene_id -> dict(symbol, biotype, product, go_terms) keyed by base reference id."""
    ann = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["chrom"] == "NA":  # no matched reference locus (extra copy)
                continue
            base = strip_copy_suffix(row["gene_id"])
            if base in ann:
                continue
            ann[base] = {
                "symbol": "" if row["symbol"] == "NA" else row["symbol"],
                "biotype": "" if row["biotype"] == "NA" else row["biotype"],
                "product": "" if row["product"] == "NA" else row["product"],
                "go_terms": "" if row["go_terms"] == "NA" else row["go_terms"],
            }
    return ann


def load_ecotype_gene_ids(path):
    """Set of base gene ids present in an ecotype's purged liftoff annotation."""
    ids = set()
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 4:
                ids.add(strip_copy_suffix(cols[3]))
    return ids


def build_annotated_gff3(out_path):
    ann = load_function_table(
        os.path.join(ANNOT_DIR, "lean.liftoff.gene_function_table.tsv")
    )
    lean_ids = load_ecotype_gene_ids(
        os.path.join(SYN_DIR, "lean.purged.liftoff.genes.bed")
    )
    sisco_ids = load_ecotype_gene_ids(
        os.path.join(SYN_DIR, "siscowet.purged.liftoff.genes.bed")
    )

    n_total = n_annot = 0
    with open(GENES_BED) as fh, open(out_path, "w") as out:
        out.write("##gff-version 3\n")
        for line in fh:
            chrom, start, end, name, _score, strand = line.rstrip("\n").split("\t")[:6]
            n_total += 1
            base = strip_copy_suffix(name)
            info = ann.get(base, {})
            if info:
                n_annot += 1

            ecos = []
            if base in lean_ids:
                ecos.append("lean")
            if base in sisco_ids:
                ecos.append("siscowet")

            symbol = info.get("symbol") or base.replace("gene-", "")
            attrs = [f"ID={gff3_escape(name)}", f"Name={gff3_escape(symbol)}"]
            attrs.append(f"gene_id={gff3_escape(base)}")
            if info.get("product"):
                attrs.append(f"product={gff3_escape(info['product'])}")
            if info.get("biotype"):
                attrs.append(f"biotype={gff3_escape(info['biotype'])}")
            if info.get("go_terms"):
                # Ontology_term is the GFF3-reserved key; values are comma-separated GO
                # ids (no reserved chars, so no escaping needed).
                go = ",".join(t for t in re.split(r"[|;,]", info["go_terms"]) if t)
                attrs.append(f"Ontology_term={go}")
            if ecos:
                attrs.append(f"ecotypes={','.join(ecos)}")

            # BED (0-based, half-open) -> GFF3 (1-based, inclusive)
            gstart = int(start) + 1
            out.write(
                "\t".join(
                    [chrom, "liftoff_annot", "gene", str(gstart), end,
                     ".", strand, ".", ";".join(attrs)]
                )
                + "\n"
            )
    print(f"  genes_annotated.gff3: {n_total} genes, {n_annot} with functional annotation")


# ---------------------------------------------------------------------------
# A2. Reference-anchored synteny-block tracks (BED, one per ecotype)
# ---------------------------------------------------------------------------

def load_gff_positions(path):
    """gene name -> (chr, start, end) from an MCScanX .gff position file."""
    pos = {}
    with open(path) as fh:
        for line in fh:
            chrom, gene, start, end = line.rstrip("\n").split("\t")[:4]
            pos[gene] = (chrom, int(start), int(end))
    return pos


def build_synteny_bed(collinearity_path, gff_path, out_path, eco_label):
    """One BED feature per synteny block = reference footprint of its anchor genes.

    MCScanX orders the two sides by input order, so the reference is not always chr_b.
    We identify the reference side by its sanitized prefix: reference chromosomes start
    with 'rf' and reference genes with 'rf_'.
    """
    positions = load_gff_positions(gff_path)
    n = 0
    rows = []
    header_re = re.compile(
        r"Alignment (\d+):.*N=(\d+)\s+(\S+)&(\S+)\s+(plus|minus)"
    )
    with open(collinearity_path) as fh:
        cur = None
        ref_genes = []

        def flush():
            nonlocal n
            if cur is None or not ref_genes:
                return
            coords = [positions[g] for g in ref_genes if g in positions]
            if not coords:
                return
            # Reference side is whichever chromosome sanitized with the 'rf' prefix.
            if cur["chr_a"].startswith("rf"):
                ref_san, eco_san = cur["chr_a"], cur["chr_b"]
            else:
                ref_san, eco_san = cur["chr_b"], cur["chr_a"]
            ref_chrom = unsanitize_ref_chr(ref_san)
            if ref_chrom is None:
                return
            starts = [c[1] for c in coords]
            ends = [c[2] for c in coords]
            eco_contig = re.sub(r"^(ln|si)", "", eco_san)  # display label only
            strand = "+" if cur["orientation"] == "plus" else "-"
            name = f"{eco_label}:{eco_contig}|block{cur['block']}|n={cur['n_anchors']}"
            score = min(1000, int(cur["n_anchors"]) * 5)
            rows.append(
                (ref_chrom, min(starts), max(ends), name, score, strand)
            )
            n += 1

        for line in fh:
            m = header_re.search(line)
            if line.startswith("## Alignment") and m:
                flush()
                cur = {
                    "block": m.group(1),
                    "n_anchors": m.group(2),
                    "chr_a": m.group(3),
                    "chr_b": m.group(4),
                    "orientation": m.group(5),
                }
                ref_genes = []
            elif re.match(r"^\s*\d+-\s*\d+:", line):
                parts = line.split("\t")
                if len(parts) >= 3:
                    for tok in (parts[1].strip(), parts[2].strip()):
                        if tok.startswith("rf_"):  # reference gene of the pair
                            ref_genes.append(tok)
                            break
        flush()

    # BED must be coordinate-sorted for tabix/BedAdapter friendliness.
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(out_path, "w") as out:
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")
    print(f"  {os.path.basename(out_path)}: {n} synteny blocks")


# ---------------------------------------------------------------------------

def main():
    for d in [OUT_DIR] + STAGE_DIRS:
        os.makedirs(d, exist_ok=True)

    print("Building reference-anchored browser tracks (Phase 1)...")
    gff3 = os.path.join(OUT_DIR, "genes_annotated.gff3")
    lean_bed = os.path.join(OUT_DIR, "synteny_lean_blocks.bed")
    sisco_bed = os.path.join(OUT_DIR, "synteny_siscowet_blocks.bed")

    build_annotated_gff3(gff3)
    build_synteny_bed(
        os.path.join(SYN_DIR, "ref_lean.collinearity"),
        os.path.join(SYN_DIR, "ref_lean.gff"),
        lean_bed, "lean",
    )
    build_synteny_bed(
        os.path.join(SYN_DIR, "ref_sisco.collinearity"),
        os.path.join(SYN_DIR, "ref_sisco.gff"),
        sisco_bed, "siscowet",
    )

    outputs = [gff3, lean_bed, sisco_bed]
    for path in outputs:
        with open(path + ".md5", "w") as fh:
            fh.write(f"{md5sum(path)}  {os.path.basename(path)}\n")

    for path in outputs:
        for dest in STAGE_DIRS:
            shutil.copy2(path, os.path.join(dest, os.path.basename(path)))
    print(f"Staged {len(outputs)} tracks into: {', '.join(STAGE_DIRS)}")


if __name__ == "__main__":
    sys.exit(main())
