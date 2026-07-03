#!/usr/bin/env python3
"""
22.2 - Build JBrowse 2 Linear Synteny View data for lean vs siscowet (Phase 2).

Converts the MCScanX lean<->siscowet collinearity (code/21.1) into the file trio a
JBrowse 2 SyntenyTrack (MCScanAnchorsAdapter) needs, so the two purged ecotype
assemblies can be shown stacked with connecting synteny ribbons:

  - lean_siscowet.anchors        jcvi/MCScan anchors (geneA<TAB>geneB<TAB>score, ### between blocks)
  - lean_purged.mcscan.bed       lean anchor-gene positions (names match the anchors file)
  - siscowet_purged.mcscan.bed   siscowet anchor-gene positions

It also stages the two purged FASTA indexes (.fai) so the assemblies can be declared with
a local faiLocation while the large FASTAs stream from Gannet.

Assembly sequence names: MCScanX sanitized contig ids by removing non-alphanumerics and
prefixing a 2-char genome code (ln/si). We invert that per-genome from the .fai seq names so
the emitted BED chromosomes match the FASTA exactly. Gene ids in the anchors keep their
`ln_`/`si_` prefixes (that is how they appear in both the collinearity and the emitted BEDs).

Usage:  python code/22.2-jbrowse-synteny-view.py
"""

import hashlib
import os
import re
import shutil
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SYN_DIR = os.path.join(ROOT, "output", "21.1-gene-anchored-synteny")
OUT_DIR = os.path.join(ROOT, "output", "22-synteny-browser-tracks")
JB_DATA = os.path.join(ROOT, "jbrowse", "data", "lake-trout")

COLLINEARITY = os.path.join(SYN_DIR, "lean_sisco.collinearity")
COMBINED_GFF = os.path.join(SYN_DIR, "lean_sisco.gff")
LEAN_FAI = os.path.join(SYN_DIR, "lean.purged.fa.fai")
SISCO_FAI = os.path.join(SYN_DIR, "siscowet.purged.fa.fai")


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_unsanitize_map(fai_path, code):
    """sanitized seq id -> original FASTA seq id, inverting 21.1's sanitization
    (drop non-alphanumerics, prefix 2-char genome code)."""
    mapping = {}
    with open(fai_path) as fh:
        for line in fh:
            name = line.split("\t", 1)[0]
            san = code + re.sub(r"[^A-Za-z0-9]", "", name)
            mapping[san] = name
    return mapping


def build_anchors_and_beds():
    lean_map = build_unsanitize_map(LEAN_FAI, "ln")
    sisco_map = build_unsanitize_map(SISCO_FAI, "si")

    # BED positions come from the combined MCScanX gff (chr, gene, start, end); its gene
    # ids exactly match the collinearity, so anchors and BEDs stay consistent.
    lean_bed_rows = {}   # gene -> (chrom, start, end)
    sisco_bed_rows = {}
    unmapped_chr = set()
    with open(COMBINED_GFF) as fh:
        for line in fh:
            chrom, gene, start, end = line.rstrip("\n").split("\t")[:4]
            if gene.startswith("ln_"):
                orig = lean_map.get(chrom)
                if orig is None:
                    unmapped_chr.add(chrom); continue
                lean_bed_rows[gene] = (orig, int(start), int(end))
            elif gene.startswith("si_"):
                orig = sisco_map.get(chrom)
                if orig is None:
                    unmapped_chr.add(chrom); continue
                sisco_bed_rows[gene] = (orig, int(start), int(end))
    if unmapped_chr:
        print(f"  WARNING: {len(unmapped_chr)} sanitized contigs had no FASTA match "
              f"(e.g. {sorted(unmapped_chr)[:3]})")

    # Parse collinearity into anchor blocks; assign each pair by prefix so orientation of
    # the header (lean-first vs siscowet-first) never matters.
    anchors_path = os.path.join(OUT_DIR, "lean_siscowet.anchors")
    header_re = re.compile(r"Alignment \d+:.*score=([0-9.]+).*N=\d+")
    n_blocks = n_pairs = 0
    with open(COLLINEARITY) as fh, open(anchors_path, "w") as out:
        score = "0"
        started = False
        for line in fh:
            m = header_re.search(line)
            if line.startswith("## Alignment") and m:
                out.write("###\n")
                started = True
                n_blocks += 1
                score = str(int(float(m.group(1))))
            elif started and re.match(r"^\s*\d+-\s*\d+:", line):
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                g1, g2 = parts[1].strip(), parts[2].strip()
                lean_g = g1 if g1.startswith("ln_") else g2 if g2.startswith("ln_") else None
                sisco_g = g1 if g1.startswith("si_") else g2 if g2.startswith("si_") else None
                if lean_g is None or sisco_g is None:
                    continue
                if lean_g not in lean_bed_rows or sisco_g not in sisco_bed_rows:
                    continue
                out.write(f"{lean_g}\t{sisco_g}\t{score}\n")
                n_pairs += 1

    # Emit the two BEDs restricted to genes that participate in anchors is unnecessary;
    # JBrowse only looks up anchor genes, extra rows are harmless. Write full sets sorted.
    def write_bed(rows, path):
        ordered = sorted(rows.items(), key=lambda kv: (kv[1][0], kv[1][1]))
        with open(path, "w") as out:
            for gene, (chrom, start, end) in ordered:
                out.write(f"{chrom}\t{start}\t{end}\t{gene}\n")

    lean_bed = os.path.join(OUT_DIR, "lean_purged.mcscan.bed")
    sisco_bed = os.path.join(OUT_DIR, "siscowet_purged.mcscan.bed")
    write_bed(lean_bed_rows, lean_bed)
    write_bed(sisco_bed_rows, sisco_bed)

    print(f"  lean_siscowet.anchors: {n_blocks} blocks, {n_pairs} anchor pairs")
    print(f"  lean_purged.mcscan.bed: {len(lean_bed_rows)} genes")
    print(f"  siscowet_purged.mcscan.bed: {len(sisco_bed_rows)} genes")
    return [anchors_path, lean_bed, sisco_bed]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(JB_DATA, exist_ok=True)

    print("Building JBrowse Linear Synteny View data (Phase 2)...")
    outputs = build_anchors_and_beds()

    for path in outputs:
        with open(path + ".md5", "w") as fh:
            fh.write(f"{md5sum(path)}  {os.path.basename(path)}\n")
        shutil.copy2(path, os.path.join(JB_DATA, os.path.basename(path)))

    # Stage the small FASTA indexes so assemblies can use a local faiLocation.
    for fai in (LEAN_FAI, SISCO_FAI):
        shutil.copy2(fai, os.path.join(JB_DATA, os.path.basename(fai)))

    print(f"Staged synteny trio + 2 .fai into {JB_DATA}")


if __name__ == "__main__":
    sys.exit(main())
