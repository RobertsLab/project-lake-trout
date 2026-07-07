#!/usr/bin/env python3
"""
23.3-sisco-sv-repair.py  (Phase 2 prerequisite of code/23-next-phase-research-plan.md)

Repair the failed siscowet-vs-reference SyRI run from Phase 1 (code/23.1-genome-sv-map.Rmd).

WHAT FAILED
-----------
sisco_vs_ref SyRI exited with "Unequal number of chromosomes in the genomes. Exiting"
(analyses/23-genome-sv/syri/sisco_vs_ref.syri.log). Cause: in 23.1's anchor step, RagTag
scaffolded a 15 kb siscowet contig onto the reference MITOCHONDRION NC_036392.1 (16.6 kb, which
is < the 10 Mb chromosome cutoff and therefore NOT part of the 42-chromosome ref.chroms.fa axis).
The anchor grep `^NC_[0-9]+\\.[0-9]+_RagTag` kept that mito scaffold, so siscowet.anchored.fa has
43 sequences vs the reference's 42 nuclear chromosomes. lean placed nothing on the mito -> clean 42.

THE FIX
-------
1. Rebuild siscowet.anchored.fa restricted to the 42 names in anchor/ref.chroms.list (drop the mito
   scaffold). Keep a .43seq backup of the original.
2. Re-run SyRI for sisco_vs_ref, reusing the existing WGA BAM (wga/sisco_vs_ref.bam, whose target
   is already the 42-chromosome ref.chroms.fa).
3. Regenerate syri/merged_sv_table.tsv and syri/sv_summary.tsv across all THREE comparisons
   (port of 23.1 Step 4), refresh md5s.
4. (caller) re-run code/23.2-sv-pav-crossvalidate.py so siscowet high-conf PAV is emitted too.

Author: Generated for project-lake-trout
Date: 2026-07-04
"""

import subprocess
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
SV_DIR = BASE / "analyses" / "23-genome-sv"
ANCHOR = SV_DIR / "anchor"
WGA = SV_DIR / "wga"
SYRI = SV_DIR / "syri"
LOGS = SV_DIR / "logs"

SYRI_BIN = Path("/home/sr320/miniconda3/envs/syri_env/bin")
SAMTOOLS = SYRI_BIN / "samtools"
SYRI_EXE = SYRI_BIN / "syri"

REF_FA = ANCHOR / "ref.chroms.fa"
REF_LIST = ANCHOR / "ref.chroms.list"
SISCO_ANCHORED = ANCHOR / "siscowet.anchored.fa"
SISCO_BAM = WGA / "sisco_vs_ref.bam"
TAG = "sisco_vs_ref"

# Coarse SyRI-type -> class map (identical to 23.1 Step 4).
SV_CLASS = {
    "SYN": "syntenic", "SYNAL": "syntenic",
    "INV": "inversion", "INVAL": "inversion",
    "TRANS": "translocation", "TRANSAL": "translocation",
    "INVTR": "translocation", "INVTRAL": "translocation",
    "DUP": "duplication", "INVDP": "duplication",
    "INS": "indel", "DEL": "indel",
    "CPG": "copy_var", "CPL": "copy_var",
    "HDR": "diverged", "TDM": "tandem_repeat",
    "SNP": "snp", "NOTAL": "unaligned",
}
COMPARISONS = ["lean_vs_ref", "sisco_vs_ref", "lean_vs_sisco"]


def run(cmd, **kw):
    print("  $", " ".join(str(c) for c in cmd))
    subprocess.run(cmd, check=True, **kw)


def sh(cmd):
    print("  $", cmd)
    subprocess.run(cmd, shell=True, check=True)


def fai_names(fa):
    return [l.split("\t")[0] for l in (Path(str(fa) + ".fai")).read_text().splitlines()]


# ------------------------------------------------------------ Step 1: rebuild fasta
def rebuild_anchored():
    print("=" * 64)
    print("STEP 1  Rebuild siscowet.anchored.fa without the mito scaffold")
    print("=" * 64)
    names = fai_names(SISCO_ANCHORED)
    keep = REF_LIST.read_text().split()
    print(f"  current anchored sequences : {len(names)}")
    print(f"  reference chromosomes      : {len(keep)}")
    extra = [n for n in names if n not in set(keep)]
    print(f"  extra (to drop)            : {extra}")
    if len(names) == len(keep) and not extra:
        print("  already 42 sequences; nothing to rebuild.")
        return
    backup = ANCHOR / "siscowet.anchored.43seq.fa"
    if not backup.exists():
        SISCO_ANCHORED.replace(backup)
        Path(str(SISCO_ANCHORED) + ".fai").replace(Path(str(backup) + ".fai"))
        print(f"  backed up original -> {backup.name}")
    # Extract the 42 reference-named scaffolds, in reference order.
    sh(f"'{SAMTOOLS}' faidx '{backup}' {' '.join(keep)} > '{SISCO_ANCHORED}'")
    run([str(SAMTOOLS), "faidx", str(SISCO_ANCHORED)])
    n = len(fai_names(SISCO_ANCHORED))
    print(f"  rebuilt siscowet.anchored.fa : {n} sequences")
    assert n == 42, f"expected 42, got {n}"


# ------------------------------------------------------------ Step 2: re-run SyRI
def filter_bam():
    """SyRI derives the query chromosome set from the BAM. The WGA BAM was built while
    siscowet.anchored.fa still carried the mito scaffold, so it contains 43 query names
    (NC_036392.1 among them). Drop that record so the BAM's query set is the 42 nuclear
    chromosomes, matching the reference. Cheaper than re-running minimap2."""
    filt = WGA / f"{TAG}.filt.bam"
    if filt.exists() and filt.stat().st_size > 0:
        print(f"  reusing existing {filt.name}")
        return filt
    sh(f"'{SAMTOOLS}' view -b -e 'qname != \"NC_036392.1\"' "
       f"-o '{filt}' '{SISCO_BAM}'")
    run([str(SAMTOOLS), "index", str(filt)])
    print(f"  wrote {filt.name} (dropped mito query contig NC_036392.1)")
    return filt


def run_syri():
    print("=" * 64)
    print("STEP 2  Re-run SyRI for sisco_vs_ref (filtered WGA BAM, 42 query chroms)")
    print("=" * 64)
    bam = filter_bam()
    log = LOGS / f"{TAG}.syri.log"
    cmd = [
        str(SYRI_EXE),
        "-c", str(bam),
        "-r", str(REF_FA),
        "-q", str(SISCO_ANCHORED),
        "-F", "B",
        "--nc", "40",
        "--dir", str(SYRI),
        "--prefix", f"{TAG}.",
    ]
    print("  $", " ".join(cmd))
    with open(log, "w") as lf:
        r = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
    out = SYRI / f"{TAG}.syri.out"
    if r.returncode != 0 or not (out.exists() and out.stat().st_size > 0):
        print(f"  SyRI FAILED (rc={r.returncode}); tail of {log.name}:")
        sh(f"tail -n 20 '{log}'")
        sys.exit(1)
    print(f"  wrote {out.name} ({out.stat().st_size:,} bytes)")
    sh(f"cut -f11 '{out}' | sort | uniq -c | sort -rn")


# ------------------------------------------------------------ Step 3: merged table + summary
def regenerate_tables():
    print("=" * 64)
    print("STEP 3  Regenerate merged_sv_table.tsv + sv_summary.tsv (all 3 comparisons)")
    print("=" * 64)
    merged = SYRI / "merged_sv_table.tsv"
    summary = SYRI / "sv_summary.tsv"
    cols = ["ref_chr", "ref_start", "ref_end", "ref_seq", "qry_seq",
            "qry_chr", "qry_start", "qry_end", "uid", "parent", "type"]

    counts = {}  # (comparison, class) -> [n, total_ref_bp]
    with open(merged, "w") as mo:
        mo.write("\t".join(cols + ["comparison", "class"]) + "\n")
        for comp in COMPARISONS:
            path = SYRI / f"{comp}.syri.out"
            if not (path.exists() and path.stat().st_size > 0):
                print(f"  WARNING missing {path.name}; skipping")
                continue
            n_rows = 0
            with open(path) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 11:
                        continue
                    typ = p[10]
                    cls = SV_CLASS.get(typ, "other")
                    mo.write("\t".join(p[:11]) + f"\t{comp}\t{cls}\n")
                    n_rows += 1
                    # summary: structural + local variants only (drop syntenic/unaligned),
                    # numeric ref coords only.
                    if cls in ("syntenic", "unaligned"):
                        continue
                    try:
                        rs, re_ = float(p[1]), float(p[2])
                    except ValueError:
                        continue
                    span = max(int(re_) - int(rs), 0)
                    key = (comp, cls)
                    agg = counts.setdefault(key, [0, 0])
                    agg[0] += 1
                    agg[1] += span
            print(f"  {comp}: {n_rows:,} rows")

    with open(summary, "w") as so:
        so.write("comparison\tclass\tn\ttotal_ref_bp\n")
        for comp in COMPARISONS:
            rows = sorted([(c, v) for (cc, c), v in counts.items() if cc == comp],
                          key=lambda kv: -kv[1][0])
            for cls, (n, bp) in rows:
                so.write(f"{comp}\t{cls}\t{n}\t{bp}\n")
    print(f"  wrote {merged.name} and {summary.name}")
    sh(f"column -t '{summary}'")

    # refresh md5s for the regenerated + repaired files
    for f in [SYRI / f"{TAG}.syri.out", merged, summary]:
        if f.exists() and f.stat().st_size > 0:
            sh(f"md5sum '{f}' | tee '{f}.md5' > /dev/null")


def main():
    rebuild_anchored()
    run_syri()
    regenerate_tables()
    print("\nDone. Next: python3 code/23.2-sv-pav-crossvalidate.py "
          "(now emits siscowet high-conf PAV too).")


if __name__ == "__main__":
    main()
