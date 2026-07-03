#!/usr/bin/env python3
"""
23.2-sv-pav-crossvalidate.py  (Phase 1 cross-validation of code/23-next-phase-research-plan.md)

Cross-validate the reference-based PAV calls (code/15-diff-pav.py) against the ecotype SVs called
by SyRI (code/23.1-genome-sv-map.Rmd). Both call sets live on reference chromosome coordinates
(NC_0523##.# ...):

  * PAV : analyses/15-diff-pav/{lean,siscowet}_specific.{insertions,deletions}.browser.bed
          — regions present in one ecotype's reads but not the other, mapped to the LEAN reference.
  * SyRI: analyses/23-genome-sv/syri/{lean_vs_ref,sisco_vs_ref}.syri.out
          — structural + large local variants of each ecotype's own genome vs. the same reference.

A reference-based PAV that falls inside a SyRI-confirmed SV of the *matching* ecotype is corroborated
by that ecotype's own assembly, so it is no longer a single-reference association: it is the
HIGH-CONFIDENCE set. Lean PAVs are checked against lean_vs_ref SVs; siscowet PAVs against
sisco_vs_ref SVs.

  lean_specific.*      x  lean_vs_ref  SVs
  siscowet_specific.*  x  sisco_vs_ref SVs

Outputs (analyses/23-genome-sv/crossvalidate/):
  {ecotype}.{vtype}.highconf.bed   PAV features overlapping a matching-ecotype SV (BED9 + SV type)
  {ecotype}.{vtype}.svtype.tsv     each high-conf PAV -> the SyRI SV type(s) it overlaps
  crossvalidate_summary.csv        counts: PAV total / overlapping / fraction, per ecotype x vtype

Author: Generated for project-lake-trout
Date: 2026-07-03
"""

import csv
import subprocess
from collections import defaultdict
from pathlib import Path

# =============================================================================
# Configuration
# =============================================================================

BASE_DIR = Path(__file__).resolve().parent.parent
PAV_DIR = BASE_DIR / "analyses" / "15-diff-pav"
SYRI_DIR = BASE_DIR / "analyses" / "23-genome-sv" / "syri"
OUTPUT_DIR = BASE_DIR / "analyses" / "23-genome-sv" / "crossvalidate"

# ecotype -> (PAV file prefix, matching SyRI comparison tag)
ECOTYPES = {
    "lean":     {"pav_prefix": "lean_specific",     "syri_tag": "lean_vs_ref"},
    "siscowet": {"pav_prefix": "siscowet_specific", "syri_tag": "sisco_vs_ref"},
}
VARIANT_TYPES = ["insertions", "deletions"]

# SyRI annotation types (column 11) that count as a structural / large-indel SV. Syntenic
# backbone (SYN/SYNAL), aligned-region markers (*AL), SNPs, and unaligned gaps (NOTAL) are
# excluded: NOTAL in particular is "no alignment", not a called variant, so overlapping it
# would corroborate nothing.
SV_TYPES = {
    "INV",    # inversion
    "TRANS",  # translocation
    "INVTR",  # inverted translocation
    "DUP",    # duplication
    "INVDP",  # inverted duplication
    "INS",    # insertion (within aligned block)
    "DEL",    # deletion
    "CPG",    # copy gain
    "CPL",    # copy loss
    "HDR",    # highly diverged region
    "TDM",    # tandem repeat
}

# Require the PAV feature to be this fraction covered by the SV region to count as corroborated.
MIN_OVERLAP_FRACTION = 0.5

BEDTOOLS = "bedtools"  # /usr/bin/bedtools

# =============================================================================
# Helpers
# =============================================================================


def nonempty(path):
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def run(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout


def read_pav_bed(path):
    """Read a 15-diff-pav browser BED (BED9 + track header). Returns list of feature dicts."""
    feats = []
    with open(path) as f:
        for line in f:
            if line.startswith(("track", "#")) or not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            feats.append({"chrom": p[0], "start": int(p[1]), "end": int(p[2]),
                          "cols": p})
    return feats


def syri_sv_to_bed(syri_out, bed_path):
    """Extract SV regions (reference coordinates) from a SyRI .syri.out into a BED4.

    Column layout of *.syri.out (11 cols, no header):
      0 ref_chr 1 ref_start 2 ref_end 3 ref_seq 4 qry_seq
      5 qry_chr 6 qry_start 7 qry_end 8 uid 9 parent 10 annotation_type
    Keep rows whose type is in SV_TYPES and whose reference coordinates are numeric.
    """
    n = 0
    with open(syri_out) as fin, open(bed_path, "w") as fout:
        for line in fin:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11:
                continue
            sv_type = p[10]
            if sv_type not in SV_TYPES:
                continue
            chrom, rstart, rend = p[0], p[1], p[2]
            try:
                start, end = int(float(rstart)), int(float(rend))
            except ValueError:
                continue  # ref coords are '-' for query-only features
            if end < start:
                start, end = end, start
            if end == start:
                end = start + 1  # BED requires end > start
            fout.write(f"{chrom}\t{start}\t{end}\t{sv_type}\n")
            n += 1
    return n


def sort_bed(path):
    run(f"sort -k1,1 -k2,2n '{path}' -o '{path}'")


# =============================================================================
# Main
# =============================================================================


def main():
    print("=" * 64)
    print("SV / PAV cross-validation (Phase 1)")
    print("=" * 64)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Build a sorted SV BED per ecotype comparison from SyRI output.
    sv_beds = {}
    for eco, cfg in ECOTYPES.items():
        syri_out = SYRI_DIR / f"{cfg['syri_tag']}.syri.out"
        if not nonempty(syri_out):
            print(f"  MISSING SyRI output for {eco}: {syri_out}")
            print(f"  Run code/23.1-genome-sv-map.Rmd first. Skipping {eco}.")
            continue
        sv_bed = OUTPUT_DIR / f"{cfg['syri_tag']}.sv.bed"
        n_sv = syri_sv_to_bed(syri_out, sv_bed)
        sort_bed(sv_bed)
        sv_beds[eco] = sv_bed
        print(f"  {eco}: {n_sv} SV regions from {cfg['syri_tag']}.syri.out -> {sv_bed.name}")

    # 2. For each ecotype x variant type, intersect PAV features with matching-ecotype SVs.
    summary = []
    for eco, cfg in ECOTYPES.items():
        if eco not in sv_beds:
            continue
        for vtype in VARIANT_TYPES:
            pav_bed = PAV_DIR / f"{cfg['pav_prefix']}.{vtype}.browser.bed"
            if not nonempty(pav_bed):
                print(f"  no PAV file: {pav_bed.name} (skipping)")
                continue

            pav_feats = read_pav_bed(pav_bed)
            n_pav = len(pav_feats)

            # Write a clean, header-less, sorted copy of the PAV features for bedtools.
            pav_clean = OUTPUT_DIR / f"{eco}.{vtype}.pav.bed"
            with open(pav_clean, "w") as f:
                for ft in pav_feats:
                    f.write("\t".join(ft["cols"]) + "\n")
            sort_bed(pav_clean)

            # PAV features overlapping a matching-ecotype SV region by >= MIN_OVERLAP_FRACTION.
            highconf = OUTPUT_DIR / f"{eco}.{vtype}.highconf.bed"
            run(f"{BEDTOOLS} intersect -a '{pav_clean}' -b '{sv_beds[eco]}' "
                f"-u -f {MIN_OVERLAP_FRACTION} > '{highconf}'")

            # Also record which SV type(s) each high-conf PAV overlaps (-wa -wb).
            svtype = OUTPUT_DIR / f"{eco}.{vtype}.svtype.tsv"
            run(f"{BEDTOOLS} intersect -a '{pav_clean}' -b '{sv_beds[eco]}' "
                f"-wa -wb -f {MIN_OVERLAP_FRACTION} > '{svtype}'")

            n_high = sum(1 for _ in open(highconf)) if nonempty(highconf) else 0

            # Tally the SV types that corroborated PAVs of this ecotype/vtype.
            type_counts = defaultdict(int)
            if nonempty(svtype):
                with open(svtype) as f:
                    for line in f:
                        cols = line.rstrip("\n").split("\t")
                        type_counts[cols[-1]] += 1

            pav_clean.unlink()  # keep only highconf + svtype outputs

            frac = (n_high / n_pav) if n_pav else 0.0
            top = ", ".join(f"{t}:{c}" for t, c in
                            sorted(type_counts.items(), key=lambda kv: -kv[1]))
            print(f"  {eco:8s} {vtype:10s}  PAV={n_pav:6d}  "
                  f"SV-confirmed={n_high:6d} ({frac:5.1%})  [{top}]")

            summary.append({
                "ecotype": eco,
                "variant_type": vtype,
                "pav_total": n_pav,
                "sv_confirmed": n_high,
                "fraction_confirmed": round(frac, 4),
                "sv_types": top,
            })

    # 3. Write summary.
    if summary:
        summary_file = OUTPUT_DIR / "crossvalidate_summary.csv"
        with open(summary_file, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
            w.writeheader()
            w.writerows(summary)
        print(f"\n  summary -> {summary_file}")

    print("\nDone. High-confidence (SV-corroborated) PAV BEDs are in:")
    print(f"  {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
