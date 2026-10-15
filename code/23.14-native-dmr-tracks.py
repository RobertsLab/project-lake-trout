#!/usr/bin/env python3
"""
23.14-native-dmr-tracks.py  (Phase 3 -> Phase 5 bridge, code/23-next-phase-research-plan.md)

Turns the 23.13 native DMRs into (1) a gene-level evidence line for the 23.7 model and (2) browser
tracks.

Status of each primary 23.13 DMR (one row per frame):
  replicated     placed on the other assembly and same sign, p < 0.05, |diff| >= 10 there
  single_frame   not placed uniquely / not testable in the other frame (FDR-controlled in its own
                 frame only, eFDR ~0.06)
  not_replicated placed and testable in the other frame but failed the replication test

1. Gene assignment (native coordinates). Every Liftoff gene of the DMR's frame assembly within
   GENE_WINDOW bp is listed with a strand-aware location class: body, promoter (<= 2 kb upstream of
   the TSS), upstream / downstream (<= 5 kb). Gene IDs are canonical (Liftoff `_N` suffix removed),
   matching 23.7. 23.7 builds its `native_dmr` line from the `replicated` rows.
2. Reference placement. Each DMR (+/- 500 bp) is placed on SaNama_1.0 with the 23.11
   `place_intervals` rule (unique primary minimap2 asm5 hit, MAPQ >= 20). DMRs from the two frames
   that land on overlapping reference intervals with the same direction are merged into one locus.
   This step needs the assemblies (raven); its output is cached in cpg_scan/native_dmrs.on_ref.tsv
   and the remaining steps run anywhere.
3. Tracks (BED9, itemRgb: blue = hypomethylated in siscowet, red = hypermethylated; dark =
   replicated, light = single_frame, grey = not_replicated -- native tracks only):
     native_dmrs.lean_purged.bed, native_dmrs.siscowet_purged.bed   JBrowse purged assemblies
     native_dmrs_hypo_siscowet.ref.bed, native_dmrs_hyper_siscowet.ref.bed
                                   reference loci (replicated + single_frame) for IGV and JBrowse
   staged into genome-browser/data/methylation/ and jbrowse/data/lake-trout/.

Outputs (analyses/23-native-methylation/cpg_scan/):
  native_dmr_gene_assignments.tsv   dmr x gene within 5 kb, with status
  native_dmrs.on_ref.tsv            reference placement per DMR (cached)
  native_dmr_loci.ref.tsv           merged reference loci
  tracks/*.bed

Run (raven, once, for step 2): /home/sr320/miniconda3/bin/python code/23.14-native-dmr-tracks.py
Then locally to stage tracks:   python3 code/23.14-native-dmr-tracks.py

Author: Generated for project-lake-trout
Date: 2026-10-15
"""

import importlib.util
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
SCAN = BASE / "analyses" / "23-native-methylation" / "cpg_scan"
TRACKS = SCAN / "tracks"
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
GENES = {"lean": SYN / "lean.purged.liftoff.genes.bed",
         "siscowet": SYN / "siscowet.purged.liftoff.genes.bed"}
OTHER = {"lean": "siscowet", "siscowet": "lean"}
STAGE = {"igv": BASE / "genome-browser" / "data" / "methylation",
         "jbrowse": BASE / "jbrowse" / "data" / "lake-trout"}

GENE_WINDOW = 5000
PROMOTER_BP = 2000
COLORS = {("hypo", "replicated"): "37,99,235", ("hypo", "single_frame"): "147,197,253",
          ("hyper", "replicated"): "220,38,38", ("hyper", "single_frame"): "252,165,165",
          ("hypo", "not_replicated"): "160,160,160", ("hyper", "not_replicated"): "160,160,160"}


def load_dmrs():
    rows = []
    for f in ("lean", "siscowet"):
        d = pd.read_csv(SCAN / f"dmrs_{f}.tsv", sep="\t")
        d["frame"] = f
        testable = d["other_frame_testable"].astype(str) == "True"
        d["status"] = np.where(d["replicated"].astype(str) == "True", "replicated",
                               np.where(testable, "not_replicated", "single_frame"))
        d["dir"] = np.where(d["meth_diff"] < 0, "hypo", "hyper")
        rows.append(d)
    return pd.concat(rows, ignore_index=True)


def gene_assignments(dmrs):
    out = []
    for f in ("lean", "siscowet"):
        g = pd.read_csv(GENES[f], sep="\t", header=None,
                        names=["contig", "gstart", "gend", "gene_raw", "score", "strand"])
        g["gene_id"] = g["gene_raw"].str.replace(r"_\d+$", "", regex=True)
        by = {c: x for c, x in g.groupby("contig")}
        for r in dmrs[dmrs["frame"] == f].itertuples(index=False):
            x = by.get(r.contig)
            if x is None:
                continue
            x = x[(x["gstart"] < r.end + GENE_WINDOW) & (x["gend"] > r.start - GENE_WINDOW)]
            for gr in x.itertuples(index=False):
                if r.start < gr.gend and r.end > gr.gstart:
                    cls, dist = "body", 0
                else:
                    # distance from the DMR to the gene; upstream relative to gene strand
                    left = r.end <= gr.gstart          # DMR lies left of the gene
                    dist = gr.gstart - r.end if left else r.start - gr.gend
                    upstream = left if gr.strand != "-" else not left
                    cls = ("promoter" if upstream and dist <= PROMOTER_BP else
                           "upstream" if upstream else "downstream")
                out.append({"dmr_id": r.dmr_id, "frame": f, "status": r.status,
                            "meth_diff": r.meth_diff, "direction": r.direction,
                            "gene_id": gr.gene_id, "gene_copy_id": gr.gene_raw,
                            "location_class": cls, "distance_bp": dist})
    a = pd.DataFrame(out)
    a.to_csv(SCAN / "native_dmr_gene_assignments.tsv", sep="\t", index=False, float_format="%.6g")
    return a


def place_on_reference(dmrs):
    cache = SCAN / "native_dmrs.on_ref.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    spec = importlib.util.spec_from_file_location("m11", BASE / "code" / "23.11-native-meth-counts.py")
    m11 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m11)
    if not m11.REF_FA.exists():
        sys.exit(f"FATAL: {m11.REF_FA} not found -- run the reference placement on raven first")
    rows = []
    for f in ("lean", "siscowet"):
        d = dmrs[dmrs["frame"] == f]
        placed = m11.place_intervals(m11.ECO_FA[f],
                                     list(zip(d["dmr_id"], d["contig"], d["start"], d["end"])),
                                     {"ref": m11.REF_FA}, f"cpgscan_{f}_ref")
        for (dmr_id, _), (chrom, s, e, mapq) in placed.items():
            rows.append({"dmr_id": dmr_id, "frame": f, "ref_chrom": chrom, "ref_start": s,
                         "ref_end": e, "ref_mapq": mapq})
    p = pd.DataFrame(rows)
    p.to_csv(cache, sep="\t", index=False)
    return p


def merge_ref_loci(dmrs, placed):
    """Overlapping same-direction reference intervals (either frame) -> one locus."""
    d = dmrs.merge(placed, on=["dmr_id", "frame"], how="inner")
    d = d[d["status"] != "not_replicated"].sort_values(["ref_chrom", "ref_start"])
    loci, cur = [], None
    for r in d.itertuples(index=False):
        if cur and r.ref_chrom == cur["chrom"] and r.ref_start < cur["end"] and r.dir == cur["dir"]:
            cur["end"] = max(cur["end"], r.ref_end)
            cur["members"].append(r)
            continue
        if cur:
            loci.append(cur)
        cur = {"chrom": r.ref_chrom, "start": r.ref_start, "end": r.ref_end, "dir": r.dir,
               "members": [r]}
    if cur:
        loci.append(cur)
    rows = []
    for i, l in enumerate(loci):
        m = l["members"]
        rows.append({
            "locus_id": f"nDMR_{i + 1}", "chrom": l["chrom"], "start": l["start"], "end": l["end"],
            "direction": f"{l['dir']}_siscowet",
            "status": "replicated" if any(x.status == "replicated" for x in m) else "single_frame",
            "frames": ",".join(sorted({x.frame for x in m})),
            "dmr_ids": ",".join(x.dmr_id for x in m),
            "meth_diff": round(float(np.mean([x.meth_diff for x in m])), 2),
            "symbol": next((x.symbol for x in m if isinstance(x.symbol, str) and x.symbol), ""),
        })
    loci = pd.DataFrame(rows)
    loci.to_csv(SCAN / "native_dmr_loci.ref.tsv", sep="\t", index=False)
    return loci


def bed9(df, chrom, start, end, name, dir_col, status_col):
    rgb = [COLORS[(d, s)] for d, s in zip(df[dir_col], df[status_col])]
    score = np.minimum(1000, (df["meth_diff"].abs() * 10).round()).astype(int)
    return pd.DataFrame({0: df[chrom], 1: df[start], 2: df[end], 3: df[name], 4: score, 5: ".",
                         6: df[start], 7: df[end], 8: rgb})


def main():
    TRACKS.mkdir(parents=True, exist_ok=True)
    dmrs = load_dmrs()
    dmrs["direction"] = dmrs["dir"] + "_siscowet"
    assign = gene_assignments(dmrs)
    placed = place_on_reference(dmrs)

    written = []
    for f, asm in (("lean", "lean_purged"), ("siscowet", "siscowet_purged")):
        d = dmrs[dmrs["frame"] == f].copy()
        d["label"] = d["dmr_id"] + "|" + d["meth_diff"].map("{:+.0f}".format) + "|" + d["status"]
        b = bed9(d.sort_values(["contig", "start"]), "contig", "start", "end", "label",
                 "dir", "status")
        path = TRACKS / f"native_dmrs.{asm}.bed"
        b.to_csv(path, sep="\t", header=False, index=False)
        written.append(path)

    loci = merge_ref_loci(dmrs, placed)
    loci["dir"] = loci["direction"].str.replace("_siscowet", "")
    loci["label"] = loci["locus_id"] + "|" + loci["meth_diff"].map("{:+.0f}".format) + "|" + \
        loci["status"] + "|" + loci["dmr_ids"]
    for d in ("hypo", "hyper"):
        x = loci[loci["dir"] == d].sort_values(["chrom", "start"])
        path = TRACKS / f"native_dmrs_{d}_siscowet.ref.bed"
        bed9(x, "chrom", "start", "end", "label", "dir", "status").to_csv(
            path, sep="\t", header=False, index=False)
        written.append(path)

    for path in written:
        targets = [STAGE["jbrowse"]] + ([STAGE["igv"]] if path.name.endswith(".ref.bed") else [])
        for t in targets:
            if t.exists():
                shutil.copy(path, t / path.name)

    rep = assign[assign["status"] == "replicated"]
    print(f"DMRs: {len(dmrs)} ({dmrs['status'].value_counts().to_dict()})")
    print(f"gene assignments within {GENE_WINDOW} bp: {len(assign)} rows, "
          f"{assign['gene_id'].nunique()} genes; replicated-DMR genes {rep['gene_id'].nunique()} "
          f"({rep.drop_duplicates('gene_id')['location_class'].value_counts().to_dict()})")
    print(f"placed on reference: {len(placed)} of {len(dmrs)}; merged loci {len(loci)} "
          f"({loci['status'].value_counts().to_dict()})")
    print("tracks:", ", ".join(p.name for p in written))


if __name__ == "__main__":
    main()
