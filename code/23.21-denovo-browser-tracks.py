#!/usr/bin/env python3
"""
23.21-denovo-browser-tracks.py  (Phase 4, step 6, of code/23-next-phase-research-plan.md)

JBrowse 2 tracks for the Phase 4 annotation on the purged ecotype assemblies (lean_purged,
siscowet_purged). IGV is reference-only, and the new genes are by definition the ones that may
have no reference counterpart, so these tracks are JBrowse-only.

Per ecotype:
  {asm}.genes_merged.bed           BED9, every merged gene (23.19 / 23.20 gene_function_table):
                                   Liftoff grey, BRAKER hq teal, BRAKER lq light teal, TE-like
                                   BRAKER light grey; name = symbol.
  {asm}.braker_new.gff3.gz(.tbi)   BRAKER-unique gene models with exon / CDS structure; the gene
                                   line carries Name, product, quality, ref_status, sp_acc,
                                   pfam_domains, te_like so they show in the feature panel.
                                   bgzip + tabix (Gff3TabixAdapter).
The full Liftoff exon-level GFF3 (375-415 MB) is not staged; Liftoff structure is already
available as the lean_purged_genes / siscowet_purged_genes tracks.

Two steps, because bgzip / tabix live on raven:
  build  (raven)  python code/23.21-denovo-browser-tracks.py build [--p4 DIR]
                  -> analyses/23-denovo-annotation/tracks/
  stage  (local)  python3 code/23.21-denovo-browser-tracks.py stage --from DIR [--config FILE]
                  copies the track files into jbrowse/data/lake-trout/ and adds the tracks to
                  jbrowse/config.json (idempotent: existing trackIds are replaced)

Author: Generated for project-lake-trout
Date: 2026-10-17
"""

import argparse
import json
import shutil
import subprocess
import urllib.parse
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
HTSBIN = "/home/sr320/miniconda3/envs/syri_env/bin"
ECOS = {"lean": "lean_purged", "siscowet": "siscowet_purged"}
RGB = {"liftoff": "107,114,128", "hq": "13,148,136", "lq": "94,234,212", "te": "209,213,219"}
JB_DATA = BASE / "jbrowse" / "data" / "lake-trout"
# JBrowse 3.6's default BED9 autoSql types column 9 as `uint reserved`, so an "R,G,B" string is
# lost and jexl colouring by itemRgb renders black. BED tracks coloured by itemRgb declare this
# schema instead (one field per line; the autoSql parser needs the newlines).
BED9_RGB_AUTOSQL = """table bed9rgb
"BED9 with a string itemRgb column"
    (
    string chrom;      "Reference sequence chromosome or scaffold"
    uint   chromStart; "Start position in chromosome"
    uint   chromEnd;   "End position in chromosome"
    string name;       "Name of item"
    uint   score;      "Score"
    char[1] strand;    "+ or - or ."
    uint   thickStart; "Start of where display should be thick"
    uint   thickEnd;   "End of where display should be thick"
    string itemRgb;    "R,G,B colour"
    )
"""
JB_CONFIG = BASE / "jbrowse" / "config.json"


def enc(v):
    return urllib.parse.quote(str(v), safe=" :|/()-.+_'")


# =============================================================================
def build(p4):
    import pandas as pd
    out = p4 / "tracks"
    out.mkdir(parents=True, exist_ok=True)
    for eco, asm in ECOS.items():
        final = p4 / "final" / eco
        ft = final / "gene_function_table.tsv"
        if not ft.exists():
            print(f"{eco}: no gene_function_table.tsv (run 23.19 + 23.20); skipped")
            continue
        t = pd.read_csv(ft, sep="\t", dtype=str, keep_default_na=False)
        t["start"], t["end"] = t["start"].astype(int), t["end"].astype(int)
        te = t.get("te_like", "") == "1"
        cls = t["provenance"].where(t["provenance"] == "liftoff", t.get("quality", "lq"))
        cls = cls.where(~(te & (t["provenance"] == "braker")), "te")
        bed = pd.DataFrame({0: t["chrom"], 1: t["start"], 2: t["end"],
                            3: t["symbol"].where(t["symbol"] != "", t["gene_id"]), 4: 0,
                            5: t["strand"], 6: t["start"], 7: t["end"], 8: cls.map(RGB)})
        bed.sort_values([0, 1]).to_csv(out / f"{asm}.genes_merged.bed", sep="\t", header=False, index=False)

        # new-gene GFF3 with annotation on the gene line
        f = t[t["provenance"] == "braker"].set_index("gene_id")
        keys = ["symbol", "product", "ref_status", "sp_acc", "pfam_domains", "te_like"]  # quality: set by 23.19
        rows = []
        for line in open(final / f"{eco}.braker_new.gff3"):
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if c[2] == "gene":
                gid = dict(kv.split("=", 1) for kv in c[8].strip(";").split(";") if "=" in kv)["ID"]
                if gid in f.index:
                    a = f.loc[gid]
                    extra = [f"Name={enc(a['symbol'])}"] + \
                            [f"{k}={enc(a[k])}" for k in keys[1:] if k in a and a[k] != ""]
                    c[8] = c[8].rstrip(";") + ";" + ";".join(extra)
            rows.append(c)
        rows.sort(key=lambda c: (c[0], int(c[3]), 0 if c[2] == "gene" else 1))
        gff = out / f"{asm}.braker_new.gff3"
        with open(gff, "w") as fh:
            fh.write("##gff-version 3\n")
            for c in rows:
                fh.write("\t".join(c) + "\n")
        subprocess.run(f"'{HTSBIN}/bgzip' -f '{gff}' && '{HTSBIN}/tabix' -f -p gff '{gff}.gz'",
                       shell=True, check=True)
        n_new = int((t["provenance"] == "braker").sum())
        print(f"{eco}: {len(t):,} genes in BED ({n_new:,} BRAKER: "
              f"{dict(cls[t['provenance'] == 'braker'].value_counts())}); new-gene GFF3 {len(rows):,} lines")


# =============================================================================
def track_defs():
    color_bed = "jexl:'rgb('+get(feature,'itemRgb')+')'"
    color_new = ("jexl:get(feature,'te_like')=='1' ? 'rgb(209,213,219)' : "
                 "(get(feature,'quality')=='hq' ? 'rgb(13,148,136)' : 'rgb(94,234,212)')")
    defs = []
    for eco, asm in ECOS.items():
        defs.append({
            "type": "FeatureTrack", "trackId": f"{asm}_genes_merged",
            "name": f"{eco.capitalize()} genes, Liftoff + BRAKER3 (merged)",
            "category": ["Annotation", "Phase 4 de novo"], "assemblyNames": [asm],
            "adapter": {"type": "BedAdapter", "bedLocation": {"uri": f"data/lake-trout/{asm}.genes_merged.bed"},
                        "autoSql": BED9_RGB_AUTOSQL},
            "displays": [{"type": "LinearBasicDisplay", "displayId": f"{asm}_genes_merged-LinearBasicDisplay",
                          "renderer": {"type": "SvgFeatureRenderer", "color1": color_bed}}]})
        defs.append({
            "type": "FeatureTrack", "trackId": f"{asm}_braker_new",
            "name": f"{eco.capitalize()} BRAKER3-only genes (not in Liftoff)",
            "category": ["Annotation", "Phase 4 de novo"], "assemblyNames": [asm],
            "adapter": {"type": "Gff3TabixAdapter",
                        "gffGzLocation": {"uri": f"data/lake-trout/{asm}.braker_new.gff3.gz"},
                        "index": {"location": {"uri": f"data/lake-trout/{asm}.braker_new.gff3.gz.tbi"}}},
            "displays": [{"type": "LinearBasicDisplay", "displayId": f"{asm}_braker_new-LinearBasicDisplay",
                          "renderer": {"type": "SvgFeatureRenderer", "color1": color_new}}]})
    return defs


def stage(src, config, data_dir):
    files = []
    for asm in ECOS.values():
        for name in (f"{asm}.genes_merged.bed", f"{asm}.braker_new.gff3.gz", f"{asm}.braker_new.gff3.gz.tbi"):
            if (src / name).exists():
                files.append(name)
    staged_asm = {asm for asm in ECOS.values() if f"{asm}.genes_merged.bed" in files}
    if not staged_asm:
        raise SystemExit(f"no track files in {src}")
    data_dir.mkdir(parents=True, exist_ok=True)
    for name in files:
        shutil.copy(src / name, data_dir / name)
    c = json.load(open(config))
    new = [d for d in track_defs() if d["assemblyNames"][0] in staged_asm]
    ids = {d["trackId"] for d in new}
    c["tracks"] = [t for t in c["tracks"] if t["trackId"] not in ids] + new
    with open(config, "w") as fh:
        json.dump(c, fh, indent=2)
        fh.write("\n")
    print(f"staged {len(files)} files for {sorted(staged_asm)} into {data_dir}; "
          f"{len(new)} tracks in {config}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--p4", type=Path, default=BASE / "analyses" / "23-denovo-annotation")
    s = sub.add_parser("stage")
    s.add_argument("--from", dest="src", type=Path, required=True)
    s.add_argument("--config", type=Path, default=JB_CONFIG)
    s.add_argument("--data-dir", type=Path, default=JB_DATA)
    a = ap.parse_args()
    if a.cmd == "build":
        build(a.p4 if a.p4.is_absolute() else BASE / a.p4)
    else:
        stage(a.src, a.config, a.data_dir)


if __name__ == "__main__":
    main()
