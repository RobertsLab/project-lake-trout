#!/usr/bin/env python3
"""
23.17-lean-reference-baseline.py  (follow-up to Phase 1, code/23-next-phase-research-plan.md)

LEAN vs PUBLISHED LEAN AS A BASELINE FOR LEAN vs SISCOWET. The reference SaNama_1.0
(GCF_016432855.1, BioSample SAMN16985724) is a Seneca-strain, gynogenetic doubled-haploid lean
female: one homozygous lean haplotype. Our lean and siscowet assemblies are hifiasm primary
contigs from four pooled fish each, purged and RagTag-anchored to the reference (23.1). Two lean
genomes should differ only by within-lean variation (strain / population) and assembly artifacts,
so lean vs reference sets the scale against which lean vs siscowet differences are read.

All three SyRI runs from 23.1 / 23.3 are reused (no new alignment):
    lean_vs_ref    ref = SaNama_1.0   qry = lean.anchored
    sisco_vs_ref   ref = SaNama_1.0   qry = siscowet.anchored
    lean_vs_sisco  ref = lean.anchored qry = siscowet.anchored

1. Genome-wide scale: SyRI summary per comparison (syntenic / inverted / translocated / duplicated
   / unaligned / highly diverged bp, SNPs, indels) and lean_vs_ref as a fraction of lean_vs_sisco.
2. Polarizing lean-siscowet SVs. Each lean_vs_sisco event (>= MIN_SV bp on either side; aligned
   children of SVs, SYN and SNPs excluded) has a lean-side interval (lean.anchored coordinates,
   shared with lean_vs_ref qry) and a siscowet-side interval (siscowet.anchored, shared with
   sisco_vs_ref qry). A side "differs from the reference" when >= MIN_COV of its interval (padded
   PAD bp) is covered by that assembly's own SVs against the reference. Classes:
       siscowet_side    lean matches the reference, siscowet does not
       lean_side        siscowet matches the reference, our lean does not: the published lean
                        agrees with siscowet, so this is not a lean-ecotype difference
       both             both differ from the reference
       neither          neither side shows an SV against the reference here
   Unaligned (NOTAL) events have one side only and are reported on that side. Null: each side's
   events shuffled within chromosome (bedtools shuffle -chrom, N_SHUF seeds); the shuffle ignores
   the clustering of SVs in repeats, so it is a lenient null.
3. Reciprocal-PAV and corroborated ecotype-only genes (23.7) are reference genes lifted by
   Liftoff, so every one is annotated in the published lean genome. Per gene: fraction of its
   reference interval aligned by each anchored assembly (aligned SyRI rows), as a check that the
   assemblies agree with the read calls.

Inputs (raven): analyses/23-genome-sv/syri/{lean_vs_ref,sisco_vs_ref,lean_vs_sisco}.syri.{out,summary},
  anchor/{lean,siscowet}.anchored.fa.fai; analyses/23-integration/refined/gene_evidence_refined.tsv
  (pass --genes if raven's copy is stale).

Outputs (analyses/23-lean-baseline/):
  genome_scale.tsv            SyRI summary rows x comparison, plus lean_vs_ref / lean_vs_sisco
  polarized_sv.tsv            counts and bp per SV class x polarization class
  polarized_sv_null.tsv       same fractions for shuffled positions
  gene_alignment.tsv          per reciprocal / corroborated gene: aligned fraction in each assembly
  gene_alignment_summary.tsv  medians per gene group
  README.md

Compute: bedtools on ~1-2 M intervals; minutes.

Author: Generated for project-lake-trout
Date: 2026-10-16
"""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

# =============================================================================
# Configuration
# =============================================================================
BASE = Path(__file__).resolve().parent.parent
SV = BASE / "analyses" / "23-genome-sv"
SYRI = SV / "syri"
ANCHOR = SV / "anchor"
GENES = BASE / "analyses" / "23-integration" / "refined" / "gene_evidence_refined.tsv"
OUT = BASE / "analyses" / "23-lean-baseline"
BEDTOOLS = "bedtools"

COMPARISONS = ["lean_vs_ref", "sisco_vs_ref", "lean_vs_sisco"]
MIN_SV = 50     # bp; drops small indels
PAD = 50        # bp each side, so point events (INS ref side, DEL qry side) can overlap
MIN_COV = 0.5   # fraction of a padded side interval covered by the other comparison's SVs
N_SHUF = 3

SV_CLASS = {
    "INV": "inversion", "TRANS": "translocation", "INVTR": "translocation",
    "DUP": "duplication", "INVDP": "duplication", "CPG": "copy_var", "CPL": "copy_var",
    "NOTAL": "unaligned", "HDR": "diverged", "DEL": "indel", "INS": "indel",
    "TDM": "tandem_repeat",
}
ALIGNED = {"SYNAL", "INVAL", "TRANSAL", "INVTRAL", "DUPAL", "INVDPAL"}
SYRI_COLS = ["ref_chr", "ref_start", "ref_end", "ref_seq", "qry_seq", "qry_chr", "qry_start",
             "qry_end", "uid", "parent", "type", "copy"]


# =============================================================================
# Helpers
# =============================================================================
def run(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"FAILED: {cmd}\n{r.stderr}")
    return r.stdout


def read_syri(comp):
    df = pd.read_csv(SYRI / f"{comp}.syri.out", sep="\t", header=None, names=SYRI_COLS,
                     usecols=["ref_chr", "ref_start", "ref_end", "qry_chr", "qry_start",
                              "qry_end", "uid", "type"], dtype=str)
    for side in ("ref", "qry"):
        has = df[f"{side}_chr"] != "-"
        for c in ("start", "end"):
            col = f"{side}_{c}"
            df[col] = pd.to_numeric(df[col].where(has), errors="coerce")
        # inverted events report query start > end
        lo = df[[f"{side}_start", f"{side}_end"]].min(axis=1)
        hi = df[[f"{side}_start", f"{side}_end"]].max(axis=1)
        df[f"{side}_start"], df[f"{side}_end"] = lo, hi
        df[f"{side}_len"] = (df[f"{side}_end"] - df[f"{side}_start"] + 1).where(has, 0)
    return df


def sv_events(df):
    ev = df[df["type"].isin(SV_CLASS)].copy()
    ev = ev[ev[["ref_len", "qry_len"]].max(axis=1) >= MIN_SV]
    ev["class"] = ev["type"].map(SV_CLASS)
    return ev


def side_bed(ev, side, path):
    s = ev[ev[f"{side}_chr"] != "-"]
    bed = pd.DataFrame({
        "chr": s[f"{side}_chr"],
        "start": (s[f"{side}_start"].astype(int) - 1 - PAD).clip(lower=0),
        "end": s[f"{side}_end"].astype(int) + PAD,
        "uid": s["uid"],
    })
    bed.sort_values(["chr", "start"]).to_csv(path, sep="\t", header=False, index=False)
    return path


def merged_mask(ev, side, path, tmp):
    raw = side_bed(ev, side, Path(tmp) / f"{path.name}.raw")
    run(f"{BEDTOOLS} merge -i {raw} > {path}")
    return path


def covered_fraction(a_bed, mask, genome):
    """uid -> fraction of the A interval covered by the merged mask."""
    out = run(f"{BEDTOOLS} coverage -sorted -g {genome} -a {a_bed} -b {mask}")
    rows = [l.split("\t") for l in out.strip().split("\n") if l]
    return {r[3]: float(r[-1]) for r in rows}


def sort_genome(fai, tmp, name):
    g = Path(tmp) / f"{name}.genome"
    run(f"cut -f1,2 {fai} | sort -k1,1 > {g}")
    return g


def sort_bed(path, genome):
    run(f"{BEDTOOLS} sort -g {genome} -i {path} > {path}.s && mv {path}.s {path}")
    return path


# =============================================================================
# 1. Genome-wide scale
# =============================================================================
def parse_summary(comp):
    rows, section = [], None
    for line in open(SYRI / f"{comp}.syri.summary"):
        line = line.rstrip("\n")
        if line.startswith("#Structural"):
            section = "structural"
        elif line.startswith("#Sequence"):
            section = "sequence"
        elif line and not line.startswith("#") and section:
            f = line.split("\t")
            rows.append({"section": section, "variation": f[0], "count": int(f[1]),
                         "len_ref": f[2], "len_qry": f[3] if len(f) > 3 else "-"})
    d = pd.DataFrame(rows)
    d.insert(0, "comparison", comp)
    return d


def genome_scale():
    d = pd.concat([parse_summary(c) for c in COMPARISONS])
    for c in ("len_ref", "len_qry"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    # one span per row: reference side, query side for query-only rows
    d["bp"] = d["len_ref"].fillna(d["len_qry"])
    wide = d.pivot_table(index=["section", "variation"], columns="comparison",
                         values=["count", "bp"], aggfunc="first")
    wide.columns = [f"{v}_{c}" for v, c in wide.columns]
    wide = wide.reset_index()
    for v in ("count", "bp"):
        wide[f"{v}_lean_ref_over_lean_sisco"] = (wide[f"{v}_lean_vs_ref"] /
                                                 wide[f"{v}_lean_vs_sisco"]).round(3)
        wide[f"{v}_sisco_ref_over_lean_sisco"] = (wide[f"{v}_sisco_vs_ref"] /
                                                  wide[f"{v}_lean_vs_sisco"]).round(3)
    wide.to_csv(OUT / "genome_scale.tsv", sep="\t", index=False)
    return wide


# =============================================================================
# 2. Polarize lean-siscowet SVs
# =============================================================================
def classify(row):
    l, s = row["lean_differs"], row["sisco_differs"]
    if pd.isna(l) or pd.isna(s):
        known = "lean" if pd.notna(l) else "siscowet"
        val = l if pd.notna(l) else s
        return f"{known}_only_side:{'differs' if val else 'matches'}"
    if l and s:
        return "both"
    if l:
        return "lean_side"
    if s:
        return "siscowet_side"
    return "neither"


def polarize(ev_ls, ev_lr, ev_sr, g_lean, g_sisco, tmp, seed=None):
    tmp = Path(tmp)
    mask_l = sort_bed(merged_mask(ev_lr, "qry", tmp / "mask_lean.bed", tmp), g_lean)
    mask_s = sort_bed(merged_mask(ev_sr, "qry", tmp / "mask_sisco.bed", tmp), g_sisco)
    a_l = side_bed(ev_ls, "ref", tmp / "ls_lean.bed")
    a_s = side_bed(ev_ls, "qry", tmp / "ls_sisco.bed")
    if seed is not None:
        for a, g in ((a_l, g_lean), (a_s, g_sisco)):
            run(f"{BEDTOOLS} shuffle -chrom -seed {seed} -i {a} -g {g} > {a}.sh && mv {a}.sh {a}")
    sort_bed(a_l, g_lean)
    sort_bed(a_s, g_sisco)
    fl = covered_fraction(a_l, mask_l, g_lean)
    fs = covered_fraction(a_s, mask_s, g_sisco)
    p = ev_ls[["uid", "type", "class", "ref_len", "qry_len"]].copy()
    p["lean_cov"] = p["uid"].map(fl)
    p["sisco_cov"] = p["uid"].map(fs)
    p["lean_differs"] = (p["lean_cov"] >= MIN_COV).where(p["lean_cov"].notna())
    p["sisco_differs"] = (p["sisco_cov"] >= MIN_COV).where(p["sisco_cov"].notna())
    p["polarity"] = p.apply(classify, axis=1)
    p["bp"] = p[["ref_len", "qry_len"]].max(axis=1)
    return p


def tabulate(p):
    t = (p.groupby(["class", "polarity"])
           .agg(n=("uid", "size"), bp=("bp", "sum")).reset_index())
    tot = p.groupby("class").agg(n_class=("uid", "size"), bp_class=("bp", "sum")).reset_index()
    t = t.merge(tot, on="class")
    t["frac_n"] = (t["n"] / t["n_class"]).round(4)
    t["frac_bp"] = (t["bp"] / t["bp_class"]).round(4)
    allc = (p.groupby("polarity").agg(n=("uid", "size"), bp=("bp", "sum")).reset_index())
    allc.insert(0, "class", "ALL")
    allc["n_class"], allc["bp_class"] = len(p), p["bp"].sum()
    allc["frac_n"] = (allc["n"] / allc["n_class"]).round(4)
    allc["frac_bp"] = (allc["bp"] / allc["bp_class"]).round(4)
    return pd.concat([t, allc], ignore_index=True)


def mask_fraction(ev, side, genome, tmp, name):
    m = sort_bed(merged_mask(ev, side, Path(tmp) / f"{name}.bed", tmp), genome)
    covered = int(run(f"awk '{{s+=$3-$2}} END{{print s}}' {m}").strip())
    total = int(run(f"awk '{{s+=$2}} END{{print s}}' {genome}").strip())
    return covered / total


# =============================================================================
# 3. Reciprocal / corroborated genes
# =============================================================================
def gene_alignment(genes_path, lr, sr, tmp):
    g = pd.read_csv(genes_path, sep="\t", dtype={"reciprocal_present_in": str,
                                                  "ecotype_only_side": str})
    g["group"] = None
    rec = g["reciprocal_pav"] == 1
    g.loc[rec, "group"] = "reciprocal_present_in_" + g.loc[rec, "reciprocal_present_in"]
    cor = (g["ecotype_only_corroborated"] == 1) & ~rec
    g.loc[cor, "group"] = "corroborated_" + g.loc[cor, "ecotype_only_side"] + "_only"
    g.loc[g["group"].isna(), "group"] = "all_other_genes"
    tmp = Path(tmp)
    ref_g = sort_genome(ANCHOR / "ref.chroms.fa.fai", tmp, "ref")
    gb = tmp / "genes.bed"
    g42 = g[g["chrom"].isin(set(pd.read_csv(ref_g, sep="\t", header=None)[0]))]
    g42[["chrom", "start", "end", "gene_id"]].assign(start=lambda d: d["start"] - 1) \
        .to_csv(gb, sep="\t", header=False, index=False)
    sort_bed(gb, ref_g)
    res = g42[["gene_id", "symbol", "product", "chrom", "start", "end", "group", "tier"]].copy()
    for name, df in (("lean", lr), ("siscowet", sr)):
        al = df[df["type"].isin(ALIGNED) & (df["ref_chr"] != "-")]
        b = tmp / f"aligned_{name}.bed"
        pd.DataFrame({"c": al["ref_chr"], "s": al["ref_start"].astype(int) - 1,
                      "e": al["ref_end"].astype(int)}).to_csv(b, sep="\t", header=False,
                                                               index=False)
        sort_bed(b, ref_g)
        run(f"{BEDTOOLS} merge -i {b} > {b}.m && mv {b}.m {b}")
        res[f"{name}_asm_aligned_frac"] = res["gene_id"].map(covered_fraction(gb, b, ref_g))
    res.to_csv(OUT / "gene_alignment.tsv", sep="\t", index=False)
    s = (res.groupby("group")
            .agg(genes=("gene_id", "size"),
                 lean_aligned_median=("lean_asm_aligned_frac", "median"),
                 siscowet_aligned_median=("siscowet_asm_aligned_frac", "median"),
                 lean_aligned_lt_0_2=("lean_asm_aligned_frac", lambda x: (x < 0.2).mean()),
                 siscowet_aligned_lt_0_2=("siscowet_asm_aligned_frac",
                                          lambda x: (x < 0.2).mean()))
            .round(3).reset_index())
    s.to_csv(OUT / "gene_alignment_summary.tsv", sep="\t", index=False)
    return s


# =============================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", default=str(GENES))
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    print("1. genome-wide scale")
    scale = genome_scale()
    print(scale.to_string(index=False))

    print("2. polarizing lean-siscowet SVs")
    syri = {c: read_syri(c) for c in COMPARISONS}
    ev = {c: sv_events(syri[c]) for c in COMPARISONS}
    with tempfile.TemporaryDirectory(dir=OUT) as tmp:
        g_lean = sort_genome(ANCHOR / "lean.anchored.fa.fai", tmp, "lean")
        g_sisco = sort_genome(ANCHOR / "siscowet.anchored.fa.fai", tmp, "sisco")
        frac = {
            "lean_genome_in_lean_vs_ref_SV": mask_fraction(ev["lean_vs_ref"], "qry", g_lean,
                                                           tmp, "fl"),
            "siscowet_genome_in_sisco_vs_ref_SV": mask_fraction(ev["sisco_vs_ref"], "qry",
                                                                g_sisco, tmp, "fs"),
        }
        p = polarize(ev["lean_vs_sisco"], ev["lean_vs_ref"], ev["sisco_vs_ref"],
                     g_lean, g_sisco, tmp)
        obs = tabulate(p)
        obs.to_csv(OUT / "polarized_sv.tsv", sep="\t", index=False)
        nulls = []
        for seed in range(1, N_SHUF + 1):
            t = tabulate(polarize(ev["lean_vs_sisco"], ev["lean_vs_ref"], ev["sisco_vs_ref"],
                                  g_lean, g_sisco, tmp, seed=seed))
            t.insert(0, "seed", seed)
            nulls.append(t)
        null = pd.concat(nulls)
        null.to_csv(OUT / "polarized_sv_null.tsv", sep="\t", index=False)

        print("3. reciprocal / corroborated genes")
        gs = gene_alignment(args.genes, syri["lean_vs_ref"], syri["sisco_vs_ref"], tmp)

    pd.Series(frac).round(4).to_csv(OUT / "mask_fractions.tsv", sep="\t", header=["fraction"])
    print(pd.Series(frac).round(4).to_string())
    print(obs[obs["class"] == "ALL"].to_string(index=False))
    print(null[null["class"] == "ALL"].groupby("polarity")[["frac_n", "frac_bp"]].mean()
              .round(4).to_string())
    print(gs.to_string(index=False))


if __name__ == "__main__":
    main()
