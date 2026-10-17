#!/usr/bin/env python3
"""
23.19-merge-annotation.py  (Phase 4, step 4, of code/23-next-phase-research-plan.md)

MERGE LIFTOFF (Track A) WITH BRAKER3 (Track B) PER ECOTYPE. Rule from plan 20 / 23 Phase 4:
prefer the lifted model where the two agree, and add BRAKER models only where no lifted gene
exists, tagging provenance.

Each BRAKER gene is classed (bedtools intersect, any overlap):
  shared         a CDS overlaps a Liftoff CDS on the same strand. The Liftoff model is kept;
                 the match (reference gene IDs, CDS Jaccard) is recorded as agreement QC.
  liftoff_locus  no same-strand CDS overlap, but it overlaps a Liftoff gene of any biotype on
                 either strand (intronic, antisense, ncRNA host). Reported, not added.
  new            overlaps no Liftoff gene. Candidate ecotype-specific / reference-absent gene;
                 added to the merged set with provenance=braker.

Evidence and quality for every BRAKER transcript (the best transcript represents the gene):
  introns matching an RNA-seq (b2h) or ProtHint (src P / C / M) intron hint exactly;
  fraction of CDS bases covered by ProtHint CDSpart hints; protein length; fraction of CDS
  bases soft-masked (repeats) in the RepeatMasker output.
  new genes are hq when aa >= MIN_AA, repeat fraction < MAX_REPEAT, and (>= 1 supported intron
  or ProtHint CDS coverage >= MIN_PROT_COV); otherwise lq. Both are written; quality is an
  attribute.

Reference homology of new genes: DIAMOND blastp (sensitive, e <= 1e-5, query cover >= 50 %) of the
new proteins against the SaNama_1.0 RefSeq proteins. ref_homolog = best hit >= 90 % identity (most
likely a reference gene Liftoff failed to place, not an ecotype-specific gene); ref_paralog = a
weaker hit (diverged copy or family member); no_ref_hit = no hit (the novel candidates).

Outputs (analyses/23-denovo-annotation/final/{eco}/):
  {eco}.merged.gff3              Liftoff features (gene lines tagged provenance=liftoff) + new
                                 BRAKER genes (IDs {eco}_braker_<gN>, provenance=braker, quality)
  {eco}.merged.genes.bed         gene BED (chrom, start, end, ID, provenance, strand)
  {eco}.braker_new.gff3 / .aa / .cds.fa   new BRAKER genes only
  {eco}.merged.proteins.fa       Liftoff proteins (gffread -y on the purged assembly) + new BRAKER
                                 proteins, for BUSCO / functional annotation
  braker_gene_classes.tsv        per BRAKER gene: class, overlaps, evidence, quality
  merge_summary.txt

Usage (raven):
  /home/sr320/miniconda3/bin/python code/23.19-merge-annotation.py [lean siscowet]
  smoke test: ... --p4 analyses/23-denovo-annotation/smoke --genome-dir analyses/23-denovo-annotation/smoke/genome

Author: Generated for project-lake-trout
Date: 2026-10-17
"""

import argparse
import io
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent.parent
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
BEDTOOLS = "/usr/bin/bedtools"
SAMTOOLS = "/usr/bin/samtools"
GFFREAD = "/home/shared/gffread-0.12.7.Linux_x86_64/gffread"
DIAMOND = "/home/shared/diamond-2.1.8"
REF_PROTEINS = "analyses/23-denovo-annotation/evidence/proteins/SaNama_1.0.refseq.protein.faa"
REF_HOMOLOG_PID = 90.0

MIN_AA = 50
MAX_REPEAT = 0.5
MIN_PROT_COV = 0.5


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=True, check=True, executable="/bin/bash", **kw)


def attrs(s):
    return dict(kv.split("=", 1) for kv in s.strip().strip(";").split(";") if "=" in kv)


def read_gff(path, contigs=None):
    cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    df = pd.read_csv(path, sep="\t", comment="#", header=None, names=cols, dtype={"seqid": str},
                     quoting=3)
    if contigs is not None:
        df = df[df["seqid"].isin(contigs)]
    df["start0"] = df["start"] - 1
    a = df["attributes"].map(attrs)
    df["ID"] = a.map(lambda d: d.get("ID", ""))
    df["Parent"] = a.map(lambda d: d.get("Parent", ""))
    return df.reset_index(drop=True)


def to_bed(df, path, name_col, extra=()):
    out = df[["seqid", "start0", "end", name_col, *extra, "strand"]] if extra else \
        df.assign(_s=0)[["seqid", "start0", "end", name_col, "_s", "strand"]]
    out.sort_values(["seqid", "start0"]).to_csv(path, sep="\t", header=False, index=False)


# =============================================================================
def liftoff_tables(gff, work):
    """Liftoff gene spans (all biotypes) and CDS intervals keyed to reference gene IDs."""
    genes = gff[gff["type"] == "gene"].copy()
    tx_gene = {}
    for r in gff[gff["Parent"].isin(set(genes["ID"]))].itertuples():
        tx_gene[r.ID] = r.Parent
    cds = gff[gff["type"] == "CDS"].copy()
    cds["gene"] = cds["Parent"].map(lambda p: tx_gene.get(p.split(",")[0], ""))
    genes_bed, cds_bed = work / "liftoff.genes.bed", work / "liftoff.cds.bed"
    to_bed(genes, genes_bed, "ID")
    to_bed(cds, cds_bed, "gene")
    return genes, cds, genes_bed, cds_bed


def braker_tables(gff):
    genes = gff[gff["type"] == "gene"].copy()
    mrna = gff[gff["type"].isin(["mRNA", "transcript"])].copy()
    mrna["gene"] = mrna["Parent"]
    cds = gff[gff["type"] == "CDS"].copy()
    cds["tx"] = cds["Parent"]
    cds["gene"] = cds["tx"].map(dict(zip(mrna["ID"], mrna["gene"])))
    intr = gff[gff["type"] == "intron"].copy()
    intr["tx"] = intr["Parent"]
    return genes, mrna, cds, intr


def hint_support(hints_path, cds, intr):
    """Per transcript: introns matching RNA / protein intron hints; ProtHint CDSpart coverage."""
    h = pd.read_csv(hints_path, sep="\t", header=None, quoting=3,
                    names=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"])
    # RNA-seq hints: source b2h (src=E). Protein hints: source ProtHint, any src (P, C chained,
    # M high-confidence); CDSpart hints are src=C.
    ih = h[h["type"] == "intron"]
    r, p = ih[ih["source"] == "b2h"], ih[ih["source"] == "ProtHint"]
    rna_int = set(zip(r["seqid"], r["start"], r["end"]))
    prot_int = set(zip(p["seqid"], p["start"], p["end"]))
    t = intr.assign(key=list(zip(intr["seqid"], intr["start"], intr["end"])))
    t["rna"] = t["key"].map(lambda k: k in rna_int)
    t["prot"] = t["key"].map(lambda k: k in prot_int)
    sup = t.groupby("tx").agg(n_introns=("key", "size"), introns_rna=("rna", "sum"),
                              introns_prot=("prot", "sum"),
                              introns_any=("key", lambda k: int(sum((x in rna_int) or (x in prot_int)
                                                                   for x in k))))
    # ProtHint CDSpart coverage of each transcript's CDS
    cp = h[(h["type"] == "CDSpart") & (h["source"] == "ProtHint")]
    return sup, cp


def cds_cover(cds, intervals, work, tag):
    """Fraction of each transcript's CDS bases covered by the given intervals (merged)."""
    a, b = work / f"_{tag}.cds.bed", work / f"_{tag}.iv.bed"
    cds[["seqid", "start0", "end", "tx"]].sort_values(["seqid", "start0"]).to_csv(
        a, sep="\t", header=False, index=False)
    iv = intervals.assign(start0=intervals["start"] - 1)[["seqid", "start0", "end"]]
    iv.sort_values(["seqid", "start0"]).to_csv(b, sep="\t", header=False, index=False)
    out = sh(f"LC_ALL=C sort -k1,1 -k2,2n '{b}' | '{BEDTOOLS}' merge -i - > '{b}.m' && "
             f"'{BEDTOOLS}' coverage -a '{a}' -b '{b}.m'", capture_output=True, text=True).stdout
    c = pd.read_csv(io.StringIO(out), sep="\t", header=None)
    c.columns = ["seqid", "s", "e", "tx", "n", "covered", "len", "frac"]
    g = c.groupby("tx")[["covered", "len"]].sum()
    return (g["covered"] / g["len"]).rename(f"{tag}_cov")


def repeat_fraction(cds, masked_fa, work):
    """Fraction of CDS bases soft-masked (lowercase) per transcript."""
    reg = work / "_cds.regions.txt"
    (cds["seqid"] + ":" + cds["start"].astype(str) + "-" + cds["end"].astype(str)).to_csv(
        reg, index=False, header=False)
    out = sh(f"'{SAMTOOLS}' faidx '{masked_fa}' -r '{reg}'", capture_output=True, text=True).stdout
    seqs, cur = [], []
    for line in out.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line)
    if cur:
        seqs.append("".join(cur))
    low = np.array([sum(ch.islower() for ch in s) for s in seqs])
    ln = np.array([len(s) for s in seqs])
    d = pd.DataFrame({"tx": cds["tx"].to_numpy(), "low": low, "len": ln}).groupby("tx").sum()
    return (d["low"] / d["len"]).rename("repeat_frac")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


# =============================================================================
def merge_one(eco, p4, genome_dir):
    out = p4 / "final" / eco
    out.mkdir(parents=True, exist_ok=True)
    work = out / "_work"
    work.mkdir(exist_ok=True)
    fa = genome_dir / f"{eco}.purged.fa"
    masked = p4 / "repeats" / eco / f"{eco}.purged.fa.masked"
    bdir = p4 / "braker" / eco
    contigs = set(pd.read_csv(f"{fa}.fai", sep="\t", header=None, usecols=[0])[0]) \
        if Path(f"{fa}.fai").exists() else None
    if contigs is None:
        sh(f"'{SAMTOOLS}' faidx '{fa}'")
        contigs = set(pd.read_csv(f"{fa}.fai", sep="\t", header=None, usecols=[0])[0])

    lo = read_gff(SYN / f"{eco}.purged.liftoff.gff3_polished", contigs)
    lgenes, lcds, lgenes_bed, lcds_bed = liftoff_tables(lo, work)
    bg = read_gff(bdir / "braker.gff3")
    bgenes, bmrna, bcds, bintr = braker_tables(bg)

    # ---- overlaps --------------------------------------------------------------------------------
    b_cds_bed = work / "braker.cds.bed"
    to_bed(bcds.assign(name=bcds["gene"]), b_cds_bed, "name")
    ov = sh(f"LC_ALL=C sort -k1,1 -k2,2n '{b_cds_bed}' | '{BEDTOOLS}' intersect -a - "
            f"-b <(LC_ALL=C sort -k1,1 -k2,2n '{lcds_bed}') -s -wo", capture_output=True, text=True).stdout
    shared = defaultdict(set)
    for line in ov.splitlines():
        p = line.split("\t")
        if p[9]:
            shared[p[3]].add(p[9])
    b_gene_bed = work / "braker.genes.bed"
    to_bed(bgenes, b_gene_bed, "ID")
    ov2 = sh(f"LC_ALL=C sort -k1,1 -k2,2n '{b_gene_bed}' | '{BEDTOOLS}' intersect -a - "
             f"-b <(LC_ALL=C sort -k1,1 -k2,2n '{lgenes_bed}') -wo", capture_output=True, text=True).stdout
    locus = defaultdict(set)
    for line in ov2.splitlines():
        p = line.split("\t")
        locus[p[3]].add(p[9])

    # CDS-level Jaccard for shared genes (BRAKER gene CDS union vs matched Liftoff genes' CDS union)
    def union_len(iv):
        iv = sorted(iv)
        tot, cs, ce = 0, None, None
        for s, e in iv:
            if cs is None or s > ce:
                if cs is not None:
                    tot += ce - cs
                cs, ce = s, e
            else:
                ce = max(ce, e)
        return tot + (ce - cs if cs is not None else 0)

    lcds_by_gene = {g: list(zip(x["start0"], x["end"])) for g, x in lcds.groupby("gene")}
    bcds_by_gene = {g: list(zip(x["start0"], x["end"])) for g, x in bcds.groupby("gene")}

    def jaccard(bgene, lgs):
        b = bcds_by_gene.get(bgene, [])
        l = [iv for g in lgs for iv in lcds_by_gene.get(g, [])]
        ub, ul, uu = union_len(b), union_len(l), union_len(b + l)
        return (ub + ul - uu) / uu if uu else np.nan

    # ---- evidence ---------------------------------------------------------------------------------
    sup, cdspart = hint_support(bdir / "hintsfile.gff", bcds, bintr)
    prot_cov = cds_cover(bcds, cdspart, work, "prothint") if len(cdspart) else pd.Series(dtype=float)
    rep = repeat_fraction(bcds, masked, work)
    aa = read_fasta(bdir / "braker.aa")
    cdsseq = read_fasta(bdir / "braker.codingseq")

    tx = bmrna[["ID", "gene", "seqid", "start", "end", "strand"]].rename(columns={"ID": "tx"})
    tx = tx.merge(sup, left_on="tx", right_index=True, how="left").fillna(
        {"n_introns": 0, "introns_rna": 0, "introns_prot": 0, "introns_any": 0})
    tx["prothint_cov"] = tx["tx"].map(prot_cov).fillna(0.0)
    tx["repeat_frac"] = tx["tx"].map(rep).fillna(0.0)
    tx["aa_len"] = tx["tx"].map(lambda t: len(aa.get(t, "").rstrip("*")))
    tx["evidence_score"] = tx["introns_any"] / tx["n_introns"].clip(lower=1) + tx["prothint_cov"]
    best = tx.sort_values(["gene", "evidence_score", "aa_len"], ascending=[True, False, False]) \
             .drop_duplicates("gene").set_index("gene")

    rows = []
    for g in bgenes.itertuples():
        b = best.loc[g.ID] if g.ID in best.index else None
        cls = "shared" if g.ID in shared else ("liftoff_locus" if g.ID in locus else "new")
        row = {"braker_gene": g.ID, "seqid": g.seqid, "start": g.start, "end": g.end,
               "strand": g.strand, "class": cls,
               "liftoff_cds_genes": ",".join(sorted(shared.get(g.ID, ()))),
               "liftoff_locus_genes": ",".join(sorted(locus.get(g.ID, set()) - shared.get(g.ID, set()))),
               "cds_jaccard": round(jaccard(g.ID, shared[g.ID]), 3) if g.ID in shared else np.nan}
        if b is not None:
            row.update({"best_tx": b["tx"], "aa_len": int(b["aa_len"]), "n_introns": int(b["n_introns"]),
                        "introns_rna": int(b["introns_rna"]), "introns_prot": int(b["introns_prot"]),
                        "introns_any": int(b["introns_any"]), "prothint_cov": round(b["prothint_cov"], 3),
                        "repeat_frac": round(b["repeat_frac"], 3)})
        rows.append(row)
    cl = pd.DataFrame(rows)
    hq = (cl["aa_len"] >= MIN_AA) & (cl["repeat_frac"] < MAX_REPEAT) & \
         ((cl["introns_any"] >= 1) | (cl["prothint_cov"] >= MIN_PROT_COV))
    cl["quality"] = np.where(cl["class"] != "new", "", np.where(hq, "hq", "lq"))
    cl["mono_exonic"] = cl["n_introns"] == 0
    cl.to_csv(out / "braker_gene_classes.tsv", sep="\t", index=False)

    # ---- write merged GFF3 + FASTA ---------------------------------------------------------------
    new = cl[cl["class"] == "new"].set_index("braker_gene")
    rename = {g: f"{eco}_braker_{g}" for g in new.index}
    keep_tx = set(bmrna.loc[bmrna["gene"].isin(new.index), "ID"])

    def new_attr(r):
        a = r.attributes
        for old, nw in rename.items():
            a = re.sub(rf"(ID|Parent)={re.escape(old)}(?=[.;]|$)", rf"\1={nw}", a)
        if r.type == "gene":
            a = a.rstrip(";") + f";provenance=braker;quality={new.loc[r.ID, 'quality']}"
        return a

    bsel = bg[(bg["type"] == "gene") & bg["ID"].isin(new.index) |
              bg["Parent"].isin(new.index) | bg["Parent"].isin(keep_tx)].copy()
    bsel["attributes"] = [new_attr(r) for r in bsel.itertuples()]
    bsel["source"] = "BRAKER3_" + bsel["source"]
    cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    lo_out = lo.copy()
    lo_out.loc[lo_out["type"] == "gene", "attributes"] = \
        lo_out.loc[lo_out["type"] == "gene", "attributes"].str.rstrip(";") + ";provenance=liftoff"
    merged = pd.concat([lo_out[cols + ["start0"]], bsel[cols + ["start0"]]])
    order = {"gene": 0, "mRNA": 1, "transcript": 1}
    merged["_o"] = merged["type"].map(order).fillna(2)
    merged = merged.sort_values(["seqid", "start", "_o"], kind="stable")
    with open(out / f"{eco}.merged.gff3", "w") as f:
        f.write("##gff-version 3\n")
        merged[cols].to_csv(f, sep="\t", header=False, index=False, quoting=3)
    with open(out / f"{eco}.braker_new.gff3", "w") as f:
        f.write("##gff-version 3\n")
        bsel.sort_values(["seqid", "start"])[cols].to_csv(f, sep="\t", header=False, index=False, quoting=3)
    gl = merged[merged["type"] == "gene"].copy()
    gl["ID"] = gl["attributes"].map(lambda s: attrs(s).get("ID", ""))
    gl["prov"] = np.where(gl["source"].str.startswith("BRAKER3"), "braker", "liftoff")
    gl[["seqid", "start0", "end", "ID", "prov", "strand"]].to_csv(
        out / f"{eco}.merged.genes.bed", sep="\t", header=False, index=False)

    with open(out / f"{eco}.braker_new.aa", "w") as fa_aa, open(out / f"{eco}.braker_new.cds.fa", "w") as fa_cds:
        for g, r in new.iterrows():
            t = r.get("best_tx")
            if isinstance(t, str) and t in aa:
                nt = t.replace(g, rename[g], 1)
                fa_aa.write(f">{nt} provenance=braker quality={r['quality']}\n{aa[t]}\n")
                if t in cdsseq:
                    fa_cds.write(f">{nt} provenance=braker quality={r['quality']}\n{cdsseq[t]}\n")
    lift_aa, lift_sub = work / "liftoff.proteins.fa", work / "liftoff.subset.gff3"
    with open(lift_sub, "w") as f:          # only contigs present in the genome FASTA
        f.write("##gff-version 3\n")
        lo[cols].to_csv(f, sep="\t", header=False, index=False, quoting=3)
    sh(f"'{GFFREAD}' '{lift_sub}' -g '{fa}' -y '{lift_aa}' 2> '{work / 'gffread.log'}'")
    with open(out / f"{eco}.merged.proteins.fa", "w") as f:
        for name, seq in read_fasta(lift_aa).items():
            if set(seq) - {"."}:
                f.write(f">{name} provenance=liftoff\n{seq.replace('.', '*')}\n")
        for name, seq in read_fasta(out / f"{eco}.braker_new.aa").items():
            f.write(f">{name}\n{seq}\n")

    # ---- reference homology of new genes ---------------------------------------------------------
    ref = BASE / REF_PROTEINS
    hits = work / "new_vs_ref.tsv"
    if ref.exists() and (out / f"{eco}.braker_new.aa").stat().st_size:
        db = work / "sanama_ref"
        sh(f"'{DIAMOND}' makedb --in '{ref}' -d '{db}' --quiet")
        sh(f"'{DIAMOND}' blastp -q '{out / f'{eco}.braker_new.aa'}' -d '{db}' -o '{hits}' --sensitive "
           f"--evalue 1e-5 --query-cover 50 --max-target-seqs 1 --threads 8 --quiet "
           f"--outfmt 6 qseqid sseqid pident qcovhsp evalue")
        hh = pd.read_csv(hits, sep="\t", header=None, names=["q", "ref_hit", "ref_pident", "ref_qcov", "ref_evalue"]) \
            if hits.stat().st_size else pd.DataFrame(columns=["q", "ref_hit", "ref_pident", "ref_qcov", "ref_evalue"])
        hh = hh.drop_duplicates("q")
        inv = {v: k for k, v in rename.items()}
        hh["braker_gene"] = hh["q"].map(lambda q: inv.get(q.rsplit(".", 1)[0], ""))
        cl = cl.merge(hh.drop(columns="q"), on="braker_gene", how="left")
        cl["ref_status"] = np.where(cl["class"] != "new", "",
                                    np.where(cl["ref_pident"] >= REF_HOMOLOG_PID, "ref_homolog",
                                             np.where(cl["ref_pident"].notna(), "ref_paralog", "no_ref_hit")))
        cl.to_csv(out / "braker_gene_classes.tsv", sep="\t", index=False)

    # ---- summary ----------------------------------------------------------------------------------
    pc = lgenes[lgenes["attributes"].str.contains("gene_biotype=protein_coding")]
    matched_lift = set().union(*shared.values()) if shared else set()
    vc = cl["class"].value_counts()
    n = cl[cl["class"] == "new"]
    lines = [f"{eco}: Liftoff genes {len(lgenes):,} (protein-coding {len(pc):,}); BRAKER genes {len(cl):,}",
             f"  BRAKER classes: shared {vc.get('shared', 0):,}, liftoff_locus {vc.get('liftoff_locus', 0):,}, "
             f"new {vc.get('new', 0):,}",
             f"  agreement: {len(matched_lift & set(pc['ID'])):,} of {len(pc):,} protein-coding Liftoff genes "
             f"({len(matched_lift & set(pc['ID'])) / max(1, len(pc)):.1%}) have a same-strand BRAKER CDS "
             f"overlap; median CDS Jaccard of shared genes {cl['cds_jaccard'].median():.2f}",
             f"  new: hq {int((n['quality'] == 'hq').sum()):,}, lq {int((n['quality'] == 'lq').sum()):,}; "
             f"mono-exonic {n['mono_exonic'].mean():.1%}; median aa {n['aa_len'].median():.0f}; "
             f"repeat_frac >= {MAX_REPEAT} {(n['repeat_frac'] >= MAX_REPEAT).mean():.1%}",
             f"  BRAKER overall: mono-exonic {cl['mono_exonic'].mean():.1%}; with >= 1 RNA-supported intron "
             f"{(cl['introns_rna'] >= 1).mean():.1%}; ProtHint CDS cov >= {MIN_PROT_COV} "
             f"{(cl['prothint_cov'] >= MIN_PROT_COV).mean():.1%}",
             f"  new by reference homology (SaNama RefSeq proteins): " + (", ".join(
                 f"{k} {int(v)} (hq {int(((n2 := cl[cl['class'] == 'new'])['ref_status'] == k)[n2['quality'] == 'hq'].sum())})"
                 for k, v in cl.loc[cl["class"] == "new", "ref_status"].value_counts().items())
                 if "ref_status" in cl else "not run"),
             f"  merged gene set: {len(gl):,} genes (liftoff {int((gl['prov'] == 'liftoff').sum()):,} + "
             f"braker {int((gl['prov'] == 'braker').sum()):,})"]
    (out / "merge_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines), flush=True)
    for f in work.glob("_*"):
        f.unlink()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ecotypes", nargs="*", default=["lean", "siscowet"])
    ap.add_argument("--p4", type=Path, default=BASE / "analyses" / "23-denovo-annotation")
    ap.add_argument("--genome-dir", type=Path, default=SYN)
    a = ap.parse_args()
    p4 = a.p4 if a.p4.is_absolute() else BASE / a.p4
    gd = a.genome_dir if a.genome_dir.is_absolute() else BASE / a.genome_dir
    for eco in a.ecotypes:
        if not (p4 / "braker" / eco / "braker.gff3").exists():
            print(f"{eco}: no BRAKER output yet in {p4 / 'braker' / eco}; skipped")
            continue
        merge_one(eco, p4, gd)


if __name__ == "__main__":
    main()
