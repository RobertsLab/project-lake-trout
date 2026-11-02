#!/usr/bin/env python3
"""
23.23 -- Target-ecotype read coverage of the reciprocal-PAV "present" regions.

Why
---
23.4 calls a region of ecotype X's purged assembly "X-specific present" when all four fish of the
OTHER ecotype have zero coverage there and the flanks are covered. It never asks whether X's own
four fish cover the region. Both assemblies were built from all four fish of their ecotype pooled,
so a region in X's assembly could be (a) fixed in X (true ecotype-specific sequence), (b) carried by
only some X fish (a within-ecotype polymorphism), (c) a retained haplotig / collapsed duplicate (all
X fish at roughly half depth), or (d) an assembly artifact or contaminant (no X fish covers it).
Every tier-A gene in 23.7 requires a reciprocal-PAV call, so the distinction matters.

What
----
For each direction (siscowet-present regions on the siscowet assembly, lean-present regions on the
lean assembly) take the four target-ecotype fish aligned to their OWN purged assembly in Phase 3
(23.8, `analyses/23-native-methylation/align/{sample}.to_own.sorted.bam`) and compute, per region
and per fish, mean depth and the fraction of bases covered at >= 1x and >= 3x. Each region's 2 kb
flanks (the same flank width 23.4 required to be covered in the querying ecotype) give a local
depth baseline, so the region/flank ratio separates single-copy (~1), half-depth (~0.5) and absent
(~0) without needing a genome-wide median per fish.

Per fish a region is   present  if frac>=1x >= 0.80 and ratio >= 0.60
                       low      if frac>=1x >= 0.20 (covered, but thin or half-depth)
                       absent   otherwise (frac>=1x < 0.20).
Per region the class is
    fixed_present   all 4 target fish present
    polymorphic     >= 1 absent and >= 1 present
    low_all         no fish present, none absent (all low: haplotig / duplicate-like)
    absent_all      all 4 absent (assembly artifact / contamination / mis-assembly)
    mixed           anything else (present + low, or low + absent)

Pooled columns sum the four fish (pooled_ratio = summed region depth / summed flank depth;
pooled_frac1x = fraction of bases covered by at least one read from any fish), which is the
"does the ecotype's read set cover this at single-copy depth" question independent of per-fish noise.

Null
----
Per-fish depth is only 2-9x, so a fish can miss a 1 kb region by chance. `--null N` shuffles the
regions N times across the same assembly (bedtools shuffle, same contig, same lengths) and runs the
identical classification, giving the chance rate of each class and of per-fish "absent" calls.

Genes are joined exactly as in 23.4 (bedtools intersect with the native Liftoff gene BED) and
each gene keeps its best region class. The tier-A list from 23.7 is then annotated with the gene's
class in the ecotype where it is "present".

Usage
-----
    python3 code/23.23-reciprocal-pav-target-coverage.py [--bam-dir DIR] [--threads N] [--null K]

--bam-dir defaults to <repo>/analyses/23-native-methylation/align (the BAMs are not committed; in a
worktree point this at the main checkout).

Outputs: analyses/23-reciprocal-pav/target_coverage/
    {eco}.region_coverage.tsv       one row per region: per-fish depth, frac1x, frac3x, ratio, status; pooled; class
    {eco}.null{k}.region_coverage.tsv  the same for shuffled regions
    {eco}.gene_coverage.tsv         one row per gene overlapping a region, with its best class
    tierA_target_coverage.tsv       the 76 tier-A genes with their class
    target_coverage_summary.txt     counts by class per direction, null rates, tier-A tally
"""
import argparse
import csv
import statistics
import subprocess
import sys
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
SYN = BASE / "output" / "21.1-gene-anchored-synteny"
RPAV = BASE / "analyses" / "23-reciprocal-pav"
OUT = RPAV / "target_coverage"
REFINED = BASE / "analyses" / "23-integration" / "refined" / "gene_evidence_refined.tsv"
SAMTOOLS = "/usr/bin/samtools"
BEDTOOLS = "/usr/bin/bedtools"

FLANK_BP = 2000           # same flank width 23.4 required to be covered in the querying ecotype
PRESENT_FRAC1 = 0.80
PRESENT_RATIO = 0.60
ABSENT_FRAC1 = 0.20

ECO = {
    "siscowet": {"samples": ["bc2071", "bc2072", "bc2073", "bc2096"],
                 "regions": RPAV / "siscowet_specific.present_regions.on_siscowet.bed",
                 "genes": SYN / "siscowet.purged.liftoff.genes.bed",
                 "fai": SYN / "siscowet.purged.fa.fai"},
    "lean":     {"samples": ["bc2041", "bc2068", "bc2069", "bc2070"],
                 "regions": RPAV / "lean_specific.present_regions.on_lean.bed",
                 "genes": SYN / "lean.purged.liftoff.genes.bed",
                 "fai": SYN / "lean.purged.fa.fai"},
}
CLASS_RANK = {"fixed_present": 0, "mixed": 1, "polymorphic": 2, "low_all": 3, "absent_all": 4}


def read_bed(p):
    rows = []
    with open(p) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c, s, e = line.split("\t")[:3]
            rows.append((c, int(s), int(e)))
    return rows


def contig_sizes(fai):
    with open(fai) as f:
        return {l.split("\t")[0]: int(l.split("\t")[1]) for l in f}


def intervals_for(regions, sizes):
    """region + left flank + right flank, labelled; flanks clipped to the contig."""
    iv = []
    for c, s, e in regions:
        rid = f"{c}:{s}-{e}"
        iv.append((c, s, e, rid, "region"))
        ls = max(0, s - FLANK_BP)
        if ls < s:
            iv.append((c, ls, s, rid, "lflank"))
        re_ = min(sizes[c], e + FLANK_BP)
        if re_ > e:
            iv.append((c, e, re_, rid, "rflank"))
    return iv


def depth_for_bam(bam, bed, threads):
    """samtools depth -a over the bed; returns {(contig,pos): depth} for every base in the bed."""
    cmd = [SAMTOOLS, "depth", "-a", "-b", str(bed), "-@", str(threads), str(bam)]
    print("  $", " ".join(cmd), flush=True)
    out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout
    d = {}
    for line in out.splitlines():
        c, p, v = line.split("\t")
        d[(c, int(p))] = int(v)
    return d


def vals(d, c, s, e):
    return [d.get((c, p), 0) for p in range(s + 1, e + 1)]


def summarise(d, c, s, e):
    """mean depth, frac>=1x, frac>=3x over bed interval [s,e)."""
    n = e - s
    if n <= 0:
        return 0.0, 0.0, 0.0
    v = vals(d, c, s, e)
    return (sum(v) / n, sum(x >= 1 for x in v) / n, sum(x >= 3 for x in v) / n)


def fish_status(frac1, ratio):
    if frac1 >= PRESENT_FRAC1 and ratio >= PRESENT_RATIO:
        return "present"
    if frac1 >= ABSENT_FRAC1:
        return "low"
    return "absent"


def region_class(statuses):
    n_p = statuses.count("present"); n_a = statuses.count("absent"); n = len(statuses)
    if n_p == n:
        return "fixed_present"
    if n_a == n:
        return "absent_all"
    if n_a >= 1 and n_p >= 1:
        return "polymorphic"
    if n_p == 0 and n_a == 0:
        return "low_all"
    return "mixed"


def classify(regions, cfg, bam_dir, threads, label):
    """Depth + classification for a list of regions; writes {label}.region_coverage.tsv."""
    sizes = contig_sizes(cfg["fai"])
    iv = intervals_for(regions, sizes)
    work = OUT / "_work"; work.mkdir(parents=True, exist_ok=True)
    bed = work / f"{label}.intervals.bed"
    with open(bed, "w") as f:
        for c, s, e, rid, kind in iv:
            f.write(f"{c}\t{s}\t{e}\t{rid}\t{kind}\n")
    sbed = work / f"{label}.intervals.sorted.bed"
    subprocess.run(f"sort -k1,1 -k2,2n '{bed}' | cut -f1-3 | '{BEDTOOLS}' merge -i - > '{sbed}'",
                   shell=True, check=True)

    bams = {s: bam_dir / f"{s}.to_own.sorted.bam" for s in cfg["samples"]}
    for s, b in bams.items():
        if not b.exists():
            sys.exit(f"FATAL: missing {b}")
    with ThreadPoolExecutor(max_workers=len(bams)) as ex:
        depth = dict(zip(bams, ex.map(lambda s: depth_for_bam(bams[s], sbed, threads), bams)))

    by_rid = defaultdict(dict)
    for c, s, e, rid, kind in iv:
        by_rid[rid][kind] = (c, s, e)

    rows = []
    flank_all = {s: [] for s in cfg["samples"]}
    for c, s, e in regions:
        rid = f"{c}:{s}-{e}"
        row = {"chrom": c, "start": s, "end": e, "length": e - s, "region_id": rid}
        statuses = []
        pooled_mean = pooled_flank = 0.0
        union = [0] * (e - s)
        for smp in cfg["samples"]:
            d = depth[smp]
            mean, f1, f3 = summarise(d, c, s, e)
            for i, x in enumerate(vals(d, c, s, e)):
                if x >= 1:
                    union[i] = 1
            fl = []
            for k in ("lflank", "rflank"):
                if k in by_rid[rid]:
                    fc, fs, fe = by_rid[rid][k]
                    fl.append(summarise(d, fc, fs, fe)[0])
            flank = statistics.mean(fl) if fl else 0.0
            flank_all[smp].append(flank)
            pooled_mean += mean; pooled_flank += flank
            ratio = mean / flank if flank > 0 else (float("inf") if mean > 0 else 0.0)
            st = fish_status(f1, ratio)
            statuses.append(st)
            row.update({f"{smp}_mean": round(mean, 2), f"{smp}_frac1x": round(f1, 3),
                        f"{smp}_frac3x": round(f3, 3), f"{smp}_flank_mean": round(flank, 2),
                        f"{smp}_ratio": round(ratio, 3) if ratio != float("inf") else "inf",
                        f"{smp}_status": st})
        pr = pooled_mean / pooled_flank if pooled_flank > 0 else (float("inf") if pooled_mean > 0 else 0.0)
        row.update({"pooled_mean": round(pooled_mean, 2), "pooled_flank_mean": round(pooled_flank, 2),
                    "pooled_ratio": round(pr, 3) if pr != float("inf") else "inf",
                    "pooled_frac1x": round(sum(union) / max(1, e - s), 3)})
        row["n_present"] = statuses.count("present")
        row["n_low"] = statuses.count("low")
        row["n_absent"] = statuses.count("absent")
        row["region_class"] = region_class(statuses)
        rows.append(row)

    with open(OUT / f"{label}.region_coverage.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    fm = {s: round(statistics.median(v), 1) for s, v in flank_all.items()}
    return rows, fm


def shuffled(regions, cfg, seed, label):
    """bedtools shuffle on the same assembly: same contig, same lengths, no overlap with the real set."""
    work = OUT / "_work"
    src = work / f"{label}.src.bed"
    with open(src, "w") as f:
        for c, s, e in regions:
            f.write(f"{c}\t{s}\t{e}\n")
    g = work / f"{label}.genome"
    subprocess.run(f"cut -f1,2 '{cfg['fai']}' > '{g}'", shell=True, check=True)
    out = work / f"{label}.shuffled.bed"
    subprocess.run(f"'{BEDTOOLS}' shuffle -i '{src}' -g '{g}' -chrom -seed {seed} -noOverlapping "
                   f"-excl '{src}' > '{out}'", shell=True, check=True)
    return read_bed(out)


def genes_for(rows, cfg, label):
    cls = {r["region_id"]: r["region_class"] for r in rows}
    work = OUT / "_work"
    rbed = work / f"{label}.regions.named.bed"
    with open(rbed, "w") as f:
        for r in rows:
            f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['region_id']}\n")
    res = subprocess.run([BEDTOOLS, "intersect", "-a", str(cfg["genes"]), "-b", str(rbed), "-wa", "-wb"],
                         check=True, capture_output=True, text=True).stdout
    gene_regions = defaultdict(list)
    for line in res.splitlines():
        p = line.split("\t")
        gene_regions[p[3]].append(p[-1])
    grows = []
    for g, rids in sorted(gene_regions.items()):
        classes = [cls[r] for r in rids]
        best = min(classes, key=lambda k: CLASS_RANK[k])
        grows.append({"gene_id": g, "n_regions": len(rids), "best_class": best,
                      "classes": ",".join(classes), "regions": ",".join(rids)})
    with open(OUT / f"{label}.gene_coverage.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["gene_id", "n_regions", "best_class", "classes", "regions"],
                           delimiter="\t")
        w.writeheader(); w.writerows(grows)
    return grows


def class_table(rows, grows=None):
    counts = Counter(); bp = Counter()
    for r in rows:
        counts[r["region_class"]] += 1; bp[r["region_class"]] += r["length"]
    gcounts = Counter(g["best_class"] for g in grows) if grows is not None else None
    lines = []
    for k in CLASS_RANK:
        s = f"  {k:14s} {counts[k]:4d} regions ({100*counts[k]/len(rows):5.1f} %)  {bp[k]:9,} bp"
        if gcounts is not None:
            s += f"   {gcounts[k]:3d} genes"
        lines.append(s)
    return lines


def fish_absent_rates(rows, samples):
    return {s: sum(r[f"{s}_status"] == "absent" for r in rows) / len(rows) for s in samples}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam-dir", default=str(BASE / "analyses" / "23-native-methylation" / "align"))
    ap.add_argument("--threads", type=int, default=4, help="samtools threads per BAM")
    ap.add_argument("--null", type=int, default=3, help="number of shuffled-region null sets per direction")
    a = ap.parse_args()
    bam_dir = Path(a.bam_dir)
    OUT.mkdir(parents=True, exist_ok=True)

    summary = []
    gene_class = {}
    for eco, cfg in ECO.items():
        print(f"\n== {eco}-present regions, {eco} fish on own assembly ==", flush=True)
        regions = read_bed(cfg["regions"])
        rows, fm = classify(regions, cfg, bam_dir, a.threads, eco)
        grows = genes_for(rows, cfg, eco)
        for g in grows:
            gene_class[(eco, g["gene_id"])] = g["best_class"]
        summary.append(f"{eco}-present regions (n={len(rows)}, {sum(r['length'] for r in rows):,} bp); "
                       f"{eco} fish median flank depth: " + ", ".join(f"{s} {v}x" for s, v in fm.items()))
        summary += class_table(rows, grows)
        ar = fish_absent_rates(rows, cfg["samples"])
        summary.append("  per-fish absent rate: " + ", ".join(f"{s} {100*v:.0f} %" for s, v in ar.items()))
        pooled = Counter("full" if r["pooled_ratio"] != "inf" and float(r["pooled_ratio"]) >= 0.6 else
                         "half" if r["pooled_ratio"] != "inf" and float(r["pooled_ratio"]) >= 0.25 else "sparse"
                         for r in rows)
        summary.append(f"  pooled 4-fish depth vs flanks: full (>=0.6) {pooled['full']}, "
                       f"half (0.25-0.6) {pooled['half']}, sparse (<0.25) {pooled['sparse']}")
        n_absent_poly = Counter(r["n_absent"] for r in rows if r["region_class"] == "polymorphic")
        summary.append("  polymorphic regions by number of absent fish: " +
                       ", ".join(f"{k} absent: {v}" for k, v in sorted(n_absent_poly.items())))

        for k in range(1, a.null + 1):
            print(f"\n-- null {k}: shuffled {eco} regions --", flush=True)
            nr = shuffled(regions, cfg, seed=100 + k, label=f"{eco}.null{k}")
            nrows, _ = classify(nr, cfg, bam_dir, a.threads, f"{eco}.null{k}")
            summary.append(f"  null {k} (shuffled, n={len(nrows)}):")
            summary += ["  " + l for l in class_table(nrows)]
            nar = fish_absent_rates(nrows, cfg["samples"])
            summary.append("    per-fish absent rate: " + ", ".join(f"{s} {100*v:.0f} %" for s, v in nar.items()))
        summary.append("")

    tier = []
    with open(REFINED) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["tier"] != "A":
                continue
            eco = r["reciprocal_present_in"]
            c = gene_class.get((eco, r["gene_id"]))
            if c is None:
                cands = [v for (e, g), v in gene_class.items()
                         if e == eco and g.split("_")[0] == r["gene_id"]]
                c = min(cands, key=lambda k: CLASS_RANK[k]) if cands else "no_region_match"
            tier.append({"gene_id": r["gene_id"], "symbol": r["symbol"], "product": r["product"],
                         "present_in": eco, "lines_refined": r["lines_refined"],
                         "caution": r.get("caution", ""), "target_coverage_class": c})
    with open(OUT / "tierA_target_coverage.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(tier[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(tier)
    tc = Counter(t["target_coverage_class"] for t in tier)
    summary.append(f"Tier-A genes (n={len(tier)}) by target-ecotype coverage class of their reciprocal-PAV region:")
    for k in list(CLASS_RANK) + ["no_region_match"]:
        if tc[k]:
            summary.append(f"  {k:16s} {tc[k]:3d}")
    summary.append("\nDefinitions: per fish present = frac>=1x >= 0.80 and region/flank depth ratio >= 0.60; "
                   "absent = frac>=1x < 0.20; low = in between. fixed_present = all 4 present; "
                   "polymorphic = >=1 absent and >=1 present; low_all = all low; absent_all = all absent; "
                   "mixed = other. Flanks = 2 kb each side on the same assembly. Null = the same regions "
                   "shuffled to random positions on the same contigs (bedtools shuffle).")
    txt = "\n".join(summary)
    (OUT / "target_coverage_summary.txt").write_text(txt + "\n")
    print("\n" + txt)


if __name__ == "__main__":
    main()
