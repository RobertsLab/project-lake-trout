#!/usr/bin/env python3
"""
18.3-go-enrichment.py  (Step 4 of code/18-diff-annotation-phenotype-plan.md)

Self-contained GO over-representation analysis (ORA) for the DMR-gene set, the stringent
PAV-gene set, and their union, against the full GO-annotated gene background. Uses the LOCAL
GO annotations (Step 0 table, derived from NCBI's GAF) so that LOC###### genes without a
symbol are retained — an ortholog-based tool (g:Profiler/zebrafish) would drop most of them.

Method:
  * parse go-basic.obo -> id->name/namespace, is_a + part_of parents, alt_id, obsolete
  * propagate each gene's direct GO terms up the DAG (true-path rule)
  * hypergeometric SF (log-space, exact) per term; Benjamini-Hochberg FDR within each set
  * flag terms matching ecotype-phenotype keywords (lipid / buoyancy / growth / hypoxia ...)

Caveats baked into the output README: GO is mostly IEA (electronic); the PAV set is on a
LEAN reference so siscowet-specific enrichment is divergence-biased; ORA ignores gene length.
"""

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
ANN = BASE / "analyses" / "18-annotation"
OBO = ANN / "raw" / "go-basic.obo"
GENE_TBL = ANN / "gene_function_table.tsv"

MIN_TERM = 5          # ignore tiny background terms
MAX_TERM = 2000       # ignore huge unspecific terms
MIN_STUDY = 3         # need >=3 study genes hitting the term

PHENO = re.compile(
    r"lipid|fatty.?acid|triglyceride|triacylglycerol|sterol|cholesterol|lipo|"
    r"buoyan|swim.?bladder|gas bladder|adipos|fat\b|"
    r"growth|myogen|muscle|skeletal|bone|ossif|"
    r"hypoxia|oxygen|response to oxygen|"
    r"osmo|ion transport|temperature|cold|"
    r"immune|immunoglobulin|complement|inflammat", re.I)


# --------------------------------------------------------------- OBO parser
def parse_obo(path):
    name, ns, parents, alt = {}, {}, defaultdict(set), {}
    cur, obsolete = None, False
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur, obsolete = None, False
                continue
            if line.startswith("[") and line != "[Term]":
                cur = None
                continue
            if line.startswith("id: GO:"):
                cur = line[4:]
            elif cur and line.startswith("name: "):
                name[cur] = line[6:]
            elif cur and line.startswith("namespace: "):
                ns[cur] = line[11:]
            elif cur and line.startswith("alt_id: "):
                alt[line[8:]] = cur
            elif cur and line.startswith("is_a: "):
                parents[cur].add(line[6:].split(" ! ")[0].strip())
            elif cur and line.startswith("relationship: part_of "):
                parents[cur].add(line.split()[2])
            elif cur and line.startswith("is_obsolete: true"):
                obsolete = True
                name.pop(cur, None)
                ns.pop(cur, None)
                parents.pop(cur, None)
    return name, ns, parents, alt


def build_ancestors(parents):
    cache = {}

    def anc(t):
        if t in cache:
            return cache[t]
        cache[t] = set()           # guard against cycles
        acc = set()
        for p in parents.get(t, ()):
            acc.add(p)
            acc |= anc(p)
        cache[t] = acc
        return acc

    return {t: anc(t) for t in parents}, anc


# --------------------------------------------------------------- gene GO
def load_gene_go(name, alt, ancestor_fn):
    """gene_id -> propagated GO set (only terms that exist & are non-obsolete)."""
    gene_go = {}
    with open(GENE_TBL) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            raw = r["go_terms"].split("|") if r["go_terms"] else []
            if not raw:
                continue
            full = set()
            for g in raw:
                g = alt.get(g, g)
                if g not in name:
                    continue
                full.add(g)
                full |= ancestor_fn(g)
            full = {t for t in full if t in name}
            if full:
                gene_go[r["gene_id"]] = full
    return gene_go


def study_genes(path, col="gene_id"):
    s = set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s.add(r[col])
    return s


# --------------------------------------------------------- hypergeometric
def logcomb(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hyper_sf(k, N, K, n):
    """P(X >= k) for drawing n from N with K successes."""
    lo, hi = k, min(K, n)
    if lo <= 0:
        return 1.0
    if lo > hi:
        return 0.0
    base = logcomb(N, n)
    logp = logcomb(K, lo) + logcomb(N - K, n - lo) - base
    p = math.exp(logp)
    total = p
    for i in range(lo, hi):
        ratio = ((K - i) / (i + 1)) * ((n - i) / (N - K - n + i + 1))
        p *= ratio
        total += p
    return min(total, 1.0)


def bh(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    fdr = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        idx = m - rank + 1
        q = pvals[i] * m / idx
        prev = min(prev, q)
        fdr[i] = prev
    return fdr


# --------------------------------------------------------------- ORA
def ora(study, term2genes, universe, name, ns):
    N = len(universe)
    study = study & universe
    n = len(study)
    rows = []
    pvals = []
    for go, genes in term2genes.items():
        K = len(genes)
        if K < MIN_TERM or K > MAX_TERM:
            continue
        hit = genes & study
        k = len(hit)
        if k < MIN_STUDY:
            continue
        p = hyper_sf(k, N, K, n)
        fold = (k / n) / (K / N) if n and K else 0.0
        rows.append({
            "go_id": go, "name": name.get(go, ""), "namespace": ns.get(go, ""),
            "study_k": k, "study_n": n, "bg_K": K, "bg_N": N,
            "fold_enrichment": round(fold, 2), "p_value": p,
            "phenotype": "Y" if PHENO.search(name.get(go, "")) else "",
            "genes": ",".join(sorted(h.replace("gene-", "") for h in
                                     (term2genes[go] & study))[:25]),
        })
        pvals.append(p)
    fdr = bh(pvals)
    for r, q in zip(rows, fdr):
        r["fdr"] = q
    rows.sort(key=lambda r: (r["fdr"], r["p_value"], -r["fold_enrichment"]))
    return rows, n, N


def write(rows, path):
    cols = ["go_id", "name", "namespace", "study_k", "study_n", "bg_K", "bg_N",
            "fold_enrichment", "p_value", "fdr", "phenotype", "genes"]
    with open(path, "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows:
            r["p_value"] = f"{r['p_value']:.3e}"
            r["fdr"] = f"{r['fdr']:.3e}"
            o.write("\t".join(str(r[c]) for c in cols) + "\n")


def main():
    print("parsing OBO ...", file=sys.stderr)
    name, ns, parents, alt = parse_obo(OBO)
    print(f"  {len(name)} live terms, {len(alt)} alt_ids", file=sys.stderr)
    _, anc_fn = build_ancestors(parents)

    print("propagating gene GO ...", file=sys.stderr)
    gene_go = load_gene_go(name, alt, anc_fn)
    universe = set(gene_go)
    print(f"  {len(universe)} genes with propagated GO (background)", file=sys.stderr)

    term2genes = defaultdict(set)
    for g, terms in gene_go.items():
        for t in terms:
            term2genes[t].add(g)

    sets = {
        "dmr": study_genes(ANN / "dmr_gene_assignments.tsv"),
        "pav": study_genes(ANN / "pav_gene_assignments.tsv"),
    }
    sets["union"] = sets["dmr"] | sets["pav"]

    for key, study in sets.items():
        rows, n, N = ora(study, term2genes, universe, name, ns)
        sig = sum(1 for r in rows if r["fdr"] < 0.1)
        pflag = sum(1 for r in rows if r["phenotype"] == "Y" and r["fdr"] < 0.25)
        out = ANN / f"go_enrichment_{key}.tsv"
        write(rows, out)
        print(f"[{key}] study={n} (of {len(study)}), bg={N}, terms tested={len(rows)}, "
              f"FDR<0.1={sig}, phenotype-flagged(FDR<0.25)={pflag} -> {out.name}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
