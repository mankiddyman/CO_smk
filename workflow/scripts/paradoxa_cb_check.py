#!/usr/bin/env python3
"""paradoxa_cb_check.py -- did switching the crossover reference to C + B do what it should?

The two crossover references of D. paradoxa:
  std  Dparadoxa_std: A (chr1_hap1, L1.P) + D (chr2_hap2, L2.Q) + chr3-6 hap1
  CB   Dparadoxa_CB:  C (chr1_hap2, L1.Q) + B (chr2_hap1, L2.P) + chr3-6 hap1
On std the pollen labels flip at the A and D joins in about half the grains (the arms on
either side are inherited independently), and the caller counts each flip as a crossover.
On CB the arms that travel together sit on one scaffold, so those flips should vanish.

  1. crossovers per grain     mean and median; obligate minimum one per bivalent = 0.5 per
                              chromosome per grain = 3.0 for six chromosomes
  2. per chromosome           mean per grain, std vs CB (the chr1/chr2 scaffolds differ)
  3. at the joins             crossovers within 10 Mb of each join, against the chromosome's
                              average for any 20 Mb
  4. hotspots                 5 Mb windows with more crossovers than chance allows: Poisson
                              tail at the genome-wide mean, p < 0.05 / number of windows
  5. transmission ratio (CB)  per 10 Mb window, the share of grains carrying the ALT copy
                              (window calls as linkage_scan.py: >= 3 molecules, ALT if >= 80 %,
                              REF if <= 20 %). A two-sided binomial p < 0.001 flags a window
                              passed on unequally: the signature of biased segregation or
                              pollen selection.

Run from the CO_smk root. Usage: paradoxa_cb_check.py OUT.txt
"""
import csv
import math
import os
import sys

import numpy as np
import pandas as pd

OUT = sys.argv[1]
REFS = [("std", "Dparadoxa_std"), ("CB", "Dparadoxa_CB")]
JOINS = {"std": [("A join L1|P", "chr1_hap1", 262.9), ("D join L2|Q", "chr2_hap2", 215.0)],
         "CB": [("C join L1|Q (approx.)", "chr1_hap2", 262.0), ("B join L2|P (approx.)", "chr2_hap1", 211.0)]}
PAIRED = [("chr1", "chr1_hap1", "chr1_hap2"), ("chr2", "chr2_hap2", "chr2_hap1")] + \
         [("chr%d" % i, "chr%d_hap1" % i, "chr%d_hap1" % i) for i in range(3, 7)]
CHROMS = {"std": [a for _, a, _ in PAIRED], "CB": [b for _, _, b in PAIRED]}
HOT_W, HOT_ALPHA, NEAR, RATIO_W, MIN_MOL, MIN_CELLS, P_FLAG = 5e6, 0.05, 10e6, 10e6, 3, 20, 1e-3


def lengths(sample):
    """Chromosome lengths of the reference the pollen were mapped to."""
    p = "results/reference/%s/genome.fa.fai" % sample
    if not os.path.exists(p):
        row = next(r for r in csv.DictReader(open("config/samples.csv"))
                   if sample in (r.get("sample_id"), list(r.values())[0]))
        p = row["assembly_fasta"] + ".fai"
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(p) if l.strip()}


def binom_p(k, n):
    """Two-sided binomial test against 50:50."""
    pk = [math.comb(n, i) * 0.5 ** n for i in range(n + 1)]
    return min(1.0, sum(p for p in pk if p <= pk[k] * (1 + 1e-9)))


def poisson_tail(k, mu):
    """P(X >= k) for X ~ Poisson(mu)."""
    return max(0.0, 1.0 - sum(math.exp(-mu) * mu ** i / math.factorial(i) for i in range(k)))


D = {}
for lab, S in REFS:
    bed = "results/crossovers/%s/co_intervals.bed" % S
    pc = "results/crossovers/%s/co_per_cell.tsv" % S
    if not (os.path.exists(bed) and os.path.exists(pc)):
        sys.exit("no crossover calls for %s yet (%s)" % (S, bed))
    co = pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "barcode"])
    co["mid"] = (co.start + co.end) / 2
    per = pd.read_csv(pc, sep="\t")
    L = lengths(S)
    D[lab] = dict(S=S, co=co, per=per, L={c: L[c] for c in CHROMS[lab] if c in L})

o = ["PARADOXA CROSSOVER REFERENCES: std (A + D) against CB (C + B)", ""]

# 1. per grain
o += ["1. CROSSOVERS PER POLLEN GRAIN (obligate minimum 3.0: one per bivalent)",
      "   %-10s %7s %11s %8s %8s" % ("reference", "grains", "crossovers", "mean", "median")]
for lab in ("std", "CB"):
    n = D[lab]["per"]["n_cos"]
    o.append("   %-10s %7d %11d %8.2f %8.0f" % (lab, len(n), len(D[lab]["co"]), n.mean(), n.median()))

# 2. per chromosome
o += ["", "2. MEAN CROSSOVERS PER GRAIN, PER CHROMOSOME (below 0.5 = fewer than one per bivalent)",
      "   %-6s %-12s %9s   %-12s %9s" % ("", "std", "per grain", "CB", "per grain")]
for name, a, b in PAIRED:
    v = [(D[lab]["co"].chrom == c).sum() / len(D[lab]["per"]) for lab, c in (("std", a), ("CB", b))]
    o.append("   %-6s %-12s %9s   %-12s %9s" % (name, a, "%.2f%s" % (v[0], " *" if v[0] < 0.5 else ""),
                                               b, "%.2f%s" % (v[1], " *" if v[1] < 0.5 else "")))
o.append("   * below 0.5")

# 3. at the joins
o += ["", "3. CROSSOVERS AT THE JOINS",
      "   %-4s %-24s %-10s %9s %15s %22s" % ("ref", "join", "scaffold", "position", "within 10 Mb",
                                           "any 20 Mb, on average")]
for lab in ("std", "CB"):
    co, L = D[lab]["co"], D[lab]["L"]
    for jn, c, mb in JOINS[lab]:
        near = ((co.chrom == c) & ((co.mid - mb * 1e6).abs() <= NEAR)).sum()
        avg = (co.chrom == c).sum() * 2 * NEAR / L[c] if c in L else float("nan")
        o.append("   %-4s %-24s %-10s %6.1f Mb %15d %22.1f" % (lab, jn, c, mb, near, avg))

# 4. hotspots
o += ["", "4. HOTSPOTS: 5 Mb windows with more crossovers than chance allows (Poisson at the genome-wide mean)"]
for lab in ("std", "CB"):
    co, L = D[lab]["co"], D[lab]["L"]
    rows = []
    for c, n in L.items():
        k = int(math.ceil(n / HOT_W))
        cnt = np.bincount(np.minimum((co.mid[co.chrom == c] // HOT_W).astype(int), k - 1), minlength=k)
        rows += [(c, i, int(x)) for i, x in enumerate(cnt)]
    mu = np.mean([r[2] for r in rows])
    cut = HOT_ALPHA / len(rows)
    hot = sorted([r for r in rows if poisson_tail(r[2], mu) < cut], key=lambda r: -r[2])
    top = max(rows, key=lambda r: r[2])
    o.append("   %s: %d windows, %.1f crossovers per window on average; a hotspot needs %d or more; busiest window %s %d-%d Mb (%d)"
             % (lab, len(rows), mu, next(k for k in range(1, 1000) if poisson_tail(k, mu) < cut),
                top[0], top[1] * 5, top[1] * 5 + 5, top[2]))
    for c, i, cnt in hot:
        o.append("      %-10s %4d-%4d Mb  %3d crossovers" % (c, i * 5, i * 5 + 5, cnt))
    if not hot:
        o.append("      no hotspot")

# 5. transmission ratio on CB
S, L = D["CB"]["S"], D["CB"]["L"]
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
nwin = {c: max(1, int(round(n / RATIO_W))) for c, n in L.items()}
wins = [(c, i) for c in L for i in range(nwin[c])]
widx = {w: k for k, w in enumerate(wins)}
G = np.full((len(cells), len(wins)), np.nan)
for ci, bc in enumerate(cells):
    p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
    if not os.path.exists(p):
        continue
    m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"])
    m = m[(m.rc != m.ac) & m.chrom.isin(L)]
    if not len(m):
        continue
    m["alt"] = (m.ac > m.rc).astype(int)
    m["w"] = [min(int(q // RATIO_W), nwin[c] - 1) for c, q in zip(m.chrom, m.pos)]
    g = m.groupby(["chrom", "w"]).alt.agg(["mean", "size"])
    for (c, w), r in g.iterrows():
        if r["size"] >= MIN_MOL and (r["mean"] >= 0.8 or r["mean"] <= 0.2):
            G[ci, widx[(c, w)]] = 1 if r["mean"] >= 0.8 else 0
res = []
for k, (c, i) in enumerate(wins):
    v = G[:, k][~np.isnan(G[:, k])]
    if len(v) >= MIN_CELLS:
        alt = int(v.sum())
        res.append((binom_p(alt, len(v)), c, i, len(v) - alt, alt))
share = np.array([r[4] / (r[3] + r[4]) for r in res])
flag = sorted(r for r in res if r[0] < P_FLAG)
o += ["", "5. TRANSMISSION RATIO ON CB: share of grains carrying the ALT copy, 10 Mb windows with >= %d grains called"
      % MIN_CELLS,
      "   windows tested: %d of %d; ALT share median %.2f (middle 90%%: %.2f-%.2f)"
      % (len(res), len(wins), np.median(share), *np.percentile(share, [5, 95])),
      "   windows passed on unequally (two-sided binomial p < 0.001; about %.1f expected by chance): %d"
      % (len(res) * P_FLAG, len(flag))]
show = flag if flag else sorted(res)[:5]
if not flag:
    o.append("   none; the five most uneven windows:")
for p, c, i, ref, alt in sorted(show, key=lambda r: (r[1], r[2])):
    o.append("      %-10s %4d-%4d Mb  REF %3d  ALT %3d  ALT share %.2f  p %.1e"
             % (c, i * 10, i * 10 + 10, ref, alt, alt / (ref + alt), p))

open(OUT, "w").write("\n".join(o) + "\n")
print("\n".join(o))
