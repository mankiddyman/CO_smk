#!/usr/bin/env python3
"""pollen_pairs_check.py -- is the pollen linkage between particular pieces solid?

Builds the per-cell window genotypes exactly as linkage_scan.py does (one row per
informative molecule; window = ALT if >= 80% of its molecules are ALT, REF if <= 20%,
at least 3 molecules), then:
  1. cells: are the good cells independent meiotic products? Pairs of cells that agree
     in >= 98% of >= 100 shared windows are the same genotype (nuclei of one pollen grain, or a
     duplicated barcode); those clusters are collapsed to one cell for step 2.
  2. pairs: for each chosen pair of windows, the 2 x 2 table of cell genotypes and r
     (share of cells where the two windows disagree), on all cells and on unique genotypes.
  3. calls: how clean the window calls are in those windows (share of cells whose ALT
     fraction is pure, 0 or 1, versus near the 20 % / 80 % cut-offs).

Run from the CO_smk root. Usage: pollen_pairs_check.py SAMPLE OUT.txt [--window_mb 10]
"""
import argparse
import csv
import os

import numpy as np
import pandas as pd

PAIRS = [  # label, (chrom, window start Mb), (chrom, window start Mb)
    ("L1 end    vs P start  (both on A = chr1_hap1)", ("chr1_hap1", 240), ("chr1_hap1", 270)),
    ("L1 end    vs Q start  (A vs D)", ("chr1_hap1", 240), ("chr2_hap2", 230)),
    ("L2 end    vs Q start  (both on D = chr2_hap2)", ("chr2_hap2", 200), ("chr2_hap2", 230)),
    ("L2 end    vs P start  (D vs A)", ("chr2_hap2", 200), ("chr1_hap1", 270)),
    ("L1 middle vs L1 end   (control: same arm)", ("chr1_hap1", 150), ("chr1_hap1", 240)),
    ("chr3 middle vs chr4 piece 0-85 (one-way move)", ("chr3_hap1", 140), ("chr4_hap1", 50)),
    ("chr3 middle vs chr4 rest", ("chr3_hap1", 140), ("chr4_hap1", 250)),
    ("chr5 middle vs chr6 middle", ("chr5_hap1", 130), ("chr6_hap1", 180)),
    ("chr5 near 61 Mb (copy of the junction repeat) vs chr6", ("chr5_hap1", 60), ("chr6_hap1", 180)),
    ("chr1 L1 middle vs chr6 middle (control: unrelated)", ("chr1_hap1", 150), ("chr6_hap1", 180)),
]

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("out")
ap.add_argument("--window_mb", type=float, default=10)
ap.add_argument("--min_mol", type=int, default=3)
A = ap.parse_args()
S, W = A.sample, int(A.window_mb * 1e6)
row = next(r for r in csv.DictReader(open("config/samples.csv")) if S in (r.get("sample_id"), list(r.values())[0]))
lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(row["assembly_fasta"] + ".fai")}
nwin = {c: max(1, int(round(lens[c] / W))) for c in lens}
wins = [(c, i) for c in lens for i in range(nwin[c])]
widx = {w: k for k, w in enumerate(wins)}
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
G = np.full((len(cells), len(wins)), np.nan)
F = np.full((len(cells), len(wins)), np.nan)          # ALT fraction per cell and window
for ci, bc in enumerate(cells):
    p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
    if not os.path.exists(p):
        continue
    m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"])
    m = m[(m.rc != m.ac) & m.chrom.isin(lens)]
    m["alt"] = (m.ac > m.rc).astype(int)
    m["w"] = [min(q // W, nwin[c] - 1) for c, q in zip(m.chrom, m.pos)]
    g = m.groupby(["chrom", "w"]).alt.agg(["mean", "size"])
    for (c, w), r in g.iterrows():
        if r["size"] >= A.min_mol:
            k = widx[(c, w)]
            F[ci, k] = r["mean"]
            if r["mean"] >= 0.8:
                G[ci, k] = 1
            elif r["mean"] <= 0.2:
                G[ci, k] = 0
out = ["POLLEN PAIRS CHECK  %s  (%d good cells, %g Mb windows)" % (S, len(cells), A.window_mb), ""]

# 1. are the cells independent?
n = len(cells)
same = np.full((n, n), np.nan)
for i in range(n):
    for j in range(i + 1, n):
        ok = ~np.isnan(G[i]) & ~np.isnan(G[j])
        if ok.sum() >= 100:
            same[i, j] = same[j, i] = float((G[i, ok] == G[j, ok]).mean())
iu = np.triu_indices(n, 1)
v = same[iu][~np.isnan(same[iu])]
parent = list(range(n))


def root(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]
        a = parent[a]
    return a


for i, j in zip(*iu):
    if same[i, j] >= 0.98:
        parent[root(i)] = root(j)
groups = {}
for i in range(n):
    groups.setdefault(root(i), []).append(i)
rep = sorted(min(g) for g in groups.values())
out += ["1. ARE THE CELLS INDEPENDENT MEIOTIC PRODUCTS?",
        "   cell pairs compared (>= 100 shared windows): %d; share of windows with the same genotype: median %.2f "
        "(independent products: ~0.5)" % (len(v), np.median(v) if len(v) else np.nan),
        "   pairs >= 0.98 identical (same pollen grain or duplicate barcode): %d; pairs <= 0.02 (exact opposites): %d"
        % ((v >= 0.98).sum(), (v <= 0.02).sum()),
        "   distinct genotypes after merging identical cells: %d of %d cells; groups larger than one: %s"
        % (len(rep), n, ", ".join(str(len(g)) for g in groups.values() if len(g) > 1) or "none"), ""]


# 2. the 2 x 2 tables
def win(c, mb):
    return widx.get((c, min(int(mb * 1e6 // W), nwin.get(c, 1) - 1)))


out += ["2. WINDOW PAIRS: cell genotypes (REF/ALT relative to the reference) and r = share of cells that disagree",
        "   %-50s %8s %8s %8s %8s %7s %7s %14s" % ("pair", "REF-REF", "REF-ALT", "ALT-REF", "ALT-ALT", "cells", "r",
                                                  "r unique cells")]
for lab, (c1, a1), (c2, a2) in PAIRS:
    k1, k2 = win(c1, a1), win(c2, a2)
    if k1 is None or k2 is None:
        out.append("   %-50s window not in this reference" % lab)
        continue
    res = []
    for rows in (range(n), rep):
        x, y = G[list(rows), k1], G[list(rows), k2]
        ok = ~np.isnan(x) & ~np.isnan(y)
        t = [int(((x == a) & (y == b) & ok).sum()) for a, b in ((0, 0), (0, 1), (1, 0), (1, 1))]
        res.append((t, ok.sum(), (t[1] + t[2]) / ok.sum() if ok.sum() else np.nan))
    (t, m_, r), (_, mu, ru) = res
    out.append("   %-50s %8d %8d %8d %8d %7d %7.2f %8.2f (n=%d)" % (lab, t[0], t[1], t[2], t[3], m_, r, ru, mu))

# 3. how clean are the calls
out += ["", "3. HOW CLEAN ARE THE CALLS IN THOSE WINDOWS (ALT fraction per cell, windows with >= %d molecules)" % A.min_mol,
        "   %-28s %7s %10s %13s %12s" % ("window", "cells", "pure 0/1", "near cut-off", "no call")]
seen = set()
for _, (c1, a1), (c2, a2) in PAIRS:
    for c, a in ((c1, a1), (c2, a2)):
        k = win(c, a)
        if k is None or k in seen:
            continue
        seen.add(k)
        f = F[:, k][~np.isnan(F[:, k])]
        if not len(f):
            continue
        pure = ((f <= 0.02) | (f >= 0.98)).mean()
        edge = (((f > 0.02) & (f <= 0.2)) | ((f >= 0.8) & (f < 0.98))).mean()
        nocall = ((f > 0.2) & (f < 0.8)).mean()
        out.append("   %-28s %7d %9.0f%% %12.0f%% %11.0f%%" % ("%s %g-%g Mb" % (c, a, a + A.window_mb), len(f),
                                                            100 * pure, 100 * edge, 100 * nocall))
open(A.out, "w").write("\n".join(out) + "\n")
print("\n".join(out))
