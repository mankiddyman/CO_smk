#!/usr/bin/env python3
"""linkage_scan.py -- do the pollen confirm the reference's joins?

Two stretches of sequence on one physical chromosome are inherited together
unless a crossover falls between them: across gametes their alleles agree far
more often than not. Stretches on different chromosomes assort independently:
they agree in half the gametes (recombination fraction r = 0.5). The gametes
therefore test every join in the reference, with no prior knowledge of where
rearrangements or assembly errors might be.

Per cell and window (default 10 Mb), from the molecule tables (one row per
informative molecule): genotype = ALT if >= 80% of its molecules are ALT, REF
if <= 20%, otherwise unknown (too few molecules, or both alleles present).
For every pair of windows, r = share of cells informative for both whose
genotypes disagree.

  joins     r between neighbouring windows along each reference chromosome.
            Real crossovers make r small (a few %); r >= MAX_R means the two
            sides are not joined in the plant, whatever the assembly says.
  partners  for each window beside an unsupported join, the windows elsewhere
            it is most tightly linked to: where that piece actually belongs.
  map       genome-wide r matrix (dark = tightly linked).

Usage: linkage_scan.py SAMPLE OUTDIR [--window_mb 10] [--min_cells 20] [--max_r 0.3]
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("sample"); ap.add_argument("outdir")
ap.add_argument("--window_mb", type=float, default=10)
ap.add_argument("--min_mol", type=int, default=3)
ap.add_argument("--min_cells", type=int, default=20)
ap.add_argument("--max_r", type=float, default=0.3)
A = ap.parse_args()
S, W = A.sample, int(A.window_mb * 1e6)
os.makedirs(A.outdir, exist_ok=True)

fai = "results/reference/%s/genome.fa.fai" % S
chroms = [l.split("\t")[0] for l in open(fai)]
lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
wins = [(c, i) for c in chroms for i in range(max(1, int(round(lens[c] / W))))]
widx = {w: k for k, w in enumerate(wins)}
nwin = {c: max(1, int(round(lens[c] / W))) for c in chroms}

cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
G = np.full((len(cells), len(wins)), np.nan)
for ci, bc in enumerate(cells):
    p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
    if not os.path.exists(p):
        continue
    m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"])
    m = m[(m.rc != m.ac) & m.chrom.isin(lens)]
    m["alt"] = (m.ac > m.rc).astype(int)
    m["w"] = [min(p_ // W, nwin[c] - 1) for c, p_ in zip(m.chrom, m.pos)]
    g = m.groupby(["chrom", "w"]).alt.agg(["mean", "size"])
    for (c, w), row in g.iterrows():
        if row["size"] >= A.min_mol:
            if row["mean"] >= 0.8:
                G[ci, widx[(c, w)]] = 1
            elif row["mean"] <= 0.2:
                G[ci, widx[(c, w)]] = 0
M = (~np.isnan(G)).astype(float)
X = np.nan_to_num(G, nan=0.0) * M
Y = (1 - np.nan_to_num(G, nan=0.0)) * M
n = M.T @ M
agree = X.T @ X + Y.T @ Y
with np.errstate(invalid="ignore", divide="ignore"):
    R = np.where(n >= A.min_cells, 1 - agree / n, np.nan)
np.save(os.path.join(A.outdir, "r_matrix.npy"), R)
pd.DataFrame(wins, columns=["chrom", "window"]).to_csv(os.path.join(A.outdir, "windows.tsv"), sep="\t", index=False)


def label(k):
    c, i = wins[k]
    return "%s:%d-%d" % (c, i * A.window_mb, min((i + 1) * A.window_mb, lens[c] / 1e6))


L = ["LINKAGE SCAN  %s  (%d cells, %g Mb windows, r = share of cells whose genotypes disagree)" % (S, len(cells), A.window_mb),
     "  every join between neighbouring windows; r >= %.2f = not joined in the plant" % A.max_r]
bad = []
for c in chroms:
    rs = []
    for i in range(nwin[c] - 1):
        k = widx[(c, i)]
        r = R[k, k + 1]
        rs.append(r)
        if not np.isnan(r) and r >= A.max_r:
            bad.append((c, i, r, int(n[k, k + 1])))
    v = np.array([x for x in rs if not np.isnan(x)])
    L.append("    %-11s %3d joins   median r %.3f   max r %.3f   untestable %d" % (
        c, len(rs), np.median(v) if len(v) else np.nan, v.max() if len(v) else np.nan, sum(np.isnan(rs))))
L.append("  unsupported joins (and where each side is most tightly linked instead):")
if not bad:
    L.append("    none -- every testable join is supported by the pollen")
for c, i, r, cnt in bad:
    L.append("    %s at %g Mb: r = %.2f over %d cells" % (c, (i + 1) * A.window_mb, r, cnt))
    for side, k in (("left ", widx[(c, i)]), ("right", widx[(c, i + 1)])):
        row = R[k].copy()
        for j, (cj, ij) in enumerate(wins):
            if cj == c and abs(ij - wins[k][1]) <= 1:
                row[j] = np.nan               # skip itself and its immediate neighbours
        order = [j for j in np.argsort(np.nan_to_num(row, nan=9)) if not np.isnan(row[j])][:3]
        L.append("      %s %-22s best linked to: %s" % (side, label(k), ", ".join(
            "%s (r %.2f)" % (label(j), row[j]) for j in order)))
open(os.path.join(A.outdir, "linkage_summary.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))

fig, ax = plt.subplots(figsize=(9, 8))
im = ax.imshow(R, cmap="magma", vmin=0, vmax=0.5, interpolation="nearest")
edges = np.cumsum([nwin[c] for c in chroms])
for e in edges[:-1]:
    ax.axhline(e - 0.5, color="white", lw=0.6); ax.axvline(e - 0.5, color="white", lw=0.6)
mids = edges - np.array([nwin[c] for c in chroms]) / 2.0
ax.set_xticks(mids); ax.set_xticklabels(chroms, rotation=90, fontsize=8)
ax.set_yticks(mids); ax.set_yticklabels(chroms, fontsize=8)
for c, i, r, cnt in bad:
    k = widx[(c, i)] + 0.5
    ax.plot([k, k], [k - 3, k + 3], color="cyan", lw=1.2); ax.plot([k - 3, k + 3], [k, k], color="cyan", lw=1.2)
plt.colorbar(im, ax=ax, fraction=0.046, label="r (share of pollen whose alleles disagree; 0.5 = unlinked)")
ax.set_title("%s: genetic linkage between %g Mb windows (cyan: unsupported joins)" % (S, A.window_mb), fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(A.outdir, "linkage_map.png"), dpi=130)
print("\nwrote %s/{linkage_summary.txt,linkage_map.png,r_matrix.npy,windows.tsv}" % A.outdir)
