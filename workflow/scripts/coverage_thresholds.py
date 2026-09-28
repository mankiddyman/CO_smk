#!/usr/bin/env python3
"""coverage_thresholds.py -- how many cells survive a total-molecule AND a per-chromosome floor?

A cell's evidence for crossovers on a chromosome is the molecules it has THERE.
A deep cell can still be nearly empty on one chromosome, and every chromosome
counts in the landscape's denominator, so (like Spondias's min_per_chrom_markers)
a cell is kept only if its weakest chromosome clears the floor.

Reads the per-chromosome molecule counts written by chrom_haploidness.py
(informative markers, one molecule per read). Prints kept-cell counts over a
grid of (total, weakest-chromosome) floors, where the review panel's cells
fall, and plots total vs weakest chromosome.

Usage: coverage_thresholds.py SAMPLE [PANEL_SEED]
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S = sys.argv[1]
SEED = sys.argv[2] if len(sys.argv) > 2 else "1"
t = pd.read_csv("qc/haplotypes/%s/%s_chrom_haploidness.tsv.gz" % (S, S), sep="\t")
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
chroms = sorted(t.chrom.unique(), key=lambda c: int("".join(ch for ch in c.split("_")[0] if ch.isdigit()) or 0))
m = t.pivot_table(index="barcode", columns="chrom", values="molecules", aggfunc="sum").reindex(
    index=cells, columns=chroms).fillna(0)
total, weakest = m.sum(axis=1), m.min(axis=1)
hap = pd.read_csv("qc/haplotypes/%s/haplotype_tracks.tsv.gz" % S, sep="\t",
                  usecols=["barcode", "haploidness"]).set_index("barcode").haploidness.reindex(cells)

TOT = [0, 800, 1000, 1500, 2000, 3000]
WK = [0, 20, 30, 40, 60, 80]
lines = ["CELLS KEPT  %s  (%d currently selected, %d chromosomes)" % (S, len(cells), len(chroms)),
         "  rows: total molecules >=   columns: weakest chromosome >=",
         "  %10s" % "" + "".join("%8d" % w for w in WK)]
for a in TOT:
    lines.append("  %10d" % a + "".join("%8d" % int(((total >= a) & (weakest >= w)).sum()) for w in WK))
lines.append("  median molecules per chromosome, cells with >= 1500 total: "
             + ", ".join("%s %d" % (c.split("_")[0], m.loc[total >= 1500, c].median())
                         for c in chroms if (total >= 1500).any()))
pf = "qc/review/%s/panel_seed%s.txt" % (S, SEED)
if os.path.exists(pf):
    panel = [l.strip() for l in open(pf) if l.strip()]
    lines.append("  review panel (seed %s): number, barcode, total, weakest chromosome (which), haploidness" % SEED)
    for i, bc in enumerate(panel, 1):
        if bc in m.index:
            lines.append("    %02d  %s  %6d  %5d (%s)  %.2f" % (i, bc, total[bc], weakest[bc],
                                                            m.loc[bc].idxmin().split("_")[0], hap[bc]))
out = "qc/review/%s" % S
os.makedirs(out, exist_ok=True)
open(os.path.join(out, "coverage_thresholds.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

fig, ax = plt.subplots(figsize=(7.5, 6))
sc = ax.scatter(total, weakest.clip(lower=0.8), c=hap, cmap="viridis", s=8, alpha=.7, vmin=0.6, vmax=1)
fig.colorbar(sc, ax=ax, label="haploidness")
for a in (1000, 1500, 2000):
    ax.axvline(a, color="grey", ls=":", lw=1)
for w in (20, 40, 60):
    ax.axhline(w, color="grey", ls=":", lw=1)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("molecules in the cell (informative markers)")
ax.set_ylabel("molecules on its WEAKEST chromosome")
ax.set_title("%s: %d selected cells" % (S, len(cells)))
fig.tight_layout()
fig.savefig(os.path.join(out, "coverage_thresholds.png"), dpi=130)
print("wrote %s/coverage_thresholds.{txt,png}" % out)
