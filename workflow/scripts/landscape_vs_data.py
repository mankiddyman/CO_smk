#!/usr/bin/env python3
"""landscape_vs_data.py -- is the landscape's shape biology, or where the caller could see?

Per chromosome, in windows:
  CO rate, midpoints        each crossover counted at the middle of its interval
  CO rate, spread           each crossover spread evenly over its interval --
                            the honest version where localisation is poor
  marker calls per Mb       where the data is (cells_called summed per window,
                            from marker_classes), scaled to the CO axis
  genes per Mb              from the annotation, scaled to the CO axis
  CO interval width         how precisely crossovers are placed, by window

If the CO rate simply mirrors marker density, the shape is at least partly
detection. If the middles are cold with NARROW intervals, they are really
cold; if the middle crossovers have WIDE intervals, the middles are poorly
observed and only the spread version should be read there.

Usage: landscape_vs_data.py OUTDIR WINDOW_MB SAMPLE FAI GFF3
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT, WMB, S, FAI, GFF = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
W = int(WMB * 1e6)
os.makedirs(OUT, exist_ok=True)
L = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"]).set_index("chrom")["len"]
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]

rows = []
for bc in cells:
    p = "results/crossovers/%s/per_cell/%s_co_pred.txt" % (S, bc)
    if not os.path.exists(p):
        continue
    for l in open(p):
        f = l.split()
        if len(f) >= 3 and f[0] in L.index:
            try:
                rows.append((f[0], int(f[1]), int(f[2])))
            except ValueError:
                pass
co = pd.DataFrame(rows, columns=["chrom", "start", "end"])
co["mid"] = (co.start + co.end) / 2.0
co["width"] = (co.end - co.start).clip(lower=1)
chroms = [c for c in L.index if c in set(co.chrom)]
n = float(len(cells))

mk = pd.read_csv("qc/markers/%s/marker_classes.tsv.gz" % S, sep="\t", usecols=["chrom", "pos", "cells_called"])
genes = []
with open(GFF) as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        f = l.split("\t", 5)
        if len(f) > 4 and f[2] == "gene":
            genes.append((f[0], int(f[3])))
genes = pd.DataFrame(genes, columns=["chrom", "pos"])

summary, allw = [], []
ncol = 4 if len(chroms) > 6 else 3
nrow = int(np.ceil(len(chroms) / float(ncol)))
fig, ax = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 2.9 * nrow), squeeze=False)
for k, ch in enumerate(chroms):
    nb = int(np.ceil(L[ch] / float(W)))
    x = (np.arange(nb) + 0.5) * W / 1e6
    c = co[co.chrom == ch]
    mid = np.bincount((c.mid // W).astype(int).clip(0, nb - 1), minlength=nb) / n * 100 / WMB
    spread = np.zeros(nb)
    for s0, e0 in zip(c.start.values, c.end.values):
        a, b = int(s0 // W), int(min(e0, L[ch] - 1) // W)
        for w in range(a, b + 1):
            lo, hi = max(s0, w * W), min(e0, (w + 1) * W)
            spread[w] += max(hi - lo, 0) / float(max(e0 - s0, 1))
    spread = spread / n * 100 / WMB
    m = mk[mk.chrom == ch]
    md = np.bincount((m.pos // W).astype(int).clip(0, nb - 1), weights=m.cells_called, minlength=nb) / WMB
    g = genes[genes.chrom == ch]
    gd = np.bincount((g.pos // W).astype(int).clip(0, nb - 1), minlength=nb) / WMB
    wid = c.groupby((c.mid // W).astype(int).clip(0, nb - 1))["width"].median().reindex(range(nb)).values / 1e6
    third = np.array_split(np.arange(nb), 3)
    summary.append((ch, len(c) / n,
                    np.nanmedian(np.concatenate([wid[third[0]], wid[third[2]]])), np.nanmedian(wid[third[1]]),
                    100 * (mid[third[0]].sum() + mid[third[2]].sum()) / max(mid.sum(), 1e-9)))
    allw.append(pd.DataFrame({"co": mid, "spread": spread, "markers": md, "genes": gd}))
    a = ax[k // ncol][k % ncol]
    top = max(mid.max(), spread.max(), 1e-9)
    a.plot(x, mid, color="black", lw=1.4, label="CO rate (midpoints)")
    a.plot(x, spread, color="#888780", lw=1.4, label="CO rate (spread over interval)")
    if md.max() > 0:
        a.plot(x, md / md.max() * top, color="#534AB7", ls="--", lw=1, label="marker calls (scaled)")
    if gd.max() > 0:
        a.plot(x, gd / gd.max() * top, color="#2E7D32", ls=":", lw=1.2, label="genes (scaled)")
    a.set_title("%s   %.2f COs/cell" % (ch, len(c) / n), fontsize=9)
    a.tick_params(labelsize=7)
    if k == 0:
        a.legend(frameon=False, fontsize=6)
    if k % ncol == 0:
        a.set_ylabel("cM / Mb", fontsize=8)
for k in range(len(chroms), nrow * ncol):
    ax[k // ncol][k % ncol].set_axis_off()
fig.suptitle("%s: crossover rate against where the data and genes are (%.0f Mb windows)" % (S, WMB))
fig.tight_layout()
fig.savefig(os.path.join(OUT, "%s_landscape_vs_data.png" % S), dpi=120)
fig.savefig(os.path.join(OUT, "%s_landscape_vs_data.pdf" % S))

A = pd.concat(allw)
lines = ["LANDSCAPE vs DATA  %s  (%d cells, %.0f Mb windows)" % (S, len(cells), WMB),
         "  Spearman across windows: CO rate vs marker calls %.2f | vs genes %.2f | markers vs genes %.2f"
         % (A.co.corr(A.markers, method="spearman"), A.co.corr(A.genes, method="spearman"),
            A.markers.corr(A.genes, method="spearman")),
         "  %-12s %9s %22s %22s %18s" % ("chrom", "COs/cell", "interval, outer 2/3rds", "interval, middle 3rd",
                                         "COs in outer 2/3")]
for ch, cpc, wo, wm, outer in summary:
    lines.append("  %-12s %9.2f %19.1f Mb %19.1f Mb %17.0f%%" % (ch, cpc, wo, wm, outer))
open(os.path.join(OUT, "%s_landscape_vs_data.txt" % S), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("wrote %s/%s_landscape_vs_data.{png,pdf,txt}" % (OUT, S))
