#!/usr/bin/env python3
"""landscape_figure.py -- crossover landscape with a y-axis scaled to the real signal.

The pipeline's landscape takes its y-axis from the tallest peak, so one or two
narrow peaks squash everything else to what looks like zero. Here the axis is
capped from the data (1.3 x the 99th percentile of the smoothed rate, rounded
up to 0.25 cM/Mb); anything taller is drawn to the cap and labelled with its
true height, so nothing is hidden. Each crossover is shown as a tick along the
bottom, and each panel states its crossovers per pollen, so low stretches read
as low, not as absent.

Rate: every crossover's interval (start-end from hapCO) is spread over the
1 Mb bins it overlaps, summed over cells, smoothed with a 5 Mb running mean:
cM/Mb = 100 x crossovers / (cells x Mb). Band: 5-95% of 200 bootstrap
resamples of the cells.

Usage: landscape_figure.py SAMPLE OUT_PREFIX [--title TEXT] [--note TEXT]
Writes OUT_PREFIX.png and OUT_PREFIX.pdf
"""
import argparse
import math
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("sample"); ap.add_argument("out")
ap.add_argument("--title", default="")
ap.add_argument("--note", default="")
ap.add_argument("--smooth_mb", type=int, default=5)
A = ap.parse_args()
S = A.sample
os.makedirs(os.path.dirname(A.out) or ".", exist_ok=True)

import csv
_rows = list(csv.DictReader(open("config/samples.csv")))
_row = next((r for r in _rows if S in (r.get("sample_id"), r.get("sample"), list(r.values())[0])), None)
if _row is None:
    raise SystemExit("%s not found in config/samples.csv" % S)
fai = _row["assembly_fasta"] + ".fai"
chroms = [l.split("\t")[0] for l in open(fai)]
lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
nb = {c: int(math.ceil(lens[c] / 1e6)) for c in chroms}
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
per = {c: np.zeros((len(cells), nb[c])) for c in chroms}
mids = {c: [] for c in chroms}
for ci, bc in enumerate(cells):
    p = "results/crossovers/%s/per_cell/%s_co_pred.txt" % (S, bc)
    if not os.path.exists(p):
        continue
    for l in open(p):
        f = l.split()
        if len(f) < 3 or f[0] not in lens or not f[1].isdigit():
            continue
        c, s, e = f[0], int(f[1]), int(f[2])
        s, e = max(0, min(s, e)), min(lens[c], max(s, e))
        mids[c].append((s + e) / 2.0)
        b0, b1 = int(s // 1e6), int(min(e, lens[c] - 1) // 1e6)
        if b1 <= b0 or e - s < 1e3:
            per[c][ci, int(((s + e) / 2.0) // 1e6)] += 1.0
        else:
            for b in range(b0, b1 + 1):
                ov = min(e, (b + 1) * 1e6) - max(s, b * 1e6)
                per[c][ci, b] += ov / float(e - s)
n = len(cells)
k = np.ones(A.smooth_mb) / A.smooth_mb


def rate(w):     # per-bin crossovers summed over cells -> smoothed cM/Mb
    return 100.0 * np.convolve(w, k, mode="same") / n


rng = np.random.default_rng(1)
boot = [rng.integers(0, n, n) for _ in range(200)]
curve, lo, hi = {}, {}, {}
for c in chroms:
    curve[c] = rate(per[c].sum(0))
    bs = np.array([rate(per[c][i].sum(0)) for i in boot])
    lo[c], hi[c] = np.percentile(bs, 5, axis=0), np.percentile(bs, 95, axis=0)
allv = np.concatenate([curve[c] for c in chroms])
cap = max(0.5, math.ceil(1.3 * np.percentile(allv, 99) / 0.25) * 0.25)
tot = sum(len(mids[c]) for c in chroms)
gmean = 100.0 * tot / (n * sum(lens.values()) / 1e6)

fig, axes = plt.subplots(len(chroms), 1, figsize=(11, 1.75 * len(chroms) + 1.6), sharex=True, sharey=True)
xmax = max(lens.values()) / 1e6
for ax, c in zip(axes, chroms):
    x = np.arange(nb[c]) + 0.5
    ax.fill_between(x, np.minimum(lo[c], cap), np.minimum(hi[c], cap), color="#B4B2A9", alpha=0.45, lw=0)
    ax.plot(x, np.minimum(curve[c], cap), color="black", lw=1.1)
    cm = 100.0 * len(mids[c]) / (n * lens[c] / 1e6)
    ax.axhline(cm, color="#2E7D32", lw=0.9, ls="--")
    ax.axhline(gmean, color="#C2185B", lw=0.8, ls=":")
    ax.axvline(lens[c] / 1e6, color="#534AB7", lw=1.2)
    ax.plot(np.array(mids[c]) / 1e6, np.full(len(mids[c]), -0.06 * cap), "|", color="#D85A30", ms=6, mew=0.6, clip_on=False)
    over = curve[c] > cap                        # clipped peaks: labelled with their true height
    i = 0
    while i < len(over):
        if over[i]:
            j = i
            while j + 1 < len(over) and over[j + 1]:
                j += 1
            top = i + int(np.argmax(curve[c][i:j + 1]))
            ax.plot(top + 0.5, cap * 0.97, marker="^", color="#D85A30", ms=6)
            right = top > 0.75 * xmax
            ax.text(top + (-3 if right else 3), cap * 0.93, "%.1f cM/Mb (clipped)" % curve[c][top], fontsize=7.5,
                    color="#D85A30", va="top", ha="right" if right else "left")
            i = j + 1
        else:
            i += 1
    ax.set_title("%s    %d crossovers  |  %.2f per pollen  |  %.2f cM/Mb" % (c, len(mids[c]), len(mids[c]) / float(n), cm),
                 loc="left", fontsize=8.5, pad=2)
    ax.set_ylim(-0.1 * cap, cap)
    ax.set_xlim(0, xmax * 1.01)
    ax.spines[["top", "right"]].set_visible(False)
axes[-1].set_xlabel("position (Mb)")
fig.text(0.015, 0.5, "crossover rate (cM/Mb)", rotation=90, va="center")
title = A.title or "Crossover landscape -- %s" % S
fig.suptitle("%s\n%d pollen, %d crossovers (%.2f per pollen), genome-wide %.2f cM/Mb; y-axis capped at %.2f cM/Mb"
             % (title, n, tot, tot / float(n), gmean, cap), fontsize=10, y=0.995)
handles = [plt.Line2D([], [], color="black", lw=1.1, label="rate (5 Mb smoothed)"),
           plt.Rectangle((0, 0), 1, 1, color="#B4B2A9", alpha=0.45, label="5-95% bootstrap over pollen"),
           plt.Line2D([], [], color="#2E7D32", ls="--", label="chromosome mean"),
           plt.Line2D([], [], color="#C2185B", ls=":", label="genome-wide mean"),
           plt.Line2D([], [], color="#D85A30", marker="|", ls="", ms=8, label="one crossover"),
           plt.Line2D([], [], color="#534AB7", lw=1.2, label="chromosome end")]
fig.legend(handles=handles, loc="upper center", ncol=6, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, 0.955))
if A.note:
    fig.text(0.5, 0.005, A.note, ha="center", va="bottom", fontsize=8, style="italic", wrap=True)
fig.tight_layout(rect=(0.03, 0.03 if A.note else 0.0, 1, 0.945))
for ext in ("png", "pdf"):
    fig.savefig("%s.%s" % (A.out, ext), dpi=200 if ext == "png" else None)
print("cap %.2f cM/Mb; wrote %s.png and %s.pdf" % (cap, A.out, A.out))
