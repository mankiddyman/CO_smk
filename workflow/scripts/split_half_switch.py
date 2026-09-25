#!/usr/bin/env python3
"""split_half_switch.py -- is a low switch rate a property of the cell, or luck?

Each cell's chromosomes are split into two independent halves (odd- and even-
numbered). Select cells on half A (rate <= threshold), then read their rate on
half B. A genuinely clean nucleus stays clean on B. A cell that only looked
clean by sampling noise regresses to the population mean on B.

Uses the per-chromosome columns of switches.tsv, so it reads no per-cell data.
Usage: split_half_switch.py OUTDIR SAMPLE [SAMPLE ...]
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1]
SAMPLES = sys.argv[2:] or ["Dbinata_hap1", "Dparadoxa_hap1"]
THR = 0.10
os.makedirs(OUT, exist_ok=True)


def chrom_no(name):
    m = re.search(r"(\d+)", name)
    return int(m.group(1)) if m else 0


def load(s):
    out = []
    with open(os.path.join("results/cell_qc", s, "switches.tsv")) as f:
        h = f.readline().rstrip("\n").split("\t")
        chroms = [c[len("markers_"):] for c in h if c.startswith("markers_")]
        im = {c: h.index("markers_" + c) for c in chroms}
        isw = {c: h.index("switches_" + c) for c in chroms}
        ir = h.index("switch_rate")
        for l in f:
            t = l.rstrip("\n").split("\t")
            try:
                half = {0: [0, 0], 1: [0, 0]}
                for c in chroms:
                    m, sw = int(float(t[im[c]])), int(float(t[isw[c]]))
                    if m > 1:
                        k = chrom_no(c) % 2
                        half[k][0] += sw
                        half[k][1] += m - 1
                if half[0][1] >= 100 and half[1][1] >= 100:
                    out.append((float(t[ir]), half[1][0] / float(half[1][1]),
                                half[0][0] / float(half[0][1]), half[0][1] + half[1][1]))
            except (ValueError, IndexError):
                continue
    return out


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else float("nan")


fig, ax = plt.subplots(1, len(SAMPLES), figsize=(6 * len(SAMPLES), 5.5), squeeze=False)
print("select on half A (odd chromosomes) at <= %.2f, then read half B (even):" % THR)
print("  %-16s %8s %12s %14s %16s %18s" % ("sample", "cells", "median B", "sel. on A", "their median B",
                                           "B also <= %.2f" % THR))
for k, s in enumerate(SAMPLES):
    d = load(s)
    if not d:
        continue
    selA = [x for x in d if x[1] <= THR]
    print("  %-16s %8s %12.3f %14s %16.3f %17.1f%%"
          % (s, format(len(d), ","), med([x[2] for x in d]), format(len(selA), ","),
             med([x[2] for x in selA]) if selA else float("nan"),
             100.0 * sum(1 for x in selA if x[2] <= THR) / len(selA) if selA else 0))
    a = ax[0][k]
    a.scatter([x[1] for x in d], [x[2] for x in d], s=2, alpha=.15, color="#888780", rasterized=True)
    a.scatter([x[1] for x in selA], [x[2] for x in selA], s=3, alpha=.4, color="#2E7D32",
              rasterized=True, label="selected on A (<=%.2f)" % THR)
    a.axvline(THR, color="black", ls="--", lw=.8)
    a.axhline(THR, color="black", ls="--", lw=.8)
    a.plot([0, .3], [0, .3], color="#D85A30", ls=":", lw=1, label="perfect agreement")
    mb = med([x[2] for x in d])
    a.axhline(mb, color="#534AB7", lw=1, label="population median on B (%.3f)" % mb)
    a.set_xlim(0, .25); a.set_ylim(0, .25)
    a.set_xlabel("switch rate, odd chromosomes (half A)")
    a.set_ylabel("switch rate, even chromosomes (half B)")
    a.set_title("%s: does a clean half predict a clean other half?" % s.replace("_hap1", ""))
    a.legend(frameon=False, fontsize=8, loc="upper left")
fig.tight_layout()
png = os.path.join(OUT, "split_half.png")
fig.savefig(png, dpi=120); fig.savefig(png.replace(".png", ".pdf"))
print("\nwrote %s" % png)
