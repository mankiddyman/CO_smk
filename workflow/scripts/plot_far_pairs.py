#!/usr/bin/env python3
"""plot_far_pairs.py -- how many cells behave like ONE haploid genome?

The all-marker switch rate is mostly a measure of how many of a cell's
neighbouring markers sit on the same read (those agree automatically). The
honest test uses only pairs of informative markers far enough apart to lie on
different molecules: a haploid nucleus keeps them in agreement except at its
few crossovers; a barcode holding both haplotypes disagrees ~40-50% of the
time. Reads marker_segregation.py's cell_switches.tsv.gz per sample.

Usage: plot_far_pairs.py OUTDIR SAMPLE [SAMPLE ...]
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT, SAMPLES = sys.argv[1], sys.argv[2:]
COL = ["#D85A30", "#534AB7", "#2E7D32", "#888780"]
THR = [0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
MINP = [50, 200]
os.makedirs(OUT, exist_ok=True)


def load(s):
    d = pd.read_csv(os.path.join("qc/markers", s, "cell_switches.tsv.gz"), sep="\t")
    old = os.path.join("results/cell_qc", s, "switches.tsv")
    if os.path.exists(old):
        o = pd.read_csv(old, sep="\t", usecols=["barcode", "switch_rate"])
        d = d.merge(o.rename(columns={"switch_rate": "pipeline_rate"}), on="barcode", how="left")
    else:
        d["pipeline_rate"] = d["all_rate"]
    return d[d["good_far_rate"].notna()]


def valley(v, lo=0.05, hi=0.35, step=0.01):
    h, e = np.histogram(v, bins=np.arange(0, 0.6 + step, step))
    sm = np.convolve(h, np.ones(3) / 3, mode="same")
    mids = (e[:-1] + e[1:]) / 2
    cand = np.where((mids >= lo) & (mids <= hi))[0]
    return mids[cand[np.argmin(sm[cand])]] if len(cand) else np.nan


D = {s: load(s) for s in SAMPLES}
lines = []
say = lambda x="": (lines.append(x), print(x))
say("CELLS THAT BEHAVE LIKE ONE HAPLOID GENOME (informative markers, pairs on different molecules)")
say("  %-38s" % "" + "".join("%16s" % s.replace("_hap1", "") for s in SAMPLES))
row = lambda lab, f: say("  %-38s" % lab + "".join("%16s" % f(D[s]) for s in SAMPLES))
row("cells with >= 50 far pairs", lambda d: format(int((d.good_far_pairs >= 50).sum()), ","))
row("median far-pair rate", lambda d: "%.3f" % d.good_far_rate.median())
row("valley in the far-pair histogram", lambda d: "%.2f" % valley(d.good_far_rate))
for mp in MINP:
    say("  kept, needing >= %d far pairs:" % mp)
    for t in THR:
        row("    far-pair rate <= %.2f" % t,
            lambda d, t=t, mp=mp: format(int(((d.good_far_rate <= t) & (d.good_far_pairs >= mp)).sum()), ","))
say("  the old pipeline set (switch rate <= 0.10):")
row("    cells", lambda d: format(int((d.pipeline_rate <= 0.10).sum()), ","))
row("    ...of which far-pair rate <= 0.15", lambda d: format(
    int(((d.pipeline_rate <= 0.10) & (d.good_far_rate <= 0.15)).sum()), ","))
with open(os.path.join(OUT, "far_pairs.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")

fig, ax = plt.subplots(2, 2, figsize=(12, 9))
a = ax[0, 0]
for k, s in enumerate(SAMPLES):
    d = D[s]
    a.hist(d.good_far_rate, bins=np.arange(0, 0.61, 0.01), histtype="step", lw=2, color=COL[k],
           label="%s (n=%s)" % (s.replace("_hap1", ""), format(len(d), ",")))
a.set_yscale("log"); a.set_xlabel("switch rate, informative markers, pairs on different molecules")
a.set_ylabel("cells (log)"); a.legend(frameon=False, fontsize=8)
a.set_title("one genome (left) vs both haplotypes (~0.4-0.5)")
for t in (0.10, 0.20):
    a.axvline(t, color="grey", ls=":", lw=1)

a = ax[0, 1]
for k, s in enumerate(SAMPLES):
    d = D[s]
    a.scatter(d.good_far_pairs, d.good_far_rate, s=2, alpha=.25, color=COL[k], rasterized=True,
              label=s.replace("_hap1", ""))
    q = pd.qcut(np.log10(d.good_far_pairs), 20, duplicates="drop")
    g = d.groupby(q, observed=True)
    a.plot(g.good_far_pairs.median(), g.good_far_rate.median(), color=COL[k], lw=2.5)
a.set_xscale("log"); a.set_xlabel("far pairs in the cell (depth)"); a.set_ylabel("far-pair switch rate")
a.set_title("does depth make cells look cleaner? (lines: median)"); a.legend(frameon=False, fontsize=8)

a = ax[1, 0]
for k, s in enumerate(SAMPLES):
    d = D[s]
    a.scatter(d.pipeline_rate, d.good_far_rate, s=2, alpha=.25, color=COL[k], rasterized=True,
              label=s.replace("_hap1", ""))
a.axvline(0.10, color="black", ls="--", lw=.8); a.axhline(0.15, color="grey", ls=":", lw=1)
a.set_xlabel("pipeline switch rate (all markers)"); a.set_ylabel("far-pair switch rate")
a.set_title("where the old 'clean' cells really sit"); a.legend(frameon=False, fontsize=8)

a = ax[1, 1]
xs = np.arange(0, 0.41, 0.005)
for k, s in enumerate(SAMPLES):
    v = np.sort(D[s].good_far_rate[D[s].good_far_pairs >= 50].values)
    a.plot(xs, np.searchsorted(v, xs, side="right"), color=COL[k], lw=2, label=s.replace("_hap1", ""))
a.set_yscale("log"); a.set_xlabel("max far-pair switch rate"); a.set_ylabel("cells kept (>= 50 far pairs)")
for t in (0.10, 0.15, 0.20):
    a.axvline(t, color="grey", ls=":", lw=1)
a.legend(frameon=False, fontsize=8); a.set_title("cells kept as the honest threshold moves")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "far_pairs.png"), dpi=120); fig.savefig(os.path.join(OUT, "far_pairs.pdf"))
print("\nwrote %s/far_pairs.{png,pdf,txt}" % OUT)
