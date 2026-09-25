#!/usr/bin/env python3
"""chrom_haploidness.py -- haploidness per chromosome, for every selected cell.

Same definition as rule cell_haploidness (good markers only; calls < 150 bp
apart = one molecule; consecutive 15-molecule windows; mean |2f - 1|), but
computed separately for each chromosome. A cell can be clean genome-wide and
still carry both haplotypes on ONE chromosome (an extra copy from an
unbalanced translocation gamete, aneuploidy, a local contaminant); the
genome-wide score barely notices, the per-chromosome one does.

Usage: chrom_haploidness.py SAMPLE OUTDIR [BARCODE ...]   (barcodes to print in full)
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

S, OUT = sys.argv[1], sys.argv[2]
SHOW = sys.argv[3:]
W, GAP, MINW = 15, 150, 2
os.makedirs(OUT, exist_ok=True)

cls = pd.read_csv("qc/markers/%s/marker_classes.tsv.gz" % S, sep="\t", usecols=["chrom", "pos", "class"])
g = cls[cls["class"] == "good"]
good = set(zip(g.chrom, g.pos.astype(int)))
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]

rows = []
for bc in cells:
    d = pd.read_csv("results/cell_data/%s/%s.tsv" % (S, bc), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "rc", "alt", "ac"])
    d = d[[(c, int(p)) in good for c, p in zip(d.chrom, d.pos)]]
    dp = d.rc + d.ac
    fr = d.ac / dp.clip(lower=1)
    d = d.assign(call=np.where(fr <= 0.2, 0, np.where(fr >= 0.8, 1, -1)))
    d = d[(dp > 0) & (d.call >= 0)].sort_values(["chrom", "pos"])
    for ch, x in d.groupby("chrom", sort=False):
        p, c = x.pos.values, x.call.values.astype(float)
        new = np.ones(len(p), dtype=bool)
        new[1:] = (p[1:] - p[:-1]) >= GAP
        mid = np.cumsum(new) - 1
        mf = np.bincount(mid, weights=c) / np.bincount(mid)
        mg = (mf[mf != 0.5] > 0.5).astype(float)
        nwin = len(mg) // W
        h = np.mean(np.abs(2 * mg[:nwin * W].reshape(nwin, W).mean(1) - 1)) if nwin >= MINW else np.nan
        rows.append((bc, ch, len(mg), nwin, h))
t = pd.DataFrame(rows, columns=["barcode", "chrom", "molecules", "windows", "haploidness"])
t.to_csv(os.path.join(OUT, "%s_chrom_haploidness.tsv.gz" % S), sep="\t", index=False)

scored = t.dropna(subset=["haploidness"])
worst = scored.groupby("barcode").haploidness.min()
lines = ["PER-CHROMOSOME HAPLOIDNESS  %s  (%d selected cells; chromosomes need >= %d windows to be scored)"
         % (S, len(cells), MINW),
         "  cell x chromosome scored: %d of %d" % (len(scored), len(t)),
         "  cells whose WORST chromosome is below:"]
for thr in (0.3, 0.4, 0.5, 0.6):
    lines.append("    %.1f   %4d cells (%.0f%%)" % (thr, (worst < thr).sum(), 100.0 * (worst < thr).mean()))
lines.append("  per chromosome: cells scored, median, share below 0.5")
for ch, x in scored.groupby("chrom", sort=False):
    lines.append("    %-12s %5d  %.2f  %5.1f%%" % (ch, len(x), x.haploidness.median(), 100 * (x.haploidness < 0.5).mean()))
for bc in SHOW:
    x = t[t.barcode == bc]
    lines.append("  %s:" % bc)
    for _, r in x.iterrows():
        lines.append("    %-12s %5d molecules  %3d windows  haploidness %s"
                     % (r.chrom, r.molecules, r.windows, "%.2f" % r.haploidness if r.haploidness == r.haploidness else "n/a"))
open(os.path.join(OUT, "%s_chrom_haploidness.txt" % S), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ax[0].hist(scored.haploidness, bins=np.linspace(0, 1, 51), color="#534AB7")
ax[0].set_xlabel("haploidness of one chromosome in one cell"); ax[0].set_ylabel("cell x chromosome")
ax[0].set_title("%s: per-chromosome haploidness, selected cells" % S)
ax[1].hist(worst, bins=np.linspace(0, 1, 51), color="#D85A30")
for v in (0.5, 0.6):
    ax[1].axvline(v, color="black", ls=":", lw=1)
ax[1].set_xlabel("each cell's WORST chromosome"); ax[1].set_ylabel("cells")
ax[1].set_title("cells hiding a both-haplotype chromosome sit on the left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "%s_chrom_haploidness.png" % S), dpi=120)
print("wrote %s/%s_chrom_haploidness.{tsv.gz,txt,png}" % (OUT, S))
