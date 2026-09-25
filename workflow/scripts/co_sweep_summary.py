#!/usr/bin/env python3
"""co_sweep_summary.py -- which hapCO settings recover crossovers without inventing them?

For each sample and each grid point (marker_num x block_size), from the
per-cell *_co_pred.txt files (one crossover per line: chrom start end ...):

  COs per cell                      rises as the caller gets more sensitive
  chromosomes below 0.5 COs/cell    the obligate-crossover floor: one CO per
                                    bivalent = 0.5 per chromosome per gamete
  map length (cM)                   100 x sum over chromosomes of COs/cell
  clean vs dirty cells              COs/cell in the cleanest cells
                                    (haploidness >= 0.85) against the noisiest
                                    admitted (0.6-0.7). Real crossovers are the
                                    same in both; noise inflates the dirty ones.
  close double crossovers           share of neighbouring COs on one chromosome
                                    closer than 10% of its length. Interference
                                    makes these rare; noise makes them common.

Usage: co_sweep_summary.py SWEEP_ROOT SAMPLE [SAMPLE ...]
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT, SAMPLES = sys.argv[1], sys.argv[2:]
rows = []
for s in SAMPLES:
    with open(os.path.join("results/cell_qc", s, "switches.tsv")) as f:
        chroms = [c[len("markers_"):] for c in f.readline().rstrip("\n").split("\t") if c.startswith("markers_")]
    hap = pd.read_csv(os.path.join("qc/haplotypes", s, "haplotype_tracks.tsv.gz"), sep="\t",
                      usecols=["barcode", "haploidness"]).set_index("barcode")["haploidness"]
    for d in sorted(glob.glob(os.path.join(ROOT, s, "mn*_bs*"))):
        m = re.search(r"mn(\d+)_bs(\d+)$", d)
        if not m or not os.path.exists(os.path.join(d, ".done")):
            continue
        mn, bs = int(m.group(1)), int(m.group(2))
        per_cell, cos = {}, []
        for p in glob.glob(os.path.join(d, "per_cell", "*_co_pred.txt")):
            bc = os.path.basename(p)[:-len("_co_pred.txt")]
            n = 0
            for l in open(p):
                f = l.split()
                if len(f) >= 3 and f[0] in chroms:
                    try:
                        cos.append((bc, f[0], (int(f[1]) + int(f[2])) / 2.0))
                        n += 1
                    except ValueError:
                        pass
            per_cell[bc] = n
        if not per_cell:
            continue
        ncell = len(per_cell)
        co = pd.DataFrame(cos, columns=["bc", "chrom", "mid"])
        per_chrom = co.groupby("chrom").size().reindex(chroms, fill_value=0) / float(ncell)
        # chromosome length ~ furthest CO midpoint seen across the grid would be biased; use max mid per chrom here
        clen = co.groupby("chrom")["mid"].max().reindex(chroms).fillna(1.0)
        close = total_pairs = 0
        for (bc, ch), g in co.groupby(["bc", "chrom"]):
            v = np.sort(g["mid"].values)
            if len(v) > 1:
                dd = np.diff(v)
                total_pairs += len(dd)
                close += int((dd < 0.1 * clen[ch]).sum())
        n = pd.Series(per_cell)
        h = hap.reindex(n.index)
        clean, dirty = n[h >= 0.85], n[(h >= 0.6) & (h < 0.7)]
        rows.append({"sample": s, "marker_num": mn, "block_size": bs, "cells": ncell,
                     "cos_per_cell": n.mean(), "median": n.median(),
                     "chrom_below_0.5": int((per_chrom < 0.5).sum()), "weakest_chrom": per_chrom.min(),
                     "map_cM": 100 * per_chrom.sum(),
                     "clean_cells": len(clean), "cos_clean": clean.mean(),
                     "dirty_cells": len(dirty), "cos_dirty": dirty.mean(),
                     "dirty_over_clean": dirty.mean() / clean.mean() if clean.mean() > 0 else np.nan,
                     "close_dco_share": close / float(total_pairs) if total_pairs else 0.0,
                     "floor": 0.5 * len(chroms)})
        per_chrom.rename("cos_per_cell").to_csv(os.path.join(d, "per_chromosome.tsv"), sep="\t")

t = pd.DataFrame(rows)
if t.empty:
    sys.exit("no finished grid points under %s" % ROOT)
t = t.sort_values(["sample", "block_size", "marker_num"])
pd.set_option("display.width", 200)
out = t[["sample", "marker_num", "block_size", "cells", "cos_per_cell", "floor", "chrom_below_0.5",
         "weakest_chrom", "map_cM", "cos_clean", "cos_dirty", "dirty_over_clean", "close_dco_share"]]
txt = out.to_string(index=False, float_format=lambda x: "%.2f" % x)
print(txt)
with open(os.path.join(ROOT, "sweep_summary.txt"), "w") as f:
    f.write(txt + "\n")
t.to_csv(os.path.join(ROOT, "sweep_summary.tsv"), sep="\t", index=False)

for s in t["sample"].unique():
    d = t[t["sample"] == s]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    styles = {bs: ls for bs, ls in zip(sorted(d.block_size.unique()), ["-", "--", ":", "-."])}
    for bs, g in d.groupby("block_size"):
        lab = "block_size %.1f Mb" % (bs / 1e6)
        ax[0, 0].plot(g.marker_num, g.cos_per_cell, "o" + styles[bs], color="#534AB7", label=lab)
        ax[0, 1].plot(g.marker_num, g["chrom_below_0.5"], "o" + styles[bs], color="#D85A30", label=lab)
        ax[1, 0].plot(g.marker_num, g.cos_clean, "o" + styles[bs], color="#2E7D32",
                      label="clean cells, " + lab)
        ax[1, 0].plot(g.marker_num, g.cos_dirty, "s" + styles[bs], color="#D85A30",
                      label="dirty cells, " + lab)
        ax[1, 1].plot(g.marker_num, g.close_dco_share, "o" + styles[bs], color="#888780", label=lab)
    ax[0, 0].axhline(d.floor.iloc[0], color="black", ls=":", lw=1, label="obligate-CO floor")
    ax[0, 0].set_ylabel("COs per cell"); ax[0, 0].set_title("%s: crossovers recovered" % s)
    ax[0, 1].set_ylabel("chromosomes below 0.5 COs/cell"); ax[0, 1].set_title("chromosomes under-called")
    ax[1, 0].set_ylabel("COs per cell"); ax[1, 0].set_title("noise check: clean vs dirty cells should match")
    ax[1, 1].set_ylabel("share of neighbouring COs < 10% of chromosome apart")
    ax[1, 1].set_title("noise check: close double crossovers")
    for a in ax.ravel():
        a.set_xlabel("marker_num (min markers per block)")
        a.legend(frameon=False, fontsize=7)
        a.invert_xaxis()
    fig.suptitle("%s: hapCO sensitivity (more sensitive to the right)" % s)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, s, "sweep.png"), dpi=120)
    fig.savefig(os.path.join(ROOT, s, "sweep.pdf"))
print("\nwrote %s/{sweep_summary.txt,.tsv} and %s/<sample>/sweep.{png,pdf}" % (ROOT, ROOT))
