#!/usr/bin/env python3
"""compare_co_runs.py -- what did a caller change add, and is it real?

Compares two sets of per-cell *_co_pred.txt calls for the same cells:
  COs per cell, and COs within EDGE_MB of a chromosome end (the end zone)
  chromosomes below the obligate 0.5 COs/cell
  noise check: COs GAINED per cell by haploidness tier within depth tertile --
  real crossovers are gained equally in clean and dirty cells at the same
  depth; noise is gained more in dirty cells.

Usage: compare_co_runs.py SAMPLE FAI OLD_PER_CELL_DIR NEW_PER_CELL_DIR EDGE_MB
"""
import os
import sys

import numpy as np
import pandas as pd

S, FAI, OLD, NEW, EDGE = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]) * 1e6
L = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"]).set_index("chrom")["len"]
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
ht = pd.read_csv("qc/haplotypes/%s/haplotype_tracks.tsv.gz" % S, sep="\t",
                 usecols=["barcode", "molecules", "haploidness"]).set_index("barcode")


def load(d):
    rows = []
    for bc in cells:
        p = os.path.join(d, "%s_co_pred.txt" % bc)
        if not os.path.exists(p):
            continue
        for l in open(p):
            f = l.split()
            if len(f) >= 3 and f[0] in L.index:
                try:
                    m = (int(f[1]) + int(f[2])) / 2.0
                except ValueError:
                    continue
                rows.append((bc, f[0], m, m < EDGE or m > L[f[0]] - EDGE))
    return pd.DataFrame(rows, columns=["bc", "chrom", "mid", "edge"])


o, n = load(OLD), load(NEW)
N = float(len(cells))
print("COMPARE  %s  (%d cells; end zone = %.1f Mb from each end)" % (S, len(cells), EDGE / 1e6))
print("  COs per cell           %6.2f -> %6.2f" % (len(o) / N, len(n) / N))
print("  in the end zone        %6.2f -> %6.2f" % (o.edge.sum() / N, n.edge.sum() / N))
print("  away from the ends     %6.2f -> %6.2f" % ((~o.edge).sum() / N, (~n.edge).sum() / N))
po = o.groupby("chrom").size().reindex(L.index).fillna(0) / N
pn = n.groupby("chrom").size().reindex(L.index).fillna(0) / N
keep = pn.index[(pn > 0) | (po > 0)]
print("  chromosomes < 0.5      %6d -> %6d" % ((po[keep] < 0.5).sum(), (pn[keep] < 0.5).sum()))
print("  %-12s %8s %8s" % ("chrom", "before", "after"))
for ch in keep:
    print("  %-12s %8.2f %8.2f%s" % (ch, po[ch], pn[ch], "   < 0.5" if pn[ch] < 0.5 else ""))
g = (n.groupby("bc").size().reindex(cells).fillna(0) - o.groupby("bc").size().reindex(cells).fillna(0))
t = ht.reindex(cells).assign(gain=g.values)
t["tier"] = pd.cut(t.haploidness, [0.6, 0.7, 0.8, 0.9, 1.0], include_lowest=True)
t["depth"] = pd.qcut(t.molecules, 3, labels=["shallow", "middle", "deep"])
piv = t.pivot_table(index="tier", columns="depth", values="gain", aggfunc="mean", observed=False)
cnt = t.pivot_table(index="tier", columns="depth", values="gain", aggfunc="size", observed=False)
print("  COs GAINED per cell, haploidness tier (rows) x depth tertile (columns); cells in brackets")
for tier in piv.index:
    print("    %-12s " % str(tier) + "  ".join("%6.2f (%4d)" % (piv.loc[tier, c], cnt.loc[tier, c]) for c in piv.columns))
