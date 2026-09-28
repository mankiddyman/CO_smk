#!/usr/bin/env python3
"""make_panel_inputs.py -- alternative hapCO inputs for the review panel's cells.

  good       informative markers only (class 'good' in marker_classes): drops the
             markers that read one allele in almost every cell, which can paint
             false blocks of that allele inside the other haplotype
  good_mol   informative markers, ONE row per molecule: markers < 150 bp apart
             usually ride on one read, so they are one observation, not several.
             The row kept is the run's deepest marker, with its own counts, so no
             read is counted twice. hapCO's marker_num then counts molecules.

Usage: make_panel_inputs.py SAMPLE SEED VARIANT
Output: qc/review/SAMPLE/inputs/VARIANT/<barcode>.tsv (hapCO's per-cell format)
"""
import os
import sys

import numpy as np
import pandas as pd

S, SEED, V = sys.argv[1], sys.argv[2], sys.argv[3]
assert V in ("good", "good_mol"), "VARIANT must be good or good_mol"
GAP = 150
out = os.path.join("qc/review", S, "inputs", V)
os.makedirs(out, exist_ok=True)
cls = pd.read_csv("qc/markers/%s/marker_classes.tsv.gz" % S, sep="\t", usecols=["chrom", "pos", "class"])
g = cls[cls["class"] == "good"]
good = set(zip(g.chrom, g.pos.astype(int)))
cells = [l.strip() for l in open(os.path.join("qc/review", S, "panel_seed%s.txt" % SEED)) if l.strip()]
before = after = 0
for bc in cells:
    d = pd.read_csv("results/cell_data/%s/%s.tsv" % (S, bc), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "rc", "alt", "ac"])
    before += len(d)
    d = d[[(c, int(p)) in good for c, p in zip(d.chrom, d.pos)]]
    if V == "good_mol":
        keep = []
        for ch, x in d.sort_values(["chrom", "pos"]).groupby("chrom", sort=False):
            p = x.pos.values
            new = np.ones(len(p), dtype=bool)
            new[1:] = (p[1:] - p[:-1]) >= GAP
            run = np.cumsum(new) - 1
            depth = (x.rc + x.ac).values
            best = {}
            for i, r in enumerate(run):
                if r not in best or depth[i] > depth[best[r]]:
                    best[r] = i
            keep.append(x.iloc[sorted(best.values())])
        d = pd.concat(keep) if keep else d.iloc[:0]
    d = d.sort_values(["chrom", "pos"], kind="stable")
    after += len(d)
    d.to_csv(os.path.join(out, "%s.tsv" % bc), sep="\t", header=False, index=False)
print("%s %s: %d panel cells, rows %s -> %s (%.0f%% kept) in %s"
      % (S, V, len(cells), format(before, ","), format(after, ","), 100.0 * after / max(before, 1), out))
