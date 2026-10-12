#!/usr/bin/env python3
"""molecule_tables.py -- hapCO input with ONE row per molecule, informative markers only.

The pipeline version of make_panel_inputs.py's good_mol input, which was scored
on the 50-cell review panel (block_support.py) before being adopted; the rule
that runs this checks the two agree byte for byte.
  - informative ('good') markers only, from marker_classes
  - markers < 150 bp apart usually ride on one read, so they are one
    observation: the run's deepest marker is kept with its own ref/alt counts,
    so no read is counted twice and hapCO's marker_num counts molecules
  - --exclude CSV (config/exclude_regions/SAMPLE.csv: chrom,start,end,reason;
    '#' lines are comments): molecules inside a listed region (start <= pos
    <= end, 1-based) are then dropped -- a stretch of the reference whose
    alleles follow another place in the genome paints false blocks. An empty
    value means none. score_verdicts.py --exclude applies the same rule.

Usage: molecule_tables.py SAMPLE CELLS_FILE OUTDIR [--exclude CSV]
"""
import argparse
import csv
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("cells")
ap.add_argument("outdir")
ap.add_argument("--exclude", default="")
A = ap.parse_args()
S, CELLS, out = A.sample, A.cells, A.outdir
GAP = 150
REGIONS = []
if A.exclude.strip():
    with open(A.exclude.strip()) as fh:
        for r in csv.DictReader(l for l in fh if not l.startswith("#")):
            REGIONS.append((r["chrom"], int(r["start"]), int(r["end"])))
os.makedirs(out, exist_ok=True)
cls = pd.read_csv("qc/markers/%s/marker_classes.tsv.gz" % S, sep="\t", usecols=["chrom", "pos", "class"])
g = cls[cls["class"] == "good"]
good = set(zip(g.chrom, g.pos.astype(int)))
cells = [l.strip() for l in open(CELLS) if l.strip()]
before = after = dropped = 0
for bc in cells:
    d = pd.read_csv("results/cell_data/%s/%s.tsv" % (S, bc), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "rc", "alt", "ac"])
    before += len(d)
    d = d[[(c, int(p)) in good for c, p in zip(d.chrom, d.pos)]]
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
    if REGIONS and len(d):
        drop = np.zeros(len(d), dtype=bool)
        for c, a, b in REGIONS:
            drop |= ((d.chrom == c) & (d.pos >= a) & (d.pos <= b)).values
        dropped += int(drop.sum())
        d = d[~drop]
    after += len(d)
    d.to_csv(os.path.join(out, "%s.tsv" % bc), sep="\t", header=False, index=False)
print("%s: %d cells, rows %s -> %s (%.0f%% kept) in %s"
      % (S, len(cells), format(before, ","), format(after, ","), 100.0 * after / max(before, 1), out))
if REGIONS:
    print("excluded regions (%s): %s -- %s molecules dropped"
          % (A.exclude.strip(), ", ".join("%s:%d-%d" % r for r in REGIONS), format(dropped, ",")))
