#!/usr/bin/env python3
"""gene_marker_coverage.py -- how many genes actually carry a heterozygous marker?

Single-cell reads come from genes, so marker coverage that matters is per gene.
For genes on the given chromosomes: the share with >= 1 good marker inside the
gene, overall and split by whether the gene's 1 Mb window has an aligned
partner in the all-vs-all (100 kb pieces, >= 50 kb aligned). If genes in
windows WITHOUT a partner still mostly carry markers, piece-level 'unmatched'
overstates the gap; if they don't, mapping cells to both haplotypes would add them.

Usage: gene_marker_coverage.py ALLVSALL_DIR GFF3 MARKER_CLASSES.tsv.gz CHROM[,CHROM...]
"""
import collections
import sys

import numpy as np
import pandas as pd

AVA, GFF, MK, CHROMS = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(",")
P, W = 100000, 1000000
hits = collections.defaultdict(set)
for l in open(AVA + "/pieces_vs_both.paf"):
    f = l.split("\t")
    q, s0 = f[0].rsplit("__", 1)
    if f[5] != q and int(f[3]) - int(f[2]) >= P // 2:
        hits[(q, int(s0) // W)].add(f[5])
mk = pd.read_csv(MK, sep="\t", usecols=["chrom", "pos", "class"])
mk = mk[mk["class"] == "good"]
rows = []
for l in open(GFF):
    if l.startswith("#"):
        continue
    f = l.split("\t")
    if len(f) > 4 and f[2] == "gene" and f[0] in CHROMS:
        rows.append((f[0], int(f[3]), int(f[4])))
g = pd.DataFrame(rows, columns=["chrom", "start", "end"])
out = []
for c, x in g.groupby("chrom"):
    pos = np.sort(mk.pos[mk.chrom == c].values)
    n = np.searchsorted(pos, x.end.values, side="right") - np.searchsorted(pos, x.start.values, side="left")
    part = [any(h != c for h in hits.get((c, s // W), ())) for s in x.start]
    out.append(pd.DataFrame({"partner": part, "markers": n}))
r = pd.concat(out)
print("genes on %s: %d; share with >= 1 good marker inside the gene:" % (",".join(CHROMS), len(r)))
for lab, d in (("  all genes", r), ("  window has an aligned partner (100 kb pieces)", r[r.partner]),
               ("  window has NO aligned partner ('unmatched')", r[~r.partner])):
    print("%-50s %6d genes  %3.0f%% with a marker" % (lab, len(d), 100 * (d.markers > 0).mean() if len(d) else 0))
