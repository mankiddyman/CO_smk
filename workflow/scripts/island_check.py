#!/usr/bin/env python3
"""island_check.py -- do the islands a smaller block_size adds behave like crossovers?

Reads the benchmark's untouched-cell islands that appear only below 2 Mb
(qc/benchmark/SAMPLE/new_islands.tsv) and asks four things per block size:
  WHERE    share of them in the windows that hold 80% of the pipeline's
           crossovers (the hot ends). Real crossovers sit there by
           definition; noise lands wherever there are molecules -- and the
           molecules are gene-rich-end heavy too, so the null is the share of
           informative MOLECULES in those windows, not of the genome.
  SUPPORT  informative molecules inside, and the share of those carrying the
           island's allele. Well supported = >= 8 molecules, >= 80% agreeing.
  FLICKER  share whose neighbouring block is also < 1 Mb: the call flickering
           around ONE crossover, which adds two false crossovers there.
  RECUR    0.5 Mb windows where islands recur in >= 3 cells -- structure, not
           meiosis. Known hotspots are marked, with the extent of their
           islands across cells, to size a mask.

Usage: island_check.py SAMPLE --fai GENOME.fa.fai [--hotspots CHROM:START-END,...]
Output: qc/benchmark/SAMPLE/island_check.txt
"""
import argparse
import collections
import glob
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("--hotspots", default="")
ap.add_argument("--fai", required=True)
A = ap.parse_args()
S, W, IW = A.sample, 1000000, 500000
B = "qc/benchmark/%s" % S
HOT = []
for h in filter(None, A.hotspots.split(",")):
    c, r = h.split(":")
    HOT.append((c, int(r.split("-")[0]), int(r.split("-")[1])))
nt = pd.read_csv(os.path.join(B, "new_islands.tsv"), sep="\t")

# where the pipeline's crossovers are: 1 Mb windows holding 80% of them
cnt = collections.Counter()
for p in glob.glob("results/crossovers/%s/per_cell/*_co_pred.txt" % S):
    for l in open(p):
        f = l.split()
        if len(f) >= 3 and f[1].isdigit():
            cnt[(f[0], int((int(f[1]) + int(f[2])) / 2 // W))] += 1
tot = float(sum(cnt.values()))
hot, acc = set(), 0.0
for k, v in sorted(cnt.items(), key=lambda kv: -kv[1]):
    if acc >= 0.8 * tot:
        break
    hot.add(k)
    acc += v
chroms_seen = {c for c, _ in cnt}
n_win = sum(int(np.ceil(int(l.split()[1]) / float(W))) for l in open(A.fai) if l.split()[0] in chroms_seen)


def molecules(bc):
    d = pd.read_csv("results/cell_data_mol/%s/%s.tsv" % (S, bc), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "rc", "alt", "ac"])
    fr = d.ac / (d.rc + d.ac).clip(lower=1)
    return d.assign(call=np.where(fr <= 0.2, 0, np.where(fr >= 0.8, 1, -1)))


def blocks(p):
    rows = []
    if os.path.exists(p):
        for l in open(p):
            f = l.split()
            if len(f) >= 4 and f[1].isdigit() and f[2].isdigit():
                rows.append((f[0], int(f[1]), int(f[2]), f[3]))
    return rows


mol_cache = {}
rows = []
for r in nt.itertuples():
    if r.barcode not in mol_cache:
        mol_cache[r.barcode] = molecules(r.barcode)
    m = mol_cache[r.barcode]
    x = m.call[(m.chrom == r.chrom) & (m.pos >= r.start) & (m.pos <= r.end) & (m.call >= 0)]
    want = 1 if r.allele == "ALT" else 0
    agree = float((x == want).mean()) if len(x) else 0.0
    bl = sorted([b for b in blocks(os.path.join(B, "runs", "base_bs%d" % r.block_size, r.barcode + "_co_block_pred.txt"))
                 if b[0] == r.chrom], key=lambda b: b[1])
    idx = [i for i, b in enumerate(bl) if b[1] == r.start and b[2] == r.end]
    flick = False
    if idx:
        i = idx[0]
        for j in (i - 1, i + 1):
            if 0 < j < len(bl) - 1 and bl[j][2] - bl[j][1] < 1000000:
                flick = True
    mid = (r.start + r.end) / 2.0
    rows.append((r.block_size, r.barcode, r.chrom, r.start, r.end, len(x), agree,
                 len(x) >= 8 and agree >= 0.8, flick, (r.chrom, int(mid // W)) in hot, r.hotspot))
t = pd.DataFrame(rows, columns=["block_size", "barcode", "chrom", "start", "end", "n", "agree", "well", "flicker",
                                "in_hot", "hotspot"])
bench_cells = sorted({os.path.basename(p).split("_co_")[0]
                      for p in glob.glob(os.path.join(B, "runs", "base_bs2000000", "*_co_block_pred.txt"))})
ncell = len(bench_cells) or 100
m_in = m_all = 0
for bc in bench_cells:
    if bc not in mol_cache:
        mol_cache[bc] = molecules(bc)
    m = mol_cache[bc]
    m = m[m.call >= 0]
    m_all += len(m)
    m_in += sum(1 for c, p in zip(m.chrom, m.pos) if (c, int(p // W)) in hot)
mol_share = m_in / float(max(m_all, 1))

lines = ["ISLAND CHECK  %s  (islands that appear only below 2 Mb, untouched cells; %d cells)" % (S, ncell),
         "  hot windows: %d of ~%d 1 Mb windows (%.0f%% of the genome) hold 80%% of the pipeline's crossovers"
         % (len(hot), n_win, 100.0 * len(hot) / max(n_win, 1)),
         "    %-10s %9s %9s %12s %14s %10s %9s" % ("block_size", "new/cell", "in hot", "median mol",
                                                   "well supported", "flicker", "hotspot")]
for bs, x in t.groupby("block_size", sort=False):
    lines.append("    %-10s %9.2f %8.0f%% %12.0f %13.0f%% %9.0f%% %8.0f%%"
                 % ("%.2g Mb" % (bs / 1e6) if bs > 1 else "off", len(x) / float(ncell), 100 * x.in_hot.mean(),
                    x.n.median(), 100 * x.well.mean(), 100 * x.flicker.mean(), 100 * x.hotspot.mean()))
lines.append("  READ: crossovers are 80%% in hot windows by construction; informative molecules %.0f%%. New islands"
             % (100 * mol_share))
lines.append("  near 80% follow crossovers; near the molecule share they follow the data, i.e. noise. Well supported")
lines.append("  high, flicker low, nothing recurring outside the masked hotspots.")

lines.append("  RECURRING positions (0.5 Mb windows with islands in >= 3 cells), per block size:")
for bs, x in t.groupby("block_size", sort=False):
    rec = collections.defaultdict(set)
    for r in x.itertuples():
        for k in range(int(r.start // IW), int(r.end // IW) + 1):
            rec[(r.chrom, k)].add(r.barcode)
    top = sorted(((len(v), c, k) for (c, k), v in rec.items() if len(v) >= 3), reverse=True)[:8]
    lines.append("    %-8s %s" % ("%.2g Mb" % (bs / 1e6) if bs > 1 else "off",
                                  "; ".join("%s %.1f-%.1f Mb x%d%s" % (c, k * IW / 1e6, (k + 1) * IW / 1e6, n,
                                            " [hotspot]" if any(hc == c and min(he, (k + 1) * IW) > max(hs, k * IW)
                                                                for hc, hs, he in HOT) else "")
                                            for n, c, k in top) or "none"))
lines.append("  HOTSPOT extents (islands overlapping each, all block sizes pooled): start / end, 10th-90th pct")
for hc, hs, he in HOT:
    x = t[(t.chrom == hc) & (t.end > hs) & (t.start < he)]
    if len(x):
        lines.append("    %s:%.2f-%.2f Mb  %d islands in %d cells  start %.2f-%.2f  end %.2f-%.2f Mb"
                     % (hc, hs / 1e6, he / 1e6, len(x), x.barcode.nunique(), x.start.quantile(.1) / 1e6,
                        x.start.quantile(.9) / 1e6, x.end.quantile(.1) / 1e6, x.end.quantile(.9) / 1e6))
    else:
        lines.append("    %s:%.2f-%.2f Mb  no islands" % (hc, hs / 1e6, he / 1e6))
t.to_csv(os.path.join(B, "island_check.tsv"), sep="\t", index=False)
open(os.path.join(B, "island_check.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
