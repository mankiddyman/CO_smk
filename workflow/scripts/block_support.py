#!/usr/bin/env python3
"""block_support.py -- does each called block have informative molecules behind it?

A called haplotype block is only as good as the reads under it. For every block
hapCO calls (<cell>_co_block_pred.txt), count the cell's informative molecules
inside it (good markers only; calls < 150 bp apart = one molecule) and ask
whether they carry the block's allele:
  supported     >= 3 molecules, >= 60% carrying the block's allele
  unsupported   fewer than 3 molecules, or most carry the other allele
An unsupported ISLAND (an interior block whose two neighbours are the other
haplotype) is a false double crossover: two false COs. Real double crossovers
go both ways (REF-in-ALT as often as ALT-in-REF); false islands from markers
that read one allele regardless do not -- so islands are split by orientation.
At each chromosome END the outermost 8 molecules are compared with the end
block: if >= 6 carry the other allele, a terminal crossover was missed.

Molecules always come from the pipeline's full per-cell tables, whatever input
the caller was given, so every calling variant is judged against the same data.

Usage: block_support.py SAMPLE PER_CELL_DIR OUTDIR [--cells FILE] [--show BARCODE ...] [--label NAME]
"""
import argparse
import collections
import glob
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("sample"); ap.add_argument("per_cell"); ap.add_argument("outdir")
ap.add_argument("--cells"); ap.add_argument("--show", nargs="*", default=[]); ap.add_argument("--label", default="")
a = ap.parse_args()
S, GAP, MIN_MOL, AGREE, END_N, END_K = a.sample, 150, 3, 0.6, 8, 6
os.makedirs(a.outdir, exist_ok=True)

cls = pd.read_csv("qc/markers/%s/marker_classes.tsv.gz" % S, sep="\t", usecols=["chrom", "pos", "class"])
g = cls[cls["class"] == "good"]
good = set(zip(g.chrom, g.pos.astype(int)))
if a.cells:
    cells = [l.strip() for l in open(a.cells) if l.strip()]
else:
    cells = [os.path.basename(p)[:-len("_co_block_pred.txt")]
             for p in glob.glob(os.path.join(a.per_cell, "*_co_block_pred.txt"))]


def molecules(bc):
    d = pd.read_csv("results/cell_data/%s/%s.tsv" % (S, bc), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "rc", "alt", "ac"])
    d = d[[(c, int(p)) in good for c, p in zip(d.chrom, d.pos)]]
    dp = d.rc + d.ac
    fr = d.ac / dp.clip(lower=1)
    d = d.assign(call=np.where(fr <= 0.2, 0, np.where(fr >= 0.8, 1, -1)))
    d = d[(dp > 0) & (d.call >= 0)].sort_values(["chrom", "pos"])
    out = {}
    for ch, x in d.groupby("chrom", sort=False):
        p, c = x.pos.values, x.call.values.astype(float)
        new = np.ones(len(p), dtype=bool)
        new[1:] = (p[1:] - p[:-1]) >= GAP
        mid = np.cumsum(new) - 1
        mf = np.bincount(mid, weights=c) / np.bincount(mid)
        keep = mf != 0.5
        out[ch] = (p[new][keep], (mf[keep] > 0.5).astype(int))
    return out


def blocks(bc):
    rows = []
    p = os.path.join(a.per_cell, "%s_co_block_pred.txt" % bc)
    if not os.path.exists(p):
        return None
    for l in open(p):
        f = l.split()
        if len(f) >= 4 and f[1].isdigit() and f[2].isdigit():
            rows.append((f[0], int(f[1]), int(f[2]), f[3]))
    return rows


# pass 1: which hapCO genotype code is ALT? learned from big blocks
data, fits = {}, collections.defaultdict(list)
for bc in cells:
    b = blocks(bc)
    if b is None:
        continue
    m = molecules(bc)
    data[bc] = (b, m)
    for ch, s0, e0, gt in b:
        if ch in m:
            pos, al = m[ch]
            sel = (pos >= s0) & (pos <= e0)
            if sel.sum() >= 10:
                fits[gt].append(al[sel].mean())
alt_code = max(fits, key=lambda k: np.mean(fits[k]))
means = {k: np.mean(v) for k, v in fits.items()}

rows, show_lines = [], []
for bc, (b, m) in data.items():
    per = collections.defaultdict(list)
    for r in b:
        per[r[0]].append(r)
    isl = uns_isl = uns_int = ends_missed = ends_checked = 0
    orient = collections.Counter()
    n_cos = 0
    for ch, bl in per.items():
        bl.sort(key=lambda r: r[1])
        n_cos += sum(1 for i in range(1, len(bl)) if bl[i][3] != bl[i - 1][3])
        pos, al = m.get(ch, (np.array([], dtype=int), np.array([], dtype=int)))
        for i, (_, s0, e0, gt) in enumerate(bl):
            want = 1 if gt == alt_code else 0
            sel = (pos >= s0) & (pos <= e0)
            n = int(sel.sum())
            agree = float((al[sel] == want).mean()) if n else 0.0
            ok = n >= MIN_MOL and agree >= AGREE
            interior = 0 < i < len(bl) - 1
            island = interior and bl[i - 1][3] == bl[i + 1][3] != gt
            if interior and not ok:
                uns_int += 1
            if island:
                isl += 1
                if not ok:
                    uns_isl += 1
                    orient["REF-in-ALT" if want == 0 else "ALT-in-REF"] += 1
            if bc in a.show:
                show_lines.append("    %s %-11s %6.1f-%6.1f Mb  %s  %4d mol  %3.0f%% agree  %s%s"
                                  % (bc, ch, s0 / 1e6, e0 / 1e6, "ALT" if want else "REF", n, 100 * agree,
                                     "ok" if ok else "UNSUPPORTED", "  (island)" if island else ""))
        for end_block, idx in ((bl[0], slice(0, END_N)), (bl[-1], slice(-END_N, None))):
            if len(pos) >= END_N:
                ends_checked += 1
                want = 1 if end_block[3] == alt_code else 0
                if int((al[idx] != want).sum()) >= END_K:
                    ends_missed += 1
                    if bc in a.show:
                        show_lines.append("    %s %-11s end at %.1f Mb: %d of the outermost %d molecules carry the "
                                          "other allele -> MISSED terminal switch"
                                          % (bc, ch, (pos[idx][0] if idx.start == 0 else pos[idx][-1]) / 1e6,
                                             int((al[idx] != want).sum()), END_N))
    rows.append((bc, n_cos, isl, uns_isl, uns_int, ends_missed, ends_checked,
                 orient["REF-in-ALT"], orient["ALT-in-REF"]))

t = pd.DataFrame(rows, columns=["barcode", "cos", "islands", "unsupported_islands", "unsupported_interior",
                                "ends_missed", "ends_checked", "uns_ref_in_alt", "uns_alt_in_ref"])
t.to_csv(os.path.join(a.outdir, "block_support_cells.tsv"), sep="\t", index=False)
N = float(max(len(t), 1))
lines = ["BLOCK SUPPORT  %s  %s  (%d cells; hapCO genotype %s = ALT, mean ALT share %s)"
         % (S, a.label or a.per_cell, len(t), alt_code,
            ", ".join("%s %.2f" % (k, v) for k, v in sorted(means.items()))),
         "  per cell: COs %.2f | islands (double COs) %.2f | UNSUPPORTED islands %.2f "
         "(REF-in-ALT %.2f, ALT-in-REF %.2f) | missed ends %.2f of %.1f checked"
         % (t.cos.mean(), t.islands.mean(), t.unsupported_islands.mean(), t.uns_ref_in_alt.mean(),
            t.uns_alt_in_ref.mean(), t.ends_missed.mean(), t.ends_checked.mean()),
         "  estimated false COs per cell (2 per unsupported island): %.2f -> supported COs ~%.2f per cell"
         % (2 * t.unsupported_islands.mean(), t.cos.mean() - 2 * t.unsupported_islands.mean())]
if show_lines:
    lines.append("  blocks in the requested cells:")
    lines += show_lines
lines.append("ROW\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (a.label or a.per_cell, t.cos.mean(), t.islands.mean(),
             t.unsupported_islands.mean(), t.uns_ref_in_alt.mean(), t.uns_alt_in_ref.mean(), t.ends_missed.mean()))
open(os.path.join(a.outdir, "block_support_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
