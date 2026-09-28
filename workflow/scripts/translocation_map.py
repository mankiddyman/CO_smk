#!/usr/bin/env python3
"""translocation_map.py -- where does each piece of hap2 sit on hap1, and what does that do to the markers?

From the hap2-on-hap1 alignment (every hap2 chromosome aligned to ALL of hap1):
  1. SYNTENY BLOCKS, hap2 coordinates: alignments from one hap2 chromosome to
     one hap1 chromosome on one strand, merged while collinear (gaps <= 5 Mb in
     both). Blocks >= 1 Mb aligned are listed. A block on a differently
     numbered hap1 chromosome is a translocation; a block on the chromosome's
     minority strand is an inversion.
  2. HOMOLOG SOURCE, hap1 coordinates (the ones crossovers are called in), per
     1 Mb window: share covered by the same-numbered hap2 chromosome, by
     another (which one), and by nothing. Runs >= 2 Mb dominated by another
     chromosome are FOREIGN: their alleles travel with the partner chromosome,
     not with their neighbours. Windows covered by nothing have no homolog to
     carry a heterozygous marker.
  3. MARKERS: informative ('good') markers per Mb in each class of window, and
     the marker deserts (>= 5 Mb below 10% of the chromosome's median) with
     what they are made of.
Figures: a dot plot per hap2 chromosome against the whole hap1 genome; hap1
tracks of homolog source with marker density on top.

Usage: translocation_map.py SAMPLE HAP1_FAI OUTDIR [--min_mapq 5] [--min_len 20000]
Reads: results/hap_align/SAMPLE/hap2_on_hap1.paf, qc/markers/SAMPLE/marker_classes.tsv.gz
"""
import argparse
import collections
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ap = argparse.ArgumentParser()
ap.add_argument("sample"); ap.add_argument("fai"); ap.add_argument("outdir")
ap.add_argument("--min_mapq", type=int, default=5)
ap.add_argument("--min_len", type=int, default=20000)
ap.add_argument("--paf")
A = ap.parse_args()
S, OUT = A.sample, A.outdir
PAF = A.paf or "results/hap_align/%s/hap2_on_hap1.paf" % S
W, GAP, MINBLOCK = 1000000, 5000000, 1000000
os.makedirs(OUT, exist_ok=True)


def base(c):
    return re.sub(r"_hap[12]$", "", c)


def cnum(c):
    m = re.search(r"(\d+)", base(c))
    return int(m.group(1)) if m else 10 ** 6


L = pd.read_csv(A.fai, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"]).set_index("chrom")["len"]
hap1 = sorted([c for c in L.index if cnum(c) < 10 ** 6], key=cnum)
nb = {c: int(np.ceil(L[c] / float(W))) for c in hap1}

# ---- alignments
alns, qlen = [], {}
with open(PAF) as fh:
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if len(f) < 12 or f[5] not in nb or "tp:A:S" in l:
            continue
        if int(f[11]) < A.min_mapq or int(f[8]) - int(f[7]) < A.min_len:
            continue
        qlen[f[0]] = int(f[1])
        alns.append((f[0], int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8])))
print("alignments kept: %s (MAPQ >= %d, >= %d bp on hap1), from %d hap2 sequences"
      % (format(len(alns), ","), A.min_mapq, A.min_len, len(qlen)), flush=True)

# ---- 1. synteny blocks (hap2 coordinates)
groups = collections.defaultdict(list)
for q, qs, qe, st, t, ts, te in alns:
    groups[(q, t, st)].append((qs, qe, ts, te))
blocks = []
for (q, t, st), lst in groups.items():
    lst.sort()
    cur = None
    for qs, qe, ts, te in lst:
        if cur is not None and qs - cur[1] <= GAP and (
                (st == "+" and -GAP <= ts - cur[3] <= GAP) or (st == "-" and -GAP <= cur[2] - te <= GAP)):
            cur = [cur[0], max(cur[1], qe), min(cur[2], ts), max(cur[3], te), cur[4] + (qe - qs)]
        else:
            if cur is not None:
                blocks.append((q, t, st) + tuple(cur))
            cur = [qs, qe, ts, te, qe - qs]
    if cur is not None:
        blocks.append((q, t, st) + tuple(cur))
bt = pd.DataFrame(blocks, columns=["hap2", "hap1", "strand", "q_start", "q_end", "t_start", "t_end", "aligned"])
bt = bt[bt.aligned >= MINBLOCK].sort_values(["hap2", "q_start"])
main_strand = {}
for q, x in bt[bt.apply(lambda r: base(r.hap1) == base(r.hap2), axis=1)].groupby("hap2"):
    main_strand[q] = x.groupby("strand").aligned.sum().idxmax()
bt["kind"] = [("TRANSLOCATED" if base(r.hap1) != base(r.hap2) else
               "INVERTED" if main_strand.get(r.hap2, r.strand) != r.strand else "")
              for r in bt.itertuples()]
bt.to_csv(os.path.join(OUT, "%s_synteny_blocks.tsv" % S), sep="\t", index=False)

lines = ["TRANSLOCATION MAP  %s  (hap2 on hap1; alignments MAPQ >= %d, >= %d bp; blocks >= %d Mb aligned, "
         "merged across <= %d Mb)" % (S, A.min_mapq, A.min_len, MINBLOCK // 10 ** 6, GAP // 10 ** 6),
         "  1. SYNTENY BLOCKS in hap2 coordinates -> where they land on hap1"]
for q in sorted(bt.hap2.unique(), key=cnum):
    x = bt[bt.hap2 == q]
    lines.append("   %s (%.1f Mb)" % (q, qlen.get(q, 0) / 1e6))
    for r in x.itertuples():
        lines.append("     %7.1f-%7.1f Mb  ->  %-11s %7.1f-%7.1f Mb  (%s)  %6.1f Mb aligned  %s"
                     % (r.q_start / 1e6, r.q_end / 1e6, r.hap1, r.t_start / 1e6, r.t_end / 1e6, r.strand,
                        r.aligned / 1e6, r.kind))

# ---- 2. homolog source per hap1 window
cov = {c: collections.defaultdict(lambda n=nb[c]: np.zeros(n)) for c in hap1}
iv = collections.defaultdict(list)
for q, qs, qe, st, t, ts, te in alns:
    iv[t].append((ts, te))
    src = base(q)
    for w in range(ts // W, min(te // W, nb[t] - 1) + 1):
        cov[t][src][w] += max(0, min(te, (w + 1) * W) - max(ts, w * W))
union = {c: np.zeros(nb[c]) for c in hap1}
for t, lst in iv.items():
    lst.sort()
    merged = []
    for s0, e0 in lst:
        if merged and s0 <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e0)
        else:
            merged.append([s0, e0])
    for s0, e0 in merged:
        for w in range(s0 // W, min(e0 // W, nb[t] - 1) + 1):
            union[t][w] += max(0, min(e0, (w + 1) * W) - max(s0, w * W))

mk = pd.read_csv("qc/markers/%s/marker_classes.tsv.gz" % S, sep="\t", usecols=["chrom", "pos", "class"])
mk = mk[mk["class"] == "good"]
gm = {c: np.bincount((mk.pos[mk.chrom == c] // W).clip(upper=nb[c] - 1), minlength=nb[c]).astype(float)
      for c in hap1}

rows = []
for c in hap1:
    for w in range(nb[c]):
        same = cov[c][base(c)][w] / W if base(c) in cov[c] else 0.0
        others = {s: v[w] / W for s, v in cov[c].items() if s != base(c)}
        osrc, oval = max(others.items(), key=lambda kv: kv[1]) if others else ("", 0.0)
        u = union[c][w] / W
        cls = "foreign" if (oval > same and oval >= 0.2) else ("none" if u < 0.2 else "same")
        rows.append((c, w * W, min((w + 1) * W, L[c]), same, osrc, oval, u, gm[c][w], cls))
wt = pd.DataFrame(rows, columns=["chrom", "start", "end", "same", "other_src", "other", "aligned",
                                 "good_markers", "class"])
wt.to_csv(os.path.join(OUT, "%s_hap1_windows.tsv" % S), sep="\t", index=False)

lines.append("  2. FOREIGN segments in hap1 coordinates (runs >= 2 Mb where another hap2 chromosome dominates)")
seg_rows = []
for c in hap1:
    x = wt[wt.chrom == c].reset_index(drop=True)
    i = 0
    while i < len(x):
        if x["class"][i] == "foreign":
            src, j, gap = x.other_src[i], i, 0
            k = i
            while k + 1 < len(x):
                if x["class"][k + 1] == "foreign" and x.other_src[k + 1] == src:
                    k += 1
                    j, gap = k, 0
                elif gap == 0 and k + 2 < len(x) and x["class"][k + 2] == "foreign" and x.other_src[k + 2] == src:
                    k += 1
                    gap = 1
                else:
                    break
            if j - i + 1 >= 2:
                y = x.iloc[i:j + 1]
                seg_rows.append((c, int(y.start.iloc[0]), int(y.end.iloc[-1]), src, y.other.mean(), y.same.mean(),
                                 y.good_markers.sum() / ((y.end.iloc[-1] - y.start.iloc[0]) / 1e6)))
            i = j + 1
        else:
            i += 1
st = pd.DataFrame(seg_rows, columns=["chrom", "start", "end", "partner", "partner_cov", "own_cov", "good_per_mb"])
st.to_csv(os.path.join(OUT, "%s_foreign_segments.tsv" % S), sep="\t", index=False)
med = {c: np.median(gm[c][gm[c] > 0]) if (gm[c] > 0).any() else 0 for c in hap1}
if len(st):
    for r in st.itertuples():
        lines.append("     %-11s %7.1f-%7.1f Mb  homolog on %s_hap2 (covers %.0f%%, own %.0f%%)  "
                     "good markers %.0f/Mb (chromosome median %.0f/Mb)"
                     % (r.chrom, r.start / 1e6, r.end / 1e6, r.partner, 100 * r.partner_cov, 100 * r.own_cov,
                        r.good_per_mb, med[r.chrom]))
else:
    lines.append("     none")

lines.append("  3. MARKERS: informative markers per Mb by window class (median over windows)")
lines.append("     %-11s %10s %10s %10s   %s" % ("chrom", "same", "foreign", "no homolog", "windows same/foreign/none"))
for c in hap1:
    x = wt[wt.chrom == c]
    f = lambda k: x[x["class"] == k].good_markers.median() if (x["class"] == k).any() else float("nan")
    lines.append("     %-11s %10.0f %10.0f %10.0f   %d / %d / %d" % (c, f("same"), f("foreign"), f("none"),
                 (x["class"] == "same").sum(), (x["class"] == "foreign").sum(), (x["class"] == "none").sum()))
lines.append("  DESERTS (>= 5 Mb below 10% of the chromosome's median marker density) and what they are made of")
nd = 0
for c in hap1:
    x = wt[wt.chrom == c].reset_index(drop=True)
    low = (x.good_markers < 0.1 * med[c]).values
    i = 0
    while i < len(x):
        if low[i]:
            j = i
            while j + 1 < len(x) and low[j + 1]:
                j += 1
            if j - i + 1 >= 5:
                y = x.iloc[i:j + 1]
                comp = y["class"].value_counts()
                lines.append("     %-11s %7.1f-%7.1f Mb  %3d windows: same %d, foreign %d, no homolog %d; "
                             "mean hap2 alignable %.0f%%"
                             % (c, y.start.iloc[0] / 1e6, y.end.iloc[-1] / 1e6, len(y), comp.get("same", 0),
                                comp.get("foreign", 0), comp.get("none", 0), 100 * y.aligned.mean()))
                nd += 1
            i = j + 1
        else:
            i += 1
if nd == 0:
    lines.append("     none")
open(os.path.join(OUT, "%s_translocation_map.txt" % S), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

# ---- figures
off, acc = {}, 0
for c in hap1:
    off[c] = acc
    acc += L[c]
qs_ = sorted(qlen, key=cnum)
qs_ = [q for q in qs_ if cnum(q) < 10 ** 6]
ncol = 3
nrow = int(np.ceil(len(qs_) / float(ncol)))
fig, ax = plt.subplots(nrow, ncol, figsize=(5.2 * ncol, 5.0 * nrow), squeeze=False)
for i, q in enumerate(qs_):
    a = ax[i // ncol][i % ncol]
    seg, col = [], []
    for qn, qs, qe, sd, t, ts, te in alns:
        if qn != q:
            continue
        y0, y1 = (off[t] + ts, off[t] + te) if sd == "+" else (off[t] + te, off[t] + ts)
        seg.append([(qs / 1e6, y0 / 1e6), (qe / 1e6, y1 / 1e6)])
        col.append("#185FA5" if sd == "+" else "#D85A30")
    a.add_collection(LineCollection(seg, colors=col, linewidths=1.2))
    for c in hap1:
        a.axhline(off[c] / 1e6, color="#B4B2A9", lw=.6)
        a.text(qlen[q] / 1e6 * 1.01, (off[c] + L[c] / 2) / 1e6, base(c), fontsize=7, va="center")
    a.set_xlim(0, qlen[q] / 1e6); a.set_ylim(acc / 1e6, 0)
    a.set_xlabel("%s (Mb)" % q, fontsize=8); a.set_ylabel("hap1 genome (Mb)", fontsize=8); a.tick_params(labelsize=7)
    a.set_title(q, fontsize=9)
for i in range(len(qs_), nrow * ncol):
    ax[i // ncol][i % ncol].set_axis_off()
fig.suptitle("%s: hap2 chromosomes (x) against the hap1 genome (y); blue + strand, orange - strand" % S)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "%s_dotplot.png" % S), dpi=140)

cmap = {base(c): plt.cm.tab10(k % 10) for k, c in enumerate(hap1)}
fig, ax = plt.subplots(len(hap1), 1, figsize=(13, 1.9 * len(hap1)), squeeze=False)
for i, c in enumerate(hap1):
    a = ax[i][0]
    x = wt[wt.chrom == c]
    xm = (x.start + x.end) / 2e6
    a.fill_between(xm, 0, x.same, color="#B4B2A9", lw=0, label="homolog on the same hap2 chromosome")
    for src in sorted(set(x.other_src) - {""}, key=lambda s: cnum(s)):
        v = np.where(x.other_src == src, x.other, 0)
        if v.max() >= 0.2:
            a.fill_between(xm, 0, v, color=cmap.get(src, "#534AB7"), alpha=.85, lw=0,
                           label="homolog on %s_hap2" % src)
    b = a.twinx()
    b.plot(xm, x.good_markers, color="black", lw=.9, label="informative markers / Mb")
    b.set_ylabel("markers/Mb", fontsize=7); b.tick_params(labelsize=6)
    a.set_ylim(0, 1.05); a.set_ylabel("share", fontsize=7); a.tick_params(labelsize=7)
    a.set_title(c, fontsize=9, loc="left")
    h1, l1 = a.get_legend_handles_labels(); h2, l2 = b.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, fontsize=6, frameon=False, loc="upper right", ncol=4)
ax[-1][0].set_xlabel("hap1 position (Mb)")
fig.suptitle("%s: where each hap1 window's homolog sits in hap2 (fill) and informative markers (line)" % S)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "%s_hap1_homolog_source.png" % S), dpi=130)
print("wrote %s/%s_{translocation_map.txt,synteny_blocks.tsv,foreign_segments.tsv,hap1_windows.tsv,"
      "dotplot.png,hap1_homolog_source.png}" % (OUT, S))
