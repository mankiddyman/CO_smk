#!/usr/bin/env python3
"""chromosome_context.py -- why is a chromosome end hot, cold, or one-sided?

Per chromosome, in windows, against the crossover rate from the pipeline's calls:
  genes per Mb          from the annotation
  repeat fraction       soft-masked (lowercase) share of the assembly, if it is
                        soft-masked; skipped otherwise
  informative markers   'good' markers per Mb -- where crossovers can be SEEN
  hap2 alignable        share of the window covered by the hap2-on-hap1
                        alignment: where the two haplotypes are homologous
                        enough to pair; a non-alignable end cannot recombine
  rDNA                  the rRNA blacklist intervals (NOR ends are CO-poor)
And a second figure, recurrent islands: per window, the share of cells whose
informative molecules there carry the OTHER allele from both flanks while the
flanks agree. A real close double crossover is rare and scattered; the same
island in many cells at one position is structural (e.g. a phase switch in the
assembly), not meiosis.

Islands are scanned in 0.5 Mb windows against the nearest informative windows
1-3 windows away on each side, so an island up to ~1.5 Mb (shorter than
hapCO's 2 Mb block minimum, i.e. one hapCO cannot call) is seen even next to a
chromosome end.

Usage: chromosome_context.py SAMPLE FAI GFF3 FASTA PAF RRNA_BED OUTDIR
         [--window_mb 1] [--region CHROM:START-END] [--cell BARCODE]
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
for k in ("sample", "fai", "gff", "fasta", "paf", "bed", "outdir"):
    ap.add_argument(k)
ap.add_argument("--window_mb", default="1"); ap.add_argument("--region"); ap.add_argument("--cell")
A = ap.parse_args()
S, FAI, GFF, FA, PAF, BED, OUT, WMB = A.sample, A.fai, A.gff, A.fasta, A.paf, A.bed, A.outdir, A.window_mb
W = int(float(WMB) * 1e6)
GAP, MIN_MOL, IW = 150, 3, 500000
os.makedirs(OUT, exist_ok=True)
L = pd.read_csv(FAI, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"]).set_index("chrom")["len"]
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
nb = {c: int(np.ceil(L[c] / float(W))) for c in L.index}


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


# ---- crossover rate from the pipeline's calls
co = {c: np.zeros(nb[c]) for c in L.index}
for bc in cells:
    p = "results/crossovers/%s/per_cell/%s_co_pred.txt" % (S, bc)
    if os.path.exists(p):
        for l in open(p):
            f = l.split()
            if len(f) >= 3 and f[0] in co and f[1].isdigit():
                m = (int(f[1]) + int(f[2])) / 2.0
                co[f[0]][min(int(m // W), nb[f[0]] - 1)] += 1
chroms = [c for c in L.index if co[c].sum() > 0]
rate = {c: co[c] / len(cells) * 100.0 / float(WMB) for c in chroms}

# ---- genes
genes = {c: np.zeros(nb[c]) for c in chroms}
with opener(GFF) as fh:
    for l in fh:
        if l.startswith("#"):
            continue
        f = l.split("\t", 5)
        if len(f) > 4 and f[2] == "gene" and f[0] in genes:
            genes[f[0]][min(int(int(f[3]) // W), nb[f[0]] - 1)] += 1

# ---- repeats from soft-masking
low = {c: np.zeros(nb[c]) for c in chroms}
tot = {c: np.zeros(nb[c]) for c in chroms}
cur, pos = None, 0
with (gzip.open(FA, "rb") if FA.endswith(".gz") else open(FA, "rb")) as fh:
    for line in fh:
        if line.startswith(b">"):
            cur, pos = line[1:].split()[0].decode(), 0
            continue
        if cur not in low:
            continue
        s = line.rstrip()
        n = len(s)
        nl = n - len(s.translate(None, b"acgtn"))
        w = min(pos // W, nb[cur] - 1)
        low[cur][w] += nl
        tot[cur][w] += n
        pos += n
masked = sum(v.sum() for v in low.values()) / max(sum(v.sum() for v in tot.values()), 1)
soft = masked >= 0.01

# ---- informative markers
mk = pd.read_csv("qc/markers/%s/marker_classes.tsv.gz" % S, sep="\t", usecols=["chrom", "pos", "class"])
mk = mk[mk["class"] == "good"]
markers = {c: np.bincount((mk.pos[mk.chrom == c] // W).clip(upper=nb[c] - 1), minlength=nb[c]).astype(float)
           for c in chroms}

# ---- hap2 alignable (union of primary alignments >= 10 kb, MAPQ >= 5)
iv = {c: [] for c in chroms}
with opener(PAF) as fh:
    for l in fh:
        f = l.split("\t")
        if len(f) < 12 or f[5] not in iv:
            continue
        if int(f[11]) < 5 or int(f[8]) - int(f[7]) < 10000 or "tp:A:S" in l:
            continue
        iv[f[5]].append((int(f[7]), int(f[8])))
alig = {c: np.zeros(nb[c]) for c in chroms}
for c, lst in iv.items():
    lst.sort()
    merged = []
    for s0, e0 in lst:
        if merged and s0 <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e0)
        else:
            merged.append([s0, e0])
    for s0, e0 in merged:
        for w in range(s0 // W, min(e0 // W, nb[c] - 1) + 1):
            alig[c][w] += max(0, min(e0, (w + 1) * W) - max(s0, w * W))
    alig[c] /= float(W)

# ---- rDNA
rdna = {c: [] for c in chroms}
if os.path.exists(BED):
    for l in open(BED):
        f = l.split()
        if len(f) >= 3 and f[0] in rdna:
            rdna[f[0]].append((int(f[1]), int(f[2])))

# ---- recurrent islands, from each cell's informative molecules (0.5 Mb windows)
good = set(zip(mk.chrom, mk.pos.astype(int)))
ni = {c: int(np.ceil(L[c] / float(IW))) for c in chroms}
isl = {c: np.zeros(ni[c]) for c in chroms}
elig = {c: np.zeros(ni[c]) for c in chroms}
cell_isl = {}


def ref_window(n, k, step):
    for j in (k + 2 * step, k + 3 * step, k + 4 * step):
        if 0 <= j < len(n) and n[j] >= MIN_MOL:
            return j
    return None


for bc in cells:
    d = pd.read_csv("results/cell_data/%s/%s.tsv" % (S, bc), sep="\t", header=None,
                    names=["chrom", "pos", "ref", "rc", "alt", "ac"])
    d = d[[(c, int(p)) in good for c, p in zip(d.chrom, d.pos)]]
    dp = d.rc + d.ac
    fr = d.ac / dp.clip(lower=1)
    d = d.assign(call=np.where(fr <= 0.2, 0, np.where(fr >= 0.8, 1, -1)))
    d = d[(dp > 0) & (d.call >= 0)].sort_values(["chrom", "pos"])
    for c, x in d.groupby("chrom", sort=False):
        if c not in isl:
            continue
        p, g = x.pos.values, x.call.values.astype(float)
        new = np.ones(len(p), dtype=bool)
        new[1:] = (p[1:] - p[:-1]) >= GAP
        mid = np.cumsum(new) - 1
        mf = np.bincount(mid, weights=g) / np.bincount(mid)
        keep = mf != 0.5
        mp, ma = p[new][keep], (mf[keep] > 0.5).astype(float)
        w = (mp // IW).clip(max=ni[c] - 1)
        n = np.bincount(w, minlength=ni[c]).astype(float)
        a = np.bincount(w, weights=ma, minlength=ni[c])
        for k in range(ni[c]):
            if n[k] < MIN_MOL:
                continue
            jl, jr = ref_window(n, k, -1), ref_window(n, k, 1)
            if jl is None or jr is None:
                continue
            fl, frr, fk = a[jl] / n[jl], a[jr] / n[jr], a[k] / n[k]
            if (fl >= 0.8 and frr >= 0.8) or (fl <= 0.2 and frr <= 0.2):
                elig[c][k] += 1
                side = 1 if fl >= 0.8 else 0
                if (side == 1 and fk <= 0.2) or (side == 0 and fk >= 0.8):
                    isl[c][k] += 1
                    if bc == A.cell:
                        cell_isl.setdefault(c, []).append(k)
freq = {c: np.where(elig[c] >= 20, isl[c] / np.maximum(elig[c], 1), np.nan) for c in chroms}

# ---- summary
lines = ["CHROMOSOME CONTEXT  %s  (%d cells, %s Mb windows; assembly %.0f%% soft-masked%s)"
         % (S, len(cells), WMB, 100 * masked, "" if soft else " -- repeats track skipped")]
lines.append("  %-12s %7s %22s %22s %24s" % ("chrom", "COs", "left 5 Mb: CO / align", "right 5 Mb: CO / align",
                                             "genes left / mid / right"))
for c in chroms:
    k = max(1, int(round(5e6 / W)))
    t = np.array_split(np.arange(nb[c]), 3)
    lines.append("  %-12s %7.2f %11.1f / %4.0f%% %13.1f / %4.0f%% %10.1f / %4.1f / %4.1f"
                 % (c, co[c].sum() / len(cells), rate[c][:k].mean(), 100 * alig[c][:k].mean(),
                    rate[c][-k:].mean(), 100 * alig[c][-k:].mean(),
                    genes[c][t[0]].mean() / float(WMB), genes[c][t[1]].mean() / float(WMB),
                    genes[c][t[2]].mean() / float(WMB)))
rec = sorted(((freq[c][k], c, k) for c in chroms for k in range(ni[c]) if freq[c][k] == freq[c][k]), reverse=True)
lines.append("  most frequent islands (share of eligible cells whose window carries the other allele from both flanks):")
for f0, c, k in rec[:12]:
    lines.append("    %-12s %6.1f-%6.1f Mb  %5.1f%% of %d cells" % (c, k * IW / 1e6, (k + 1) * IW / 1e6, 100 * f0,
                                                                 int(elig[c][k])))
if A.region:
    rc_, rng = A.region.split(":")
    r0, r1 = [float(v) for v in rng.split("-")]
    lines.append("  region %s:" % A.region)
    for k in range(int(r0 // IW), min(int(r1 // IW) + 1, ni.get(rc_, 0))):
        lines.append("    %.1f-%.1f Mb  island in %d of %d eligible cells (%.1f%%)"
                     % (k * IW / 1e6, (k + 1) * IW / 1e6, isl[rc_][k], elig[rc_][k],
                        100 * isl[rc_][k] / max(elig[rc_][k], 1)))
if A.cell:
    lines.append("  cell %s islands: %s" % (A.cell, "; ".join("%s %s" % (c, ", ".join(
        "%.1f-%.1f Mb" % (k * IW / 1e6, (k + 1) * IW / 1e6) for k in ks)) for c, ks in sorted(cell_isl.items()))
        or "none found"))
open(os.path.join(OUT, "%s_chromosome_context.txt" % S), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

# ---- figures
ncol = 4 if len(chroms) > 6 else 3
nrow = int(np.ceil(len(chroms) / float(ncol)))
fig, ax = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.0 * nrow), squeeze=False)
for i, c in enumerate(chroms):
    a = ax[i // ncol][i % ncol]
    x = (np.arange(nb[c]) + 0.5) * W / 1e6
    b = a.twinx()
    if soft:
        b.fill_between(x, low[c] / np.maximum(tot[c], 1), color="#B4B2A9", alpha=.45, lw=0, label="repeats (fraction)")
    b.plot(x, alig[c], color="#534AB7", lw=1.1, label="hap2 alignable (fraction)")
    if genes[c].max() > 0:
        b.plot(x, genes[c] / genes[c].max(), color="#2E7D32", lw=1.1, label="genes (scaled)")
    if markers[c].max() > 0:
        b.plot(x, markers[c] / markers[c].max(), color="#1D9E75", ls="--", lw=.9, label="informative markers (scaled)")
    for s0, e0 in rdna[c]:
        b.axvspan(s0 / 1e6, max(e0, s0 + W / 5) / 1e6, ymin=.93, ymax=1, color="#D85A30")
    b.set_ylim(0, 1.05); b.tick_params(labelsize=6)
    a.plot(x, rate[c], color="black", lw=1.6, label="CO rate (cM/Mb)")
    a.set_zorder(b.get_zorder() + 1); a.patch.set_visible(False)
    a.set_title("%s   %.2f COs/cell" % (c, co[c].sum() / len(cells)), fontsize=9)
    a.tick_params(labelsize=7)
    if i == 0:
        h1, l1 = a.get_legend_handles_labels(); h2, l2 = b.get_legend_handles_labels()
        a.legend(h1 + h2, l1 + l2, frameon=False, fontsize=6, loc="upper center")
for i in range(len(chroms), nrow * ncol):
    ax[i // ncol][i % ncol].set_axis_off()
fig.suptitle("%s: crossover rate (black, left axis) against genome context (right axis); rDNA in orange at top" % S)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "%s_chromosome_context.png" % S), dpi=130)
fig.savefig(os.path.join(OUT, "%s_chromosome_context.pdf" % S))

fig, ax = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 2.6 * nrow), squeeze=False)
for i, c in enumerate(chroms):
    a = ax[i // ncol][i % ncol]
    x = (np.arange(ni[c]) + 0.5) * IW / 1e6
    a.bar(x, np.nan_to_num(freq[c]), width=IW / 1e6 * .9, color="#D85A30")
    a.set_ylim(0, max(0.2, np.nanmax(freq[c]) * 1.1 if np.isfinite(np.nanmax(freq[c])) else 0.2))
    a.set_title(c, fontsize=9); a.tick_params(labelsize=7)
    if i % ncol == 0:
        a.set_ylabel("share of cells", fontsize=8)
for i in range(len(chroms), nrow * ncol):
    ax[i // ncol][i % ncol].set_axis_off()
fig.suptitle("%s: recurrent islands -- scattered and low = real close doubles; one tall bar = structural" % S)
fig.tight_layout()
a0 = ax[0][0]
a0.set_xlabel("Mb (0.5 Mb windows)", fontsize=7)
fig.savefig(os.path.join(OUT, "%s_recurrent_islands.png" % S), dpi=130)
print("wrote %s/%s_{chromosome_context,recurrent_islands}.png" % (OUT, S))
