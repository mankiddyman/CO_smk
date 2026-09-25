#!/usr/bin/env python3
"""haplotype_tracks.py -- can a sliding window recover ONE haplotype per cell?

Two steps, both on informative ("good") markers only:

1. Collapse to molecules. Markers closer than GAP bp usually sit on one read,
   so they are one observation, not several: each run of close markers in a
   cell becomes one molecule carrying the majority allele.

2. Slide along each chromosome in windows of W molecules and take the ALT
   fraction per window.
     one haploid nucleus + some ambient: windows near 0 or near 1, flipping
       only at crossovers. Ambient pulls them inward a little, but a majority
       haplotype is always there to recover.
     a barcode holding both haplotypes (diploid, tetrad, pooled nuclei):
       every window near 0.5. There is no majority haplotype, so no window
       size can recover one.
     two haploid nuclei together: windows at 0, 0.5 and 1.

   haploidness = mean |window fraction - 0.5| x 2   (1 = clean haploid, ~0.2 = flat 50:50)

Usage: haplotype_tracks.py SAMPLE OUTDIR
Needs results/snps/SAMPLE/* and qc/markers/SAMPLE/marker_classes.tsv.gz.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLE, OUT = sys.argv[1], sys.argv[2]
SNP = os.path.join("results/snps", SAMPLE)
W, GAP, MIN_WIN = 15, 150, 5
os.makedirs(OUT, exist_ok=True)
lines = []


def say(s=""):
    lines.append(s)
    print(s, flush=True)


def read_mtx(path):
    df = pd.read_csv(path, sep=r"\s+", comment="%", header=None, dtype=np.int64)
    d = df.iloc[0].values
    b = df.iloc[1:].values
    return int(d[0]), int(d[1]), b[:, 0] - 1, b[:, 1] - 1, b[:, 2]


sites = pd.read_csv(os.path.join(SNP, "cellSNP.base.vcf.gz"), sep="\t", comment="#", header=None,
                    usecols=[0, 1], names=["chrom", "pos"], dtype={0: str})
cls = pd.read_csv(os.path.join("qc/markers", SAMPLE, "marker_classes.tsv.gz"), sep="\t", usecols=["class"])
good = (cls["class"].values == "good")
assert len(good) == len(sites), "marker_classes.tsv.gz does not match this sample's sites"
nsite, ncell, s_dp, c_dp, dp = read_mtx(os.path.join(SNP, "cellSNP.tag.DP.mtx"))
_, _, s_ad, c_ad, ad = read_mtx(os.path.join(SNP, "cellSNP.tag.AD.mtx"))
key = s_dp * ncell + c_dp
o = np.argsort(key, kind="stable")
k_ad = s_ad * ncell + c_ad
at = np.searchsorted(key[o], k_ad)
alt = np.zeros(len(dp), dtype=np.int64)
ok = at < len(o)
ok[ok] &= key[o][at[ok]] == k_ad[ok]
alt[o[at[ok]]] = ad[ok]
del key, o, at, k_ad
fr = alt / np.maximum(dp, 1)
call = np.full(len(dp), -1, dtype=np.int8)
call[fr <= 0.2] = 0
call[fr >= 0.8] = 1

chrom_code, chrom_names = pd.factorize(sites["chrom"])
chrom_code = chrom_code.astype(np.int64)
spos = sites["pos"].values.astype(np.int64)
i = np.nonzero(good[s_dp] & (call >= 0))[0]
c, ch, p, g = c_dp[i], chrom_code[s_dp[i]], spos[s_dp[i]], call[i].astype(np.float64)
o = np.lexsort((p, ch, c))
c, ch, p, g = c[o], ch[o], p[o], g[o]

# ---- 1. molecules
new = np.ones(len(c), dtype=bool)
new[1:] = (c[1:] != c[:-1]) | (ch[1:] != ch[:-1]) | ((p[1:] - p[:-1]) >= GAP)
mid = np.cumsum(new) - 1
nm = int(mid[-1]) + 1
mfrac = np.bincount(mid, weights=g, minlength=nm) / np.bincount(mid, minlength=nm)
first = np.nonzero(new)[0]
keep = mfrac != 0.5
mc, mch, mp = c[first][keep], ch[first][keep], p[first][keep]
mg = (mfrac[keep] > 0.5).astype(np.int8)
nm = len(mc)
say("HAPLOTYPE TRACKS  %s   (window = %d molecules; markers < %d bp apart = one molecule)" % (SAMPLE, W, GAP))
say("  %s informative observations -> %s molecules" % (format(len(i), ","), format(nm, ",")))

run_new = np.ones(nm, dtype=bool)
run_new[1:] = (mc[1:] != mc[:-1]) | (mch[1:] != mch[:-1])
same_run = ~run_new[1:]
mol_sw = np.bincount(mc[1:][same_run & (mg[1:] != mg[:-1])], minlength=ncell)
mol_pairs = np.bincount(mc[1:][same_run], minlength=ncell)
n_mol = np.bincount(mc, minlength=ncell)

# ---- 2. windows
run_id = np.cumsum(run_new) - 1
run_start = np.maximum.accumulate(np.where(run_new, np.arange(nm), 0))
rank = np.arange(nm) - run_start
wkey = run_id.astype(np.int64) * (int(rank.max()) // W + 2) + rank // W
_, w_first, w_inv = np.unique(wkey, return_index=True, return_inverse=True)
w_n = np.bincount(w_inv)
w_f = np.bincount(w_inv, weights=mg) / w_n
full = w_n == W
w_cell = mc[w_first][full]
w_f = w_f[full]
n_win = np.bincount(w_cell, minlength=ncell)
hap = np.bincount(w_cell, weights=np.abs(w_f - 0.5) * 2, minlength=ncell) / np.maximum(n_win, 1)
dec = np.bincount(w_cell, weights=((w_f <= 0.2) | (w_f >= 0.8)).astype(float), minlength=ncell) / np.maximum(n_win, 1)
mid_share = np.bincount(w_cell, weights=((w_f > 0.35) & (w_f < 0.65)).astype(float),
                        minlength=ncell) / np.maximum(n_win, 1)

bcf = os.path.join(SNP, "cellSNP.samples.tsv")
bc = [l.strip() for l in open(bcf)] if os.path.exists(bcf) else ["cell%d" % k for k in range(ncell)]
cell = pd.DataFrame({"barcode": bc[:ncell], "molecules": n_mol, "windows": n_win,
                     "molecule_switch_rate": np.where(mol_pairs >= 20, mol_sw / np.maximum(mol_pairs, 1), np.nan),
                     "haploidness": np.where(n_win >= MIN_WIN, hap, np.nan),
                     "decisive_windows": np.where(n_win >= MIN_WIN, dec, np.nan),
                     "middle_windows": np.where(n_win >= MIN_WIN, mid_share, np.nan)})
cell.to_csv(os.path.join(OUT, "haplotype_tracks.tsv.gz"), sep="\t", index=False)

v = cell.dropna(subset=["haploidness"])
h, e = np.histogram(v.haploidness, bins=np.arange(0, 1.02, 0.02))
sm = np.convolve(h, np.ones(3) / 3, mode="same")
mids = (e[:-1] + e[1:]) / 2
cand = np.where((mids >= 0.3) & (mids <= 0.8))[0]
valley = mids[cand[np.argmin(sm[cand])]] if len(cand) else np.nan
say("  cells with >= %d full windows: %s of %s" % (MIN_WIN, format(len(v), ","), format(ncell, ",")))
say("  median haploidness %.2f; valley between the modes at %.2f" % (v.haploidness.median(), valley))
say()
say("  %-34s %10s %12s %12s" % ("cells with haploidness >=", "cells", "median mol.", "median win."))
for t in (0.4, 0.5, 0.6, 0.7, 0.8):
    s = v[v.haploidness >= t]
    say("  %-34s %10s %12s %12s" % ("  %.1f" % t, format(len(s), ","),
                                    format(int(s.molecules.median()) if len(s) else 0, ","),
                                    format(int(s.windows.median()) if len(s) else 0, ",")))
for mw in (10, 20):
    s = v[(v.haploidness >= valley) & (v.windows >= mw)]
    say("  at the valley, with >= %-2d windows      %10s" % (mw, format(len(s), ",")))

# ---- plots: the classifier
fig, ax = plt.subplots(2, 2, figsize=(12, 9))
a = ax[0, 0]
a.hist(v.haploidness, bins=e, color="#534AB7", alpha=.85)
a.axvline(valley, color="black", ls="--", lw=.8, label="valley %.2f" % valley)
a.set_xlabel("haploidness (1 = one clean haplotype, ~0.2 = flat 50:50)"); a.set_ylabel("cells")
a.set_title("%s: one haplotype, or both?" % SAMPLE); a.legend(frameon=False, fontsize=8)

a = ax[0, 1]
a.scatter(v.molecules, v.haploidness, s=3, alpha=.3, color="#534AB7", rasterized=True)
a.axhline(valley, color="black", ls="--", lw=.8)
a.set_xscale("log"); a.set_xlabel("molecules in the cell"); a.set_ylabel("haploidness")
a.set_title("does it depend on depth?")

a = ax[1, 0]
a.scatter(v.molecule_switch_rate, v.haploidness, s=3, alpha=.3, color="#534AB7", rasterized=True)
a.axhline(valley, color="black", ls="--", lw=.8)
a.set_xlabel("switch rate between consecutive molecules"); a.set_ylabel("haploidness")
a.set_title("haploidness vs the far-pair switch rate")

a = ax[1, 1]
hi = set(v.index[v.haploidness >= max(valley, 0.6)])
lo = set(v.index[v.haploidness < min(valley, 0.35)])
wh = w_f[np.isin(w_cell, list(hi))]
wl = w_f[np.isin(w_cell, list(lo))]
bb = np.linspace(0, 1, W + 2)
if len(wh):
    a.hist(wh, bins=bb, histtype="step", lw=2, density=True, color="#2E7D32",
           label="haploid-like cells (%s)" % format(len(hi), ","))
if len(wl):
    a.hist(wl, bins=bb, histtype="step", lw=2, density=True, color="#D85A30",
           label="both-haplotype cells (%s)" % format(len(lo), ","))
a.set_xlabel("ALT fraction in a %d-molecule window" % W); a.set_ylabel("density")
a.set_title("what the windows look like"); a.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "haploidness.png"), dpi=120); fig.savefig(os.path.join(OUT, "haploidness.pdf"))

# ---- plots: example cells along the genome
clen = np.array([spos[chrom_code == k].max() for k in range(len(chrom_names))], dtype=np.int64)
offs = np.concatenate([[0], np.cumsum(clen)[:-1]])
pick = []
hv = v[v.haploidness >= max(valley, 0.6)]
pick += [("haploid, deep", x) for x in hv.sort_values("molecules", ascending=False).index[:3]]
mv = hv[hv.molecules <= hv.molecules.median()].sort_values("molecules", ascending=False)
pick += [("haploid, shallow", x) for x in mv.index[:2]]
lv = v[v.haploidness < min(valley, 0.35)]
pick += [("both haplotypes", x) for x in lv.sort_values("molecules", ascending=False).index[:3]]
if pick:
    fig, ax = plt.subplots(len(pick), 1, figsize=(14, 1.9 * len(pick)), squeeze=False)
    rng = np.random.default_rng(1)
    for r, (lab, k) in enumerate(pick):
        a = ax[r][0]
        m = mc == k
        x = offs[mch[m]] + mp[m]
        y = mg[m].astype(float)
        a.scatter(x / 1e6, y + rng.uniform(-.08, .08, len(y)), s=2, alpha=.35, color="#888780", rasterized=True)
        for kc in np.unique(mch[m]):
            sel = mch[m] == kc
            if sel.sum() >= W:
                rm = np.convolve(y[sel], np.ones(W) / W, mode="valid")
                xs = x[sel][W // 2: W // 2 + len(rm)]
                a.plot(xs / 1e6, rm, color="#2E7D32" if lab.startswith("haploid") else "#D85A30", lw=1.5)
        for b0 in offs[1:]:
            a.axvline(b0 / 1e6, color="black", lw=.4, alpha=.5)
        a.axhline(0.5, color="grey", ls=":", lw=.8)
        a.set_ylim(-.15, 1.15); a.set_yticks([0, .5, 1]); a.tick_params(labelsize=7)
        a.set_ylabel("ALT", fontsize=8)
        a.set_title("%s  |  %s  |  %s molecules, haploidness %.2f" % (lab, v.loc[k, "barcode"],
                    format(int(v.loc[k, "molecules"]), ","), v.loc[k, "haploidness"]), fontsize=9, loc="left")
    ax[-1][0].set_xlabel("genome position (Mb; vertical lines = chromosome boundaries)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "example_tracks.png"), dpi=110)
with open(os.path.join(OUT, "haplotype_tracks.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
say("\nwrote %s/{haploidness.png,.pdf, example_tracks.png, haplotype_tracks.tsv.gz,.txt}" % OUT)
