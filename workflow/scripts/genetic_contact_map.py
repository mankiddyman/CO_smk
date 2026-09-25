#!/usr/bin/env python3
"""genetic_contact_map.py -- which parts of the genome travel together?

The genetic twin of a Hi-C map, built from haploid cells. Each cell is one
meiotic product, so for every pair of genome bins we can ask: across cells,
does knowing the haplotype in bin i predict the haplotype in bin j?

  same chromosome, close together   r near 1 (few crossovers between them)
  same chromosome, far apart        r decays toward 0 (crossovers accumulate)
  different chromosomes             r near 0 (independent assortment)
  ...EXCEPT a reciprocal translocation: in balanced (alternate) segregation
  hap2's two rearranged chromosomes travel together, so the regions flanking
  the breakpoints on BOTH hap1 chromosomes carry the same parental state.
  That shows as an off-diagonal peak whose position marks the breakpoint on
  each chromosome.

Per cell and bin: molecules as in haplotype_tracks.py (good markers, calls
closer than 150 bp = one molecule); bin state REF/ALT when >= MIN_MOL
molecules agree at >= 80%, otherwise missing. r = phi correlation over the
cells where both bins are called.

Usage: genetic_contact_map.py SAMPLE CELL_LIST OUTDIR BIN_MB
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLE, CELLS, OUT, BIN_MB = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
SNP = os.path.join("results/snps", SAMPLE)
GAP, MIN_MOL, MIN_SHARED = 150, 3, 20
B = int(BIN_MB * 1e6)
os.makedirs(OUT, exist_ok=True)


def read_mtx(path):
    df = pd.read_csv(path, sep=r"\s+", comment="%", header=None, dtype=np.int64)
    d = df.iloc[0].values
    b = df.iloc[1:].values
    return int(d[0]), int(d[1]), b[:, 0] - 1, b[:, 1] - 1, b[:, 2]


sites = pd.read_csv(os.path.join(SNP, "cellSNP.base.vcf.gz"), sep="\t", comment="#", header=None,
                    usecols=[0, 1], names=["chrom", "pos"], dtype={0: str})
good = pd.read_csv(os.path.join("qc/markers", SAMPLE, "marker_classes.tsv.gz"), sep="\t",
                   usecols=["class"])["class"].values == "good"
bcf = os.path.join(SNP, "cellSNP.samples.tsv")
barcodes = [l.strip() for l in open(bcf)]
want = set(l.strip() for l in open(CELLS) if l.strip())
keep_cell = np.array([b in want for b in barcodes])
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
fr = alt / np.maximum(dp, 1)
call = np.full(len(dp), -1, dtype=np.int8)
call[fr <= 0.2] = 0
call[fr >= 0.8] = 1

chrom_code, chrom_names = pd.factorize(sites["chrom"])
chrom_code = chrom_code.astype(np.int64)
spos = sites["pos"].values.astype(np.int64)
i = np.nonzero(good[s_dp] & (call >= 0) & keep_cell[c_dp])[0]
c, ch, p, g = c_dp[i], chrom_code[s_dp[i]], spos[s_dp[i]], call[i].astype(np.float64)
o = np.lexsort((p, ch, c))
c, ch, p, g = c[o], ch[o], p[o], g[o]
new = np.ones(len(c), dtype=bool)
new[1:] = (c[1:] != c[:-1]) | (ch[1:] != ch[:-1]) | ((p[1:] - p[:-1]) >= GAP)
mid = np.cumsum(new) - 1
mf = np.bincount(mid, weights=g) / np.bincount(mid)
first = np.nonzero(new)[0]
k = mf != 0.5
mc, mch, mp, mg = c[first][k], ch[first][k], p[first][k], (mf[k] > 0.5).astype(np.int64)

# genome bins
clen = np.array([spos[chrom_code == j].max() + 1 for j in range(len(chrom_names))])
nb = (clen + B - 1) // B
boff = np.concatenate([[0], np.cumsum(nb)[:-1]])
nbin = int(nb.sum())
gbin = boff[mch] + mp // B
cells = np.unique(mc)
cidx = np.searchsorted(cells, mc)
n = np.zeros((len(cells), nbin)); a = np.zeros((len(cells), nbin))
np.add.at(n, (cidx, gbin), 1)
np.add.at(a, (cidx, gbin), mg)
f = np.where(n > 0, a / np.maximum(n, 1), np.nan)
S = np.full(f.shape, np.nan)
S[(n >= MIN_MOL) & (f >= 0.8)] = 1
S[(n >= MIN_MOL) & (f <= 0.2)] = 0
called = ~np.isnan(S)
print("GENETIC CONTACT MAP  %s: %d cells, %d bins of %.0f Mb; %.0f%% of cell-bins called"
      % (SAMPLE, len(cells), nbin, BIN_MB, 100 * called.mean()), flush=True)

R = np.full((nbin, nbin), np.nan)
NS = np.zeros((nbin, nbin), dtype=int)
for x in range(nbin):
    for y in range(x, nbin):
        both = called[:, x] & called[:, y]
        m = int(both.sum())
        NS[x, y] = NS[y, x] = m
        if m >= MIN_SHARED:
            u, v = S[both, x], S[both, y]
            if u.std() > 0 and v.std() > 0:
                R[x, y] = R[y, x] = np.corrcoef(u, v)[0, 1]

bin_chrom = np.repeat(np.arange(len(chrom_names)), nb)
bin_start = np.concatenate([np.arange(k) * B for k in nb])
lines = ["GENETIC CONTACT MAP  %s  (%d haploid cells, %.0f Mb bins)" % (SAMPLE, len(cells), BIN_MB),
         "  inter-chromosomal linkage (independent assortment gives |r| ~ 0):",
         "  %-14s %-14s %8s %14s %14s %8s" % ("chrom A", "chrom B", "max |r|", "at A (Mb)", "at B (Mb)", "cells")]
flag = []
for x in range(len(chrom_names)):
    for y in range(x + 1, len(chrom_names)):
        blk = np.abs(R[np.ix_(bin_chrom == x, bin_chrom == y)])
        if np.all(np.isnan(blk)):
            continue
        ij = np.unravel_index(np.nanargmax(blk), blk.shape)
        bx, by = np.where(bin_chrom == x)[0][ij[0]], np.where(bin_chrom == y)[0][ij[1]]
        flag.append((blk[ij], x, y, bx, by))
flag.sort(reverse=True)
for r, x, y, bx, by in flag[:8]:
    lines.append("  %-14s %-14s %8.2f %14s %14s %8d" % (chrom_names[x], chrom_names[y], r,
                 "%.0f-%.0f" % (bin_start[bx] / 1e6, (bin_start[bx] + B) / 1e6),
                 "%.0f-%.0f" % (bin_start[by] / 1e6, (bin_start[by] + B) / 1e6), NS[bx, by]))
off = np.abs(R[bin_chrom[:, None] != bin_chrom[None, :]])
lines.append("  background: median inter-chromosomal |r| %.3f, 99th percentile %.3f"
             % (np.nanmedian(off), np.nanpercentile(off, 99)))
print("\n".join(lines[1:]))
with open(os.path.join(OUT, "genetic_contact_map.txt"), "w") as fh:
    fh.write("\n".join(lines) + "\n")
pd.DataFrame(R).to_csv(os.path.join(OUT, "genetic_contact_map.tsv.gz"), sep="\t", index=False)

fig, ax = plt.subplots(1, 2, figsize=(15, 7), gridspec_kw={"width_ratios": [1.15, 1]})
a0 = ax[0]
im = a0.imshow(R, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
for b0 in boff[1:]:
    a0.axhline(b0 - .5, color="black", lw=.5); a0.axvline(b0 - .5, color="black", lw=.5)
ticks = boff + nb / 2.0 - .5
a0.set_xticks(ticks); a0.set_yticks(ticks)
lab = [str(x).replace("_hap1", "").replace("chr", "") for x in chrom_names]
a0.set_xticklabels(lab, fontsize=8); a0.set_yticklabels(lab, fontsize=8)
a0.set_title("%s: do these bins travel together? (r across %d haploid cells)" % (SAMPLE, len(cells)), fontsize=10)
fig.colorbar(im, ax=a0, fraction=.046, pad=.02, label="r  (+1 same parent, 0 independent)")
a1 = ax[1]
a1.hist(off[~np.isnan(off)].ravel(), bins=np.linspace(0, 1, 51), color="#534AB7")
a1.set_yscale("log"); a1.set_xlabel("|r| between bins on DIFFERENT chromosomes")
a1.set_ylabel("bin pairs (log)")
a1.set_title("independent assortment piles up near 0; a translocation adds a tail")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "genetic_contact_map.png"), dpi=120)
fig.savefig(os.path.join(OUT, "genetic_contact_map.pdf"))
print("wrote %s/genetic_contact_map.{png,pdf,txt,tsv.gz}" % OUT)
