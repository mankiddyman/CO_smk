#!/usr/bin/env python3
"""marker_segregation.py -- which markers carry information, and does using
only those rescue clean cells?

Built on one fact: in a pool of haploid meiotic products from a heterozygous
parent, half the nuclei inherit hap2 at any locus. So across cells:

  good marker      reads ALT in ~50% of cells, and one allele per cell
  sticky marker    reads REF (or ALT) in nearly every cell -- hap2 reads fail
                   to map there, so it carries no haplotype information
  paralog-like     both alleles inside single haploid cells, because two
                   copies collapsed onto one locus

Then per cell, switch rates recomputed on:
  all called markers | good markers only | good markers, pairs >=150 bp apart
and, separately, neighbouring pairs <100 bp (usually one read, one molecule)
against pairs >=150 bp (usually two molecules). If near pairs agree far more
than far pairs, a dense marker set flatters the switch rate.

Usage: marker_segregation.py SAMPLE ASM_VCF OUTDIR
Reads results/snps/SAMPLE/cellSNP.tag.{DP,AD}.mtx and cellSNP.base.vcf.gz.
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLE, ASM, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
SNP = os.path.join("results/snps", SAMPLE)
MIN_CELLS = 20           # a marker needs this many called cells to be judged
STICKY = 0.10            # ALT fraction below this (or above 1 - this): sticky
GOOD = (0.25, 0.75)      # ALT fraction window for a segregating marker
PARALOG = 0.40           # share of DP>=2 observations showing both alleles
NEAR, FAR = 100, 150     # bp: same read vs different molecules
os.makedirs(OUT, exist_ok=True)
say_lines = []


def say(s=""):
    say_lines.append(s)
    print(s, flush=True)


def read_mtx(path):
    df = pd.read_csv(path, sep=r"\s+", comment="%", header=None, dtype=np.int64)
    dims = df.iloc[0].values
    body = df.iloc[1:].values
    return int(dims[0]), int(dims[1]), body[:, 0] - 1, body[:, 1] - 1, body[:, 2]


# ---------------------------------------------------------------- load
say("MARKER SEGREGATION  %s" % SAMPLE)
sites = pd.read_csv(os.path.join(SNP, "cellSNP.base.vcf.gz"), sep="\t", comment="#", header=None,
                    usecols=[0, 1, 3, 4], names=["chrom", "pos", "ref", "alt"], dtype={0: str})
sites["ref"] = sites["ref"].str.upper(); sites["alt"] = sites["alt"].str.upper()
nsite, ncell, s_dp, c_dp, dp = read_mtx(os.path.join(SNP, "cellSNP.tag.DP.mtx"))
_, _, s_ad, c_ad, ad = read_mtx(os.path.join(SNP, "cellSNP.tag.AD.mtx"))
say("  %s sites x %s cells; %s observations" % (format(nsite, ","), format(ncell, ","), format(len(dp), ",")))

key_dp = s_dp * ncell + c_dp
order = np.argsort(key_dp, kind="stable")
pos_in = np.searchsorted(key_dp[order], s_ad * ncell + c_ad)
alt = np.zeros(len(dp), dtype=np.int64)
ok = (pos_in < len(order))
ok[ok] &= key_dp[order][pos_in[ok]] == (s_ad * ncell + c_ad)[ok]
alt[order[pos_in[ok]]] = ad[ok]
del key_dp, order, pos_in

frac = alt / np.maximum(dp, 1)
call = np.full(len(dp), -1, dtype=np.int8)
call[frac <= 0.2] = 0
call[frac >= 0.8] = 1

# ---------------------------------------------------------------- per marker
n_ref = np.bincount(s_dp[call == 0], minlength=nsite)
n_alt = np.bincount(s_dp[call == 1], minlength=nsite)
d2 = dp >= 2
n_d2 = np.bincount(s_dp[d2], minlength=nsite)
n_both = np.bincount(s_dp[d2 & (alt > 0) & (alt < dp)], minlength=nsite)
called = n_ref + n_alt
altf = np.where(called > 0, n_alt / np.maximum(called, 1), np.nan)
bothf = np.where(n_d2 >= 5, n_both / np.maximum(n_d2, 1), np.nan)

asm = pd.read_csv(ASM, sep="\t", comment="#", header=None, usecols=[0, 1, 3, 4],
                  names=["chrom", "pos", "ref", "alt"], dtype={0: str})
asm = asm[(asm["ref"].str.len() == 1) & (asm["alt"].str.len() == 1)]
asm["ref"] = asm["ref"].str.upper(); asm["alt"] = asm["alt"].str.upper()
asm["conf"] = True
conf = sites.merge(asm.drop_duplicates(["chrom", "pos"]), how="left",
                   on=["chrom", "pos", "ref", "alt"])["conf"].fillna(False).values.astype(bool)

judged = called >= MIN_CELLS
sticky = judged & ((altf < STICKY) | (altf > 1 - STICKY))
paralog = judged & ~sticky & (bothf >= PARALOG)
good = judged & (altf >= GOOD[0]) & (altf <= GOOD[1]) & ~(bothf >= PARALOG)
other = judged & ~sticky & ~paralog & ~good

say()
say("MARKERS seen in >= %d cells: %s of %s observed sites" % (MIN_CELLS, format(int(judged.sum()), ","),
                                                             format(nsite, ",")))
say("  %-28s %14s %14s %14s" % ("class", "all", "confirmed", "unconfirmed"))
for name, m in (("good (segregating ~50%)", good), ("sticky (ALT <10% or >90%)", sticky),
                ("paralog-like (both alleles)", paralog), ("in between", other)):
    say("  %-28s %8s %4.1f%% %8s %4.1f%% %8s %4.1f%%"
        % (name, format(int(m.sum()), ","), 100 * m.sum() / max(1, judged.sum()),
           format(int((m & conf).sum()), ","), 100 * (m & conf).sum() / max(1, (judged & conf).sum()),
           format(int((m & ~conf).sum()), ","), 100 * (m & ~conf).sum() / max(1, (judged & ~conf).sum())))
say("  median ALT fraction: confirmed %.3f, unconfirmed %.3f"
    % (np.nanmedian(altf[judged & conf]), np.nanmedian(altf[judged & ~conf])))

# ---------------------------------------------------------------- per cell
chrom_code = pd.factorize(sites["chrom"])[0].astype(np.int64)
spos = sites["pos"].values.astype(np.int64)


def pairs(mask, lo=0, hi=None):
    """Per-cell adjacent-pair and switch counts over called observations in mask."""
    i = np.nonzero(mask & (call >= 0))[0]
    c, s, g = c_dp[i], s_dp[i], call[i]
    ch, p = chrom_code[s], spos[s]
    o = np.lexsort((p, ch, c))
    c, ch, p, g = c[o], ch[o], p[o], g[o]
    same = (c[1:] == c[:-1]) & (ch[1:] == ch[:-1])
    gap = p[1:] - p[:-1]
    sel = same & (gap >= lo)
    if hi is not None:
        sel &= gap < hi
    npair = np.bincount(c[1:][sel], minlength=ncell)
    nsw = np.bincount(c[1:][sel & (g[1:] != g[:-1])], minlength=ncell)
    return npair, nsw


def rates(npair, nsw, min_pairs=50):
    ok = npair >= min_pairs
    return np.where(ok, nsw / np.maximum(npair, 1), np.nan)


everything = np.ones(len(dp), dtype=bool)
on_good = good[s_dp]
r_all = rates(*pairs(everything))
r_good = rates(*pairs(on_good))
r_good_far = rates(*pairs(on_good, lo=FAR))
r_near = rates(*pairs(everything, hi=NEAR))
r_far = rates(*pairs(everything, lo=FAR))

say()
say("SAME READ vs DIFFERENT MOLECULES (all markers, per-cell medians)")
say("  neighbours < %d bp apart:  %.3f   (usually one read)" % (NEAR, np.nanmedian(r_near)))
say("  neighbours >= %d bp apart: %.3f   (usually two molecules)" % (FAR, np.nanmedian(r_far)))

say()
say("CELLS: switch rate recomputed on")
say("  %-32s %10s %10s %12s %12s" % ("marker set", "cells", "median", "<= 0.10", "<= 0.05"))
for name, r in (("all called markers", r_all), ("good markers only", r_good),
                ("good markers, pairs >= %d bp" % FAR, r_good_far)):
    v = r[~np.isnan(r)]
    say("  %-32s %10s %10.3f %12s %12s" % (name, format(len(v), ","), np.median(v) if len(v) else np.nan,
                                           format(int((v <= 0.10).sum()), ","), format(int((v <= 0.05).sum()), ",")))

# ---------------------------------------------------------------- plots
fig, ax = plt.subplots(2, 2, figsize=(12, 9))
b = np.linspace(0, 1, 51)
a = ax[0, 0]
for m, col, lab in ((judged & conf, "#2E7D32", "confirmed"), (judged & ~conf, "#534AB7", "unconfirmed")):
    a.hist(altf[m], bins=b, histtype="step", lw=2, color=col, density=True,
           label="%s (n=%s)" % (lab, format(int(m.sum()), ",")))
for x in (STICKY, 1 - STICKY):
    a.axvline(x, color="grey", ls=":", lw=1)
a.axvline(0.5, color="black", ls="--", lw=.8)
a.set_xlabel("fraction of cells reading ALT at the marker"); a.set_ylabel("density")
a.set_title("%s: do markers segregate 1:1?" % SAMPLE); a.legend(frameon=False, fontsize=8)

a = ax[0, 1]
for m, col, lab in ((judged & conf, "#2E7D32", "confirmed"), (judged & ~conf, "#534AB7", "unconfirmed")):
    v = bothf[m & ~np.isnan(bothf)]
    if len(v):
        a.hist(v, bins=b, histtype="step", lw=2, color=col, density=True, label=lab)
a.axvline(PARALOG, color="grey", ls=":", lw=1)
a.set_xlabel("share of DP>=2 observations showing BOTH alleles in one cell")
a.set_ylabel("density"); a.set_title("paralog-like markers (right of dotted line)")
a.legend(frameon=False, fontsize=8)

a = ax[1, 0]
bs = np.linspace(0, 0.4, 81)
a.hist(r_near[~np.isnan(r_near)], bins=bs, histtype="step", lw=2, color="#D85A30",
       label="pairs < %d bp (same read)" % NEAR)
a.hist(r_far[~np.isnan(r_far)], bins=bs, histtype="step", lw=2, color="#534AB7",
       label="pairs >= %d bp (two molecules)" % FAR)
a.set_xlabel("switch rate per cell"); a.set_ylabel("cells")
a.set_title("does one read hide the mixing?"); a.legend(frameon=False, fontsize=8)

a = ax[1, 1]
for r, col, lab in ((r_all, "#888780", "all markers"), (r_good, "#2E7D32", "good markers"),
                    (r_good_far, "#534AB7", "good, pairs >= %d bp" % FAR)):
    v = r[~np.isnan(r)]
    a.hist(v, bins=bs, histtype="step", lw=2, color=col, label="%s (n=%s)" % (lab, format(len(v), ",")))
a.axvline(0.10, color="black", ls="--", lw=.8)
a.set_xlabel("switch rate per cell"); a.set_ylabel("cells")
a.set_title("does a clean population appear on good markers?"); a.legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "marker_segregation.png"), dpi=120)
fig.savefig(os.path.join(OUT, "marker_segregation.pdf"))

# genome tracks: where do sticky markers sit?
chroms = list(dict.fromkeys(sites["chrom"]))
nc = len(chroms); ncol = 4 if nc > 6 else 3
nrow = int(np.ceil(nc / float(ncol)))
fig, ax = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.6 * nrow), squeeze=False)
W = 2_000_000
for k, ch in enumerate(chroms):
    a = ax[k // ncol][k % ncol]
    m = judged & (sites["chrom"].values == ch)
    if not m.any():
        a.set_axis_off(); continue
    w = spos[m] // W
    nw = int(w.max()) + 1
    cnt = np.bincount(w, minlength=nw)
    mean_alt = np.bincount(w, weights=np.nan_to_num(altf[m]), minlength=nw) / np.maximum(cnt, 1)
    f_st = np.bincount(w, weights=sticky[m].astype(float), minlength=nw) / np.maximum(cnt, 1)
    f_cf = np.bincount(w, weights=conf[m].astype(float), minlength=nw) / np.maximum(cnt, 1)
    x = (np.arange(nw) + 0.5) * W / 1e6
    keep = cnt >= 20
    a.plot(x[keep], mean_alt[keep], color="#2E7D32", lw=1.2, label="mean ALT fraction")
    a.plot(x[keep], f_st[keep], color="#D85A30", lw=1.2, label="share sticky")
    a.plot(x[keep], f_cf[keep], color="#534AB7", lw=.8, alpha=.7, label="share confirmed")
    a.axhline(0.5, color="grey", ls=":", lw=.8)
    a.set_ylim(0, 1); a.set_title(ch, fontsize=9); a.tick_params(labelsize=7)
    if k == 0:
        a.legend(frameon=False, fontsize=6)
for k in range(nc, nrow * ncol):
    ax[k // ncol][k % ncol].set_axis_off()
fig.suptitle("%s: marker behaviour along the genome (2 Mb windows; x in Mb)" % SAMPLE)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "marker_tracks.png"), dpi=120)

pd.DataFrame({"chrom": sites["chrom"], "pos": spos, "confirmed": conf, "cells_called": called,
              "alt_frac": altf, "both_frac_dp2": bothf,
              "class": np.select([good, sticky, paralog, other], ["good", "sticky", "paralog", "other"],
                                 default="unjudged")}).to_csv(
    os.path.join(OUT, "marker_classes.tsv.gz"), sep="\t", index=False)
with open(os.path.join(OUT, "marker_segregation.txt"), "w") as f:
    f.write("\n".join(say_lines) + "\n")
say("\nwrote %s/{marker_segregation.png,.pdf,.txt, marker_tracks.png, marker_classes.tsv.gz}" % OUT)
