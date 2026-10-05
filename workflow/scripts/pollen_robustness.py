#!/usr/bin/env python3
"""pollen_robustness.py -- could cell QC, SNP calling or marker choices be producing the paradoxa linkage?

Recomputes r between whole arms under deliberately different choices. A real
inheritance pattern stays put; an artefact of one choice moves when that choice
changes. Windows are 10 Mb; each cell-window is called ALT/REF from the molecule
tables exactly as linkage_scan.py does unless a variant changes it.

  variants   baseline   80/20 calls, >= 3 molecules, all good cells (= linkage_scan)
             strict     90/10 calls
             deep       >= 8 molecules per window
             pure       only cell-windows whose molecules all agree (ALT share 0 or 1)
             unique     identical cells (same pollen grain) merged
             clean half / dirty half   cells split by their genome-wide share of pure windows
  blocks     median r over all window pairs between two arms (>= 10 cells each)

Also compares the cells that break L1-P with those that keep it (depth, purity).
Run from the CO_smk root. Usage: pollen_robustness.py SAMPLE OUT.txt
"""
import csv
import os
import sys

import numpy as np
import pandas as pd

S, OUT = sys.argv[1:3]
W = 10_000_000
ARMS = {  # composite (Dparadoxa_std) coordinates, Mb; 10 Mb on either side of each join left out
    "L1 (A)": ("chr1_hap1", 0, 250), "P (A)": ("chr1_hap1", 275, 410),
    "L2 (D)": ("chr2_hap2", 0, 205), "Q (D)": ("chr2_hap2", 225, 330),
    "chr3": ("chr3_hap1", 0, 290), "chr4 piece 0-85": ("chr4_hap1", 0, 80), "chr4 rest": ("chr4_hap1", 95, 480),
    "chr5": ("chr5_hap1", 70, 300), "chr6": ("chr6_hap1", 60, 330),
}
BLOCKS = [  # (arm, arm, what the tissue's chromosomes predict, what the pollen showed)
    ("L1 (A)", "P (A)", "linked", "0.5"), ("L2 (D)", "Q (D)", "linked", "0.5"),
    ("L1 (A)", "Q (D)", "linked", "0"), ("L2 (D)", "P (A)", "linked", "0"),
    ("L1 (A)", "L2 (D)", "linked", "0.5"), ("P (A)", "Q (D)", "linked", "0.5"),
    ("chr3", "chr4 piece 0-85", "linked", "0"), ("chr5", "chr6", "0.5", "low"),
    ("L1 (A)", "chr6", "0.5", "0.5"), ("chr3", "chr5", "0.5", "0.5"),
]

row = next(r for r in csv.DictReader(open("config/samples.csv")) if r["sample_id"] == S)
lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(row["assembly_fasta"] + ".fai")}
nwin = {c: max(1, int(round(lens[c] / W))) for c in lens}
wins = [(c, i) for c in lens for i in range(nwin[c])]
widx = {w: k for k, w in enumerate(wins)}
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
n, m = len(cells), len(wins)
F = np.full((n, m), np.nan)       # ALT share of molecules per cell and window
N = np.zeros((n, m))              # molecules per cell and window
for ci, bc in enumerate(cells):
    p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
    if not os.path.exists(p):
        continue
    t = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"])
    t = t[(t.rc != t.ac) & t.chrom.isin(lens)]
    t["alt"] = (t.ac > t.rc).astype(int)
    t["w"] = [widx[(c, min(q // W, nwin[c] - 1))] for c, q in zip(t.chrom, t.pos)]
    g = t.groupby("w").alt.agg(["mean", "size"])
    F[ci, g.index.to_numpy()] = g["mean"].to_numpy()
    N[ci, g.index.to_numpy()] = g["size"].to_numpy()


def calls(lo=0.2, hi=0.8, min_mol=3, pure=False):
    G = np.full((n, m), np.nan)
    ok = N >= min_mol
    if pure:
        G[ok & (F == 0)] = 0
        G[ok & (F == 1)] = 1
    else:
        G[ok & (F <= lo)] = 0
        G[ok & (F >= hi)] = 1
    return G


def rmat(G, rows):
    X = G[rows]
    Mk = (~np.isnan(X)).astype(float)
    A = np.nan_to_num(X) * Mk
    B = (1 - np.nan_to_num(X)) * Mk
    nn = Mk.T @ Mk
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(nn >= 10, 1 - (A.T @ A + B.T @ B) / nn, np.nan)


def arm_idx(a):
    c, s, e = ARMS[a]
    return [widx[(c, i)] for i in range(nwin.get(c, 0)) if s * 1e6 <= i * W and (i + 1) * W <= e * 1e6 + W / 2]


def block(R, a, b):
    v = R[np.ix_(arm_idx(a), arm_idx(b))]
    v = v[np.isfinite(v)]
    return (np.median(v), len(v)) if len(v) else (np.nan, 0)


# unique genotypes: merge cells that agree in >= 98% of >= 100 shared windows
G0 = calls()
parent = list(range(n))


def root(a):
    while parent[a] != a:
        parent[a] = parent[parent[a]]
        a = parent[a]
    return a


for i in range(n):
    for j in range(i + 1, n):
        ok = ~np.isnan(G0[i]) & ~np.isnan(G0[j])
        if ok.sum() >= 100 and (G0[i, ok] == G0[j, ok]).mean() >= 0.98:
            parent[root(i)] = root(j)
uniq = sorted({root(i) for i in range(n)})
called = ~np.isnan(G0)
purity = np.array([((F[i] == 0) | (F[i] == 1))[called[i]].mean() if called[i].any() else np.nan for i in range(n)])
order = np.argsort(-np.nan_to_num(purity, nan=-1))
clean, dirty = sorted(order[: n // 2]), sorted(order[n // 2:])

VARIANTS = [("baseline", calls(), list(range(n))), ("strict 90/10", calls(0.1, 0.9), list(range(n))),
            ("deep >=8 mol", calls(min_mol=8), list(range(n))), ("pure only", calls(pure=True), list(range(n))),
            ("unique cells", G0, uniq), ("clean half", G0, clean), ("dirty half", G0, dirty)]
Rs = [(name, rmat(G, rows)) for name, G, rows in VARIANTS]
out = ["POLLEN ROBUSTNESS  %s  (%d good cells; %d distinct genotypes; clean half purity >= %.2f)" % (
           S, n, len(uniq), purity[order[n // 2 - 1]]),
       "  median r over all window pairs between two arms; a real pattern stays put across columns", "",
       "  %-32s %-8s %-7s | %s" % ("arms", "tissue", "pollen", " ".join("%12s" % v for v, _ in Rs))]
for a, b, tissue, pol in BLOCKS:
    vals = [block(R, a, b) for _, R in Rs]
    out.append("  %-32s %-8s %-7s | %s" % ("%s x %s" % (a, b), tissue, pol,
                                          " ".join("%12s" % ("%.2f (%d)" % v if v[1] else "-") for v in vals)))
out += ["  (in brackets: window pairs with >= 10 cells)", ""]

# cells that break L1-P versus cells that keep it (arm-level calls: majority of called windows)
def arm_call(G, i, a):
    v = G[i, arm_idx(a)]
    v = v[~np.isnan(v)]
    return np.nan if len(v) < 2 else (1.0 if v.mean() >= 0.75 else (0.0 if v.mean() <= 0.25 else np.nan))


l1 = np.array([arm_call(G0, i, "L1 (A)") for i in range(n)])
pp = np.array([arm_call(G0, i, "P (A)") for i in range(n)])
qq = np.array([arm_call(G0, i, "Q (D)") for i in range(n)])
both = ~np.isnan(l1) & ~np.isnan(pp)
keep, brk = np.where(both & (l1 == pp))[0], np.where(both & (l1 != pp))[0]
tot = N.sum(axis=1)
out += ["CELLS THAT BREAK L1-P VS CELLS THAT KEEP IT  (arm-level call = >= 75% of an arm's called windows agree)",
        "  %-34s %6s %18s %16s %22s" % ("", "cells", "median molecules", "median purity", "same call at L1 and Q"),
        ]
for lab, idx in (("keep L1 = P", keep), ("break L1 != P", brk)):
    lq = [(l1[i] == qq[i]) for i in idx if not np.isnan(qq[i])]
    out.append("  %-34s %6d %18.0f %16.2f %15d of %d" % (lab, len(idx), np.median(tot[idx]) if len(idx) else np.nan,
                                                         np.nanmedian(purity[idx]) if len(idx) else np.nan,
                                                         sum(lq), len(lq)))
out.append("  (if the break were noise, breaking cells would be the shallow or impure ones, and L1 would not follow Q)")
open(OUT, "w").write("\n".join(out) + "\n")
print("\n".join(out))
