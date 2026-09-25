#!/usr/bin/env python3
"""co_sweep_where.py -- where do the suspicious crossovers come from?

For chosen sweep settings, from the per-cell *_co_pred.txt files:
  1. per chromosome: COs per cell and the share of neighbouring CO pairs
     closer than 10% of the chromosome (real length, from the .fai)
  2. COs per cell by haploidness tier AND depth tertile, so noise (more COs
     in dirtier cells at the same depth) is separated from detection power
     (more COs in deeper cells at the same cleanliness)
  3. for samples with known translocation breakpoints: are close-double COs
     piled up near them (unbalanced gametes) or spread everywhere (sparse-data
     artefacts)? Plotted as CO position histograms per chromosome.

Usage: co_sweep_where.py SWEEP_ROOT OUTDIR SAMPLE:FAI:SETTING[,SETTING] ...
       SETTING like mn8_bs2000000
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT, OUT = sys.argv[1], sys.argv[2]
SPECS = [a.split(":") for a in sys.argv[3:]]
BREAKS = {  # genetic contact map, 25 Mb bins: (chrom, start, end) pairs that travel together
    "Dparadoxa_hap1": [("chr5_hap1", 100e6, 125e6), ("chr6_hap1", 125e6, 150e6),
                       ("chr3_hap1", 275e6, 300e6), ("chr4_hap1", 50e6, 75e6)],
}
os.makedirs(OUT, exist_ok=True)
lines = []


def say(s=""):
    lines.append(s)
    print(s)


def load(s, setting, chroms):
    d = os.path.join(ROOT, s, setting, "per_cell")
    rows, per_cell = [], {}
    for f in os.listdir(d):
        if not f.endswith("_co_pred.txt"):
            continue
        bc = f[:-len("_co_pred.txt")]
        per_cell[bc] = 0
        for l in open(os.path.join(d, f)):
            x = l.split()
            if len(x) >= 3 and x[0] in chroms:
                try:
                    rows.append((bc, x[0], (int(x[1]) + int(x[2])) / 2.0))
                    per_cell[bc] += 1
                except ValueError:
                    pass
    return pd.DataFrame(rows, columns=["bc", "chrom", "mid"]), pd.Series(per_cell)


for s, fai, settings in SPECS:
    L = pd.read_csv(fai, sep="\t", header=None, usecols=[0, 1], names=["chrom", "len"]).set_index("chrom")["len"]
    with open(os.path.join("results/cell_qc", s, "switches.tsv")) as f:
        chroms = [c[len("markers_"):] for c in f.readline().rstrip("\n").split("\t") if c.startswith("markers_")]
    ht = pd.read_csv(os.path.join("qc/haplotypes", s, "haplotype_tracks.tsv.gz"), sep="\t",
                     usecols=["barcode", "molecules", "haploidness"]).set_index("barcode")
    for setting in settings.split(","):
        co, n = load(s, setting, chroms)
        co["close"] = False
        for (bc, ch), g in co.groupby(["bc", "chrom"]):
            if len(g) > 1:
                v = g.sort_values("mid")
                dd = np.diff(v["mid"].values)
                near = np.zeros(len(v), dtype=bool)
                near[:-1] |= dd < 0.1 * L[ch]
                near[1:] |= dd < 0.1 * L[ch]
                co.loc[v.index, "close"] = near
        say("=== %s  %s  (%d cells, %.2f COs/cell)" % (s, setting, len(n), n.mean()))
        say("  %-12s %8s %10s %14s" % ("chrom", "Mb", "COs/cell", "in close pair"))
        for ch in chroms:
            g = co[co.chrom == ch]
            say("  %-12s %8.0f %10.2f %13.0f%%" % (ch, L[ch] / 1e6, len(g) / float(len(n)),
                                                   100.0 * g.close.mean() if len(g) else 0))
        # cleanliness x depth
        t = ht.reindex(n.index).assign(cos=n.values)
        t["tier"] = pd.cut(t.haploidness, [0.6, 0.7, 0.8, 0.9, 1.0], include_lowest=True)
        t["depth"] = pd.qcut(t.molecules, 3, labels=["shallow", "middle", "deep"])
        piv = t.pivot_table(index="tier", columns="depth", values="cos", aggfunc="mean", observed=False)
        cnt = t.pivot_table(index="tier", columns="depth", values="cos", aggfunc="size", observed=False)
        say("  COs/cell by haploidness tier (rows) and depth tertile (columns); cells in brackets")
        for tier in piv.index:
            say("    %-12s " % str(tier) + "  ".join("%6.2f (%4d)" % (piv.loc[tier, c], cnt.loc[tier, c])
                                                    for c in piv.columns))
        if s in BREAKS:
            nearbp = np.zeros(len(co), dtype=bool)
            frac_genome = 0.0
            for ch, a, b in BREAKS[s]:
                nearbp |= ((co.chrom == ch) & (co.mid >= a - 12.5e6) & (co.mid <= b + 12.5e6)).values
                frac_genome += (b - a + 25e6)
            frac_genome /= float(L[chroms].sum())
            say("  near a translocation breakpoint (+/- 12.5 Mb): %.0f%% of the genome"
                % (100 * frac_genome))
            say("    all COs:            %.0f%% of calls" % (100 * nearbp.mean()))
            say("    close-double COs:   %.0f%% of calls" % (100 * nearbp[co.close.values].mean()
                                                         if co.close.any() else 0))
            fig, ax = plt.subplots(2, (len(chroms) + 1) // 2, figsize=(4 * ((len(chroms) + 1) // 2), 6.5),
                                   squeeze=False)
            for k, ch in enumerate(chroms):
                a = ax[k // ax.shape[1]][k % ax.shape[1]]
                g = co[co.chrom == ch]
                bins = np.linspace(0, L[ch] / 1e6, 40)
                a.hist(g.mid / 1e6, bins=bins, color="#888780", label="all COs")
                a.hist(g[g.close].mid / 1e6, bins=bins, color="#D85A30", alpha=.85, label="in a close pair")
                for c2, lo, hi in BREAKS[s]:
                    if c2 == ch:
                        a.axvspan(lo / 1e6, hi / 1e6, color="#534AB7", alpha=.2, label="breakpoint bin")
                a.set_title(ch, fontsize=9); a.set_xlabel("Mb", fontsize=8); a.tick_params(labelsize=7)
                if k == 0:
                    a.legend(frameon=False, fontsize=7)
            fig.suptitle("%s %s: where the crossovers are called" % (s, setting))
            fig.tight_layout()
            fig.savefig(os.path.join(OUT, "%s_%s_positions.png" % (s, setting)), dpi=120)
        say()

with open(os.path.join(OUT, "co_sweep_where.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
print("wrote %s/co_sweep_where.txt and *_positions.png" % OUT)
