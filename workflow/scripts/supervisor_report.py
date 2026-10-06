#!/usr/bin/env python3
"""supervisor_report.py -- D. binata and D. paradoxa crossover landscapes, start to finish.

One self-contained HTML page for the supervisor. Every number, table and figure is computed here
from the pipeline's own outputs (each file read is listed in the appendix of the page); nothing
is typed in by hand. Sentences that state a number build it from the computed value.

  0 summary  1 material  2 pipeline at a glance  3 markers  4 cells  5 crossovers
  6 landscapes  7 is each reference right?  8 limitations  9 next steps  appendix

Usage (from the CO_smk root):
  supervisor_report.py [--binata Dbinata_hap1] [--paradoxa Dparadoxa_std] [--out reports/supervisor]
Writes OUT/report.html (figures embedded), OUT/fig/*.png and OUT/numbers.txt (every number used).
"""
import argparse
import base64
import csv
import datetime
import html
import math
import os
import re
import subprocess

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("--binata", default="Dbinata_hap1")
ap.add_argument("--paradoxa", default="Dparadoxa_std")
ap.add_argument("--out", default="reports/supervisor")
ap.add_argument("--hic_dir", default="/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/"
                                     "results/Drosera_paradoxa/qc/hic_remap")
A = ap.parse_args()
OUT, FIG = A.out, os.path.join(A.out, "fig")
os.makedirs(FIG, exist_ok=True)
SPP = [("binata", A.binata), ("paradoxa", A.paradoxa)]
COL = {"binata": "#0F6E56", "paradoxa": "#534AB7"}
NAME = {"binata": "D. binata", "paradoxa": "D. paradoxa"}
# translocation joins on the paradoxa crossover references, as in paradoxa_cb_check.py
JOINS = {"Dparadoxa_std": [("A join L1|P", "chr1_hap1", 262.9), ("D join L2|Q", "chr2_hap2", 215.0)],
         "Dparadoxa_CB": [("C join L1|Q", "chr1_hap2", 262.0), ("B join L2|P", "chr2_hap1", 211.0)]}
WIN_MB, MIN_MOL, BINS, NBOOT = 10, 3, 20, 500
USED, MISSING, NUM = [], [], {}
rng = np.random.default_rng(1)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


def have(p):
    if p and os.path.exists(p):
        USED.append(p)
        return True
    MISSING.append(p)
    return False


def natural(c):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", c)]


def put(sp, key, val):
    NUM["%s.%s" % (sp, key)] = val
    return val


def ok(x):
    return x is not None and not (isinstance(x, float) and math.isnan(x))


def fi(x):
    return "{:,}".format(int(round(x))) if ok(x) else "–"


def ff(x, d=1):
    return "%.*f" % (d, x) if ok(x) else "–"


def kv_file(p):
    out = {}
    for l in open(p):
        if ":" in l:
            k, v = l.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def lines_of(p):
    v = [l.strip() for l in open(p) if l.strip()]
    return v[1:] if v and not re.match(r"^[ACGTN]+(-\d+)?$", v[0]) else v


SAMPLES = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}
CFG = yaml.safe_load(open("config/config.yaml")) or {}

# =========================================================================== load
D = {}
for sp, S in SPP:
    row = SAMPLES[S]
    d = {"S": S, "row": row}
    fai = row["assembly_fasta"] + ".fai"
    if not have(fai):
        raise SystemExit("no reference index for %s: %s" % (S, fai))
    L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai) if l.strip()}
    n = int(row["chr_number_2n"]) // 2
    main = sorted(sorted(L, key=lambda c: -L[c])[:n], key=natural)
    d.update(n=n, main=main, L={c: L[c] for c in main}, G=sum(L[c] for c in main),
             G_all=sum(L.values()), n_seq=len(L))
    put(sp, "chromosomes", n)
    put(sp, "genome_mb_chromosomes", d["G"] / 1e6)

    # ---- markers
    p = "results/markers/%s/filter_summary.txt" % S
    if have(p):
        t = open(p).read()
        m = re.search(r"Raw variants:\s*(\d+)", t)
        d["raw"] = int(m.group(1)) if m else None
        m = re.search(r"Markers after filter:\s*(\d+)", t)
        d["filtered"] = int(m.group(1)) if m else None
        m = re.search(r"Filter params:\s*(.*)", t)
        d["filter_params"] = m.group(1).strip() if m else ""
    p = "results/blacklist/%s/rrna_exclude.bed" % S
    if have(p):
        b = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2], comment="#")
        d["rdna_bp"], d["rdna_n"] = int((b[2] - b[1]).sum()), len(b)
    p = "qc/markers/%s/marker_classes.tsv.gz" % S
    if have(p):
        mc = pd.read_csv(p, sep="\t", usecols=["chrom", "pos", "class"], dtype={"chrom": str})
        d["classes"] = mc["class"].value_counts().to_dict()
        d["seen"] = len(mc)
        d["good"] = mc[(mc["class"] == "good") & mc.chrom.isin(main)][["chrom", "pos"]].copy()
    p = row["annotation_gff"]
    if have(p):
        g = pd.read_csv(p, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4],
                        names=["chrom", "type", "start", "end"], dtype={0: str}, low_memory=False)
        d["genes"] = g[(g.type == "gene") & g.chrom.isin(main)][["chrom", "start", "end"]].copy()

    # ---- scRNA, cells
    p = "results/starsolo/%s/Log.final.out" % S
    if have(p):
        for l in open(p):
            if "Number of input reads" in l:
                d["reads"] = int(l.split("|")[1])
            if "Uniquely mapped reads %" in l:
                d["uniq_pct"] = float(l.split("|")[1].strip().rstrip("%"))
    p = "results/cells/%s/barcode_stats.tsv" % S
    if have(p):
        d["bs"] = pd.read_csv(p, sep="\t", usecols=["barcode", "total_umi", "n_genes"])
    p = "results/cells/%s/barcodes_called.tsv" % S
    if have(p):
        d["called"] = lines_of(p)
    p = "qc/haplotypes/%s/haplotype_tracks.tsv.gz" % S
    if have(p):
        d["ht"] = pd.read_csv(p, sep="\t")
    p = "results/cell_qc/%s/selection_summary.txt" % S
    if have(p):
        d["sel"] = kv_file(p)
    p = "results/cell_qc/%s/good_cells.tsv" % S
    d["cells"] = lines_of(p) if have(p) else []

    # ---- crossovers
    p = "results/crossovers/%s/co_per_cell.tsv" % S
    if have(p):
        d["per"] = pd.read_csv(p, sep="\t")
    p = "results/crossovers/%s/co_intervals.bed" % S
    if have(p):
        co = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "barcode"],
                         dtype={0: str})
        co = co[co.chrom.isin(main)].copy()
        co["mid"] = (co.start + co.end) / 2
        co["u"] = co.mid / co.chrom.map(d["L"])
        d["co"] = co
    p = "results/crossovers/%s/co_calling_params.txt" % S
    if have(p):
        d["co_params_text"] = open(p).read()
    p = "results/landscape/%s/gamma_interference_summary.tsv" % S
    if have(p):
        gm = pd.read_csv(p, sep="\t")
        d["gamma"] = gm[gm.iloc[:, 0].astype(str) == "POOLED"].to_dict("records")
    D[sp] = d

# ====================================================== window genotypes per cell
for sp, d in D.items():
    S, W = d["S"], int(WIN_MB * 1e6)
    nwin = {c: max(1, int(round(d["L"][c] / W))) for c in d["main"]}
    wins = [(c, i) for c in d["main"] for i in range(nwin[c])]
    widx = {w: k for k, w in enumerate(wins)}
    G = np.full((len(d["cells"]), len(wins)), np.nan)
    mol_cache = {}
    for ci, bc in enumerate(d["cells"]):
        p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
        if not os.path.exists(p):
            continue
        m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"],
                        dtype={0: str})
        m = m[(m.rc != m.ac) & m.chrom.isin(d["L"])].copy()
        m["alt"] = (m.ac > m.rc).astype(int)
        mol_cache[bc] = len(m)
        m["w"] = [min(int(q // W), nwin[c] - 1) for c, q in zip(m.chrom, m.pos)]
        g = m.groupby(["chrom", "w"]).alt.agg(["mean", "size"])
        for (c, w), r in g.iterrows():
            if r["size"] >= MIN_MOL and (r["mean"] >= 0.8 or r["mean"] <= 0.2):
                G[ci, widx[(c, w)]] = 1 if r["mean"] >= 0.8 else 0
    USED.append("results/cell_data_mol/%s/<barcode>.tsv (%d cells)" % (S, len(mol_cache)))
    d["mol_per_cell"] = pd.Series(mol_cache)
    M = (~np.isnan(G)).astype(float)
    X = np.nan_to_num(G) * M
    Y = (1 - np.nan_to_num(G)) * M
    shared = M @ M.T
    agree = X @ X.T + Y @ Y.T
    iu = np.triu_indices(len(d["cells"]), 1)
    sh, ag = shared[iu], agree[iu]
    keep = sh >= 20
    sim = ag[keep] / sh[keep]
    d["pair_sim_median"] = float(np.median(sim)) if len(sim) else float("nan")
    dup = keep.copy()
    dup[keep] = sim >= 0.95
    d["dup_pairs"] = int(dup.sum())
    d["dup_cells"] = len(set(iu[0][dup]) | set(iu[1][dup]))
    d["called_share"] = float(M.mean())

# ================================================================== computations
for sp, d in D.items():
    if "good" in d:
        gm = d["good"]
        gaps, ing, tot = [], 0, 0
        gi = {}
        if "genes" in d:
            for c, g in d["genes"].groupby("chrom"):
                s = g.sort_values("start")
                st, en = [], []
                for a, b in zip(s.start, s.end):
                    if st and a <= en[-1]:
                        en[-1] = max(en[-1], b)
                    else:
                        st.append(a)
                        en.append(b)
                gi[c] = (np.array(st), np.array(en))
        for c, g in gm.groupby("chrom"):
            pos = np.sort(g.pos.values)
            gaps += list(np.diff(np.concatenate([[0], pos, [d["L"][c]]])))
            if c in gi:
                st, en = gi[c]
                k = np.searchsorted(st, pos, side="right") - 1
                ing += int(((k >= 0) & (pos <= en[np.maximum(k, 0)])).sum())
            tot += len(pos)
        d["good_n"] = len(gm)
        d["good_per_mb"] = len(gm) / (d["G"] / 1e6)
        d["gap_median_kb"] = float(np.median(gaps)) / 1e3
        d["gap_max_mb"] = float(np.max(gaps)) / 1e6
        d["good_in_genes"] = ing / tot if tot and gi else float("nan")
    if "per" in d and "co" in d:
        per, co = d["per"], d["co"]
        ncell = len(per)
        d["ncell_co"] = ncell
        d["co_total"] = len(co)
        d["co_mean"], d["co_sd"], d["co_median"] = per.n_cos.mean(), per.n_cos.std(), per.n_cos.median()
        d["per_chrom"] = {c: (co.chrom == c).sum() / ncell for c in d["main"]}
        d["below_half"] = [c for c, v in d["per_chrom"].items() if v < 0.5]
        d["width_median_kb"] = float((co.end - co.start).median()) / 1e3
        # relative-position profile: per-cell bin counts, bootstrapped over cells
        bix = np.minimum((co.u.values * BINS).astype(int), BINS - 1)
        cidx = {b: k for k, b in enumerate(per.barcode)}
        C = np.zeros((ncell, BINS))
        rows = co.barcode.map(cidx)
        okr = rows.notna().values
        np.add.at(C, (rows[okr].astype(int).values, bix[okr]), 1)
        outer = [0, 1, BINS - 2, BINS - 1]
        inner = list(range(4, BINS - 4))

        def uidx(v):
            return (v[outer].sum() / 0.2) / (v[inner].sum() / 0.6) if v[inner].sum() > 0 else float("inf")

        prof, us = [], []
        for _ in range(NBOOT):
            v = C[rng.integers(0, ncell, ncell)].sum(0)
            prof.append(v / max(v.mean(), 1e-9))
            us.append(uidx(v))
        v = C.sum(0)
        d["prof"], d["prof_lo"], d["prof_hi"] = v / v.mean(), np.percentile(prof, 2.5, 0), np.percentile(prof, 97.5, 0)
        d["U"], d["U_lo"], d["U_hi"] = uidx(v), np.percentile(us, 2.5), np.percentile(us, 97.5)
        if "good" in d:
            mu = (d["good"].pos / d["good"].chrom.map(d["L"])).values
            h = np.histogram(mu, bins=BINS, range=(0, 1))[0].astype(float)
            d["mprof"], d["U_markers"] = h / h.mean(), uidx(h)
        if "genes" in d:
            g = d["genes"]
            gu = (((g.start + g.end) / 2) / g.chrom.map(d["L"])).values
            h = np.histogram(gu, bins=BINS, range=(0, 1))[0].astype(float)
            d["gprof"], d["U_genes"] = h / h.mean(), uidx(h)
        # paradoxa std: crossovers at the translocation joins (reference artefact)
        if d["S"] in JOINS:
            near = np.zeros(len(co), bool)
            for _, c, mb in JOINS[d["S"]]:
                near |= ((co.chrom == c) & ((co.mid - mb * 1e6).abs() <= 10e6)).values
            d["join_cos"] = int(near.sum())
            v2 = np.histogram(co.u.values[~near], bins=BINS, range=(0, 1))[0].astype(float)
            d["U_nojoin"] = uidx(v2)
        # window-level: crossovers against good markers (detectability)
        if "good" in d:
            a, b = [], []
            for c in d["main"]:
                k = int(math.ceil(d["L"][c] / 5e6))
                a += list(np.bincount(np.minimum((co.mid[co.chrom == c] // 5e6).astype(int), k - 1), minlength=k))
                b += list(np.bincount(np.minimum((d["good"].pos[d["good"].chrom == c] // 5e6).astype(int), k - 1),
                                      minlength=k))
            d["rho_markers"] = pd.Series(a).corr(pd.Series(b), method="spearman")

# ================================================================== numbers
for sp, d in D.items():
    for k in ("raw", "filtered", "seen", "good_n", "good_per_mb", "gap_median_kb", "gap_max_mb", "good_in_genes",
              "rdna_bp", "reads", "uniq_pct", "pair_sim_median", "dup_pairs", "dup_cells", "called_share",
              "ncell_co", "co_total", "co_mean", "co_sd", "co_median", "width_median_kb", "U", "U_lo", "U_hi",
              "U_markers", "U_genes", "join_cos", "U_nojoin", "rho_markers"):
        if k in d:
            put(sp, k, d[k])
    if "classes" in d:
        for k, v in d["classes"].items():
            put(sp, "marker_class_" + k, v)
    if "bs" in d:
        put(sp, "barcodes_tested", len(d["bs"]))
    if "called" in d:
        put(sp, "cells_called", len(d["called"]))
    put(sp, "cells_selected", len(d["cells"]))
    if "sel" in d:
        for k, v in d["sel"].items():
            put(sp, "selection." + k, v)
    if "ht" in d and d["cells"]:
        h = d["ht"][d["ht"].barcode.isin(d["cells"])]
        put(sp, "selected_molecules_median", h.molecules.median())
        put(sp, "selected_molecules_q25", h.molecules.quantile(0.25))
        put(sp, "selected_molecules_q75", h.molecules.quantile(0.75))
        put(sp, "selected_weakest_chrom_molecules_median", h.min_chrom_molecules.median())
        put(sp, "selected_haploidness_median", h.haploidness.median())
    if "bs" in d and d["cells"]:
        b = d["bs"][d["bs"].barcode.isin(d["cells"])]
        put(sp, "selected_umi_median", b.total_umi.median())
        put(sp, "selected_genes_median", b.n_genes.median())
    if "per_chrom" in d:
        for c, v in d["per_chrom"].items():
            put(sp, "cos_per_grain." + c, v)
        put(sp, "chromosomes_below_half", len(d["below_half"]))
    if d.get("gamma"):
        for k, v in d["gamma"][0].items():
            put(sp, "gamma_pooled." + str(k), v)


def shape(d):
    if not ok(d.get("U")):
        return "–"
    if math.isinf(d["U"]):
        return "U-shaped: no crossovers in the middle 60%"
    if d["U_lo"] > 1.5:
        return "U-shaped: chromosome ends %s× the middle" % ff(d["U"])
    if d["U_lo"] <= 1.2 and d["U_hi"] >= 0.83:
        return "flat within error: ends %s× the middle" % ff(d["U"], 2)
    return "intermediate: ends %s× the middle" % ff(d["U"])


for sp, d in D.items():
    d["shape"] = put(sp, "shape", shape(d))

# ================================================================== figures


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return p


FIGS = {}
# fig 1: good markers along the genome
fig, axs = plt.subplots(2, 1, figsize=(10, 4.2))
for ax, (sp, d) in zip(axs, D.items()):
    if "good" not in d:
        ax.text(0.5, 0.5, "marker classes missing", ha="center", transform=ax.transAxes)
        continue
    off = 0
    for i, c in enumerate(d["main"]):
        pos = d["good"].pos[d["good"].chrom == c].values
        k = int(math.ceil(d["L"][c] / 2e6))
        h = np.bincount(np.minimum((pos // 2e6).astype(int), k - 1), minlength=k) / 2.0
        x = off / 1e6 + (np.arange(k) + 0.5) * 2
        if i % 2:
            ax.axvspan(off / 1e6, (off + d["L"][c]) / 1e6, color="#F1EFE8", lw=0)
        ax.plot(x, h, color=COL[sp], lw=0.8)
        ax.text((off + d["L"][c] / 2) / 1e6, -0.08, re.sub(r"_hap\d", "", c).replace("chr", ""),
                transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7)
        off += d["L"][c]
    ax.set_xlim(0, off / 1e6)
    ax.set_xticks([])
    ax.set_ylabel("good markers / Mb")
    ax.set_title("%s (%s): %s good markers, %s per Mb" % (NAME[sp], d["S"], fi(d.get("good_n")),
                                                        ff(d.get("good_per_mb"))), loc="left", fontsize=9)
fig.tight_layout()
FIGS["markers"] = save(fig, "fig1_markers_along_genome.png")

# fig 2: cells -- knee plot and molecules vs haploidness
fig, axs = plt.subplots(2, 2, figsize=(10, 7))
for j, (sp, d) in enumerate(D.items()):
    ax = axs[0, j]
    if "bs" in d:
        u = np.sort(d["bs"].total_umi.values)[::-1]
        ax.loglog(np.arange(1, len(u) + 1), u, color="#888780", lw=1)
        if "called" in d:
            cu = d["bs"].set_index("barcode").total_umi.reindex(d["called"]).dropna()
            ax.axhline(cu.min(), color=COL[sp], lw=0.8, ls="--")
            ax.text(1.5, cu.min() * 1.15, "called: %s barcodes" % fi(len(d["called"])), color=COL[sp], fontsize=8)
    ax.set_xlabel("barcode rank")
    ax.set_ylabel("UMIs per barcode")
    ax.set_title("%s: cell calling (knee plot)" % NAME[sp], loc="left")
    ax = axs[1, j]
    if "ht" in d:
        h = d["ht"].dropna(subset=["haploidness"])
        sel = h.barcode.isin(d["cells"])
        ax.scatter(h.molecules[~sel], h.haploidness[~sel], s=4, color="#B4B2A9", alpha=0.6, lw=0,
                   label="not selected (%s)" % fi((~sel).sum()))
        ax.scatter(h.molecules[sel], h.haploidness[sel], s=6, color=COL[sp], alpha=0.8, lw=0,
                   label="selected (%s)" % fi(sel.sum()))
        ax.set_xscale("log")
        s = d.get("sel", {})
        for key, axis in (("min_haploidness", "h"), ("min_molecules", "v")):
            try:
                val = float(s.get(key, "nan"))
            except ValueError:
                val = float("nan")
            if ok(val) and val > 0:
                (ax.axhline if axis == "h" else ax.axvline)(val, color="#444441", lw=0.8, ls=":")
        ax.legend(frameon=False, fontsize=7, loc="lower right")
    ax.set_xlabel("informative molecules per barcode")
    ax.set_ylabel("haploidness (1 = one clean haplotype)")
    ax.set_title("%s: which barcodes are single haploid nuclei" % NAME[sp], loc="left")
fig.tight_layout()
FIGS["cells"] = save(fig, "fig2_cells.png")

# fig 3: one typical cell per species
fig, axs = plt.subplots(1, 2, figsize=(11, 6.5), gridspec_kw={"width_ratios": [1, 1]})
for ax, (sp, d) in zip(axs, D.items()):
    if "ht" not in d or not d["cells"]:
        continue
    h = d["ht"][d["ht"].barcode.isin(d["cells"])]
    bc = h.iloc[(h.molecules - h.molecules.median()).abs().argsort().iloc[0]].barcode
    d["example"] = bc
    m = pd.read_csv("results/cell_data_mol/%s/%s.tsv" % (d["S"], bc), sep="\t", header=None, usecols=[0, 1, 3, 5],
                    names=["chrom", "pos", "rc", "ac"], dtype={0: str})
    m = m[(m.rc != m.ac) & m.chrom.isin(d["L"])]
    for i, c in enumerate(d["main"]):
        y0 = len(d["main"]) - 1 - i
        mm = m[m.chrom == c].sort_values("pos")
        alt = (mm.ac > mm.rc).astype(float).values
        ax.scatter(mm.pos / 1e6, y0 + 0.08 + alt * 0.64, s=2, color="#B4B2A9", lw=0)
        if len(alt) >= 15:
            sm = pd.Series(alt).rolling(15, center=True, min_periods=8).mean().values
            ax.plot(mm.pos / 1e6, y0 + 0.08 + sm * 0.64, color=COL[sp], lw=1)
        if "co" in d:
            for _, r in d["co"][(d["co"].barcode == bc) & (d["co"].chrom == c)].iterrows():
                ax.fill_betweenx([y0, y0 + 0.8], r.start / 1e6, r.end / 1e6, color="#D85A30", alpha=0.3, lw=0)
                ax.plot([r.mid / 1e6] * 2, [y0, y0 + 0.8], color="#D85A30", lw=1)
        ax.plot([0, d["L"][c] / 1e6], [y0, y0], color="#D3D1C7", lw=0.5)
    ax.set_yticks([len(d["main"]) - 1 - i + 0.4 for i in range(len(d["main"]))])
    ax.set_yticklabels([re.sub(r"_hap\d", "", c) for c in d["main"]], fontsize=7)
    ax.set_ylim(-0.1, len(d["main"]))
    ax.set_xlabel("Mb")
    nco = int((d["co"].barcode == bc).sum()) if "co" in d else 0
    ax.set_title("%s: one typical cell (%s molecules, %d crossovers)" % (
        NAME[sp], fi(d["mol_per_cell"].get(bc, float("nan"))), nco), loc="left")
fig.tight_layout()
FIGS["example"] = save(fig, "fig3_example_cells.png")

# fig 4: crossovers per grain
fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
for ax, (sp, d) in zip(axs, D.items()):
    if "per" not in d:
        continue
    v = d["per"].n_cos.values
    ax.hist(v, bins=np.arange(-0.5, v.max() + 1.5, 1), color=COL[sp], alpha=0.85)
    ax.axvline(d["n"] / 2, color="#444441", ls="--", lw=1)
    ax.text(d["n"] / 2, ax.get_ylim()[1] * 0.95, " obligate minimum %s" % ff(d["n"] / 2), fontsize=7, va="top")
    ax.axvline(v.mean(), color="#D85A30", lw=1)
    ax.text(v.mean(), ax.get_ylim()[1] * 0.82, " mean %s" % ff(v.mean(), 2), fontsize=7, va="top", color="#D85A30")
    ax.set_xlabel("crossovers per pollen grain")
    ax.set_ylabel("grains")
    ax.set_title("%s (%s grains)" % (NAME[sp], fi(len(v))), loc="left")
fig.tight_layout()
FIGS["per_grain"] = save(fig, "fig4_crossovers_per_grain.png")

# fig 5: landscape along the chromosome, relative position
fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
x = (np.arange(BINS) + 0.5) / BINS * 100
for ax, (sp, d) in zip(axs, D.items()):
    if "prof" not in d:
        continue
    ax.fill_between(x, d["prof_lo"], d["prof_hi"], color=COL[sp], alpha=0.18, lw=0)
    ax.plot(x, d["prof"], color=COL[sp], lw=2, label="crossovers (95% CI over grains)")
    if "mprof" in d:
        ax.plot(x, d["mprof"], color="#888780", lw=1.2, label="good markers")
    if "gprof" in d:
        ax.plot(x, d["gprof"], color="#888780", lw=1, ls="--", label="genes")
    ax.axhline(1, color="#D3D1C7", lw=0.8)
    ax.set_xlabel("position along chromosome (% of length, all chromosomes pooled)")
    ax.set_title(NAME[sp], loc="left")
    ax.text(0.02, 0.97, d["shape"], transform=ax.transAxes, va="top", fontsize=8, color=COL[sp])
    ax.legend(frameon=False, fontsize=7, loc="upper right")
axs[0].set_ylabel("density relative to the chromosome mean")
fig.tight_layout()
FIGS["landscape"] = save(fig, "fig5_landscape_relative.png")

# fig 6: crossovers per chromosome
fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), gridspec_kw={"width_ratios": [max(len(D[s]["main"]), 1)
                                                                                 for s in D]})
for ax, (sp, d) in zip(axs, D.items()):
    if "per_chrom" not in d:
        continue
    cs = sorted(d["main"], key=lambda c: d["L"][c])
    ax.bar(range(len(cs)), [d["per_chrom"][c] for c in cs], color=COL[sp])
    ax.axhline(0.5, color="#444441", ls="--", lw=1)
    ax.set_xticks(range(len(cs)))
    ax.set_xticklabels(["%s\n%.0f" % (re.sub(r"_hap\d", "", c).replace("chr", ""), d["L"][c] / 1e6) for c in cs],
                       fontsize=6)
    ax.set_xlabel("chromosome (length, Mb), shortest to longest")
    ax.set_title("%s: crossovers per grain per chromosome" % NAME[sp], loc="left")
axs[0].set_ylabel("mean crossovers per grain")
fig.tight_layout()
FIGS["per_chrom"] = save(fig, "fig6_crossovers_per_chromosome.png")

# fig 7 (appendix): per-chromosome landscapes
ncol = 4
nrow = int(math.ceil(len(D["binata"]["main"]) / ncol)) + int(math.ceil(len(D["paradoxa"]["main"]) / ncol))
fig, axs = plt.subplots(nrow, ncol, figsize=(12, 1.9 * nrow))
axs = np.atleast_2d(axs)
r0 = 0
for sp, d in D.items():
    k = 0
    for c in d["main"]:
        ax = axs[r0 + k // ncol, k % ncol]
        if "co" in d:
            nb = int(math.ceil(d["L"][c] / 5e6))
            cnt = np.bincount(np.minimum((d["co"].mid[d["co"].chrom == c] // 5e6).astype(int), nb - 1), minlength=nb)
            cmmb = 100.0 * cnt / max(d.get("ncell_co", 1), 1) / 5
            ax.bar((np.arange(nb) + 0.5) * 5, cmmb, width=5, color=COL[sp], alpha=0.85)
        if "good" in d:
            pos = d["good"].pos[d["good"].chrom == c].values
            nb = int(math.ceil(d["L"][c] / 5e6))
            mk = np.bincount(np.minimum((pos // 5e6).astype(int), nb - 1), minlength=nb)
            ax2 = ax.twinx()
            ax2.plot((np.arange(nb) + 0.5) * 5, mk / 5, color="#888780", lw=0.8)
            ax2.set_yticks([])
            ax2.spines["right"].set_visible(False)
        ax.set_title("%s %s" % (NAME[sp].split()[1], c), fontsize=7, loc="left")
        ax.tick_params(labelsize=6)
        k += 1
    for j in range(k, int(math.ceil(len(d["main"]) / ncol)) * ncol):
        axs[r0 + j // ncol, j % ncol].axis("off")
    r0 += int(math.ceil(len(d["main"]) / ncol))
fig.text(0.0, 0.5, "cM/Mb (bars); good markers per Mb (grey line, own scale)", rotation=90, va="center", fontsize=8)
fig.tight_layout()
FIGS["per_chrom_landscapes"] = save(fig, "fig7_per_chromosome_landscapes.png")

# ================================================================== html helpers


def img(p, width=100):
    if not p or not os.path.exists(p):
        MISSING.append(p)
        return "<p class='missing'>figure missing: %s</p>" % html.escape(str(p))
    if p not in USED:
        USED.append(p)
    b = base64.b64encode(open(p, "rb").read()).decode()
    return "<img src='data:image/png;base64,%s' style='width:%d%%'>" % (b, width)


def table(head, rows, cls=""):
    h = "".join("<th>%s</th>" % html.escape(str(x)) for x in head)
    b = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % html.escape(str(x)) for x in r) for r in rows)
    return "<table class='%s'><tr>%s</tr>%s</table>" % (cls, h, b)


def details(title, p):
    if not p or not os.path.exists(p):
        MISSING.append(p)
        return "<p class='missing'>missing: %s</p>" % html.escape(str(p))
    if p not in USED:
        USED.append(p)
    return "<details><summary>%s <span class='path'>%s</span></summary><pre>%s</pre></details>" % (
        html.escape(title), html.escape(p), html.escape(open(p, errors="replace").read()))


def P(text):
    return "<p>%s</p>" % text


B, Pd = D["binata"], D["paradoxa"]


# ================================================================== paradoxa evidence
pairs = []
p = "qc/linkage/%s/pollen_pairs_check.txt" % Pd["S"]
if have(p):
    for l in open(p):
        if len(l) > 60 and l.startswith("   ") and not l.startswith("    "):
            lab, rest = l[3:53].strip(), l[53:].split()
            if len(rest) >= 6 and all(t.isdigit() for t in rest[:5]):
                pairs.append((lab, rest[0], rest[1], rest[2], rest[3], rest[4], rest[5]))
for lab, *vals in pairs:
    NUM["paradoxa.pair." + lab] = " ".join(vals)

# ================================================================== html
H = []
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
try:
    commit = subprocess.run("git rev-parse --short HEAD", shell=True, stdout=subprocess.PIPE,
                            universal_newlines=True).stdout.strip()
except OSError:
    commit = "?"

TAG = subprocess.run("git tag -l 'binata-landscape*' | tail -n 1", shell=True, stdout=subprocess.PIPE,
                     universal_newlines=True).stdout.strip()
H.append("<h1>Meiotic crossover landscapes of <i>Drosera binata</i> and <i>D. paradoxa</i> from single pollen nuclei</h1>")
H.append("<p class='meta'>Generated %s from CO_smk commit %s by workflow/scripts/supervisor_report.py. "
         "Every number, table and figure below is computed from the pipeline's outputs (files listed in the appendix).</p>"
         % (now, commit))

# 0 summary
H.append("<h2>0. Summary</h2>")
H.append(P("Single pollen nuclei of <i>D. binata</i> and <i>D. paradoxa</i> were sequenced (scRNA-seq) and genotyped at "
           "heterozygous SNPs called from HiFi reads; a crossover is called where a nucleus switches between the "
           "two haplotypes. "
           "<b>%s</b> (%s): %s grains, %s crossovers per grain (obligate minimum %s); landscape %s. "
           "<b>%s</b> (%s): %s grains, %s crossovers per grain (obligate minimum %s); landscape %s. "
           "The <i>D. paradoxa</i> landscape is provisional: pollen and tissue disagree on how chromosomes 1 and 2 "
           "are built (section 7)."
           % (NAME["binata"], B["row"]["centromere"], fi(B.get("ncell_co")), ff(B.get("co_mean"), 2), ff(B["n"] / 2),
              B["shape"], NAME["paradoxa"], Pd["row"]["centromere"], fi(Pd.get("ncell_co")), ff(Pd.get("co_mean"), 2),
              ff(Pd["n"] / 2), Pd["shape"])))
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["centromere type (config)", B["row"]["centromere"], Pd["row"]["centromere"]],
    ["chromosomes (n)", B["n"], Pd["n"]],
    ["good markers", fi(B.get("good_n")), fi(Pd.get("good_n"))],
    ["pollen nuclei used", fi(len(B["cells"])), fi(len(Pd["cells"]))],
    ["informative molecules per nucleus, median", fi(NUM.get("binata.selected_molecules_median")),
     fi(NUM.get("paradoxa.selected_molecules_median"))],
    ["crossovers per grain, mean (SD)", "%s (%s)" % (ff(B.get("co_mean"), 2), ff(B.get("co_sd"), 2)),
     "%s (%s)" % (ff(Pd.get("co_mean"), 2), ff(Pd.get("co_sd"), 2))],
    ["obligate minimum (one per bivalent)", ff(B["n"] / 2), ff(Pd["n"] / 2)],
    ["genetic map length, cM", fi(100 * B["co_mean"]) if ok(B.get("co_mean")) else "–",
     fi(100 * Pd["co_mean"]) if ok(Pd.get("co_mean")) else "–"],
    ["ends vs middle (95% CI)", "%s (%s–%s)" % (ff(B.get("U"), 2), ff(B.get("U_lo"), 2), ff(B.get("U_hi"), 2)),
     "%s (%s–%s)" % (ff(Pd.get("U"), 2), ff(Pd.get("U_lo"), 2), ff(Pd.get("U_hi"), 2))],
    ["status", "final" + (" (git tag %s)" % TAG if TAG else ""), "provisional (reference under review)"]]))

# 1 material
H.append("<h2>1. Material and data</h2>")
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["sample in the pipeline", B["S"], Pd["S"]],
    ["reference", B["row"]["haplotype"], Pd["row"]["haplotype"]],
    ["reference notes", B["row"]["notes"][:160], Pd["row"]["notes"][:160]],
    ["2n (config)", B["row"]["chr_number_2n"], Pd["row"]["chr_number_2n"]],
    ["chromosome sequence, Mb", ff(B["G"] / 1e6, 0), ff(Pd["G"] / 1e6, 0)],
    ["genes on chromosomes (annotation)", fi(len(B["genes"])) if "genes" in B else "–",
     fi(len(Pd["genes"])) if "genes" in Pd else "–"],
    ["scRNA chemistry", B["row"]["chemistry"], Pd["row"]["chemistry"]],
    ["scRNA reads", fi(B.get("reads")), fi(Pd.get("reads"))],
    ["uniquely mapped, %", ff(B.get("uniq_pct")), ff(Pd.get("uniq_pct"))]]))

# 2 pipeline
H.append("<h2>2. Pipeline at a glance</h2>")
H.append(P("Each step, what it does, and what survives it."))
H.append(table(["step", "what it does", NAME["binata"], NAME["paradoxa"]], [
    ["1 markers", "heterozygous SNPs called from HiFi reads (bcftools), filtered on quality, depth, "
                  "allele balance and rDNA", "%s raw → %s" % (fi(B.get("raw")), fi(B.get("filtered"))),
     "%s raw → %s" % (fi(Pd.get("raw")), fi(Pd.get("filtered")))],
    ["2 markers seen in pollen", "markers covered by the pollen reads, classed by how they segregate; "
                                 "'good' = ALT in about half the nuclei, one allele per nucleus",
     "%s seen → %s good" % (fi(B.get("seen")), fi(B.get("good_n"))),
     "%s seen → %s good" % (fi(Pd.get("seen")), fi(Pd.get("good_n")))],
    ["3 align pollen RNA", "STARsolo to the reference", "%s reads, %s%% unique" % (fi(B.get("reads")),
                                                                                   ff(B.get("uniq_pct"))),
     "%s reads, %s%% unique" % (fi(Pd.get("reads")), ff(Pd.get("uniq_pct")))],
    ["4 call cells", "EmptyDrops on UMI counts", "%s tested → %s called" % (
        fi(NUM.get("binata.barcodes_tested")), fi(NUM.get("binata.cells_called"))),
     "%s tested → %s called" % (fi(NUM.get("paradoxa.barcodes_tested")), fi(NUM.get("paradoxa.cells_called")))],
    ["5 keep single haploid nuclei", "one haplotype per window (haploidness) and enough molecules on every "
                                     "chromosome", "%s → %s" % (B.get("sel", {}).get("barcodes counted by cellsnp", "–"),
                                                               fi(len(B["cells"]))),
     "%s → %s" % (Pd.get("sel", {}).get("barcodes counted by cellsnp", "–"), fi(len(Pd["cells"])))],
    ["6 call crossovers", "hapCO: blocks of informative molecules, a crossover where the haplotype switches",
     "%s crossovers" % fi(B.get("co_total")), "%s crossovers" % fi(Pd.get("co_total"))]]))

# 3 markers
H.append("<h2>3. Markers: what defines the two haplotypes</h2>")
H.append(img(FIGS["markers"]))
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["raw variants (HiFi)", fi(B.get("raw")), fi(Pd.get("raw"))],
    ["after filters", fi(B.get("filtered")), fi(Pd.get("filtered"))],
    ["filters", B.get("filter_params", "–"), Pd.get("filter_params", "–")],
    ["rDNA excluded, kb", fi(B.get("rdna_bp", float("nan")) / 1e3) if "rdna_bp" in B else "–",
     fi(Pd.get("rdna_bp", float("nan")) / 1e3) if "rdna_bp" in Pd else "–"],
    ["seen in pollen", fi(B.get("seen")), fi(Pd.get("seen"))]] + [
    ["class: %s" % k, fi(B.get("classes", {}).get(k)), fi(Pd.get("classes", {}).get(k))]
    for k in ("good", "sticky", "paralog", "other", "unjudged")] + [
    ["good markers per Mb", ff(B.get("good_per_mb")), ff(Pd.get("good_per_mb"))],
    ["gap between good markers, median kb", ff(B.get("gap_median_kb")), ff(Pd.get("gap_median_kb"))],
    ["largest gap, Mb", ff(B.get("gap_max_mb")), ff(Pd.get("gap_max_mb"))],
    ["good markers inside genes, %", ff(100 * B.get("good_in_genes", float("nan"))),
     ff(100 * Pd.get("good_in_genes", float("nan")))]]))
H.append(P("Markers come from RNA, so they sit in expressed genes; their density sets where crossovers can be seen "
           "(section 6 checks the landscape against it)."))

# 4 cells
H.append("<h2>4. Cells: how many, how deep, how clean</h2>")
H.append(img(FIGS["cells"]))
sb, sp_ = B.get("sel", {}), Pd.get("sel", {})
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["barcodes tested / called", "%s / %s" % (fi(NUM.get("binata.barcodes_tested")), fi(NUM.get("binata.cells_called"))),
     "%s / %s" % (fi(NUM.get("paradoxa.barcodes_tested")), fi(NUM.get("paradoxa.cells_called")))]] + [
    [k, sb.get(k, "–"), sp_.get(k, "–")] for k in
    ("barcodes counted by cellsnp", "scored (>= 5 full windows)", "pass haploidness and windows",
     "min_haploidness", "min_windows", "min_molecules", "min_chrom_molecules")] + [
    ["selected nuclei", fi(len(B["cells"])), fi(len(Pd["cells"]))],
    ["UMIs per selected nucleus, median", fi(NUM.get("binata.selected_umi_median")),
     fi(NUM.get("paradoxa.selected_umi_median"))],
    ["genes per selected nucleus, median", fi(NUM.get("binata.selected_genes_median")),
     fi(NUM.get("paradoxa.selected_genes_median"))],
    ["informative molecules per nucleus, median (IQR)",
     "%s (%s–%s)" % (fi(NUM.get("binata.selected_molecules_median")), fi(NUM.get("binata.selected_molecules_q25")),
                     fi(NUM.get("binata.selected_molecules_q75"))),
     "%s (%s–%s)" % (fi(NUM.get("paradoxa.selected_molecules_median")), fi(NUM.get("paradoxa.selected_molecules_q25")),
                     fi(NUM.get("paradoxa.selected_molecules_q75")))],
    ["weakest chromosome, molecules, median", fi(NUM.get("binata.selected_weakest_chrom_molecules_median")),
     fi(NUM.get("paradoxa.selected_weakest_chrom_molecules_median"))],
    ["haploidness, median", ff(NUM.get("binata.selected_haploidness_median"), 2),
     ff(NUM.get("paradoxa.selected_haploidness_median"), 2)],
    ["%g Mb windows with a clean call, %%" % WIN_MB, ff(100 * B["called_share"]), ff(100 * Pd["called_share"])],
    ["genotype similarity between nuclei, median (unrelated = 0.5)", ff(B["pair_sim_median"], 2),
     ff(Pd["pair_sim_median"], 2)],
    ["near-identical pairs (>= 95% of shared windows)", "%d pairs, %d nuclei" % (B["dup_pairs"], B["dup_cells"]),
     "%d pairs, %d nuclei" % (Pd["dup_pairs"], Pd["dup_cells"])]]))
H.append(img(FIGS["example"]))
H.append(P("One nucleus of median depth per species: grey dots are molecules (top = ALT, bottom = REF), the line is "
           "a 15-molecule running ALT share, orange spans are the called crossover intervals."))

# 5 crossovers
H.append("<h2>5. Crossover calling: settings and checks</h2>")
cc = CFG.get("co_calling", {})
rowsp = []
for k in ("input_rows", "block_size", "marker_num", "terminal_marker_num", "base_af", "window_af", "genotype",
          "_approved"):
    rowsp.append([k, (cc.get(B["row"]["species"]) or {}).get(k, "–"), (cc.get(Pd["row"]["species"]) or {}).get(k, "–")])
H.append(table(["hapCO setting (config co_calling)", NAME["binata"], NAME["paradoxa"]], rowsp))
H.append(img(FIGS["per_grain"]))
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["grains", fi(B.get("ncell_co")), fi(Pd.get("ncell_co"))],
    ["crossovers", fi(B.get("co_total")), fi(Pd.get("co_total"))],
    ["per grain, mean (SD); median", "%s (%s); %s" % (ff(B.get("co_mean"), 2), ff(B.get("co_sd"), 2),
                                                       ff(B.get("co_median"), 0)),
     "%s (%s); %s" % (ff(Pd.get("co_mean"), 2), ff(Pd.get("co_sd"), 2), ff(Pd.get("co_median"), 0))],
    ["obligate minimum per grain (n / 2)", ff(B["n"] / 2), ff(Pd["n"] / 2)],
    ["chromosomes below 0.5 per grain", "%s of %d" % (fi(len(B.get("below_half", []))), B["n"]),
     "%s of %d" % (fi(len(Pd.get("below_half", []))), Pd["n"])],
    ["crossover interval, median kb (resolution)", fi(B.get("width_median_kb")), fi(Pd.get("width_median_kb"))]] +
    ([["crossovers within 10 Mb of the translocation joins", "–", fi(Pd["join_cos"])]] if "join_cos" in Pd else [])))

# 6 landscapes
H.append("<h2>6. Landscapes</h2>")
H.append(img(FIGS["landscape"]))
rows6 = [["crossovers: ends vs middle (95% CI)", "%s (%s–%s)" % (ff(B.get("U"), 2), ff(B.get("U_lo"), 2),
                                                                ff(B.get("U_hi"), 2)),
          "%s (%s–%s)" % (ff(Pd.get("U"), 2), ff(Pd.get("U_lo"), 2), ff(Pd.get("U_hi"), 2))],
         ["good markers: ends vs middle", ff(B.get("U_markers"), 2), ff(Pd.get("U_markers"), 2)],
         ["genes: ends vs middle", ff(B.get("U_genes"), 2), ff(Pd.get("U_genes"), 2)],
         ["crossovers vs markers per 5 Mb window, Spearman", ff(B.get("rho_markers"), 2), ff(Pd.get("rho_markers"), 2)]]
if "U_nojoin" in Pd:
    rows6.append(["crossovers: ends vs middle, without the join artefact", "–", ff(Pd["U_nojoin"], 2)])
H.append(table(["", NAME["binata"], NAME["paradoxa"]], rows6))
H.append(P("'Ends vs middle' = crossover density in the outer 20% of each chromosome (10% at each end) over the "
           "inner 60%; 1 = flat. The same ratio for markers and genes shows how much end bias detection alone could "
           "produce."))
H.append(img(FIGS["per_chrom"]))
if B.get("gamma") or Pd.get("gamma"):
    H.append(table(["crossover interference, gamma model, pooled (recombination_landscape.R)", NAME["binata"],
                    NAME["paradoxa"]],
                   [[k, B.get("gamma", [{}])[0].get(k, "–") if B.get("gamma") else "–",
                     Pd.get("gamma", [{}])[0].get(k, "–") if Pd.get("gamma") else "–"]
                    for k in (B.get("gamma") or Pd.get("gamma"))[0].keys()]))

# 7 references
H.append("<h2>7. Is each reference right?</h2>")
H.append(P("Pollen linkage scan: for every pair of 10 Mb windows, the share of nuclei whose genotypes disagree "
           "(r; near 0 = inherited together, near 0.5 = independent). A join in the reference that the pollen do not "
           "support shows up as r near 0.5 between neighbouring windows."))
for sp, d in D.items():
    H.append("<h3>%s (%s)</h3>" % (NAME[sp], d["S"]))
    H.append(details("linkage scan summary", "qc/linkage/%s/linkage_summary.txt" % d["S"]))
    H.append(img("qc/linkage/%s/linkage_map.png" % d["S"], 70))
if pairs:
    H.append("<h3><i>D. paradoxa</i>: window pairs across the translocations</h3>")
    H.append(P("Names: L1 and L2 are the left arms of chr1 and chr2, P and Q their right arms. "
               "A = chr1_hap1 (L1·P), B = chr2_hap1 (L2·P), C = chr1_hap2 (L1·Q), D = chr2_hap2 (L2·Q). "
               "Genotypes are REF/ALT relative to the reference used (%s)." % Pd["S"]))
    H.append(table(["window pair", "REF-REF", "REF-ALT", "ALT-REF", "ALT-ALT", "nuclei", "r"], pairs))
H.append(img(os.path.join(A.hic_dir, "hicmap_chr1_chr2_zoom_unique.png"), 80))
H.append(img(os.path.join(A.hic_dir, "hicmap_chr5_chr6.png"), 60))
H.append(img("qc/rdna/rdna_overview.png", 90))
for t, p in (("Hi-C, every library at the joins", os.path.join(A.hic_dir, "library_joins.txt")),
             ("chr5 against chr6", "qc/linkage/%s/chr5_chr6_check.txt" % Pd["S"]),
             ("rDNA arrays, crossover hotspots and joins", "qc/rdna/rdna_overview.txt"),
             ("crossover reference C + B against A + D", "qc/linkage/Dparadoxa_CB/cb_check.txt")):
    H.append(details(t, p))

# 8, 9
H.append("<h2>8. Limitations</h2>")
ref_line = ("The <i>D. paradoxa</i> reference used here joins chr1 and chr2 the way the tissue (Hi-C) does, not "
            "the way the pollen inherit them; crossovers at those joins are an artefact of that choice (section 7)."
            if Pd["S"] == "Dparadoxa_std" else
            "The <i>D. paradoxa</i> reference used here joins chr1 and chr2 the way the pollen inherit them; the "
            "tissue (Hi-C) shows the other arrangement, unresolved until cytology (section 7).")
H.append("<ul><li>Markers come from pollen RNA, so they cluster in expressed genes (%s%% and %s%% of good markers "
         "fall inside annotated genes); the largest marker gap is %s Mb in <i>D. binata</i> and %s Mb in "
         "<i>D. paradoxa</i>, and the median crossover interval is %s kb and %s kb.</li>"
         "<li><i>D. paradoxa</i> rests on %s nuclei against %s for <i>D. binata</i>.</li><li>%s</li></ul>"
         % (ff(100 * B.get("good_in_genes", float("nan")), 0), ff(100 * Pd.get("good_in_genes", float("nan")), 0),
            ff(B.get("gap_max_mb")), ff(Pd.get("gap_max_mb")), fi(B.get("width_median_kb")),
            fi(Pd.get("width_median_kb")), fi(len(Pd["cells"])), fi(len(B["cells"])), ref_line))
H.append("<h2>9. Next steps and asks</h2>")
H.append("<ul><li>Rerun <i>D. paradoxa</i> on the reference joined as the pollen inherit it "
         "(chr1_hap2 + chr2_hap1): %s.</li>"
         "<li>Cytology: meiotic chromosome spreads with oligo-FISH paints for the chr1/chr2 arms, chr5 and chr6, "
         "and 45S rDNA FISH.</li>"
         "<li>Then fix the crossover-calling sensitivity for <i>D. paradoxa</i> on the review panel, as done for "
         "<i>D. binata</i>.</li></ul>"
         % ("results included in section 7" if os.path.exists("qc/linkage/Dparadoxa_CB/cb_check.txt")
            else "running"))

# appendix
H.append("<h2>Appendix</h2>")
H.append("<h3>A. Per-chromosome landscapes</h3>")
H.append(img(FIGS["per_chrom_landscapes"]))
H.append("<h3>B. Crossover-calling parameters as run</h3>")
for sp, d in D.items():
    H.append(details("%s: co_calling_params.txt" % NAME[sp], "results/crossovers/%s/co_calling_params.txt" % d["S"]))
H.append("<h3>C. Files read</h3><pre>%s</pre>" % html.escape("\n".join(sorted(set(USED)))))
if MISSING:
    H.append("<h3>D. Files not found</h3><pre>%s</pre>" % html.escape("\n".join(sorted(set(map(str, MISSING))))))

CSS = """body{font-family:Helvetica,Arial,sans-serif;max-width:1050px;margin:24px auto;padding:0 16px;color:#2C2C2A;
font-size:13px;line-height:1.45}h1{font-size:21px}h2{font-size:16px;border-bottom:1px solid #D3D1C7;margin-top:30px}
h3{font-size:13px}table{border-collapse:collapse;margin:8px 0 14px;font-size:12px}td,th{border:1px solid #D3D1C7;
padding:3px 8px;text-align:left;vertical-align:top}th{background:#F1EFE8}.meta{color:#5F5E5A;font-size:11px}
.missing{color:#A32D2D}.path{color:#888780;font-size:11px}pre{font-size:10.5px;background:#F7F6F2;padding:8px;
overflow-x:auto}img{display:block;margin:8px 0}details{margin:4px 0}@media print{h2{page-break-before:auto}
img{page-break-inside:avoid}table{page-break-inside:avoid}}"""
open(os.path.join(OUT, "report.html"), "w").write(
    "<!doctype html><html><head><meta charset='utf-8'><title>Drosera crossover landscapes</title><style>%s</style>"
    "</head><body>%s</body></html>" % (CSS, "\n".join(H)))
with open(os.path.join(OUT, "numbers.txt"), "w") as f:
    for k in sorted(NUM):
        v = NUM[k]
        f.write("%s\t%s\n" % (k, ("%.4g" % v) if isinstance(v, (float, np.floating)) else v))
print("\n".join("%s\t%s" % (k, ("%.4g" % NUM[k]) if isinstance(NUM[k], (float, np.floating)) else NUM[k])
                for k in sorted(NUM)))
print("\nmissing files: %d%s" % (len(set(map(str, MISSING))), "".join("\n  " + str(m) for m in sorted(set(map(str, MISSING))))))
print("wrote %s/report.html, %s/numbers.txt, %d figures in %s" % (OUT, OUT, len(FIGS), FIG))
