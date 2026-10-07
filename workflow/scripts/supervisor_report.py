#!/usr/bin/env python3
"""supervisor_report.py -- D. binata and D. paradoxa crossover landscapes, start to finish.

One self-contained HTML page for the supervisor, plus one PDF per species with every nucleus's
crossover plot. Every number, table and figure is computed here from the pipeline's own outputs
(each file read is listed in the appendix of the page); nothing is typed in by hand. Sentences
that state a number build it from the computed value. A section that cannot be built is shown
in red with the error instead of stopping the report.

  0 summary  1 material  2 pipeline at a glance  3 markers  4 cells  5 crossovers
  6 landscapes  7 is each reference right?  8 limitations  9 next steps  appendix

Usage (from the CO_smk root):
  supervisor_report.py [--binata Dbinata_hap1] [--paradoxa Dparadoxa_std] [--out reports/supervisor]
                       [--no_cell_pdfs]
Writes OUT/report.html, OUT/cells_<sample>.pdf, OUT/fig/*.png and OUT/numbers.txt.

D. paradoxa references: Dparadoxa_std (A + D, chr1/chr2 joined as assembled) or Dparadoxa_CB (C + B, joined as the
pollen inherit them). Crossovers near joins the pollen do not support (ARTEFACT_JOINS) are shown with and without;
on CB they are kept and only counted against the chromosome's average. All text names the arms of the reference used.
"""
import argparse
import base64
import csv
import datetime
import glob
import html
import math
import os
import re
import shlex
import shutil
import subprocess
import traceback

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ap = argparse.ArgumentParser()
ap.add_argument("--binata", default="Dbinata_hap1")
ap.add_argument("--paradoxa", default="Dparadoxa_std")
ap.add_argument("--out", default="reports/supervisor")
ap.add_argument("--refs_root", default="/netscratch/dep_mercier/grp_marques/Aaryan/refs")
ap.add_argument("--phd_root", default="/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd")
ap.add_argument("--no_cell_pdfs", action="store_true")
A = ap.parse_args()
OUT, FIG = A.out, os.path.join(A.out, "fig")
os.makedirs(FIG, exist_ok=True)
SPP = [("binata", A.binata), ("paradoxa", A.paradoxa)]
COL = {"binata": "#0F6E56", "paradoxa": "#534AB7"}
NAME = {"binata": "D. binata", "paradoxa": "D. paradoxa"}
ARM = {"L1": "#5DCAA5", "L2": "#F0997B", "P": "#AFA9EC", "Q": "#FAC775"}
# translocation joins on the paradoxa crossover references, as in paradoxa_cb_check.py
JOINS = {"Dparadoxa_std": [("A join L1|P", "chr1_hap1", 262.9), ("D join L2|Q", "chr2_hap2", 215.0)],
         "Dparadoxa_CB": [("C join L1|Q", "chr1_hap2", 262.0), ("B join L2|P", "chr2_hap1", 211.0)]}
ARMSPAN = {"chr1_hap1": (262.9, "L1", "P"), "chr2_hap2": (215.0, "L2", "Q"),   # A, D: the joins of JOINS
           "chr1_hap2": (262.0, "L1", "Q"), "chr2_hap1": (211.0, "L2", "P")}   # C, B
ARMNAME = {"chr1_hap1": "A", "chr1_hap2": "C", "chr2_hap1": "B", "chr2_hap2": "D"}
PAIRS_FILE = "qc/linkage/Dparadoxa_std/pollen_pairs_check.txt"   # the conflict is measured on A + D
# References whose chr1/chr2 joins the pollen do NOT support (linkage r ~ 0.5 across them): crossovers within 10 Mb
# of these joins are an artefact and are shown with and without. On Dparadoxa_CB the joins are as the pollen inherit
# them, so crossovers there are kept and only counted against the chromosome's average.
ARTEFACT_JOINS = {"Dparadoxa_std"}


def join_names(S):
    """('chr1_hap2 and chr2_hap1', 'C and B', 'Q and P') for the reference's two chr1/chr2 joins."""
    js = JOINS.get(S, [])
    return (" and ".join(c for _, c, _ in js), " and ".join(ARMNAME.get(c, c) for _, c, _ in js),
            " and ".join(ARMSPAN[c][2] for _, c, _ in js if c in ARMSPAN))
WIN_MB, MIN_MOL, BINS, NBOOT = 10, 3, 20, 500
USED, MISSING, NUM, PROBLEMS = [], [], {}, []
rng = np.random.default_rng(1)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


# ============================================================================ helpers
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
    return "{:,}".format(int(round(x))) if ok(x) and not (isinstance(x, float) and math.isinf(x)) else "–"


def ff(x, d=1):
    return "%.*f" % (d, x) if ok(x) else "–"


def kv_file(p, sep=":"):
    out = {}
    for l in open(p):
        if sep in l:
            k, v = l.split(sep, 1)
            out[k.strip()] = v.strip()
    return out


def lines_of(p):
    v = [l.strip() for l in open(p) if l.strip()]
    return v[1:] if v and not re.match(r"^[ACGTN]+(-\d+)?$", v[0]) else v


def tool(name):
    c = sorted(glob.glob(".snakemake/conda/*/bin/%s" % name))
    return os.path.abspath(c[0]) if c else shutil.which(name)


def fai_lens(p):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(p) if l.strip()}


def guarded(name, fn):
    try:
        return fn()
    except Exception as e:
        PROBLEMS.append((name, "%s: %s" % (type(e).__name__, e)))
        traceback.print_exc()
        return None


def poisson_tail(k, mu):
    """P(X >= k), X ~ Poisson(mu), in log space."""
    if k <= 0:
        return 1.0
    s = sum(math.exp(-mu + i * math.log(mu) - math.lgamma(i + 1)) for i in range(int(k)))
    return max(0.0, 1.0 - s)


def uidx(v):
    """Crossover density in the outer 20% of the chromosome over the inner 60% (BINS bins)."""
    o, i = v[[0, 1, BINS - 2, BINS - 1]].sum(), v[4:BINS - 4].sum()
    return (o / 0.2) / (i / 0.6) if i > 0 else float("inf")


def profile(C):
    """Relative-position profile and ends-vs-middle ratio, bootstrapped over cells (rows of C)."""
    v = C.sum(0)
    prof, us = [], []
    for _ in range(NBOOT):
        vv = C[rng.integers(0, len(C), len(C))].sum(0)
        prof.append(vv / max(vv.mean(), 1e-9))
        us.append(uidx(vv))
    return dict(prof=v / max(v.mean(), 1e-9), lo=np.percentile(prof, 2.5, 0), hi=np.percentile(prof, 97.5, 0),
                U=uidx(v), U_lo=np.percentile(us, 2.5), U_hi=np.percentile(us, 97.5))


def shape(U, lo, hi):
    if not ok(U):
        return "–"
    if math.isinf(U):
        return "U-shaped: no crossovers in the middle 60%"
    ci = "95%% CI %s–%s" % (ff(lo, 2), ff(hi, 2))
    if lo > 1.5:
        return "U-shaped: ends %s× the middle (%s)" % (ff(U), ci)
    if lo <= 1 <= hi:
        return "flat: ends %s× the middle (%s)" % (ff(U, 2), ci)
    if hi < 1:
        return "ends depleted: %s× the middle (%s)" % (ff(U, 2), ci)
    return "mild end bias: ends %s× the middle (%s)" % (ff(U, 2), ci)


def code_const(path, pattern, n):
    """Read constants from a pipeline script, so definitions in the text match the code."""
    m = re.search(pattern, open(path).read()) if have(path) else None
    return [float(x) for x in m.groups()][:n] if m else [float("nan")] * n


SAMPLES = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}
CFG = yaml.safe_load(open("config/config.yaml")) or {}
MIN_CELLS, = code_const("workflow/scripts/marker_segregation.py", r"MIN_CELLS\s*=\s*(\d+)", 1)
STICKY, = code_const("workflow/scripts/marker_segregation.py", r"STICKY\s*=\s*([\d.]+)", 1)
GOOD_LO, GOOD_HI = code_const("workflow/scripts/marker_segregation.py", r"GOOD\s*=\s*\(([\d.]+),\s*([\d.]+)\)", 2)
PARALOG, = code_const("workflow/scripts/marker_segregation.py", r"PARALOG\s*=\s*([\d.]+)", 1)
HT_W, HT_GAP, HT_MINWIN = code_const("workflow/scripts/haplotype_tracks.py", r"W, GAP, MIN_WIN = (\d+), (\d+), (\d+)", 3)

# =========================================================================== load
D = {}
for sp, S in SPP:
    row = SAMPLES[S]
    d = {"S": S, "row": row, "species": row["species"]}
    fai = row["assembly_fasta"] + ".fai"
    if not have(fai):
        raise SystemExit("no reference index for %s: %s" % (S, fai))
    L = fai_lens(fai)
    n = int(row["chr_number_2n"]) // 2
    main = sorted(sorted(L, key=lambda c: -L[c])[:n], key=natural)
    d.update(n=n, main=main, L={c: L[c] for c in main}, G=sum(L[c] for c in main))
    put(sp, "chromosomes", n)
    put(sp, "reference_1C_mb", d["G"] / 1e6)
    # the phased two-haplotype assembly the reference came from (MANIFEST of the published hap1)
    for man in (os.path.join(os.path.dirname(row["assembly_fasta"]), "MANIFEST.txt"),
                os.path.join(A.refs_root, "%s_hap1" % row["species"], "MANIFEST.txt")):
        if os.path.exists(man):
            USED.append(man)
            mf = kv_file(man)
            d["asm2n"] = mf.get("source_fasta")
            if d["asm2n"] and have(d["asm2n"] + ".fai"):
                L2 = fai_lens(d["asm2n"] + ".fai")
                d["L2n"] = L2
                d["asm2n_chr_mb"] = sum(v for c, v in L2.items() if re.match(r"^chr\d+_hap[12]$", c)) / 1e6
                d["asm2n_all_mb"] = sum(L2.values()) / 1e6
                put(sp, "assembly_2C_chromosomes_mb", d["asm2n_chr_mb"])
                put(sp, "assembly_2C_all_mb", d["asm2n_all_mb"])
            break

    # ---- markers
    p = "results/markers/%s/filter_summary.txt" % S
    if have(p):
        t = open(p).read()
        m = re.search(r"Raw variants:\s*(\d+)", t)
        d["raw"] = int(m.group(1)) if m else None
        m = re.search(r"Markers after filter:\s*(\d+)", t)
        d["filtered"] = int(m.group(1)) if m else None
        m = re.search(r"min_dp=(\d+), max_dp=(\d+), min_qual=([\d.]+), ratio=\[([\d.]+),([\d.]+)\]", t)
        if m:
            d["fp"] = dict(min_dp=int(m.group(1)), max_dp=int(m.group(2)), min_qual=float(m.group(3)),
                           lo=float(m.group(4)), hi=float(m.group(5)))
    p = "results/blacklist/%s/rrna_exclude.bed" % S
    if have(p):
        b = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2], comment="#", dtype={0: str})
        d["rdna_bed"] = b
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
        d["star"] = {k.strip(): v.strip() for k, v in (l.split("|", 1) for l in open(p) if "|" in l)}
    for sub in ("GeneFull", "Gene"):
        p = "results/starsolo/%s/Solo.out/%s/Summary.csv" % (S, sub)
        if os.path.exists(p):
            USED.append(p)
            d["solo"] = [l.strip().split(",", 1) for l in open(p) if "," in l]
            d["solo_feature"] = sub
            break
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
        co["at_join"] = False
        for _, c, mb in JOINS.get(S, []):
            co.loc[(co.chrom == c) & ((co.mid - mb * 1e6).abs() <= 10e6), "at_join"] = True
        co["near_join"] = co.at_join & (S in ARTEFACT_JOINS)     # left out of the "without joins" versions
        if JOINS.get(S):
            # observed within 10 Mb of the joins, against the chromosome's own average for any 20 Mb
            d["at_join_cos"] = int(co.at_join.sum())
            d["at_join_exp"] = float(sum((co.chrom == c).sum() * 20e6 / d["L"][c]
                                         for _, c, _ in JOINS[S] if c in d["L"]))
        d["co"] = co
    p = "results/landscape/%s/coc_table.tsv" % S
    if have(p):
        d["coc"] = pd.read_csv(p, sep="\t")
    p = "results/landscape/%s/gamma_interference_summary.tsv" % S
    if have(p):
        gm = pd.read_csv(p, sep="\t")
        pooled = gm[gm.iloc[:, 0].astype(str) == "POOLED"]
        for k, v in (pooled.iloc[0].to_dict().items() if len(pooled) else []):
            put(sp, "gamma_pooled." + str(k), v)
    D[sp] = d
B, Pd = D["binata"], D["paradoxa"]


# ============================================================ window genotypes per cell
def window_genotypes(d):
    S, W = d["S"], int(WIN_MB * 1e6)
    nwin = {c: max(1, int(round(d["L"][c] / W))) for c in d["main"]}
    wins = [(c, i) for c in d["main"] for i in range(nwin[c])]
    widx = {w: k for k, w in enumerate(wins)}
    G = np.full((len(d["cells"]), len(wins)), np.nan)
    nmol = {}
    for ci, bc in enumerate(d["cells"]):
        p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
        if not os.path.exists(p):
            continue
        m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"],
                        dtype={0: str})
        m = m[(m.rc != m.ac) & m.chrom.isin(d["L"])].copy()
        m["alt"] = (m.ac > m.rc).astype(int)
        nmol[bc] = len(m)
        m["w"] = [min(int(q // W), nwin[c] - 1) for c, q in zip(m.chrom, m.pos)]
        g = m.groupby(["chrom", "w"]).alt.agg(["mean", "size"])
        for (c, w), r in g.iterrows():
            if r["size"] >= MIN_MOL and (r["mean"] >= 0.8 or r["mean"] <= 0.2):
                G[ci, widx[(c, w)]] = 1 if r["mean"] >= 0.8 else 0
    USED.append("results/cell_data_mol/%s/<barcode>.tsv (%d cells)" % (S, len(nmol)))
    d["mol_per_cell"] = pd.Series(nmol)
    M = (~np.isnan(G)).astype(float)
    X, Y = np.nan_to_num(G) * M, (1 - np.nan_to_num(G)) * M
    shared, agree = M @ M.T, X @ X.T + Y @ Y.T
    iu = np.triu_indices(len(d["cells"]), 1)
    sh, ag = shared[iu], agree[iu]
    keep = sh >= 20
    sim = np.full(len(sh), np.nan)
    sim[keep] = ag[keep] / sh[keep]
    dup = keep & (sim >= 0.95)
    d["pair_sim_median"] = float(np.nanmedian(sim)) if keep.any() else float("nan")
    d["dup_pairs"] = int(dup.sum())
    d["dup_cells"] = len(set(iu[0][dup]) | set(iu[1][dup]))
    # nuclei left if each group of near-identical nuclei is counted once (connected components)
    parent = list(range(len(d["cells"])))

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a, b in zip(iu[0][dup], iu[1][dup]):
        parent[root(a)] = root(b)
    d["distinct_genotypes"] = len({root(a) for a in range(len(d["cells"]))})
    d["called_share"] = float(M.mean())


for sp, d in D.items():
    guarded("window genotypes %s" % sp, lambda d=d: window_genotypes(d))


# ======================================================================= markers
def marker_stats(d):
    gm = d["good"]
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
    gaps, ing, tot, ncl = [], 0, 0, 0
    for c, g in gm.groupby("chrom"):
        pos = np.sort(g.pos.values)
        br = np.where(np.diff(pos) > HT_GAP)[0]
        cs, ce = np.concatenate([[pos[0]], pos[br + 1]]), np.concatenate([pos[br], [pos[-1]]])
        ncl += len(cs)
        gaps += list(np.concatenate([[cs[0]], cs[1:] - ce[:-1], [d["L"][c] - ce[-1]]]))
        if c in gi:
            st, en = gi[c]
            k = np.searchsorted(st, pos, side="right") - 1
            ing += int(((k >= 0) & (pos <= en[np.maximum(k, 0)])).sum())
        tot += len(pos)
    gaps = np.array(gaps, float)
    d.update(good_n=len(gm), good_per_mb=len(gm) / (d["G"] / 1e6), clusters=ncl,
             gap_median_kb=float(np.median(gaps)) / 1e3, gap_p95_kb=float(np.percentile(gaps, 95)) / 1e3,
             gap_max_mb=float(gaps.max()) / 1e6, gap_1mb_share=float(gaps[gaps > 1e6].sum()) / d["G"],
             good_in_genes=ing / tot if tot and gi else float("nan"))


def raw_variants(d):
    """A regular sample of the raw HiFi variant calls: what each marker filter removes."""
    p = "results/markers/%s/hifi_raw.vcf.gz" % d["S"]
    if not have(p) or "fp" not in d:
        return
    k = max(1, int((d.get("raw") or 3e5) // 300000))
    cmd = "zcat %s | awk -v k=%d '!/^#/ { if (++n %% k == 0) print }'" % (shlex.quote(p), k)
    out = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True).stdout
    rows = []
    for l in out.splitlines():
        f = l.split("\t")
        if len(f) < 10:
            continue
        fmt = dict(zip(f[8].split(":"), f[9].split(":")))
        snp = len(f[3]) == 1 and len(f[4]) == 1 and f[4] != "."
        het = fmt.get("GT", "") in ("0/1", "0|1", "1|0", "1/0")
        dp = int(fmt["DP"]) if fmt.get("DP", ".") not in (".", "") else np.nan
        ab = np.nan
        ad = fmt.get("AD", "")
        if "," in ad and ok(dp) and dp > 0 and ad.split(",")[1] != ".":
            ab = int(ad.split(",")[1]) / dp
        rows.append((f[0], int(f[1]), snp and het, float(f[5]) if f[5] != "." else np.nan, dp, ab))
    v = pd.DataFrame(rows, columns=["chrom", "pos", "cand", "qual", "dp", "ab"])
    fp = d["fp"]
    c = v[v.cand].copy()
    c["q_ok"] = c.qual >= fp["min_qual"]
    c["dp_low"], c["dp_high"] = c.dp < fp["min_dp"], c.dp > fp["max_dp"]
    c["ab_ok"] = (c.ab >= fp["lo"]) & (c.ab <= fp["hi"])
    c["pass"] = c.q_ok & ~c.dp_low & ~c.dp_high & c.ab_ok
    inr = np.zeros(len(c), bool)
    if "rdna_bed" in d:
        for ch, b in d["rdna_bed"].groupby(0):
            s, e = np.sort(b[1].values), b.sort_values(1)[2].values
            m = (c.chrom == ch).values
            kk = np.searchsorted(s, c.pos.values[m], side="right") - 1
            inr[m] = (kk >= 0) & (c.pos.values[m] <= e[np.maximum(kk, 0)])
    c["in_rdna"] = inr
    d["rawv"], d["raw_k"] = c, k
    n = len(v)
    d["raw_share_cand"] = len(c) / max(n, 1)
    for key in ("dp_low", "dp_high"):
        d["raw_share_" + key] = float(c[key].mean())
    d["raw_share_qual_fail"] = float((~c.q_ok).mean())
    d["raw_share_ab_fail"] = float((~c.ab_ok).mean())
    d["raw_share_pass"] = float(c["pass"].mean())
    d["rdna_markers_est"] = int((c["pass"] & c.in_rdna).sum() * k)


for sp, d in D.items():
    if "good" in d:
        guarded("marker statistics %s" % sp, lambda d=d: marker_stats(d))
    guarded("raw variants %s" % sp, lambda d=d: raw_variants(d))


# ===================================================================== crossovers
def co_stats(d):
    per, co = d["per"], d["co"]
    ncell = len(per)
    d.update(ncell_co=ncell, co_total=len(co), co_mean=per.n_cos.mean(), co_sd=per.n_cos.std(),
             co_median=per.n_cos.median())
    d["per_chrom"] = {c: (co.chrom == c).sum() / ncell for c in d["main"]}
    d["below_half"] = [c for c, v in d["per_chrom"].items() if v < 0.5]
    d["below_obligate"] = int((per.n_cos < d["n"] / 2).sum())
    d["width_median_kb"] = float((co.end - co.start).median()) / 1e3
    cidx = {b: k for k, b in enumerate(per.barcode)}
    rows = co.barcode.map(cidx)
    okr = rows.notna().values
    bix = np.minimum((co.u.values * BINS).astype(int), BINS - 1)

    def counts(mask):
        C = np.zeros((ncell, BINS))
        np.add.at(C, (rows[okr & mask].astype(int).values, bix[okr & mask]), 1)
        return C
    d["land"] = profile(counts(np.ones(len(co), bool)))
    d["shape"] = shape(d["land"]["U"], d["land"]["U_lo"], d["land"]["U_hi"])
    if co.near_join.any():
        d["join_cos"] = int(co.near_join.sum())
        d["land_nj"] = profile(counts(~co.near_join.values))
        d["shape_nj"] = shape(d["land_nj"]["U"], d["land_nj"]["U_lo"], d["land_nj"]["U_hi"])
        nj = co[co.near_join].groupby("barcode").size()
        d["per_nj"] = per.n_cos - per.barcode.map(nj).fillna(0).values
        d["co_mean_nj"] = float(d["per_nj"].mean())
    if "good" in d:
        mu = (d["good"].pos / d["good"].chrom.map(d["L"])).values
        h = np.histogram(mu, bins=BINS, range=(0, 1))[0].astype(float)
        d["mprof"], d["U_markers"] = h / h.mean(), uidx(h)
    if "genes" in d:
        g = d["genes"]
        h = np.histogram((((g.start + g.end) / 2) / g.chrom.map(d["L"])).values, bins=BINS, range=(0, 1))[0]
        d["gprof"], d["U_genes"] = h / h.mean(), uidx(h.astype(float))
    # 5 Mb windows: detectability and the busiest windows
    a, b, wl = [], [], []
    for c in d["main"]:
        k = int(math.ceil(d["L"][c] / 5e6))
        a += list(np.bincount(np.minimum((co.mid[co.chrom == c] // 5e6).astype(int), k - 1), minlength=k))
        if "good" in d:
            b += list(np.bincount(np.minimum((d["good"].pos[d["good"].chrom == c] // 5e6).astype(int), k - 1),
                                  minlength=k))
        wl += [(c, i) for i in range(k)]
    if b:
        d["rho_markers"] = pd.Series(a).corr(pd.Series(b), method="spearman")
    mu = float(np.mean(a))
    top = sorted(range(len(a)), key=lambda i: -a[i])[:6]
    d["busiest"] = []
    for i in top:
        c, w = wl[i]
        near = any(c == jc and abs((w + 0.5) * 5 - mb) <= 10 for _, jc, mb in JOINS.get(d["S"], []))
        d["busiest"].append([c, "%d–%d" % (w * 5, w * 5 + 5), a[i], ff(mu), "%.1e" % poisson_tail(a[i], mu),
                             "yes" if near else ""])
    d["window_mean"] = mu


for sp, d in D.items():
    if "per" in d and "co" in d:
        guarded("crossover statistics %s" % sp, lambda d=d: co_stats(d))
    d.setdefault("shape", "–")
HEAD = {sp: d.get("shape_nj", d["shape"]) for sp, d in D.items()}


# ========================================================================== rDNA
def rdna_depth(d):
    """HiFi reads over each rDNA array against a typical 20 kb window: collapsed copies."""
    gff = "results/blacklist/%s/rrna.gff" % d["S"]
    bam = "results/markers/%s/hifi.bam" % d["S"]
    sam = tool("samtools")
    if not (have(gff) and have(bam)) or not sam:
        d["rdna_note"] = "skipped: %s" % ("samtools not found" if not sam else "rrna.gff or hifi.bam missing")
        return
    g = pd.read_csv(gff, sep="\t", comment="#", header=None, usecols=[0, 3, 4, 8], names=["chrom", "start", "end", "att"],
                    dtype={0: str})
    g = g[g.chrom.isin(d["main"])].sort_values(["chrom", "start"])
    g["gene"] = g.att.str.extract(r"Name=([0-9_.]+S)", expand=False).fillna("?")
    arrays = []
    for c, x in g.groupby("chrom"):
        cur = None
        for _, r in x.iterrows():
            if cur and r.start - cur["end"] < 50000:
                cur["end"] = max(cur["end"], r.end)
                cur["genes"].append(r.gene)
            else:
                cur = {"chrom": c, "start": r.start, "end": r.end, "genes": [r.gene]}
                arrays.append(cur)

    def count(c, s, e):
        mid = (s + e) // 2
        s, e = min(s, mid - 10000), max(e, mid + 10000)
        out = subprocess.run([sam, "view", "-c", "-F", "0x904", bam, "%s:%d-%d" % (c, max(1, s), e)],
                             stdout=subprocess.PIPE, universal_newlines=True).stdout.strip()
        return int(out or 0), e - max(1, s)
    base = []
    for _ in range(120):
        c = d["main"][rng.integers(0, len(d["main"]))]
        s = int(rng.integers(20000, max(20001, d["L"][c] - 20000)))
        base.append(count(c, s, s + 20000)[0])
    b0 = max(float(np.median(base)), 1.0)
    rows = []
    for a in arrays:
        n, span = count(a["chrom"], a["start"], a["end"])
        typ = "45S" if any(x in ("18S", "28S", "5_8S", "5.8S") for x in a["genes"]) else "5S"
        fold = n / (b0 * span / 20000.0)
        near = ["%s (%.1f Mb away)" % (j, abs((a["start"] + a["end"]) / 2e6 - mb)) for j, jc, mb in JOINS.get(d["S"], [])
                if jc == a["chrom"] and abs((a["start"] + a["end"]) / 2e6 - mb) <= 5]
        gs = pd.Series(a["genes"]).value_counts()
        complete = "18S" in a["genes"] and "28S" in a["genes"]
        n5s = int(gs.get("5S", 0))
        trust = complete or fold >= 5 or n5s >= 50
        rows.append([a["chrom"], "%.2f" % (a["start"] / 1e6), typ, ", ".join("%s x%d" % (k, v) for k, v in gs.items()),
                     n, fold, "; ".join(near), complete, n5s, trust, a["start"], a["end"]])
    rows.sort(key=lambda r: -r[5])
    d["rdna_n_arrays"] = len(rows)
    d["rdna_n_trust"] = sum(r[9] for r in rows)
    d["rdna_rows"], d["rdna_base"] = rows, b0
    d["rdna_45s_n"] = sum(r[2] == "45S" for r in rows)
    d["rdna_45s_complete"] = sum(r[7] for r in rows)
    d["rdna_at_join"] = []
    for j, jc, mb in JOINS.get(d["S"], []):
        near = [r for r in rows if r[0] == jc and abs((r[10] + r[11]) / 2e6 - mb) <= 5]
        if near:
            r = max(near, key=lambda r: r[5])
            d["rdna_at_join"].append("%s: %s array at %.1f Mb, %s× the reads of a typical window%s" % (
                j, r[2], (r[10] + r[11]) / 2e6, fi(r[5]), ", a complete 45S unit" if r[7] else ""))
        else:
            d["rdna_at_join"].append("%s: no rDNA within 5 Mb" % j)


for sp, d in D.items():
    guarded("rDNA read depth %s" % sp, lambda d=d: rdna_depth(d))

# ===================================================== structure between the haplotypes
KINDS = ["collinear", "INVERTED", "TRANSLOCATED"]


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def covered(iv, a, b):
    return sum(max(0, min(e, b) - max(s0, a)) for s0, e in iv)


def structure(d):
    """hap2-on-hap1 synteny blocks (translocation_map.py) in this reference's coordinates: for each chromosome,
    where its other-haplotype copy aligns collinear, inverted, from another chromosome, or not at all."""
    sid = next((k for k, r in SAMPLES.items() if r["species"] == d["species"] and r["haplotype"] == "hap1"), None)
    p = "qc/translocations/%s/%s_synteny_blocks.tsv" % (sid, sid)
    if not have(p):
        return
    bt = pd.read_csv(p, sep="\t")
    bt["kind"] = bt["kind"].fillna("").replace("", "collinear")
    st = {}
    for c in d["main"]:
        x = bt[bt.hap1 == c] if c.endswith("_hap1") else bt[bt.hap2 == c]
        s0, e0 = (x.t_start, x.t_end) if c.endswith("_hap1") else (x.q_start, x.q_end)
        st[c] = {k: merge([(int(a), int(b)) for a, b, kk in zip(s0, e0, x["kind"]) if kk == k]) for k in KINDS}
        st[c]["any"] = merge([(int(a), int(b)) for a, b in zip(s0, e0)])
    d["struct"], d["struct_src"] = st, p
    for k in KINDS:
        put("paradoxa" if d is Pd else "binata", "struct_mb." + k,
            sum(covered(st[c][k], 0, d["L"][c]) for c in d["main"]) / 1e6)


def dead_spots(d, w=5e6):
    """Runs of 5 Mb windows without a crossover (join artefact excluded) and what they overlap. Only runs long
    enough that fewer than one is expected by chance (same crossovers per chromosome, placed uniformly) are kept."""
    co = d["co"][~d["co"].near_join]
    gm = d.get("good")

    def runs(cnt, L):
        out, i = [], 0
        while i < len(cnt):
            if cnt[i] == 0:
                j = i
                while j + 1 < len(cnt) and cnt[j + 1] == 0:
                    j += 1
                out.append((i * w, min((j + 1) * w, L)))
                i = j + 1
            else:
                i += 1
        return out
    sim = []
    for _ in range(200):
        lens = []
        for c in d["main"]:
            L = d["L"][c]
            nb = int(math.ceil(L / w))
            pos = rng.uniform(0, L, int((co.chrom == c).sum()))
            lens += [b - a for a, b in runs(np.bincount(np.minimum((pos // w).astype(int), nb - 1), minlength=nb), L)]
        sim.append(np.array(lens))
    thr = next((t for t in range(20, 205, 5) if np.mean([(x >= t * 1e6).sum() for x in sim]) < 1.0), 200)
    chance = np.array([(x >= thr * 1e6).sum() for x in sim])
    rows = []
    for c in d["main"]:
        L = d["L"][c]
        nb = int(math.ceil(L / w))
        cnt = np.bincount(np.minimum((co.mid[co.chrom == c] // w).astype(int), nb - 1), minlength=nb)
        mk = np.bincount(np.minimum((gm.pos[gm.chrom == c] // w).astype(int), nb - 1), minlength=nb) if gm is not None \
            else np.zeros(nb)
        med = max(float(np.median(mk)), 1.0)
        for a, b in runs(cnt, L):
            ln = b - a
            if ln < thr * 1e6:
                continue
            st = d.get("struct", {}).get(c)
            sh = {k: covered(st[k], a, b) / ln for k in KINDS} if st else {}
            none = 1 - covered(st["any"], a, b) / ln if st else float("nan")
            m = mk[int(a // w):int(math.ceil(b / w))].mean() / med
            join = any(jc == c and a <= mb * 1e6 <= b for _, jc, mb in JOINS.get(d["S"], []))
            cause = ("inverted in the other haplotype" if sh.get("INVERTED", 0) >= 0.5 else
                     "no homolog in the other haplotype" if ok(none) and none >= 0.5 else
                     "translocated segment" if sh.get("TRANSLOCATED", 0) >= 0.5 else
                     "few markers" if m < 0.25 else "not explained by structure or markers")
            rows.append([c, "%.0f–%.0f" % (a / 1e6, b / 1e6), ff(ln / 1e6, 0),
                         ff(100 * sh.get("INVERTED", float("nan")), 0), ff(100 * sh.get("TRANSLOCATED", float("nan")), 0),
                         ff(100 * none, 0), ff(m, 2), "yes" if join else "", cause])
    d["dead"], d["dead_thr"] = rows, thr
    d["dead_chance"] = (float(chance.mean()), np.percentile(chance, 2.5), np.percentile(chance, 97.5))
    d["dead_mb"] = sum(float(r[2]) for r in rows)


for sp, d in D.items():
    if sp == "paradoxa":
        guarded("structure %s" % sp, lambda d=d: structure(d))
    if "co" in d:
        guarded("dead spots %s" % sp, lambda d=d: dead_spots(d))

# =================================================================== numbers out
for sp, d in D.items():
    for k in ("raw", "filtered", "seen", "good_n", "good_per_mb", "clusters", "gap_median_kb", "gap_p95_kb",
              "gap_max_mb", "gap_1mb_share", "good_in_genes", "rdna_bp", "rdna_markers_est", "raw_k",
              "raw_share_cand", "raw_share_dp_low", "raw_share_dp_high", "raw_share_qual_fail", "raw_share_ab_fail",
              "raw_share_pass", "pair_sim_median", "dup_pairs", "dup_cells", "distinct_genotypes", "called_share",
              "ncell_co", "co_total", "co_mean", "co_sd", "co_median", "co_mean_nj", "below_obligate",
              "width_median_kb", "U_markers", "U_genes", "join_cos", "at_join_cos", "at_join_exp", "rho_markers", "window_mean", "shape", "shape_nj",
              "asm2n", "rdna_base", "rdna_45s_n", "rdna_45s_complete", "rdna_n_arrays", "rdna_n_trust", "dead_mb"):
        if k in d:
            put(sp, k, d[k])
    for lab in ("land", "land_nj"):
        if lab in d:
            for k in ("U", "U_lo", "U_hi"):
                put(sp, "%s.%s" % (lab, k), d[lab][k])
    for k, v in d.get("classes", {}).items():
        put(sp, "marker_class_" + k, v)
    if "bs" in d:
        put(sp, "barcodes_tested", len(d["bs"]))
    if "called" in d:
        put(sp, "cells_called", len(d["called"]))
    put(sp, "cells_selected", len(d["cells"]))
    for k, v in d.get("sel", {}).items():
        put(sp, "selection." + k, v)
    if "ht" in d and d["cells"]:
        h = d["ht"][d["ht"].barcode.isin(d["cells"])]
        for k, f in (("selected_molecules_median", h.molecules.median()),
                     ("selected_molecules_q25", h.molecules.quantile(0.25)),
                     ("selected_molecules_q75", h.molecules.quantile(0.75)),
                     ("selected_weakest_chrom_molecules_median", h.min_chrom_molecules.median()),
                     ("selected_haploidness_median", h.haploidness.median())):
            d[k] = put(sp, k, f)
    if "bs" in d and d["cells"]:
        b = d["bs"][d["bs"].barcode.isin(d["cells"])]
        d["umi_med"] = put(sp, "selected_umi_median", b.total_umi.median())
        d["genes_med"] = put(sp, "selected_genes_median", b.n_genes.median())
    for c, v in d.get("per_chrom", {}).items():
        put(sp, "cos_per_grain." + c, v)
    if "dead" in d:
        put(sp, "dead_spots", len(d["dead"]))
        put(sp, "dead_spots_chance", "%.1f (95%% %.0f–%.0f)" % d["dead_chance"])
        put(sp, "dead_spots_min_mb", d["dead_thr"])
        for r in d["dead"]:
            put(sp, "dead.%s:%s" % (r[0], r[1]), "%s Mb, inv %s%%, transl %s%%, no homolog %s%%, markers %s, %s"
                % (r[2], r[3], r[4], r[5], r[6], r[8]))
    for r in d.get("busiest", []):
        put(sp, "busiest_5mb.%s:%s" % (r[0], r[1]), "%d crossovers, p %s%s" % (r[2], r[4], ", near join" if r[5] else ""))
    for r in [r for r in d.get("rdna_rows", []) if r[9] or r[6]][:12]:
        put(sp, "rdna_fold.%s:%s" % (r[0], r[1]), "%s %s reads %d fold %.0f %s" % (r[2], r[3], r[4], r[5], r[6]))

# ======================================================================== figures


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return p


FIGS = {}


def fig_markers():
    fig, axs = plt.subplots(2, 1, figsize=(10, 4.2))
    for ax, (sp, d) in zip(axs, D.items()):
        off = 0
        for i, c in enumerate(d["main"]):
            pos = d["good"].pos[d["good"].chrom == c].values
            k = int(math.ceil(d["L"][c] / 2e6))
            h = np.bincount(np.minimum((pos // 2e6).astype(int), k - 1), minlength=k) / 2.0
            if i % 2:
                ax.axvspan(off / 1e6, (off + d["L"][c]) / 1e6, color="#F1EFE8", lw=0)
            ax.plot(off / 1e6 + (np.arange(k) + 0.5) * 2, h, color=COL[sp], lw=0.8)
            ax.text((off + d["L"][c] / 2) / 1e6, -0.08, re.sub(r"_hap\d", "", c).replace("chr", ""),
                    transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7)
            off += d["L"][c]
        ax.set_xlim(0, off / 1e6)
        ax.set_xticks([])
        ax.set_ylabel("good markers / Mb")
        ax.set_title("%s (%s): %s good markers" % (NAME[sp], d["S"], fi(d.get("good_n"))), loc="left", fontsize=9)
    fig.tight_layout()
    return save(fig, "fig_markers_along_genome.png")


def fig_raw():
    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    for j, (sp, d) in enumerate(D.items()):
        if "rawv" not in d:
            continue
        c, fp = d["rawv"], d["fp"]
        ax = axs[j, 0]
        top = max(fp["max_dp"] * 2.5, float(np.nanpercentile(c.dp, 99)))
        ax.hist(c.dp.clip(upper=top), bins=np.arange(0, top + 2, 2), color=COL[sp], alpha=0.85)
        for x in (fp["min_dp"], fp["max_dp"]):
            ax.axvline(x, color="#444441", ls="--", lw=1)
        ax.axvspan(fp["min_dp"], fp["max_dp"], color="#F1EFE8", zorder=0)
        ax.set_xlabel("HiFi read depth at the SNP")
        ax.set_ylabel("SNPs (1 in %d sampled)" % d["raw_k"])
        ax.set_title("%s: depth kept %d–%d (too low %s%%, too high %s%%)" % (
            NAME[sp], fp["min_dp"], fp["max_dp"], ff(100 * d["raw_share_dp_low"]), ff(100 * d["raw_share_dp_high"])),
            loc="left", fontsize=8.5)
        ax = axs[j, 1]
        ax.hist(c.ab.dropna(), bins=np.linspace(0, 1, 51), color=COL[sp], alpha=0.85)
        for x in (fp["lo"], fp["hi"]):
            ax.axvline(x, color="#444441", ls="--", lw=1)
        ax.axvspan(fp["lo"], fp["hi"], color="#F1EFE8", zorder=0)
        ax.set_xlabel("ALT allele share of the HiFi reads")
        ax.set_title("%s: allele balance kept %s–%s (outside: %s%%)" % (
            NAME[sp], fp["lo"], fp["hi"], ff(100 * d["raw_share_ab_fail"])), loc="left", fontsize=8.5)
    fig.suptitle("Heterozygous biallelic SNPs in the raw HiFi calls (a regular sample), with the marker filters",
                 x=0.01, ha="left", fontsize=9)
    fig.tight_layout()
    return save(fig, "fig_raw_variants.png")


def fig_cells():
    fig, axs = plt.subplots(2, 2, figsize=(10, 7.2))
    for j, (sp, d) in enumerate(D.items()):
        ax = axs[0, j]
        if "bs" in d:
            u = np.sort(d["bs"].total_umi.values)[::-1]
            ax.loglog(np.arange(1, len(u) + 1), u, color="#888780", lw=1)
            if "called" in d:
                cu = d["bs"].set_index("barcode").total_umi.reindex(d["called"]).dropna()
                ax.axhline(cu.min(), color=COL[sp], lw=0.8, ls="--")
                ax.text(1.5, cu.min() * 1.25, "called: %s barcodes (≥ %s UMIs)" % (fi(len(d["called"])), fi(cu.min())),
                        color=COL[sp], fontsize=8)
                ax.text(0.98, 0.97, "cut-off kept generous: it removes only empty\ndroplets and barcode errors; "
                        "single nuclei\nare chosen below, by haploidness", transform=ax.transAxes, ha="right",
                        va="top", fontsize=7, color="#5F5E5A")
        ax.set_xlabel("barcode rank")
        ax.set_ylabel("UMIs per barcode")
        ax.set_title("%s: step 1, call barcodes (knee plot)" % NAME[sp], loc="left")
        ax = axs[1, j]
        if "ht" in d:
            h = d["ht"].dropna(subset=["haploidness"])
            sel = h.barcode.isin(d["cells"])
            ax.scatter(h.molecules[~sel], h.haploidness[~sel], s=4, color="#B4B2A9", alpha=0.6, lw=0,
                       label="called, not selected (%s)" % fi((~sel).sum()))
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
        ax.set_xlabel("informative molecules per called barcode")
        ax.set_ylabel("haploidness (1 = one clean haplotype)")
        ax.set_title("%s: step 2, called barcodes only: single haploid nuclei?" % NAME[sp], loc="left", fontsize=8.5)
    fig.tight_layout()
    return save(fig, "fig_cells.png")


def fig_per_grain():
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
    for ax, (sp, d) in zip(axs, D.items()):
        v = d["per"].n_cos.values
        bins = np.arange(-0.5, v.max() + 1.5, 1)
        ax.hist(v, bins=bins, color=COL[sp], alpha=0.85, label="as called")
        if "per_nj" in d:
            ax.hist(d["per_nj"], bins=bins, histtype="step", color="#D85A30", lw=1.5,
                    label="without crossovers at the joins")
            ax.legend(frameon=False, fontsize=7, loc="upper right")
        ax.axvline(d["n"] / 2, color="#444441", ls="--", lw=1)
        ax.text(d["n"] / 2, ax.get_ylim()[1] * 0.97, " obligate minimum %s" % ff(d["n"] / 2), fontsize=7, va="top")
        ax.set_xlabel("crossovers per pollen grain")
        ax.set_ylabel("grains")
        ax.set_title("%s: %s grains, mean %s" % (NAME[sp], fi(len(v)), ff(v.mean(), 2)), loc="left")
    fig.tight_layout()
    return save(fig, "fig_crossovers_per_grain.png")


def fig_landscape():
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    x = (np.arange(BINS) + 0.5) / BINS * 100
    for ax, (sp, d) in zip(axs, D.items()):
        L = d["land"]
        ax.fill_between(x, L["lo"], L["hi"], color=COL[sp], alpha=0.18, lw=0)
        ax.plot(x, L["prof"], color=COL[sp], lw=2, label="crossovers (95% CI over grains)")
        if "land_nj" in d:
            ax.plot(x, d["land_nj"]["prof"], color=COL[sp], lw=1.2, ls=":", label="without the join crossovers")
        if "mprof" in d:
            ax.plot(x, d["mprof"], color="#888780", lw=1.2, label="good markers")
        if "gprof" in d:
            ax.plot(x, d["gprof"], color="#888780", lw=1, ls="--", label="genes")
        ax.axhline(1, color="#D3D1C7", lw=0.8)
        ax.set_xlabel("position along chromosome (% of length, all chromosomes pooled)")
        ax.set_title(NAME[sp], loc="left")
        ax.text(0.02, 0.97, HEAD[sp], transform=ax.transAxes, va="top", fontsize=7.5, color=COL[sp])
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    axs[0].set_ylabel("density relative to the mean")
    fig.tight_layout()
    return save(fig, "fig_landscape_relative.png")


def fig_per_chrom():
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4), gridspec_kw={"width_ratios": [len(D[s]["main"]) for s in D]})
    for ax, (sp, d) in zip(axs, D.items()):
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
    return save(fig, "fig_crossovers_per_chromosome.png")


def interference(d):
    """Adjacent crossovers on one chromosome of one nucleus: their distance as a share of the chromosome,
    against the same numbers of crossovers placed independently (no interference)."""
    co = d["co"][~d["co"].near_join]
    obs, groups = [], []
    for (bc, c), x in co.groupby(["barcode", "chrom"]):
        if len(x) >= 2:
            obs += list(np.diff(np.sort(x.u.values)))
            groups.append((c, len(x)))
    pool = {c: co.u[co.chrom == c].values for c in d["main"]}
    edges = np.linspace(0, 1, 11)
    o = np.histogram(obs, bins=edges)[0].astype(float)
    nulls = []
    for _ in range(300):
        dd = []
        for c, k in groups:
            dd += list(np.diff(np.sort(rng.choice(pool[c], k, replace=True))))
        nulls.append(np.histogram(dd, bins=edges)[0])
    nulls = np.array(nulls, float)
    se = nulls[:, :2].sum(1)
    d["intf"] = dict(edges=edges, obs=o, exp=nulls.mean(0), lo=np.percentile(nulls, 2.5, 0),
                     hi=np.percentile(nulls, 97.5, 0), pairs=len(obs), short_obs=o[:2].sum(),
                     short_exp=se.mean(), ratio=o[:2].sum() / max(se.mean(), 1e-9),
                     ratio_lo=o[:2].sum() / max(np.percentile(se, 97.5), 1e-9),
                     ratio_hi=o[:2].sum() / max(np.percentile(se, 2.5), 1e-9))


def fig_interference():
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.4))
    for ax, (sp, d) in zip(axs, D.items()):
        f = d.get("intf")
        if not f:
            continue
        x = (f["edges"][:-1] + f["edges"][1:]) / 2 * 100
        ax.bar(x, f["obs"], width=9, color=COL[sp], alpha=0.8, label="observed")
        ax.fill_between(x, f["lo"], f["hi"], color="#888780", alpha=0.25, lw=0)
        ax.plot(x, f["exp"], color="#444441", lw=1.5, marker="o", ms=3, label="expected without interference")
        ax.set_xlabel("distance between adjacent crossovers, % of the chromosome")
        ax.set_ylabel("pairs of adjacent crossovers")
        ax.set_title("%s: %d pairs; closer than 20%%: %s× expected" % (NAME[sp], f["pairs"], ff(f["ratio"], 2)),
                     loc="left", fontsize=9)
        ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    return save(fig, "fig_interference.png")


def fig_diagram():
    pairs = PAIRS
    rr = {k: v for k, v in pairs}
    L2n = Pd.get("L2n", {})
    jpos = {"chr1_hap1": 262.9, "chr2_hap2": 215.0, "chr1_hap2": 262.0, "chr2_hap1": 211.0}
    chrs = [("A", "chr1_hap1", "L1", "P"), ("C", "chr1_hap2", "L1", "Q"),
            ("B", "chr2_hap1", "L2", "P"), ("D", "chr2_hap2", "L2", "Q")]
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.6), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axs[0]
    scale = 1 / 480.0
    for i, (nm, c, left, right) in enumerate(chrs):
        y = 3.2 - i * 0.95
        tot = L2n.get(c, 0) / 1e6 or (jpos[c] + 120)
        j = jpos[c]
        for x0, w, arm in ((0, j, left), (j, tot - j, right)):
            ax.add_patch(FancyBboxPatch((0.18 + x0 * scale, y), w * scale - 0.005, 0.42,
                                        boxstyle="round,pad=0,rounding_size=0.06", fc=ARM[arm], ec="none"))
            ax.text(0.18 + (x0 + w / 2) * scale, y + 0.21, arm, ha="center", va="center", fontsize=10,
                    fontweight="bold", color="#2C2C2A")
        ax.text(0.0, y + 0.21, "%s  %s" % (nm, c), ha="left", va="center", fontsize=8)
        ax.text(0.18 + tot * scale + 0.02, y + 0.21, "%.0f Mb" % tot, va="center", fontsize=7, color="#5F5E5A")
    ax.set_xlim(-0.05, 1.35)
    ax.set_ylim(0.1, 3.85)
    ax.axis("off")
    ax.set_title("In the plant's tissue (assembly + Hi-C): four chromosomes", loc="left", fontsize=9)
    ax = axs[1]
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("In the pollen: which arms are inherited together", loc="left", fontsize=9)

    def r_of(prefix):
        for k, v in rr.items():
            if k.startswith(prefix):
                return v
        return None
    lines = [("L1", "Q", r_of("L1 end    vs Q start")), ("L2", "P", r_of("L2 end    vs P start")),
             ("L1", "P", r_of("L1 end    vs P start")), ("L2", "Q", r_of("L2 end    vs Q start"))]
    for i, (a, b, r) in enumerate(lines):
        y = 0.82 - i * 0.2
        for k, arm in enumerate((a, b)):
            ax.add_patch(FancyBboxPatch((0.03 + k * 0.13, y - 0.06), 0.11, 0.12,
                                        boxstyle="round,pad=0,rounding_size=0.02", fc=ARM[arm], ec="none"))
            ax.text(0.085 + k * 0.13, y, arm, ha="center", va="center", fontsize=10, fontweight="bold")
        if r:
            together = float(r[5]) < 0.15
            ax.text(0.31, y, "r = %s over %s nuclei: %s" % (r[5], r[4], "inherited together" if together
                                                            else "inherited independently"),
                    va="center", fontsize=8.5, color="#2C2C2A" if together else "#A32D2D")
    ax.text(0.03, 0.02, "r = share of nuclei whose two windows disagree (0 together, 0.5 independent).\n"
            "The pollen behave as if L1 were joined to Q and L2 to P, in both homologs.",
            fontsize=7.5, color="#5F5E5A", va="bottom")
    fig.tight_layout()
    return save(fig, "fig_arm_diagram.png")


def hic_agg(M, chroms, bin_bp, lens, order, agg_bp):
    """Sum the bin-level contact matrix into agg_bp bins over the chromosomes in `order`."""
    f = max(1, int(round(agg_bp / bin_bp)))
    nb = {c: -(-lens[c] // bin_bp) for c in chroms}
    off, o = {}, 0
    for c in chroms:
        off[c], o = o, o + nb[c]
    sel, grp, nper, k = [], [], [], 0
    for c in order:
        na = -(-nb[c] // f)
        nper.append(na)
        for i in range(na):
            b0, b1 = i * f, min((i + 1) * f, nb[c])
            sel += list(range(off[c] + b0, off[c] + b1))
            grp += [k] * (b1 - b0)
            k += 1
    sel, grp = np.array(sel), np.array(grp)
    starts = np.r_[0, np.where(np.diff(grp) != 0)[0] + 1]
    Ms = M[np.ix_(sel, sel)].astype(np.float32)
    return np.add.reduceat(np.add.reduceat(Ms, starts, axis=0), starts, axis=1), nper, f * bin_bp


def draw_hic(ax, H, nper, order, agg_bp, short=True):
    """Contacts on a log scale, white to red: white = the typical contact between two different chromosomes,
    red = the diagonal. A chromosome is a red square; a misjoin would break a square or put red off the diagonal."""
    lab = np.repeat(np.arange(len(nper)), nper)
    inter = lab[:, None] != lab[None, :]
    X = np.log10(H.astype(float) + 1)
    lo = np.percentile(X[inter], 75) if inter.any() else np.percentile(X, 50)
    hi = np.percentile(X[~inter], 99.5)
    tot = sum(nper) * agg_bp / 1e6
    im = ax.imshow(X, cmap="Reds", vmin=lo, vmax=max(hi, lo + 1e-6), interpolation="nearest",
                   extent=(0, tot, tot, 0))
    edges = np.cumsum([0] + list(nper)) * agg_bp / 1e6
    for e in edges[1:-1]:
        ax.axhline(e, color="#888780", lw=0.4)
        ax.axvline(e, color="#888780", lw=0.4)
    names = [(c.replace("_hap", " h") if short else c) for c in order]
    ax.set_xticks((edges[:-1] + edges[1:]) / 2)
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_yticks((edges[:-1] + edges[1:]) / 2)
    ax.set_yticklabels(names, fontsize=6)
    arm_strips(ax, order, edges, tot)
    return im


def arm_strips(ax, order, edges, tot, line="#2C2C2A"):
    """Arm colours of the chr1/chr2 translocation along the top and left of a chromosome-by-chromosome map."""
    band = tot * 0.018
    for i, c in enumerate(order):
        spans = [(0, ARMSPAN[c][0], ARMSPAN[c][1]), (ARMSPAN[c][0], edges[i + 1] - edges[i], ARMSPAN[c][2])] \
            if c in ARMSPAN else [(0, edges[i + 1] - edges[i], None)]
        for a, b, arm in spans:
            col = ARM[arm] if arm else "#E9E7E0"
            ax.add_patch(plt.Rectangle((edges[i] + a, -2.2 * band), b - a, band, color=col, clip_on=False, lw=0))
            ax.add_patch(plt.Rectangle((-2.2 * band, edges[i] + a), band, b - a, color=col, clip_on=False, lw=0))
            if arm and b - a > tot * 0.03:
                ax.text(edges[i] + (a + b) / 2, -1.7 * band, arm, ha="center", va="center", fontsize=6)
        if c in ARMSPAN:
            x = edges[i] + ARMSPAN[c][0]
            ax.plot([x, x], [0, tot], color=line, lw=0.6, ls=":")
            ax.plot([0, tot], [x, x], color=line, lw=0.6, ls=":")
    ax.set_xlim(-2.4 * band, tot)
    ax.set_ylim(tot, -2.4 * band)
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)


def hic_load(sp, name):
    d = D[sp]
    p = os.path.join(A.phd_root, "results", d["species"], "qc/hic_remap", name)
    if not (have(p) and "L2n" in d):
        return None
    z = np.load(p)
    return z["M"], [str(x) for x in z["chroms"]], int(z["bin"]), p


def fig_hic_dual(sp):
    """Hi-C on both haplotypes, all reads with MAPQ >= 1: reads that fit either copy of a sequence put contacts
    between the copies, so homologous sequence shows as a red line off the diagonal."""
    d = D[sp]
    got = hic_load(sp, "contacts_all.npz") or hic_load(sp, "contacts_unique.npz")
    if not got:
        return None
    M, chroms, bin_bp, p = got
    order = sorted([c for c in chroms if re.match(r"^chr\d+_hap[12]$", c)],
                   key=lambda c: (int(re.search(r"\d+", c).group()), c[-1]))
    focus = [c for c in ("chr1_hap1", "chr1_hap2", "chr2_hap1", "chr2_hap2") if c in chroms]
    H1, n1, b1 = hic_agg(M, chroms, bin_bp, d["L2n"], order, 2e6)
    H2, n2, b2 = hic_agg(M, chroms, bin_bp, d["L2n"], focus, 1e6)
    del M
    fig, axs = plt.subplots(1, 2, figsize=(15, 7.6), gridspec_kw={"width_ratios": [1.05, 1]})
    im = draw_hic(axs[0], H1, n1, order, b1)
    axs[0].set_title("all twelve chromosomes, hap1 next to hap2 (2 Mb bins)", loc="left", fontsize=9, pad=16)
    draw_hic(axs[1], H2, n2, focus, b2, short=False)
    axs[1].set_title("chr1 and chr2 only (1 Mb bins): A, C, B, D", loc="left", fontsize=9, pad=16)
    cb = fig.colorbar(im, ax=axs, fraction=0.02, pad=0.01)
    cb.set_label("log10 contacts (white = typical between two chromosomes)", fontsize=7)
    fig.suptitle("%s Hi-C on both haplotypes (%s). Arm colours as in the diagram; dotted = joins" % (
        NAME[sp], os.path.basename(p).replace(".npz", "")), x=0.01, ha="left", fontsize=9)
    return save(fig, "fig_hic_both_haplotypes_%s.png" % sp)


def fig_link_hic(sp):
    d = D[sp]
    S = d["S"]
    rp, wp = "qc/linkage/%s/r_matrix.npy" % S, "qc/linkage/%s/windows.tsv" % S
    if not (have(rp) and have(wp)):
        return None
    R = np.load(rp)
    win = pd.read_csv(wp, sep="\t")
    idx = [i for c in d["main"] for i in win.index[win.chrom == c]]
    Rm = R[np.ix_(idx, idx)]
    wsub = win.loc[idx].reset_index(drop=True)
    got = hic_load(sp, "contacts_unique.npz")
    H = None
    if got:
        M, chroms, bin_bp, p = got
        if all(d["L2n"].get(c) == d["L"][c] for c in d["main"]) and all(c in chroms for c in d["main"]):
            H, nper, aggb = hic_agg(M, chroms, bin_bp, d["L2n"], d["main"], 2e6)
        else:
            PROBLEMS.append(("Hi-C map %s" % sp, "reference chromosomes not in the 2n assembly with the same lengths"))
        del M
    ncol = 2 if H is not None else 1
    fig, axs = plt.subplots(1, ncol, figsize=(6.6 * ncol, 6.6))
    axs = np.atleast_1d(axs)
    nwin = [int((wsub.chrom == c).sum()) for c in d["main"]]
    tot = sum(nwin) * WIN_MB
    im = axs[0].imshow(np.clip(1 - 2 * Rm, 0, 1), cmap="Greens" if sp == "binata" else "Purples", vmin=0, vmax=1,
                       interpolation="nearest", extent=(0, tot, tot, 0))
    edges = np.cumsum([0] + nwin) * WIN_MB
    for e in edges[1:-1]:
        axs[0].axhline(e, color="#888780", lw=0.4)
        axs[0].axvline(e, color="#888780", lw=0.4)
    axs[0].set_xticks((edges[:-1] + edges[1:]) / 2)
    axs[0].set_xticklabels(d["main"], rotation=90, fontsize=6)
    axs[0].set_yticks((edges[:-1] + edges[1:]) / 2)
    axs[0].set_yticklabels(d["main"], fontsize=6)
    arm_strips(axs[0], d["main"], edges, tot)
    axs[0].set_title("pollen linkage, %g Mb windows (1 − 2r; dark = inherited together)" % WIN_MB, loc="left",
                     fontsize=9, pad=16)
    if H is not None:
        draw_hic(axs[1], H, nper, d["main"], aggb, short=False)
        axs[1].set_title("Hi-C on this reference, MAPQ ≥ 30 (2 Mb bins; log, white = between chromosomes)",
                         loc="left", fontsize=9, pad=16)
    fig.suptitle("%s (%s): every chromosome%s" % (NAME[sp], S, "; arm colours as in the diagram, dotted = chr1/chr2 "
                 "joins" if S in JOINS else ""), x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return save(fig, "fig_linkage_hic_%s.png" % sp)


def fig_chrom_landscapes():
    ncol = 4
    nrows = [int(math.ceil(len(D[s]["main"]) / ncol)) for s in D]
    fig, axs = plt.subplots(sum(nrows), ncol, figsize=(12, 1.9 * sum(nrows)))
    axs = np.atleast_2d(axs)
    r0 = 0
    for (sp, d), nr in zip(D.items(), nrows):
        for k, c in enumerate(d["main"]):
            ax = axs[r0 + k // ncol, k % ncol]
            nb = int(math.ceil(d["L"][c] / 5e6))
            sub = d["co"][d["co"].chrom == c]
            cnt = np.bincount(np.minimum((sub.mid // 5e6).astype(int), nb - 1), minlength=nb)
            ax.bar((np.arange(nb) + 0.5) * 5, 100.0 * cnt / d["ncell_co"] / 5, width=5, color=COL[sp], alpha=0.85)
            pos = d["good"].pos[d["good"].chrom == c].values
            mk = np.bincount(np.minimum((pos // 5e6).astype(int), nb - 1), minlength=nb)
            ax2 = ax.twinx()
            ax2.plot((np.arange(nb) + 0.5) * 5, mk / 5, color="#888780", lw=0.8)
            ax2.set_yticks([])
            ax2.spines["right"].set_visible(False)
            for _, jc, mb in JOINS.get(d["S"], []):
                if jc == c:
                    ax.axvline(mb, color="#D85A30", lw=0.8, ls=":")
            ax.set_title("%s %s" % (NAME[sp].split()[1], c), fontsize=7, loc="left")
            ax.tick_params(labelsize=6)
        for j in range(len(d["main"]), nr * ncol):
            axs[r0 + j // ncol, j % ncol].axis("off")
        r0 += nr
    fig.text(0.0, 0.5, "cM/Mb (bars); good markers per Mb (grey, own scale); orange = translocation join",
             rotation=90, va="center", fontsize=8)
    fig.tight_layout()
    return save(fig, "fig_per_chromosome_landscapes.png")


def fig_rdna(sp):
    d = D[sp]
    rows = [r for r in d.get("rdna_rows", []) if r[9]]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(11, 0.75 * len(d["main"]) + 1.2))
    co = d.get("co")
    big = max(d["L"].values()) / 1e6
    for i, c in enumerate(d["main"]):
        y = len(d["main"]) - 1 - i
        ax.add_patch(plt.Rectangle((0, y - 0.12), d["L"][c] / 1e6, 0.24, color="#D3D1C7", lw=0))
        ax.text(-big * 0.01, y, c, ha="right", va="center", fontsize=8)
        if co is not None:
            nb = int(math.ceil(d["L"][c] / 5e6))
            cnt = np.bincount(np.minimum((co.mid[co.chrom == c] // 5e6).astype(int), nb - 1), minlength=nb)
            top = max(1.0, float(np.bincount(np.minimum((co.mid // 5e6).astype(int), 10 ** 6)).max()))
            ax.plot((np.arange(nb) + 0.5) * 5, y + 0.16 + 0.3 * cnt / top, color="#1D9E75", lw=0.9)
        for _, jc, mb in JOINS.get(d["S"], []):
            if jc == c:
                ax.plot([mb], [y - 0.3], marker="^", color="#2C2C2A", ms=6)
        for r in rows:
            if r[0] != c:
                continue
            x = (r[10] + r[11]) / 2e6
            sz = 4 + 4 * math.log10(max(r[5], 1))
            ax.plot([x], [y], marker="s" if r[2] == "45S" else "o", ms=sz, mec="none",
                    color="#D85A30" if r[2] == "45S" else "#7F77DD", alpha=0.9)
            if r[5] >= 5:
                ax.text(x, y + 0.2, "%s×" % fi(r[5]), ha="center", va="bottom", fontsize=6.5)
    ax.set_xlim(-big * 0.12, big * 1.02)
    ax.set_ylim(-0.6, len(d["main"]) - 0.3)
    ax.set_yticks([])
    ax.set_xlabel("Mb")
    ax.set_title("%s (%s): the %d of %d rDNA arrays worth trusting\ncomplete 45S units (squares), 5S arrays "
                 "(circles); size and label = HiFi reads × a typical window; triangles = joins; green = crossovers"
                 % (NAME[sp], d["S"], len(rows), d["rdna_n_arrays"]), loc="left", fontsize=8.5)
    for k in ("top", "right", "left"):
        ax.spines[k].set_visible(False)
    fig.tight_layout()
    return save(fig, "fig_rdna_%s.png" % sp)


def fig_structure(sp):
    d = D[sp]
    if "struct" not in d:
        return None
    n = len(d["main"])
    fig = plt.figure(figsize=(11, 1.9 * n + 1.0))
    gs = fig.add_gridspec(3 * n, 1, height_ratios=[3.0, 0.3, 0.42] * n, hspace=0.16)
    big = max(d["L"].values()) / 1e6
    colk = {"collinear": "#D3D1C7", "INVERTED": "#A32D2D", "TRANSLOCATED": "#185FA5"}
    co = d["co"]
    for i, c in enumerate(d["main"]):
        ax, axa, axs_ = fig.add_subplot(gs[3 * i]), fig.add_subplot(gs[3 * i + 1]), fig.add_subplot(gs[3 * i + 2])
        L = d["L"][c] / 1e6
        nb = int(math.ceil(d["L"][c] / 5e6))
        cnt = np.bincount(np.minimum((co.mid[co.chrom == c] // 5e6).astype(int), nb - 1), minlength=nb)
        y = 100.0 * cnt / d["ncell_co"] / 5
        ax.bar((np.arange(nb) + 0.5) * 5, y, width=5, color=COL[sp], alpha=0.85)
        pos = d["good"].pos[d["good"].chrom == c].values
        mk = np.bincount(np.minimum((pos // 5e6).astype(int), nb - 1), minlength=nb)
        a3 = ax.twinx()
        a3.plot((np.arange(nb) + 0.5) * 5, mk / 5, color="#888780", lw=0.8)
        a3.set_yticks([])
        a3.spines["right"].set_visible(False)
        top = max(float(y.max()), 0.1) * 1.18
        ax.set_ylim(0, top)
        for r in d.get("dead", []):
            if r[0] == c:
                a, b = [float(x) for x in r[1].split("–")]
                ax.plot([a, b], [top * 0.94] * 2, color="#2C2C2A", lw=3, solid_capstyle="butt")
        if c in ARMSPAN:
            for a_ in (ax, axa, axs_):
                a_.axvline(ARMSPAN[c][0], color="#2C2C2A", lw=0.9, ls=":")
        ax.set_xlim(0, big)
        ax.set_ylabel("cM/Mb", fontsize=6)
        ax.tick_params(labelsize=6, labelbottom=False)
        ax.set_title("%s%s" % (c, " = %s (%s·%s)" % (ARMNAME[c], ARMSPAN[c][1], ARMSPAN[c][2]) if c in ARMSPAN else ""),
                     loc="left", fontsize=8, pad=2)
        spans = [(0, ARMSPAN[c][0], ARMSPAN[c][1]), (ARMSPAN[c][0], L, ARMSPAN[c][2])] if c in ARMSPAN else [(0, L, None)]
        for a, b, arm in spans:
            axa.add_patch(plt.Rectangle((a, 0), b - a, 1, color=ARM[arm] if arm else "#F1EFE8", lw=0))
            if arm:
                axa.text((a + b) / 2, 0.5, arm, ha="center", va="center", fontsize=7, fontweight="bold")
        st = d["struct"][c]
        axs_.add_patch(plt.Rectangle((0, 0), L, 1, color="white", ec="#B4B2A9", lw=0.5))
        for k in KINDS:
            for a, b in st[k]:
                axs_.add_patch(plt.Rectangle((a / 1e6, 0), (b - a) / 1e6, 1, color=colk[k], lw=0))
        for a_ in (axa, axs_):
            a_.set_xlim(0, big)
            a_.set_ylim(0, 1)
            a_.set_yticks([])
            for k in ("top", "right", "left"):
                a_.spines[k].set_visible(False)
        axa.tick_params(labelbottom=False, bottom=False)
        axa.spines["bottom"].set_visible(False)
        axs_.tick_params(labelsize=6, labelbottom=(i == n - 1))
    handles = [plt.Rectangle((0, 0), 1, 1, color=ARM[a]) for a in ("L1", "L2", "P", "Q")] + \
        [plt.Rectangle((0, 0), 1, 1, color=colk[k]) for k in KINDS] + \
        [plt.Rectangle((0, 0), 1, 1, fc="white", ec="#B4B2A9"), plt.Line2D([0], [0], color="#2C2C2A", lw=3)]
    fig.legend(handles, ["arm L1", "arm L2", "arm P", "arm Q", "homolog collinear", "homolog inverted",
                         "homolog from another chromosome", "no homolog alignment",
                         "no crossover over >= %d Mb" % d.get("dead_thr", 20)],
               loc="lower center", ncol=5, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.035, "Mb. Per chromosome: crossovers (bars, cM/Mb; grey line = good markers), the arms of the "
             "chr1/chr2 translocation, and where the other haplotype's copy aligns (hap2 on hap1)", ha="center",
             fontsize=7)
    fig.subplots_adjust(top=0.985, bottom=0.085, left=0.07, right=0.98)
    return save(fig, "fig_structure_landscape_%s.png" % sp)


def render_pdf(pdf, png_prefix, dpi=75):
    if shutil.which("pdftoppm"):
        subprocess.run(["pdftoppm", "-r", str(dpi), "-png", "-singlefile", pdf, png_prefix], check=True)
    elif shutil.which("gs"):
        subprocess.run(["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=png16m", "-r%d" % dpi, "-dFirstPage=1",
                        "-dLastPage=1", "-sOutputFile=%s.png" % png_prefix, pdf], check=True)
    else:
        raise RuntimeError("neither pdftoppm nor gs on PATH")
    return png_prefix + ".png"


def example_cells():
    out = {}
    for sp, d in D.items():
        h = d["ht"][d["ht"].barcode.isin(d["cells"])]
        bc = h.iloc[(h.molecules - h.molecules.median()).abs().argsort().iloc[0]].barcode
        pdf = "results/crossovers/%s/per_cell/%s_co.pdf" % (d["S"], bc)
        if have(pdf):
            out[sp] = (bc, render_pdf(pdf, os.path.join(FIG, "example_cell_%s" % sp)))
    return out


def cell_browser(sp):
    d = D[sp]
    per = d["per"].set_index("barcode")
    ht = d["ht"].set_index("barcode") if "ht" in d else pd.DataFrame()
    order = per.n_cos.sort_values(ascending=False, kind="mergesort").index.tolist()
    pdfs = [(bc, "results/crossovers/%s/per_cell/%s_co.pdf" % (d["S"], bc)) for bc in order]
    pdfs = [(bc, p) for bc, p in pdfs if os.path.exists(p)]
    outp = os.path.join(OUT, "cells_%s.pdf" % d["S"])
    files = [p for _, p in pdfs]
    if shutil.which("pdfunite"):
        subprocess.run(["pdfunite"] + files + [outp], check=True)
    elif shutil.which("gs"):
        subprocess.run(["gs", "-dBATCH", "-dNOPAUSE", "-q", "-sDEVICE=pdfwrite", "-sOutputFile=%s" % outp] + files,
                       check=True)
    else:
        raise RuntimeError("neither pdfunite nor gs on PATH")
    USED.append("results/crossovers/%s/per_cell/<barcode>_co.pdf (%d cells)" % (d["S"], len(files)))
    idx = [[k + 1, bc, int(per.loc[bc, "n_cos"]), fi(ht.molecules.get(bc, float("nan")) if len(ht) else float("nan")),
            ff(ht.haploidness.get(bc, float("nan")) if len(ht) else float("nan"), 2)] for k, (bc, _) in enumerate(pdfs)]
    return outp, idx


for sp, d in D.items():
    if "co" in d:
        guarded("interference %s" % sp, lambda d=d: interference(d))
    if "intf" in d:
        f = d["intf"]
        put(sp, "interference_pairs", f["pairs"])
        put(sp, "interference_short_ratio", "%.2f (null range %.2f–%.2f)" % (f["ratio"], f["ratio_lo"], f["ratio_hi"]))

for name, fn in (("markers", fig_markers), ("raw", fig_raw), ("cells", fig_cells), ("per_grain", fig_per_grain),
                 ("landscape", fig_landscape), ("per_chrom", fig_per_chrom), ("interference", fig_interference),
                 ("chrom_landscapes", fig_chrom_landscapes)):
    FIGS[name] = guarded("figure " + name, fn)

PAIRS = []
if have(PAIRS_FILE):
    for l in open(PAIRS_FILE):
        if len(l) > 60 and l.startswith("   ") and not l.startswith("    "):
            lab, rest = l[3:53].strip(), l[53:].split()
            if len(rest) >= 6 and all(t.isdigit() for t in rest[:5]):
                PAIRS.append((lab, rest[:6]))
for lab, vals in PAIRS:
    NUM["paradoxa.pair." + lab] = " ".join(vals)
FIGS["diagram"] = guarded("figure arm diagram", fig_diagram)
FIGS["structure"] = guarded("figure structure landscape", lambda: fig_structure("paradoxa"))
FIGS["rdna"] = guarded("figure rDNA", lambda: fig_rdna("paradoxa"))
FIGS["hic_dual"] = guarded("figure Hi-C both haplotypes", lambda: fig_hic_dual("paradoxa"))
for sp in D:
    FIGS["linkhic_" + sp] = guarded("figure linkage vs Hi-C " + sp, lambda sp=sp: fig_link_hic(sp))
EXAMPLES = guarded("example cells", example_cells) or {}
BROWSER = {}
if not A.no_cell_pdfs:
    for sp in D:
        BROWSER[sp] = guarded("cell browser " + sp, lambda sp=sp: cell_browser(sp))

# =========================================================================== html


def img(p, width=100):
    if not p or not os.path.exists(p):
        return "<p class='missing'>figure not built: %s</p>" % html.escape(str(p))
    b = base64.b64encode(open(p, "rb").read()).decode()
    return "<img src='data:image/png;base64,%s' style='width:%d%%'>" % (b, width)


def table(head, rows):
    h = "".join("<th>%s</th>" % html.escape(str(x)) for x in head)
    b = "".join("<tr>%s</tr>" % "".join("<td>%s</td>" % html.escape(str(x)) for x in r) for r in rows)
    return "<table><tr>%s</tr>%s</table>" % (h, b)


def details(title, p=None, text=None):
    if text is None:
        if not p or not os.path.exists(p):
            MISSING.append(p)
            return "<p class='missing'>missing: %s</p>" % html.escape(str(p))
        if p not in USED:
            USED.append(p)
        text = open(p, errors="replace").read()
    return "<details><summary>%s <span class='path'>%s</span></summary><pre>%s</pre></details>" % (
        html.escape(title), html.escape(p or ""), html.escape(text))


def P(text):
    return "<p>%s</p>" % text


def two(f):
    return [f(B), f(Pd)]


def sel(d, k):
    return d.get("sel", {}).get(k, "–")


H = []
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
commit = subprocess.run("git rev-parse --short HEAD", shell=True, stdout=subprocess.PIPE,
                        universal_newlines=True).stdout.strip()
TAG = subprocess.run("git tag -l 'binata-landscape*' | tail -n 1", shell=True, stdout=subprocess.PIPE,
                     universal_newlines=True).stdout.strip()
GOOD_TXT = ("heterozygous SNPs that segregate cleanly in the pollen: ALT in %s–%s%% of the nuclei that cover them, "
            "and one allele per nucleus" % (ff(100 * GOOD_LO, 0), ff(100 * GOOD_HI, 0)))

H.append("<h1>Meiotic crossover landscapes of <i>Drosera binata</i> and <i>D. paradoxa</i> from single pollen nuclei</h1>")
H.append("<p class='meta'>Generated %s from CO_smk commit %s by workflow/scripts/supervisor_report.py. Every number, "
         "table and figure is computed from the pipeline's outputs (files listed in the appendix).</p>" % (now, commit))
if PROBLEMS:
    H.append("<p class='missing'>Parts that could not be built: %s</p>" % html.escape(
        "; ".join("%s (%s)" % p for p in PROBLEMS)))

# ---- 0
H.append("<h2>0. Summary</h2>")
nj_txt = (" (%s without the crossovers at the translocation joins, an artefact of the reference)" % ff(Pd["co_mean_nj"], 2)
          if "co_mean_nj" in Pd else "")
H.append(P(
    "Single pollen nuclei of <i>D. binata</i> and <i>D. paradoxa</i> were sequenced (scRNA-seq) and genotyped at "
    "heterozygous SNPs called from HiFi reads; a crossover is called where a nucleus switches between the two "
    "haplotypes. <b>%s</b> (%s): %s nuclei, %s crossovers per grain against an obligate minimum of %s; landscape %s. "
    "<b>%s</b> (%s): %s nuclei, %s crossovers per grain%s against a minimum of %s; landscape %s. "
    "The <i>D. paradoxa</i> landscape is provisional: the pollen inherit the arms of chromosomes 1 and 2 in a different "
    "combination than the plant's tissue carries them (section 7)."
    % (NAME["binata"], B["row"]["centromere"], fi(B.get("ncell_co")), ff(B.get("co_mean"), 2), ff(B["n"] / 2), HEAD["binata"],
       NAME["paradoxa"], Pd["row"]["centromere"], fi(Pd.get("ncell_co")), ff(Pd.get("co_mean"), 2), nj_txt,
       ff(Pd["n"] / 2), HEAD["paradoxa"])))
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["centromere type", B["row"]["centromere"], Pd["row"]["centromere"]],
    ["chromosomes (n)", B["n"], Pd["n"]],
    ["assembly, both haplotypes (2C), chromosomes, Mb", ff(B.get("asm2n_chr_mb"), 0), ff(Pd.get("asm2n_chr_mb"), 0)],
    ["reference used for crossovers (1C), Mb", ff(B["G"] / 1e6, 0), ff(Pd["G"] / 1e6, 0)],
    ["good markers (%s)" % GOOD_TXT, fi(B.get("good_n")), fi(Pd.get("good_n"))],
    ["pollen nuclei used", fi(len(B["cells"])), fi(len(Pd["cells"]))],
    ["informative molecules per nucleus, median", fi(B.get("selected_molecules_median")),
     fi(Pd.get("selected_molecules_median"))],
    ["crossovers per grain, mean (SD)", "%s (%s)" % (ff(B.get("co_mean"), 2), ff(B.get("co_sd"), 2)),
     "%s (%s)%s" % (ff(Pd.get("co_mean"), 2), ff(Pd.get("co_sd"), 2),
                    "; %s without the join artefact" % ff(Pd["co_mean_nj"], 2) if "co_mean_nj" in Pd else "")],
    ["obligate minimum (one per bivalent)", ff(B["n"] / 2), ff(Pd["n"] / 2)],
    ["genetic map length, cM", fi(100 * B["co_mean"]) if ok(B.get("co_mean")) else "–",
     fi(100 * Pd["co_mean_nj"]) + " (without the join artefact)" if "co_mean_nj" in Pd
     else fi(100 * Pd.get("co_mean", float("nan")))],
    ["landscape", HEAD["binata"], HEAD["paradoxa"]],
    ["status", "near-final%s: crossover-calling settings need minor optimisation" % (" (git tag %s)" % TAG if TAG else ""),
     "provisional: reference under review, crossover-calling settings not yet tuned"]]))
H.append(P("2n assemblies: <i>D. binata</i> <code>%s</code>; <i>D. paradoxa</i> <code>%s</code>."
           % (html.escape(str(B.get("asm2n", "not found"))), html.escape(str(Pd.get("asm2n", "not found"))))))

# ---- 1
H.append("<h2>1. Material and data</h2>")


def ref_desc(d):
    if d["row"]["haplotype"] == "composite":
        return "composite of the two haplotypes: %s" % ", ".join(d["main"])
    return "haplotype %s of the phased assembly: %s" % (d["row"]["haplotype"].replace("hap", ""),
                                                        "%s–%s" % (d["main"][0], d["main"][-1]))


star_keys = [("Number of input reads", "reads into STARsolo"), ("Average input read length", "read length, bp"),
             ("Uniquely mapped reads %", "mapped uniquely"), ("% of reads mapped to multiple loci", "mapped to several loci"),
             ("% of reads mapped to too many loci", "mapped to too many loci"),
             ("% of reads unmapped: too short", "unmapped: too short"), ("% of reads unmapped: other", "unmapped: other")]
rows1 = [["sample in the pipeline", B["S"], Pd["S"]],
         ["reference used", ref_desc(B), ref_desc(Pd) + " (hap1, but chr2 from hap2, so each arm of the chr1/chr2 "
                                                           "translocation occurs once; section 7)"
          if Pd["S"] == "Dparadoxa_std" else ref_desc(Pd)],
         ["2n (config)", B["row"]["chr_number_2n"], Pd["row"]["chr_number_2n"]],
         ["assembly, both haplotypes (2C): chromosomes / all sequence, Mb",
          "%s / %s" % (ff(B.get("asm2n_chr_mb"), 0), ff(B.get("asm2n_all_mb"), 0)),
          "%s / %s" % (ff(Pd.get("asm2n_chr_mb"), 0), ff(Pd.get("asm2n_all_mb"), 0))],
         ["reference used (1C), chromosomes, Mb", ff(B["G"] / 1e6, 0), ff(Pd["G"] / 1e6, 0)],
         ["genes on these chromosomes (annotation)", fi(len(B["genes"])) if "genes" in B else "–",
          (fi(len(Pd["genes"])) + " (inflated: includes many transposable-element gene models)") if "genes" in Pd else "–"],
         ["scRNA chemistry", B["row"]["chemistry"], Pd["row"]["chemistry"]]]
for k, lab in star_keys:
    rows1.append([lab, B.get("star", {}).get(k, "–"), Pd.get("star", {}).get(k, "–")])
solo_keys = []
for d in (B, Pd):
    for k, _ in d.get("solo", []):
        if k not in solo_keys and "Cell" not in k:
            solo_keys.append(k)


def solo_val(d, k):
    v = dict(d.get("solo", [])).get(k, "–")
    try:
        x = float(v)
        return ("%.1f%%" % (100 * x)) if x <= 1 else fi(x)
    except ValueError:
        return v


for k in solo_keys:
    rows1.append(["STARsolo %s: %s" % (B.get("solo_feature", Pd.get("solo_feature", "")), k), solo_val(B, k),
                  solo_val(Pd, k)])
H.append(table(["", NAME["binata"], NAME["paradoxa"]], rows1))

# ---- 2
H.append("<h2>2. Pipeline at a glance</h2>")


def sel_txt(d):
    s = d.get("sel", {})
    t = "haploidness ≥ %s over ≥ %s windows; ≥ %s molecules on the weakest chromosome" % (
        s.get("min_haploidness", "–"), s.get("min_windows", "–"), s.get("min_chrom_molecules", "–"))
    mm = s.get("min_molecules", "0")
    return t + ("; ≥ %s molecules in total" % mm if mm not in ("0", "0.0", "–") else "; no total-molecule floor")


H.append(table(["step", "what it does", NAME["binata"], NAME["paradoxa"]], [
    ["1 markers", "heterozygous SNPs called from HiFi reads (bcftools), filtered on quality, depth, allele balance "
                  "and rDNA (section 3)", "%s raw → %s" % (fi(B.get("raw")), fi(B.get("filtered"))),
     "%s raw → %s" % (fi(Pd.get("raw")), fi(Pd.get("filtered")))],
    ["2 markers seen in pollen", "markers covered by pollen reads, classed by how they segregate; 'good' = %s" % GOOD_TXT,
     "%s seen → %s good" % (fi(B.get("seen")), fi(B.get("good_n"))),
     "%s seen → %s good" % (fi(Pd.get("seen")), fi(Pd.get("good_n")))],
    ["3 align pollen RNA", "STARsolo to the reference",
     "%s reads, %s unique" % (B.get("star", {}).get("Number of input reads", "–"),
                              B.get("star", {}).get("Uniquely mapped reads %", "–")),
     "%s reads, %s unique" % (Pd.get("star", {}).get("Number of input reads", "–"),
                              Pd.get("star", {}).get("Uniquely mapped reads %", "–"))],
    ["4 call barcodes", "EmptyDrops on UMI counts; generous, removes only empty droplets and barcode errors",
     "%s → %s" % (fi(NUM.get("binata.barcodes_tested")), fi(NUM.get("binata.cells_called"))),
     "%s → %s" % (fi(NUM.get("paradoxa.barcodes_tested")), fi(NUM.get("paradoxa.cells_called")))],
    ["5 keep single haploid nuclei", "one haplotype per window and enough evidence on every chromosome: " + sel_txt(B)
     + " (binata); " + sel_txt(Pd) + " (paradoxa)",
     "%s → %s" % (sel(B, "barcodes counted by cellsnp"), fi(len(B["cells"]))),
     "%s → %s" % (sel(Pd, "barcodes counted by cellsnp"), fi(len(Pd["cells"])))],
    ["6 call crossovers", "hapCO: blocks of informative molecules; a crossover where the haplotype switches",
     "%s crossovers" % fi(B.get("co_total")), "%s crossovers" % fi(Pd.get("co_total"))]]))

# ---- 3
H.append("<h2>3. Markers: what defines the two haplotypes</h2>")
H.append(img(FIGS.get("raw")))
H.append(P("Why the filters sit where they do: the depth window brackets the main heterozygous peak (single-copy "
           "sequence) and cuts the high-depth tail (collapsed repeats), and the allele-balance window keeps SNPs "
           "near 50:50, as true heterozygous sites are. Shares are of heterozygous biallelic SNPs in a regular "
           "1-in-k sample of the raw calls; %s%% (<i>D. binata</i>) and %s%% (<i>D. paradoxa</i>) pass all filters."
           % (ff(100 * B.get("raw_share_pass", float("nan"))), ff(100 * Pd.get("raw_share_pass", float("nan"))))))
H.append(img(FIGS.get("markers")))
cls_def = {"good": GOOD_TXT,
           "sticky": "nearly always the same allele: ALT in < %s%% or > %s%% of nuclei (assembly or genotyping error)"
                     % (ff(100 * STICKY, 0), ff(100 * (1 - STICKY), 0)),
           "paralog": "both alleles inside one nucleus in ≥ %s%% of deep observations: two copies collapsed in the "
                      "assembly" % ff(100 * PARALOG, 0),
           "other": "judged, but between sticky and good (ALT in %s–%s%% or %s–%s%%)"
                    % (ff(100 * STICKY, 0), ff(100 * GOOD_LO, 0), ff(100 * GOOD_HI, 0), ff(100 * (1 - STICKY), 0)),
           "unjudged": "covered in fewer than %s nuclei: too few to judge, not used" % ff(MIN_CELLS, 0)}
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["raw variants (HiFi)", fi(B.get("raw")), fi(Pd.get("raw"))],
    ["after filters", fi(B.get("filtered")), fi(Pd.get("filtered"))],
    ["filters (depth, quality, allele balance)",
     "%(min_dp)d–%(max_dp)d reads, QUAL ≥ %(min_qual)g, ALT %(lo)g–%(hi)g" % B["fp"] if "fp" in B else "–",
     "%(min_dp)d–%(max_dp)d reads, QUAL ≥ %(min_qual)g, ALT %(lo)g–%(hi)g" % Pd["fp"] if "fp" in Pd else "–"],
    ["markers removed for lying in rDNA (estimate from the sample) / rDNA sequence masked, kb",
     "≈%s / %s" % (fi(B.get("rdna_markers_est")), fi(B.get("rdna_bp", float("nan")) / 1e3)),
     "≈%s / %s" % (fi(Pd.get("rdna_markers_est")), fi(Pd.get("rdna_bp", float("nan")) / 1e3))],
    ["seen in pollen", fi(B.get("seen")), fi(Pd.get("seen"))]] + [
    ["class %s (%s)" % (k, cls_def[k]), fi(B.get("classes", {}).get(k)), fi(Pd.get("classes", {}).get(k))]
    for k in ("good", "sticky", "paralog", "other", "unjudged")] + [
    ["good markers per Mb", ff(B.get("good_per_mb")), ff(Pd.get("good_per_mb"))],
    ["marker clusters (good markers ≤ %s bp apart merged, as on one read)" % ff(HT_GAP, 0), fi(B.get("clusters")),
     fi(Pd.get("clusters"))],
    ["gap between clusters: median kb / 95th percentile kb / largest Mb",
     "%s / %s / %s" % (ff(B.get("gap_median_kb")), ff(B.get("gap_p95_kb"), 0), ff(B.get("gap_max_mb"))),
     "%s / %s / %s" % (ff(Pd.get("gap_median_kb")), ff(Pd.get("gap_p95_kb"), 0), ff(Pd.get("gap_max_mb")))],
    ["chromosome sequence inside gaps > 1 Mb, %", ff(100 * B.get("gap_1mb_share", float("nan"))),
     ff(100 * Pd.get("gap_1mb_share", float("nan")))],
    ["good markers inside annotated genes, %", ff(100 * B.get("good_in_genes", float("nan"))),
     ff(100 * Pd.get("good_in_genes", float("nan")))]]))
H.append(P("Markers come from pollen RNA, so they cluster in expressed genes, several per read; the distance between "
           "neighbouring markers is therefore tens of bp, and the gaps that matter are those between clusters."))

# ---- 4
H.append("<h2>4. Cells: how many, how deep, how clean</h2>")
H.append(img(FIGS.get("cells")))
H.append(P("Doublets: two nuclei in one droplet carry both haplotypes, so their 15-molecule windows sit at 0, 0.5 "
           "and 1 instead of 0 or 1, giving haploidness near 0.5 or below; the haploidness floor removes them, together "
           "with diploid barcodes. This filter replaced the switch-rate filter used for <i>Cuscuta</i> and "
           "<i>Spondias</i>: in these libraries most neighbouring markers sit on one read and agree automatically, "
           "so the switch rate measured read clustering rather than doublets (rule select_cells)."))
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["barcodes tested / called", "%s / %s" % (fi(NUM.get("binata.barcodes_tested")), fi(NUM.get("binata.cells_called"))),
     "%s / %s" % (fi(NUM.get("paradoxa.barcodes_tested")), fi(NUM.get("paradoxa.cells_called")))],
    ["called barcodes with allele counts (cellsnp)", sel(B, "barcodes counted by cellsnp"),
     sel(Pd, "barcodes counted by cellsnp")],
    ["scored for haploidness (≥ %s windows of %s molecules)" % (ff(HT_MINWIN, 0), ff(HT_W, 0)), sel(B, "scored (>= 5 full windows)"),
     sel(Pd, "scored (>= 5 full windows)")],
    ["pass haploidness and window count", sel(B, "pass haploidness and windows"), sel(Pd, "pass haploidness and windows")],
    ["floors", sel_txt(B), sel_txt(Pd)],
    ["selected nuclei", fi(len(B["cells"])), fi(len(Pd["cells"]))],
    ["UMIs per selected nucleus, median", fi(B.get("umi_med")), fi(Pd.get("umi_med"))],
    ["genes per selected nucleus, median", fi(B.get("genes_med")), fi(Pd.get("genes_med"))],
    ["informative molecules per nucleus, median (IQR)",
     "%s (%s–%s)" % (fi(B.get("selected_molecules_median")), fi(B.get("selected_molecules_q25")),
                     fi(B.get("selected_molecules_q75"))),
     "%s (%s–%s)" % (fi(Pd.get("selected_molecules_median")), fi(Pd.get("selected_molecules_q25")),
                     fi(Pd.get("selected_molecules_q75")))],
    ["weakest chromosome, molecules, median", fi(B.get("selected_weakest_chrom_molecules_median")),
     fi(Pd.get("selected_weakest_chrom_molecules_median"))],
    ["haploidness, median", ff(B.get("selected_haploidness_median"), 2), ff(Pd.get("selected_haploidness_median"), 2)],
    ["%g Mb windows with a clean genotype call, %%" % WIN_MB, ff(100 * B.get("called_share", float("nan"))),
     ff(100 * Pd.get("called_share", float("nan")))],
    ["genotype similarity between nuclei, median (unrelated = 0.5)", ff(B.get("pair_sim_median"), 2),
     ff(Pd.get("pair_sim_median"), 2)],
    ["near-identical pairs (≥ 95% of ≥ 20 shared windows agree)",
     "%s pairs among %s nuclei" % (B.get("dup_pairs", "–"), B.get("dup_cells", "–")),
     "%s pairs among %s nuclei" % (Pd.get("dup_pairs", "–"), Pd.get("dup_cells", "–"))],
    ["distinct genotypes if each near-identical group counts once", fi(B.get("distinct_genotypes")),
     fi(Pd.get("distinct_genotypes"))]]))
H.append(P("Near-identical pairs are most likely two nuclei from one pollen grain (vegetative and generative or sperm "
           "nuclei come from one meiotic product and share its genotype). They are one meiosis observed twice and "
           "should be counted once; collapsing them is listed under next steps."))
if EXAMPLES:
    H.append("<h3>One nucleus of median depth per species (hapCO's own plot, as used for manual review)</h3>")
    H.append("<table class='plain'><tr>%s</tr></table>" % "".join(
        "<td><b>%s</b> %s<br>%s</td>" % (NAME[sp], html.escape(bc), img(png, 100)) for sp, (bc, png) in EXAMPLES.items()))
if BROWSER:
    links = ["<a href='%s'>%s</a> (%s nuclei, %.0f MB)" % (os.path.basename(p), os.path.basename(p), len(ix),
                                                            os.path.getsize(p) / 1e6)
             for sp, v in BROWSER.items() if v for p, ix in [v]]
    H.append(P("Every nucleus, one page each, sorted by crossovers per grain (most first): %s. Keep the PDFs in the "
               "same folder as this page." % ", ".join(links)))
    for sp, v in BROWSER.items():
        if v:
            H.append(details("%s: page index (page, barcode, crossovers, molecules, haploidness)" % NAME[sp],
                             text="\n".join("%4d  %s  %3d  %6s  %s" % tuple(r) for r in v[1])))

# ---- 5
H.append("<h2>5. Crossover calling: settings and checks</h2>")
cc = CFG.get("co_calling", {})
meaning = [("input_rows", "what a block counts: molecules = informative markers on one read, counted once"),
           ("block_size", "minimum length of a haplotype block, bp"),
           ("marker_num", "minimum molecules in a block"),
           ("terminal_marker_num", "molecules needed by a chromosome's first and last block (no length rule there)"),
           ("base_af", "smoothing threshold, pass 1"), ("window_af", "smoothing threshold, pass 2"),
           ("genotype", "ALT share that decides a block's genotype")]


def status(sp):
    a = str((cc.get(D[sp]["species"]) or {}).get("_approved", "–"))
    return a.split("|")[0].strip() + (": " + a.split("|", 1)[1].strip() if "|" in a else "")


H.append(table(["hapCO setting", "meaning", NAME["binata"], NAME["paradoxa"]],
               [[k, m, (cc.get(B["species"]) or {}).get(k, "–"), (cc.get(Pd["species"]) or {}).get(k, "–")]
                for k, m in meaning] + [["status", "when and how the settings were chosen", status("binata"),
                                         status("paradoxa")]]))
H.append(img(FIGS.get("per_grain")))
H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
    ["grains", fi(B.get("ncell_co")), fi(Pd.get("ncell_co"))],
    ["crossovers", fi(B.get("co_total")), fi(Pd.get("co_total"))],
    ["per grain, mean (SD); median", "%s (%s); %s" % (ff(B.get("co_mean"), 2), ff(B.get("co_sd"), 2), ff(B.get("co_median"), 0)),
     "%s (%s); %s" % (ff(Pd.get("co_mean"), 2), ff(Pd.get("co_sd"), 2), ff(Pd.get("co_median"), 0))],
    ["crossovers within 10 Mb of the chr1/chr2 joins (%s), observed / expected at the chromosome's average"
     % ("joins the pollen do not support: reference artefact" if Pd["S"] in ARTEFACT_JOINS
        else "joined as the pollen inherit them: no excess expected"),
     "–", "%s / %s" % (fi(Pd.get("at_join_cos")), ff(Pd.get("at_join_exp"), 1)) if "at_join_cos" in Pd else "–"],
    ["per grain without them", "–", ff(Pd.get("co_mean_nj"), 2) if "co_mean_nj" in Pd else "not applicable"],
    ["obligate minimum per grain (n / 2)", ff(B["n"] / 2), ff(Pd["n"] / 2)],
    ["grains below the minimum", "%s of %s" % (fi(B.get("below_obligate")), fi(B.get("ncell_co"))),
     "%s of %s" % (fi(Pd.get("below_obligate")), fi(Pd.get("ncell_co")))],
    ["chromosomes below 0.5 per grain", "%s of %d" % (len(B.get("below_half", [])), B["n"]),
     "%s of %d" % (len(Pd.get("below_half", [])), Pd["n"])],
    ["crossover interval, median kb (resolution)", fi(B.get("width_median_kb")), fi(Pd.get("width_median_kb"))]]))
H.append(P("<i>D. paradoxa</i> nuclei carry %s informative molecules against %s in <i>D. binata</i> (medians), so its "
           "crossovers are placed less precisely (median interval %s kb against %s kb) and a nucleus with few molecules "
           "on a chromosome can miss one; its settings are binata's, not yet tuned."
           % (fi(Pd.get("selected_molecules_median")), fi(B.get("selected_molecules_median")),
              fi(Pd.get("width_median_kb")), fi(B.get("width_median_kb")))))

# ---- 6
H.append("<h2>6. Landscapes</h2>")
H.append(img(FIGS.get("landscape")))
rows6 = [["crossovers: ends vs middle", B["shape"], Pd["shape"]]]
if "shape_nj" in Pd:
    rows6.append(["crossovers without the join artefact", "–", Pd["shape_nj"]])
rows6 += [["good markers: ends vs middle", ff(B.get("U_markers"), 2), ff(Pd.get("U_markers"), 2)],
          ["genes: ends vs middle", ff(B.get("U_genes"), 2), ff(Pd.get("U_genes"), 2)],
          ["crossovers vs good markers per 5 Mb window, Spearman", ff(B.get("rho_markers"), 2), ff(Pd.get("rho_markers"), 2)]]
H.append(table(["", NAME["binata"], NAME["paradoxa"]], rows6))
H.append(P("'Ends vs middle' = density in the outer 20%% of each chromosome (10%% at each end) over the inner 60%%; "
           "1 = flat. Markers follow genes, so the same ratio for markers and genes shows how much end bias detection "
           "alone could make: in <i>D. binata</i> crossovers are %s× end-biased against %s× for markers, so the U is "
           "not a detection artefact, although crossovers do track marker density window by window (Spearman %s)."
           % (ff(B.get("land", {}).get("U")), ff(B.get("U_markers")), ff(B.get("rho_markers"), 2))))
H.append("<h3>Busiest 5 Mb windows (Poisson p against the genome-wide mean per window)</h3>")
for sp, d in D.items():
    if d.get("busiest"):
        H.append("<p><b>%s</b> (mean %s crossovers per window)</p>" % (NAME[sp], ff(d["window_mean"])))
        H.append(table(["chromosome", "Mb", "crossovers", "expected", "p", "within 10 Mb of a join"], d["busiest"]))
H.append(img(FIGS.get("per_chrom")))
H.append("<h3>Crossover interference</h3>")
H.append(img(FIGS.get("interference")))
H.append(P("For every nucleus and chromosome with two or more crossovers, the distance between neighbouring crossovers "
           "as a share of the chromosome; the line is what the same crossovers would give if placed independently "
           "(drawn from that chromosome's crossover positions, 300 times; band = 95%%). Interference shows as too few "
           "close pairs. Closer than 20%% of the chromosome: <i>D. binata</i> %s× expected (%s pairs), "
           "<i>D. paradoxa</i> %s× (%s pairs%s). This replaces the coefficient-of-"
           "coincidence curve, which needs more double crossovers than <i>D. paradoxa</i>'s %s nuclei give."
           % (ff(B.get("intf", {}).get("ratio"), 2), B.get("intf", {}).get("pairs", "–"),
              ff(Pd.get("intf", {}).get("ratio"), 2), Pd.get("intf", {}).get("pairs", "–"),
              "; crossovers at the joins left out" if Pd["S"] in ARTEFACT_JOINS else "", fi(Pd.get("ncell_co")))))
if FIGS.get("structure"):
    H.append("<h3><i>D. paradoxa</i>: crossovers along each chromosome against the structure of its homolog</h3>")
    H.append(img(FIGS.get("structure")))
    H.append(P("Under each chromosome: where the other haplotype's copy aligns to it (hap2 on hap1, synteny blocks from "
               "translocation_map.py), collinear, inverted or from a differently numbered chromosome, and where nothing "
               "aligns. Black bars: stretches without a crossover in any of the %s nuclei, long enough to be unlikely by "
               "chance (table below). Arm colours as in the chr1/chr2 diagram. On %s the arms beyond the joins "
               "(%s) have no partner in the other haplotype by construction: their partners are the other chr1/chr2 "
               "homolog (section 7)." % (fi(Pd.get("ncell_co")), join_names(Pd["S"])[0], join_names(Pd["S"])[2])))
if Pd.get("dead") is not None:
    dc = Pd["dead_chance"]
    H.append(P("<b>Stretches without crossovers</b>%s. With %s crossovers per grain, "
               "gaps up to %d Mb arise by chance, so only stretches of at least %d Mb are listed: %d, %s Mb in total, "
               "against %s expected if the same crossovers fell uniformly along each chromosome (95%% range %s–%s). "
               "What each overlaps:"
               % (" (crossovers at the joins left out)" if Pd["S"] in ARTEFACT_JOINS else "",
                  ff(Pd.get("co_mean_nj", Pd.get("co_mean")), 2), Pd["dead_thr"] - 5, Pd["dead_thr"], len(Pd["dead"]),
                  ff(Pd["dead_mb"], 0), ff(dc[0]), ff(dc[1], 0), ff(dc[2], 0))))
    H.append(table(["chromosome", "Mb", "length, Mb", "inverted, %", "translocated, %", "no homolog, %",
                    "markers vs chromosome median", "contains a join", "most likely reason"], Pd["dead"]))

# ---- 7
H.append("<h2>7. Is each reference right?</h2>")
H.append(P("Left, pollen linkage: for every pair of %g Mb windows, r = share of nuclei whose genotypes disagree (0 = "
           "always inherited together, 0.5 = independent), shown as 1 − 2r; in a correct reference each chromosome is "
           "one dark block on the diagonal and everything off it is pale. Right, Hi-C on the same chromosomes, raw "
           "contacts on a log scale from white (the typical contact between two different chromosomes) to red: each "
           "chromosome of the plant's tissue is one red square. A misassembly would break a square or put a red block "
           "off the diagonal. Hi-C reads come from one chromosome copy at a time, so a translocation between the two "
           "copies cannot show here; it shows in the pollen (off-diagonal blocks on the left) and in the map on both "
           "haplotypes below." % WIN_MB))
for sp in D:
    H.append("<h3>%s</h3>" % NAME[sp])
    H.append(img(FIGS.get("linkhic_" + sp)))
    H.append(details("linkage scan summary", "qc/linkage/%s/linkage_summary.txt" % D[sp]["S"]))
if FIGS.get("hic_dual"):
    H.append("<h3><i>D. paradoxa</i>: Hi-C on both haplotypes</h3>")
    H.append(img(FIGS.get("hic_dual"), 80))
    H.append(P("Reads that fit either copy of a sequence put contacts between the two copies, so homologous "
               "sequence on two chromosomes shows as a red line off the diagonal: for an ordinary pair (chr3 h1 and "
               "chr3 h2, say) one line through their block. Right, chromosomes 1 and 2 with their arms coloured as in "
               "the diagram below: if the tissue chromosomes are as assembled, L1 links A with C, L2 links B with D, "
               "P links A with B and Q links C with D, so each line changes partner at a dotted join."))
H.append("<h3><i>D. paradoxa</i>: chromosomes 1 and 2 in the tissue and in the pollen</h3>")
H.append(img(FIGS.get("diagram")))
H.append(P("Names: L1 and L2 are the left arms of chr1 and chr2, P and Q the right arms; each arm exists twice. "
           "A, B, C and D are the four assembled chromosomes; the crossover reference used here (%s) contains %s "
           "plus chr3–6 hap1. Hi-C supports all four as continuous molecules, yet the pollen inherit L1 with Q and L2 "
           "with P. Window pairs behind the diagram (measured on A + D):" % (Pd["S"], join_names(Pd["S"])[1])))


def binom_p(k, n):
    pk = [math.comb(n, i) * 0.5 ** n for i in range(n + 1)]
    return min(1.0, sum(p for p in pk if p <= pk[k] * (1 + 1e-9)))


if PAIRS:
    prow = []
    for lab, v in PAIRS:
        rr, ra, ar, aa, n, r = int(v[0]), int(v[1]), int(v[2]), int(v[3]), int(v[4]), float(v[5])
        r1, r2 = rr + ra, rr + ar
        e_same = (r1 * r2 + (n - r1) * (n - r2)) / max(n, 1)
        reading = "inherited together" if r <= 0.15 else ("independent" if r >= 0.35 else "partly linked")
        prow.append([lab, rr, ra, ar, aa, n, "%.2f" % r, "%d / %.1f" % (rr + aa, e_same),
                     "%d:%d (p %.2g)" % (r1, n - r1, binom_p(r1, n)), "%d:%d (p %.2g)" % (r2, n - r2, binom_p(r2, n)),
                     reading])
    H.append(table(["window pair", "REF-REF", "REF-ALT", "ALT-REF", "ALT-ALT", "nuclei", "r",
                    "same-type pairs: observed / expected if independent", "window 1 REF:ALT (p)",
                    "window 2 REF:ALT (p)", "reading"], prow))
    H.append(P("How to read it: when two windows are inherited independently all four classes appear, including the "
               "two that look like an intact tissue chromosome (REF-REF and ALT-ALT), and they appear about as often as "
               "chance predicts (the 'observed / expected if independent' column). When they are inherited together "
               "only REF-REF and ALT-ALT appear. So nuclei that carry A's L1 with A's P exist, but no more often than "
               "if L1 and P were on different chromosomes. The REF:ALT columns test each window for segregation bias "
               "(two-sided binomial against 50:50)."))
H.append("<h3>rDNA: which arrays are real and large (HiFi read depth)</h3>")
H.append(P("barrnap finds rRNA gene fragments in many places, but an assembled copy says little about array size: a "
           "large array collapses to a few copies in the assembly while its reads pile up there. Reads over each array "
           "against a typical 20 kb window estimate how many copies collapsed into it. Trusted here: a complete 45S "
           "unit (18S and 28S together), at least 5× the reads of a typical window, or at least 50 assembled 5S genes."))
H.append(img(FIGS.get("rdna")))
if Pd.get("rdna_at_join"):
    H.append(P("At the joins: %s. A collapsed array shows as many times the typical depth; a complete unit at about "
               "1× is a single copy, not an array." % "; ".join(Pd["rdna_at_join"])))
for sp, d in D.items():
    if d.get("rdna_rows"):
        keep = [r for r in d["rdna_rows"] if r[9] or r[6]]
        H.append("<p><b>%s</b> (%s): %d arrays, %d trusted; typical 20 kb window %s reads. Trusted arrays and every "
                 "array within 5 Mb of a join:</p>" % (NAME[sp], d["S"], d["rdna_n_arrays"], d["rdna_n_trust"],
                                                       ff(d["rdna_base"], 0)))
        H.append(table(["chromosome", "Mb", "type", "genes", "complete 45S unit", "HiFi reads", "× a typical window",
                        "trusted", "join nearby"],
                       [[r[0], r[1], r[2], r[3], "yes" if r[7] else "", fi(r[4]), fi(r[5]), "yes" if r[9] else "",
                         r[6]] for r in keep[:15]]))
    elif d.get("rdna_note"):
        H.append(P("%s: rDNA depth %s" % (NAME[sp], d["rdna_note"])))
for t, p in (("Hi-C, every library at the joins", os.path.join(A.phd_root, "results/Drosera_paradoxa/qc/hic_remap/library_joins.txt")),
             ("chr5 against chr6", "qc/linkage/Dparadoxa_std/chr5_chr6_check.txt"),
             ("rDNA arrays, crossover hotspots and joins", "qc/rdna/rdna_overview.txt"),
             ("crossover reference C + B against A + D", "qc/linkage/Dparadoxa_CB/cb_check.txt")):
    H.append(details(t, p))

# ---- 8, 9
H.append("<h2>8. Limitations</h2>")
ref_line = ("The <i>D. paradoxa</i> reference used here joins chr1 and chr2 the way the tissue does, not the way the "
            "pollen inherit them; %s crossovers at those joins are an artefact of that choice (section 7)."
            % fi(Pd.get("join_cos")) if Pd["S"] == "Dparadoxa_std" else
            "The <i>D. paradoxa</i> reference used here joins chr1 and chr2 the way the pollen inherit them; the tissue "
            "(Hi-C) shows the other arrangement, unresolved until cytology (section 7).")
H.append("<ul><li>Markers come from pollen RNA and sit in expressed genes (%s%% and %s%% of good markers inside annotated "
         "genes); %s%% and %s%% of the chromosome sequence lies in marker gaps over 1 Mb, where crossovers are "
         "placed only roughly.</li>"
         "<li><i>D. paradoxa</i> rests on %s nuclei with a median of %s informative molecules, against %s nuclei with "
         "%s in <i>D. binata</i>.</li>"
         "<li>Near-identical nuclei (%s pairs in <i>D. binata</i>, %s in <i>D. paradoxa</i>) are still counted "
         "separately.</li><li>%s</li></ul>"
         % (ff(100 * B.get("good_in_genes", float("nan")), 0), ff(100 * Pd.get("good_in_genes", float("nan")), 0),
            ff(100 * B.get("gap_1mb_share", float("nan"))), ff(100 * Pd.get("gap_1mb_share", float("nan"))),
            fi(len(Pd["cells"])), fi(Pd.get("selected_molecules_median")), fi(len(B["cells"])),
            fi(B.get("selected_molecules_median")), B.get("dup_pairs", "–"), Pd.get("dup_pairs", "–"), ref_line))
H.append("<h2>9. Next steps and asks</h2>")
H.append("<ul><li>%s</li>"
         "<li>Cytology: meiotic chromosome spreads with oligo-FISH paints for the chr1/chr2 arms, chr5 and chr6, and "
         "45S rDNA FISH.</li>"
         "<li>Count near-identical nuclei once, in both species.</li>"
         "<li>Tune the crossover-calling settings: minor for <i>D. binata</i>; for <i>D. paradoxa</i> on the review panel "
         "once the reference is settled.</li></ul>"
         % ("<i>D. paradoxa</i> on the reference joined as the pollen inherit it (chr1_hap2 + chr2_hap1): %s."
            % ("results in section 7" if os.path.exists("qc/linkage/Dparadoxa_CB/cb_check.txt") else "running")
            if Pd["S"] == "Dparadoxa_std" else
            "<i>D. paradoxa</i> crossovers are now called on %s (%s); the comparison with A + D, join by join, is in "
            "section 7 (crossover reference C + B against A + D)." % (Pd["S"], join_names(Pd["S"])[1])))

# ---- appendix
H.append("<h2>Appendix</h2><h3>A. Per-chromosome landscapes</h3>")
H.append(img(FIGS.get("chrom_landscapes")))
H.append("<h3>B. Crossover-calling parameters as run</h3>")
for sp, d in D.items():
    H.append(details("%s: co_calling_params.txt" % NAME[sp], "results/crossovers/%s/co_calling_params.txt" % d["S"]))
H.append("<h3>C. Files read</h3><pre>%s</pre>" % html.escape("\n".join(sorted(set(map(str, USED))))))
if MISSING:
    H.append("<h3>D. Files not found</h3><pre>%s</pre>" % html.escape("\n".join(sorted(set(map(str, MISSING))))))

CSS = """body{font-family:Helvetica,Arial,sans-serif;max-width:1080px;margin:24px auto;padding:0 16px;color:#2C2C2A;
font-size:13px;line-height:1.45}h1{font-size:21px}h2{font-size:16px;border-bottom:1px solid #D3D1C7;margin-top:30px}
h3{font-size:13px;margin-top:18px}table{border-collapse:collapse;margin:8px 0 14px;font-size:12px}td,th{border:1px
solid #D3D1C7;padding:3px 8px;text-align:left;vertical-align:top}th{background:#F1EFE8}table.plain td{border:none}
.meta{color:#5F5E5A;font-size:11px}.missing{color:#A32D2D}.path{color:#888780;font-size:11px}code{font-size:11px}
pre{font-size:10.5px;background:#F7F6F2;padding:8px;overflow-x:auto}img{display:block;margin:8px 0}details{margin:4px 0}
@media print{img{page-break-inside:avoid}table{page-break-inside:avoid}}"""
open(os.path.join(OUT, "report.html"), "w").write(
    "<!doctype html><html><head><meta charset='utf-8'><title>Drosera crossover landscapes</title><style>%s</style>"
    "</head><body>%s</body></html>" % (CSS, "\n".join(H)))
with open(os.path.join(OUT, "numbers.txt"), "w") as f:
    for k in sorted(NUM):
        v = NUM[k]
        f.write("%s\t%s\n" % (k, ("%.4g" % v) if isinstance(v, (float, np.floating)) else v))
print("\n".join("%s\t%s" % (k, ("%.4g" % NUM[k]) if isinstance(NUM[k], (float, np.floating)) else NUM[k])
                for k in sorted(NUM)))
print("\nparts not built: %d%s" % (len(PROBLEMS), "".join("\n  %s: %s" % p for p in PROBLEMS)))
print("missing files: %d%s" % (len(set(map(str, MISSING))), "".join("\n  " + str(m) for m in sorted(set(map(str, MISSING))))))
print("wrote %s/report.html%s" % (OUT, "".join(", %s" % v[0] for v in BROWSER.values() if v)))
