#!/usr/bin/env python3
"""supervisor_report.py -- D. binata and D. paradoxa crossover landscapes, start to finish (version 5).

One self-contained HTML page for the supervisor, plus one PDF per species with every nucleus's
crossover plot. Every number, table and figure is computed here from the pipeline's own outputs
(each file read is listed in the appendix of the page); nothing is typed in by hand. Sentences
that state a number build it from the computed value. A section that cannot be built is shown
in red with the error instead of stopping the report.

  0 summary  1 material  2 pipeline at a glance  3 markers  4 cells  5 crossover calling: settings from the
  review  6 landscapes  7 is each reference right? (pollen co-segregation, Hi-C, flags)  8 limitations
  9 next steps  appendix

Usage (from the CO_smk root):
  supervisor_report.py [--binata Dbinata_hap1] [--paradoxa Dparadoxa_split] [--out reports/supervisor_v5]
                       [--previous reports/supervisor_CB/report.html,...] [--no_cell_pdfs]
Writes OUT/report.html, OUT/cells_<sample>.pdf, OUT/fig/*.png and OUT/numbers.txt.

Sequences and chromosomes. A crossover reference can be cut into pieces (results/reference/<sample>/pieces.tsv;
Dparadoxa_split is Dparadoxa_CB cut at its three translocation breakpoints, so no join between translocated
pieces is assumed). Pieces are the units of crossover calling, dosage masks and the pollen co-segregation test;
landscapes, interference and the one-crossover-per-chromosome checks use the chromosomes the pieces were cut
from (positions mapped back, cut points drawn). Without pieces.tsv both are the reference's n longest sequences.
Settings: from the review of the calls (docs/co_review/<sample>_{events.tsv.gz,settings.tsv,verdicts.csv}) and its
scores (qc/co_review/<sample>/scores.tsv and scores.txt from score_verdicts.py).
Version 4 (Dparadoxa_std / Dparadoxa_CB, the chr1/chr2 joins, rDNA at the joins) is not repeated; section 7
points to it (--previous).
"""
import argparse
import base64
import collections
import csv
import datetime
import gzip
import html
import math
import os
import re
import shlex
import shutil
import subprocess
import traceback
import zlib

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ap = argparse.ArgumentParser()
ap.add_argument("--binata", default="Dbinata_hap1")
ap.add_argument("--paradoxa", default="Dparadoxa_split")
ap.add_argument("--out", default="reports/supervisor_v5")
ap.add_argument("--refs_root", default="/netscratch/dep_mercier/grp_marques/Aaryan/refs")
ap.add_argument("--phd_root", default="/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd")
ap.add_argument("--previous", default="reports/supervisor_CB/report.html,reports/supervisor/report.html")
ap.add_argument("--no_cell_pdfs", action="store_true")
A = ap.parse_args()
OUT, FIG = A.out, os.path.join(A.out, "fig")
os.makedirs(FIG, exist_ok=True)
SPP = [("binata", A.binata), ("paradoxa", A.paradoxa)]
COL = {"binata": "#0F6E56", "paradoxa": "#534AB7"}
NAME = {"binata": "D. binata", "paradoxa": "D. paradoxa"}
CMAP = {sp: LinearSegmentedColormap.from_list("phi_" + sp, ["#FFFFFF", c]) for sp, c in COL.items()}
for _cm in CMAP.values():
    _cm.set_bad("#E9E7E0")
# windows of WIN_MB per reference sequence; a nucleus's genotype in a window needs MIN_MOL molecules, >= GT_MAJ of
# them one allele. Pollen co-segregation: phi between two windows over the nuclei that have both (>= PHI_N), against
# NPERM shuffles. Near-identical nuclei: >= DUP_AGREE of >= DUP_MIN windows called in both agree; NSHUF_DUP shuffles.
WIN_MB, MIN_MOL, GT_MAJ, BINS, NBOOT, NPERM, PHI_N = 10, 3, 0.8, 20, 500, 200, 20
DUP_AGREE, DUP_MIN, NSHUF_DUP = 0.95, 20, 500
USED, MISSING, NUM, PROBLEMS, WARN = [], [], {}, [], []
def stream(*key):
    """A random stream of its own for each analysis and species, so changing one leaves the others as they were."""
    return np.random.default_rng([1] + [zlib.crc32(str(k).encode()) for k in key])


plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})


# ============================================================================ helpers
def have(p):
    if p and os.path.exists(p):
        USED.append(p)
        return True
    MISSING.append(p)
    return False


def have_opt(p):
    """An optional input: listed as read when there, not as missing when absent."""
    if p and os.path.exists(p):
        USED.append(p)
        return True
    return False


def natural(c):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", c)]


def short(c):
    """chr4_hap1 -> chr4"""
    return re.sub(r"_hap\d+", "", c)


def put(sp, key, val):
    NUM["%s.%s" % (sp, key)] = val
    return val


def ok(x):
    return x is not None and not (isinstance(x, float) and math.isnan(x))


def fi(x):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "–"
    return "{:,}".format(int(round(x))) if ok(x) and not math.isinf(x) else "–"


def ff(x, d=1):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "–"
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


def fai_lens(p):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(p) if l.strip()}


def guarded(name, fn):
    try:
        return fn()
    except Exception as e:
        PROBLEMS.append((name, "%s: %s" % (type(e).__name__, e)))
        traceback.print_exc()
        return None


def poisson_logtail(k, mu):
    """log P(X >= k), X ~ Poisson(mu): the upper terms summed in log space (no 1 - CDF cancellation, no
    underflow); stops once the terms fall 30 log units below the sum."""
    if k <= 0:
        return 0.0
    if mu <= 0:
        return float("-inf")
    i = int(k)
    ls = -mu + i * math.log(mu) - math.lgamma(i + 1)
    while True:
        i += 1
        lt = -mu + i * math.log(mu) - math.lgamma(i + 1)
        ls = max(ls, lt) + math.log1p(math.exp(-abs(ls - lt)))
        if i > mu and lt < ls - 30:
            return min(0.0, ls)


def p_txt(lp):
    """A p-value from its natural log, readable at any size: 3.1e-05, 2e-180."""
    if lp >= math.log(1e-300):
        return "%.1e" % math.exp(lp)
    e = lp / math.log(10)
    return "%.1fe%d" % (10 ** (e - math.floor(e)), math.floor(e))


def uidx(v):
    """Crossover density in the outer 20% of the chromosome over the inner 60% (BINS bins)."""
    o, i = v[[0, 1, BINS - 2, BINS - 1]].sum(), v[4:BINS - 4].sum()
    return (o / 0.2) / (i / 0.6) if i > 0 else float("inf")


def profile(C, rng):
    """Relative-position profile and ends-vs-middle ratio, bootstrapped over cells (rows of C)."""
    v = C.sum(0)
    prof, us = [], []
    for _ in range(NBOOT):
        vv = C[rng.integers(0, len(C), len(C))].sum(0)
        prof.append(vv / max(vv.mean(), 1e-9))
        us.append(uidx(vv))
    return dict(prof=v / max(v.mean(), 1e-9), lo=np.percentile(prof, 2.5, 0), hi=np.percentile(prof, 97.5, 0),
                U=uidx(v), U_lo=np.percentile(us, 2.5), U_hi=np.percentile(us, 97.5))


def fg(x):
    """A ratio with as many decimals as it needs: 0.83, 12.4, 178."""
    if not ok(x) or math.isinf(float(x)):
        return "–"
    x = float(x)
    return "%.0f" % x if abs(x) >= 100 else ("%.1f" % x if abs(x) >= 10 else "%.2f" % x)


def shape(U, lo, hi):
    if not ok(U):
        return "–"
    if math.isinf(U):
        return "U-shaped: no crossovers in the middle 60%"
    ci = "95%% CI %s–%s" % (fg(lo), fg(hi))
    if lo > 1.5:
        return "U-shaped: ends %s× the middle (%s)" % (fg(U), ci)
    if lo <= 1 <= hi:
        return "flat: ends %s× the middle (%s)" % (fg(U), ci)
    if hi < 1:
        return "ends depleted: %s× the middle (%s)" % (fg(U), ci)
    return "mild end bias: ends %s× the middle (%s)" % (fg(U), ci)


def code_const(path, pattern, n):
    """Read constants from a pipeline script, so definitions in the text match the code."""
    m = re.search(pattern, open(path).read()) if have(path) else None
    return [float(x) for x in m.groups()][:n] if m else [float("nan")] * n


def sid_of(base_af, mn, bs_bp, t):
    return "b%g_mn%d_bs%g_t%d" % (float(base_af), int(mn), float(bs_bp) / 1e6, int(t))


def sid_parse(s):
    m = re.fullmatch(r"b([\d.]+)_mn(\d+)_bs([\d.]+)_t(\d+)", str(s))
    return (float(m.group(1)), int(m.group(2)), float(m.group(3)), int(m.group(4))) if m else None


def label(d, s):
    """A reference sequence as people say it: chr7 (a chromosome), chr1 L1 or chr4 core (pieces of chr1, chr4)."""
    c = d["src"][s][0]
    if s == c:
        return short(c)
    return short(c) + " " + s[len(c) + 1:] if s.startswith(c + "_") else s


def to_chrom(d, seq, pos):
    """Reference-sequence coordinates to chromosome coordinates (a piece's offset in the chromosome it was cut from)."""
    seq = list(seq)
    return (np.array([d["src"][s][0] for s in seq], dtype=object),
            np.asarray(pos, dtype=float) + np.array([d["src"][s][1] for s in seq], dtype=float))


def pieces_of(d, c):
    """The reference sequences of chromosome c, in chromosome order."""
    return sorted([s for s in d["seqs"] if d["src"][s][0] == c], key=lambda s: d["src"][s][1])


def bin_cells(d, c, nb, w):
    """Nuclei that count for each w-bp bin of chromosome c (those of the piece holding the bin's middle)."""
    out = np.zeros(nb)
    for i in range(nb):
        x = (i + 0.5) * w
        s = next((s for s in pieces_of(d, c) if d["src"][s][1] <= x < d["src"][s][1] + d["Ls"][s]), pieces_of(d, c)[-1])
        out[i] = d.get("nseq_used", {}).get(s, d.get("ncell_co", 1)) or 1
    return out


SAMPLES = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}
CFG = yaml.safe_load(open("config/config.yaml")) or {}
CC = CFG.get("co_calling", {}) or {}
MIN_CELLS, = code_const("workflow/scripts/marker_segregation.py", r"MIN_CELLS\s*=\s*(\d+)", 1)
STICKY, = code_const("workflow/scripts/marker_segregation.py", r"STICKY\s*=\s*([\d.]+)", 1)
GOOD_LO, GOOD_HI = code_const("workflow/scripts/marker_segregation.py", r"GOOD\s*=\s*\(([\d.]+),\s*([\d.]+)\)", 2)
PARALOG, = code_const("workflow/scripts/marker_segregation.py", r"PARALOG\s*=\s*([\d.]+)", 1)
HT_W, HT_GAP, HT_MINWIN = code_const("workflow/scripts/haplotype_tracks.py", r"W, GAP, MIN_WIN = (\d+), (\d+), (\d+)", 3)
AUTO_FRAC, = code_const("workflow/scripts/review_page.py", r'"--auto_frac", type=float, default=([\d.]+)', 1)
AUTO_SUP, = code_const("workflow/scripts/review_page.py", r'"--auto_support", type=int, default=(\d+)', 1)
AUTO_PUR, = code_const("workflow/scripts/review_page.py", r'"--auto_purity", type=float, default=([\d.]+)', 1)
EXP_LOW, = code_const("workflow/scripts/co_masks.py", r"if exp < (\d+):", 1)
MISS_COV, = code_const("workflow/scripts/co_masks.py", r"elif cov < ([\d.]+):", 1)
SW_MIX, DBL_COV = code_const("workflow/scripts/co_masks.py", r"elif sw >= ([\d.]+):\s*\n\s*cl = \"doubled\" if cov >= ([\d.]+)", 2)

# =========================================================================== load
D = {}
for sp, S in SPP:
    if S not in SAMPLES:
        raise SystemExit("sample %s is not in config/samples.csv" % S)
    row = SAMPLES[S]
    d = {"S": S, "row": row, "species": row["species"]}
    fai = row["assembly_fasta"] + ".fai"
    if not have(fai):
        raise SystemExit("no reference index for %s: %s" % (S, fai))
    L = fai_lens(fai)
    n = int(row["chr_number_2n"]) // 2
    pp = os.path.join(os.path.dirname(row["assembly_fasta"]), "pieces.tsv")
    if have_opt(pp):
        pcs = list(csv.DictReader(open(pp), delimiter="\t"))
        absent = [r["piece"] for r in pcs if r["piece"] not in L]
        if absent:
            raise SystemExit("pieces of %s not in %s: %s" % (S, fai, ", ".join(absent)))
        seqs = [r["piece"] for r in pcs]
        src = {r["piece"]: (r["source"], int(r["start"]) - 1) for r in pcs}
        Lc = {}
        for r in pcs:
            Lc[r["source"]] = max(Lc.get(r["source"], 0), int(r["end"]))
            ln = int(r["end"]) - int(r["start"]) + 1
            if L[r["piece"]] != ln:
                PROBLEMS.append(("reference pieces %s" % S, "%s is %s bp in pieces.tsv but %s bp in %s"
                                 % (r["piece"], fi(ln), fi(L[r["piece"]]), fai)))
        chroms = sorted(Lc, key=natural)
        cuts = {c: sorted(int(r["start"]) - 1 for r in pcs if r["source"] == c and int(r["start"]) > 1) for c in chroms}
    else:
        seqs = sorted(sorted(L, key=lambda c: -L[c])[:n], key=natural)
        src = {c: (c, 0) for c in seqs}
        Lc = {c: L[c] for c in seqs}
        chroms, cuts = list(seqs), {c: [] for c in seqs}
    d.update(n=n, seqs=seqs, Ls={s: L[s] for s in seqs}, src=src, main=chroms, L=Lc, cuts=cuts,
             G=sum(Lc.values()), pieces=any(cuts.values()))
    d["how"] = {}
    plan = "config/split_plans/%s.csv" % S
    if d["pieces"] and have_opt(plan):
        d["how"] = {r["piece"]: r.get("how", "") for r in csv.DictReader(l for l in open(plan) if not l.startswith("#"))}
    put(sp, "chromosomes", n)
    put(sp, "reference_sequences", len(seqs))
    put(sp, "reference_1C_mb", d["G"] / 1e6)
    # the phased two-haplotype assembly the reference came from (MANIFEST of the reference, else of the published hap1)
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

    # ---- markers (reference-sequence coordinates in the files; 'chrom'/'pos' here are chromosome coordinates)
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
        g = mc[(mc["class"] == "good") & mc.chrom.isin(seqs)][["chrom", "pos"]].copy()
        g["seq"] = g.chrom.values
        g["chrom"], g["pos"] = to_chrom(d, g.seq, g.pos)
        d["good"] = g
    p = row["annotation_gff"]
    if have(p):
        g = pd.read_csv(p, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4],
                        names=["chrom", "type", "start", "end"], dtype={0: str}, low_memory=False)
        g = g[(g.type == "gene") & g.chrom.isin(seqs)][["chrom", "start", "end"]].copy()
        s0 = g.chrom.values
        g["chrom"], g["start"] = to_chrom(d, s0, g.start)
        _, g["end"] = to_chrom(d, s0, g.end)
        d["genes"] = g

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

    # ---- dosage masks, crossovers
    p = "results/crossovers/%s/masks.tsv" % S
    if have(p):
        mk = pd.read_csv(p, sep="\t", dtype={"barcode": str, "chrom": str})
        if (~mk.chrom.isin(seqs)).any():
            WARN.append("%s: masks.tsv has %d rows for sequences not in the reference (left out): %s"
                        % (S, int((~mk.chrom.isin(seqs)).sum()), ", ".join(sorted(set(mk.chrom[~mk.chrom.isin(seqs)]))[:5])))
            mk = mk[mk.chrom.isin(seqs)].copy()
        d["masks"] = mk
        d["masked_pairs"] = set(zip(mk.barcode[mk.masked == "yes"], mk.chrom[mk.masked == "yes"]))
    p = "results/crossovers/%s/cells_per_chrom.tsv" % S
    if have(p):
        cpc = pd.read_csv(p, sep="\t", dtype={"chrom": str})
        d["nseq"] = {c: v for c, v in zip(cpc.chrom, cpc.n_cells) if c in seqs}
    p = "results/crossovers/%s/co_intervals_masked.bed" % S
    if have(p):
        d["masked_calls"] = sum(1 for l in open(p) if l.strip())
    p = "results/crossovers/%s/co_per_cell.tsv" % S
    if have(p):
        d["per"] = pd.read_csv(p, sep="\t")
    p = "results/crossovers/%s/co_intervals.bed" % S
    if have(p):
        co = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["seq", "start", "end", "barcode"],
                         dtype={0: str})
        co = co[co.seq.isin(seqs)].copy()
        s0 = co.seq.values
        co["chrom"], co["start"] = to_chrom(d, s0, co.start)
        _, co["end"] = to_chrom(d, s0, co.end)
        co["mid"] = (co.start + co.end) / 2
        co["u"] = co.mid / co.chrom.map(d["L"]).astype(float)
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

    # ---- regions left out of the molecule tables (config/exclude_regions; optional)
    p = "config/exclude_regions/%s.csv" % S
    if have_opt(p):
        d["excl"] = list(csv.DictReader(l for l in open(p) if not l.startswith("#")))
        d["excl_file"] = p
        lg = "logs/cell_data_molecules/%s.log" % S
        if have(lg):
            m = re.search(r"excluded regions \((.*?)\): (.*?) -- ([\d,]+) molecules dropped", open(lg).read())
            if m:
                d["excl_dropped"] = int(m.group(3).replace(",", ""))
                want = sorted("%s:%d-%d" % (r["chrom"], int(r["start"]), int(r["end"])) for r in d["excl"])
                if sorted(x.strip() for x in m.group(2).split(",")) != want:
                    WARN.append("%s: the molecule tables left out %s, but %s lists %s: rebuild the tables"
                                % (S, m.group(2), p, ", ".join(want)))
    D[sp] = d
B, Pd = D["binata"], D["paradoxa"]


def consistency(sp, d):
    """The run on disk against the config this report quotes: hapCO settings (co_calling_params.txt), dosage masks
    (masks.tsv) and regions left out (the molecule-table log). A mismatch is listed in red at the top of the page.
    Also the Python port of hapCO against the pipeline's own calls (co_sweep.py --check log), where it was run."""
    S, c = d["S"], CC.get(d["species"]) or {}
    p = "results/crossovers/%s/co_calling_params.txt" % S
    if have(p):
        t = open(p).read()
        got = {}
        for k, pat in (("block_size", r"block_size\s+\(-s\):\s*(\S+)"), ("marker_num", r"marker_num\s+\(-n\):\s*(\S+)"),
                       ("base_af", r"baseAF:\s*(\S+)"), ("window_af", r"windowAF:\s*(\S+)"),
                       ("genotype", r"genotype:\s*(\S+)"), ("terminal_marker_num", r"terminal_marker_num:\s*(\S+)")):
            m = re.search(pat, t)
            if m:
                got[k] = m.group(1)
        d["ran"] = got
        for k, v in got.items():
            try:
                same = float(v) == float(c[k])
            except (KeyError, TypeError, ValueError):
                same = k not in c
            if not same:
                WARN.append("%s: crossovers were called with %s = %s, but the config (quoted in this report) says %s"
                            % (S, k, v, c[k]))
    if "masks" in d:
        mc = set(c.get("mask_classes", []) or [])
        mm = d["masks"]
        bad = int((mm["class"].isin(mc) != (mm.masked == "yes")).sum())
        if bad:
            WARN.append("%s: %d rows of masks.tsv are masked differently from config mask_classes (%s): masks from "
                        "an older run?" % (S, bad, ", ".join(sorted(mc)) or "none"))
    if d.get("excl") and "excl_dropped" not in d:
        WARN.append("%s: config/exclude_regions/%s.csv lists %s, but logs/cell_data_molecules/%s.log has no "
                    "'excluded regions' line: were the molecule tables rebuilt with it?"
                    % (S, S, ", ".join("%s:%s-%s" % (r["chrom"], r["start"], r["end"]) for r in d["excl"]), S))
    p = "logs/%s_sweep_check.log" % S
    if have_opt(p):
        t = open(p).read()
        m = re.search(r"identical in (\d+) of (\d+) cells", t)
        m2 = re.search(r"at the settings it ran with \(([^;)]+)", t)
        if m:
            d["port_check"] = (int(m.group(1)), int(m.group(2)), m2.group(1) if m2 else "?")
            try:
                chosen = sid_of(c["base_af"], c["marker_num"], c["block_size"], c["terminal_marker_num"])
            except (KeyError, TypeError, ValueError):
                chosen = None
            if d["port_check"][0] != d["port_check"][1] or (chosen and d["port_check"][2] != chosen):
                WARN.append("%s: the Python port matched the R's calls in %s of %s nuclei at %s (chosen: %s; %s)"
                            % (S, m.group(1), m.group(2), d["port_check"][2], chosen, p))


for sp, d in D.items():
    guarded("consistency checks %s" % sp, lambda sp=sp, d=d: consistency(sp, d))


# ============================================================ window genotypes per cell
def window_genotypes(d):
    """WIN_MB windows of every reference sequence: per nucleus 1 (ALT) or 0 (REF) where >= MIN_MOL molecules agree
    >= 80%, else NaN; NaN also where the sequence is masked for dosage in that nucleus."""
    S, W = d["S"], int(WIN_MB * 1e6)
    nwin = {c: max(1, int(round(d["Ls"][c] / W))) for c in d["seqs"]}
    wins = [(c, i) for c in d["seqs"] for i in range(nwin[c])]
    widx = {w: k for k, w in enumerate(wins)}
    masked = d.get("masked_pairs", set())
    G = np.full((len(d["cells"]), len(wins)), np.nan)
    nmol = {}
    for ci, bc in enumerate(d["cells"]):
        p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            continue
        m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"],
                        dtype={0: str})
        m = m[(m.rc != m.ac) & m.chrom.isin(d["Ls"])].copy()
        m["alt"] = (m.ac > m.rc).astype(int)
        nmol[bc] = len(m)
        if not len(m):
            continue
        m["w"] = [min(int(q // W), nwin[c] - 1) for c, q in zip(m.chrom, m.pos)]
        g = m.groupby(["chrom", "w"]).alt.agg(["mean", "size"])
        for (c, w), r in g.iterrows():
            if (bc, c) in masked:
                continue
            if r["size"] >= MIN_MOL and (r["mean"] >= GT_MAJ - 1e-9 or 1 - r["mean"] >= GT_MAJ - 1e-9):
                G[ci, widx[(c, w)]] = 1 if r["mean"] >= GT_MAJ - 1e-9 else 0
    USED.append("results/cell_data_mol/%s/<barcode>.tsv (%d cells)" % (S, len(nmol)))
    d["mol_per_cell"] = pd.Series(nmol, dtype=float)
    d["Gw"], d["wins"], d["nwin"], d["widx"] = G, wins, nwin, widx
    sim, keep, dup, iu = near_identical(G)
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
    seen, rep = set(), []
    for a in range(len(d["cells"])):
        if root(a) not in seen:
            seen.add(root(a))
            rep.append(a)
    d["rep_idx"] = np.array(rep, dtype=int)
    d["distinct_genotypes"] = len(rep)
    d["called_share"] = float((~np.isnan(G)).mean()) if G.size else float("nan")


def near_identical(G):
    """Pairs of nuclei whose window genotypes agree in >= DUP_AGREE of >= DUP_MIN windows called in both."""
    M = (~np.isnan(G)).astype(float)
    X, Y = np.nan_to_num(G) * M, (1 - np.nan_to_num(G)) * M
    shared, agree = M @ M.T, X @ X.T + Y @ Y.T
    iu = np.triu_indices(len(G), 1)
    sh, ag = shared[iu], agree[iu]
    keep = sh >= DUP_MIN
    sim = np.full(len(sh), np.nan)
    sim[keep] = ag[keep] / sh[keep]
    return sim, keep, keep & (sim >= DUP_AGREE), iu


for sp, d in D.items():
    guarded("window genotypes %s" % sp, lambda d=d: window_genotypes(d))


# ======================================================= pollen co-segregation (phi)
def phi_matrix(G):
    """phi between every two windows over the nuclei with a genotype in both (NaN below PHI_N nuclei): 2 x 2 counts
    from matrix products, phi = (n11 n00 - n10 n01) / sqrt(n1. n0. n.1 n.0)."""
    X1, X0 = (G == 1).astype(float), (G == 0).astype(float)
    n11, n00, n10, n01 = X1.T @ X1, X0.T @ X0, X1.T @ X0, X0.T @ X1
    n = n11 + n00 + n10 + n01
    den = np.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    with np.errstate(divide="ignore", invalid="ignore"):
        ph = (n11 * n00 - n10 * n01) / den
    ph[(den == 0) | (n < PHI_N)] = np.nan
    return ph, n


def relations(d):
    """What the assembly says about two reference sequences: the same chromosome (cut for calling), or a split-plan
    note that a piece's copy in the other haplotype lies on the other's chromosome. Returns {(a, b): text}."""
    rel = {}
    for a in d["seqs"]:
        for b in d["seqs"]:
            if a >= b:
                continue
            ca, cb = d["src"][a][0], d["src"][b][0]
            if ca == cb:
                rel[(a, b)] = "one chromosome in the assembly (%s), cut for crossover calling" % ca
                continue
            notes = []
            for x, y in ((a, b), (b, a)):
                m = re.search(r"homolog (chr\d+)_hap\d", d["how"].get(x, ""))
                if m and m.group(1) != short(d["src"][x][0]) and m.group(1) == short(d["src"][y][0]):
                    notes.append("the other haplotype's copy of %s is part of %s (split plan)" % (label(d, x), m.group(1)))
            if notes:
                rel[(a, b)] = "; ".join(notes)
    return rel


def cosegregation(d, one_per_group=False):
    """Which reference sequences are inherited together? phi between windows over the nuclei (one per group of
    near-identical nuclei when those exceed chance: the shuffles assume independent nuclei); per pair of sequences
    the strongest |phi| over their window pairs, against NPERM shuffles of the nuclei done independently per
    sequence (each sequence keeps its own linkage, any between sequences is broken). Chance level: 95th percentile
    of the strongest pair per shuffle (family-wise)."""
    rng = stream("cosegregation", d["S"], one_per_group)
    G = d["Gw"][d["rep_idx"]] if one_per_group else d["Gw"]
    G = G[~np.all(np.isnan(G), axis=1)]
    d["phi_one_per_group"] = one_per_group
    seqs, wins = d["seqs"], d["wins"]
    blk = {s: np.array([k for k, (c, _) in enumerate(wins) if c == s]) for s in seqs}
    ph, nn = phi_matrix(G)
    ns = len(seqs)
    iu = np.triu_indices(ns, 1)

    def pair_max(P):
        M = np.full((ns, ns), np.nan)
        for a, b in zip(*iu):
            x = np.abs(P[np.ix_(blk[seqs[a]], blk[seqs[b]])])
            if np.isfinite(x).any():
                M[a, b] = M[b, a] = np.nanmax(x)
        return M
    obs = pair_max(ph)
    perm = np.full((NPERM, ns, ns), np.nan)
    for k in range(NPERM):
        Gp = np.empty_like(G)
        for s in seqs:
            Gp[:, blk[s]] = G[np.ix_(rng.permutation(len(G)), blk[s])]
        perm[k] = pair_max(phi_matrix(Gp)[0])
    fw = np.nanmax(perm[:, iu[0], iu[1]], axis=1)
    thr = float(np.nanpercentile(fw, 95))
    pair_thr = np.nanpercentile(perm, 95, axis=0)
    rel = relations(d)
    # sequences connected through what the assembly says (a path of relations) explain each other
    adj = collections.defaultdict(set)
    for a, b in rel:
        adj[a].add(b)
        adj[b].add(a)

    def path(a, b):
        prev, todo = {a: None}, [a]
        while todo:
            x = todo.pop(0)
            for y in adj[x]:
                if y not in prev:
                    prev[y] = x
                    todo.append(y)
        if b not in prev:
            return None
        out, x = [], b
        while x is not None:
            out.append(x)
            x = prev[x]
        return out[::-1]
    rows = []
    W = WIN_MB * 1e6
    for a, b in zip(*iu):
        if not (np.isfinite(obs[a, b]) and obs[a, b] > thr):
            continue
        sa, sb = seqs[a], seqs[b]
        sub = ph[np.ix_(blk[sa], blk[sb])]
        i, j = np.unravel_index(np.nanargmax(np.abs(sub)), sub.shape)
        key = (min(sa, sb), max(sa, sb))
        pth = path(sa, sb)
        if key in rel:
            why, expl = rel[key], True
        elif pth:
            why, expl = "through %s" % " and ".join(label(d, x) for x in pth[1:-1]), True
        else:
            why, expl = "different chromosomes in the assembly: not explained", False
        span = lambda s, k: ((d["src"][s][1] + k * W) / 1e6,
                             (d["src"][s][1] + ((k + 1) * W if k < d["nwin"][s] - 1 else d["Ls"][s])) / 1e6)
        rows.append(dict(a=sa, b=sb, phi=float(sub[i, j]), absphi=float(abs(sub[i, j])),
                         n=int(nn[blk[sa][i], blk[sb][j]]), wa=span(sa, i), wb=span(sb, j),
                         chance=float(pair_thr[a, b]), as_strong=int(np.sum(fw >= obs[a, b])),
                         p=(1 + int(np.sum(fw >= obs[a, b]))) / (1.0 + NPERM), why=why, explained=expl))
    rows.sort(key=lambda r: -r["absphi"])
    # pairs the assembly says belong together (same chromosome, or a split-plan homolog note) but not above chance
    missing = []
    for (x, y), why in sorted(rel.items()):
        a, b = seqs.index(x), seqs.index(y)
        if not (np.isfinite(obs[a, b]) and obs[a, b] > thr):
            missing.append(dict(a=x, b=y, absphi=float(obs[a, b]), why=why))
    d.update(phi=ph, phi_n=nn, pair_obs=obs, pair_thr=pair_thr, fw_thr=thr, co_pairs=rows, phi_nuclei=len(G),
             phi_blk=blk, unexplained=[r for r in rows if not r["explained"]], expected_missing=missing,
             expected_pairs=len(rel))


def dup_chance(d, nshuf=NSHUF_DUP):
    """Near-identical pairs among unrelated nuclei: the nuclei shuffled independently per linkage group (sequences
    inherited together above chance form one group), which keeps how the genome is linked and breaks any
    relationship between nuclei; the same rule counted on each shuffle."""
    G, rng = d["Gw"], stream("near-identical", d["S"])
    par = {s: s for s in d["seqs"]}

    def f(x):
        while par[x] != x:
            x = par[x]
        return x
    for r in d.get("co_pairs", []):
        par[f(r["a"])] = f(r["b"])
    groups = collections.defaultdict(list)
    for k, (s, _) in enumerate(d["wins"]):
        groups[f(s)].append(k)
    cnt = []
    for _ in range(nshuf):
        Gp = np.empty_like(G)
        for cols in groups.values():
            cols = np.array(cols)
            Gp[:, cols] = G[np.ix_(rng.permutation(len(G)), cols)]
        cnt.append(int(near_identical(Gp)[2].sum()))
    cnt = np.array(cnt)
    d["dup_chance"] = (float(np.mean(cnt)), float(np.percentile(cnt, 2.5)), float(np.percentile(cnt, 97.5)))
    d["dup_p"] = (1 + int(np.sum(cnt >= d["dup_pairs"]))) / (1.0 + nshuf)
    d["link_groups"] = len(groups)


for sp, d in D.items():
    if "Gw" in d:
        guarded("pollen co-segregation %s" % sp, lambda d=d: cosegregation(d))
        guarded("near-identical nuclei by chance %s" % sp, lambda d=d: dup_chance(d))
        if d.get("dup_p", 1.0) < 0.05:          # siblings present: the shuffles assume independent nuclei
            guarded("pollen co-segregation %s, one nucleus per group" % sp,
                    lambda d=d: cosegregation(d, one_per_group=True))


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
    # at each cut, the stretch between the nearest good markers on either side: a crossover there cannot be called
    d["cut_gaps"] = []
    for c, xs in d["cuts"].items():
        pos = np.sort(gm.pos[gm.chrom == c].values)
        for x in xs:
            k = np.searchsorted(pos, x, side="right")
            a = pos[k - 1] if k > 0 else 0.0
            b = pos[k] if k < len(pos) else float(d["L"][c])
            d["cut_gaps"].append((c, x, a, b))


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
    nseq = {s: int(d.get("nseq", {}).get(s, ncell)) for s in d["seqs"]}
    d["nseq_used"] = nseq
    d.update(ncell_co=ncell, co_total=len(co), co_mean=per.n_cos.mean(), co_sd=per.n_cos.std(),
             co_median=per.n_cos.median())
    # per sequence: crossovers per nucleus that counts for it (masked nuclei do not); per chromosome: sum of its pieces
    d["co_seq"] = {s: int((co.seq == s).sum()) for s in d["seqs"]}
    d["per_seq"] = {s: d["co_seq"][s] / max(nseq[s], 1) for s in d["seqs"]}
    d["per_chrom"] = {c: sum(d["per_seq"][s] for s in pieces_of(d, c)) for c in d["main"]}
    d["map_cm"] = 100 * sum(d["per_seq"].values())
    d["below_half"] = [c for c, v in d["per_chrom"].items() if v < 0.5]
    d["below_obligate"] = int((per.n_cos < d["n"] / 2).sum())
    # with exactly one crossover per bivalent, each grain gets each crossover with p 1/2: Binomial(n, 1/2) per grain
    d["below_obligate_exp"] = ncell * sum(math.comb(d["n"], k) for k in range(d["n"] + 1) if k < d["n"] / 2) / 2.0 ** d["n"]
    d["width_median_kb"] = float((co.end - co.start).median()) / 1e3
    cidx = {b: k for k, b in enumerate(per.barcode)}
    rows = co.barcode.map(cidx)
    okr = rows.notna().values
    bix = np.minimum((co.u.values * BINS).astype(int), BINS - 1)
    C = np.zeros((ncell, BINS))
    np.add.at(C, (rows[okr].astype(int).values, bix[okr]), 1)
    d["land"] = profile(C, stream("landscape", d["S"]))
    d["shape"] = shape(d["land"]["U"], d["land"]["U_lo"], d["land"]["U_hi"])
    if "good" in d:
        mu = (d["good"].pos / d["good"].chrom.map(d["L"]).astype(float)).values
        h = np.histogram(mu, bins=BINS, range=(0, 1))[0].astype(float)
        d["mprof"], d["U_markers"] = h / h.mean(), uidx(h)
    if "genes" in d:
        g = d["genes"]
        h = np.histogram((((g.start + g.end) / 2) / g.chrom.map(d["L"]).astype(float)).values, bins=BINS,
                         range=(0, 1))[0]
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
        d["busiest"].append([c, "%d–%d" % (w * 5, w * 5 + 5), a[i], ff(mu), p_txt(poisson_logtail(a[i], mu))])
    d["window_mean"] = mu


for sp, d in D.items():
    if "per" in d and "co" in d:
        guarded("crossover statistics %s" % sp, lambda d=d: co_stats(d))
    d.setdefault("shape", "–")
HEAD = {sp: d["shape"] for sp, d in D.items()}

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


def structure(sp, d):
    """hap2-on-hap1 synteny blocks (translocation_map.py) on each chromosome of this reference: where its
    other-haplotype copy aligns collinear, inverted, from another chromosome, or not at all."""
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
        put(sp, "struct_mb." + k, sum(covered(st[c][k], 0, d["L"][c]) for c in d["main"]) / 1e6)


def dead_spots(d, w=5e6):
    """Runs of 5 Mb windows without a crossover and what they overlap. Only runs long enough that fewer than one is
    expected by chance (same crossovers per chromosome, placed uniformly) are kept."""
    co, rng = d["co"], stream("dead spots", d["S"])
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
            cut = any(a < x < b for x in d["cuts"].get(c, []))
            cause = ("inverted in the other haplotype" if sh.get("INVERTED", 0) >= 0.5 else
                     "no homolog in the other haplotype" if ok(none) and none >= 0.5 else
                     "translocated segment" if sh.get("TRANSLOCATED", 0) >= 0.5 else
                     "few markers" if m < 0.25 else "not explained by structure or markers")
            rows.append([c, "%.0f–%.0f" % (a / 1e6, b / 1e6), ff(ln / 1e6, 0),
                         ff(100 * sh.get("INVERTED", float("nan")), 0), ff(100 * sh.get("TRANSLOCATED", float("nan")), 0),
                         ff(100 * none, 0), ff(m, 2), "yes" if cut else "", cause])
    d["dead"], d["dead_thr"] = rows, thr
    d["dead_chance"] = (float(chance.mean()), np.percentile(chance, 2.5), np.percentile(chance, 97.5))
    d["dead_mb"] = sum(float(r[2]) for r in rows)


for sp, d in D.items():
    if sp == "paradoxa":
        guarded("structure %s" % sp, lambda sp=sp, d=d: structure(sp, d))
    if "co" in d:
        guarded("dead spots %s" % sp, lambda d=d: dead_spots(d))


def interference(d):
    """Adjacent crossovers on one chromosome of one nucleus: their distance as a share of the chromosome,
    against the same numbers of crossovers placed independently (no interference)."""
    co, rng = d["co"], stream("interference", d["S"])
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


for sp, d in D.items():
    if "co" in d:
        guarded("interference %s" % sp, lambda d=d: interference(d))


# ============================================================ the review behind the settings
def review(sp, d):
    """The review of the calls (docs/co_review) and the scores of every setting against it (qc/co_review)."""
    S, rv = d["S"], {}
    c = CC.get(d["species"]) or {}
    try:
        rv["chosen"] = sid_of(c["base_af"], c["marker_num"], c["block_size"], c["terminal_marker_num"])
    except (KeyError, TypeError, ValueError):
        rv["chosen"] = None
    a = str(c.get("_approved", ""))
    rv["stamp_date"] = a.split("|")[0].strip()
    p = "docs/co_review/%s_settings.tsv" % S
    if have(p):
        st = list(csv.DictReader(open(p), delimiter="\t"))
        rv["n_settings"] = len(st)
        rv["previous"] = next((r["id"] for r in st if r.get("label") == "current"), None)
        rv["grid"] = {k: sorted({float(r[k]) for r in st}) for k in ("base_af", "marker_num", "block_size",
                                                                      "terminal_marker_num")}
    p = "docs/co_review/%s_events.tsv.gz" % S
    if have(p):
        with gzip.open(p, "rt") as fh:
            ev = list(csv.DictReader(fh, delimiter="\t"))
        rv["page_events"] = len(ev)
        rv["page_cells"] = len({e["barcode"] for e in ev})
    p = "docs/co_review/%s_verdicts.csv" % S
    if have(p):
        lines = open(p).read().splitlines()
        head = [l for l in lines if l.startswith("#")]
        m = re.search(r"page (\w+)", head[0]) if head else None
        rv["page"] = m.group(1) if m else "?"
        rv["notes"] = head[1:]
        rows = list(csv.DictReader(l for l in lines if not l.startswith("#")))
        ev = [r for r in rows if r["type"] == "event"]
        masked = {(r["barcode"], r["chrom"]) for r in rows if r["type"] == "mask" and r["verdict"] == "masked"}
        bad = {r["barcode"] for r in rows if r["type"] == "cell" and r["verdict"] == "bad"}
        out = [r for r in ev if (r["barcode"], r["chrom"]) in masked or r["barcode"] in bad]
        sc = [r for r in ev if not ((r["barcode"], r["chrom"]) in masked or r["barcode"] in bad)]
        rv["events"] = len(ev)
        rv["cells"] = len({r["barcode"] for r in ev} | {r["barcode"] for r in rows if r["type"] == "missed"})
        rv["eye"] = collections.Counter(r["verdict"] for r in sc if r["source"] == "user" and r["verdict"])
        rv["auto"] = collections.Counter(r["verdict"] for r in sc if r["source"] == "auto" and r["verdict"])
        rv["none"] = sum(1 for r in sc if not r["verdict"])
        rv["on_masked"] = len(out)
        rv["masked_seqs"] = len(masked)
        rv["missed"] = sum(1 for r in rows if r["type"] == "missed")
        rv["bad_cells"] = len(bad)
    ps = "qc/co_review/%s/scores.tsv" % S
    if have(ps):
        t = pd.read_csv(ps, sep="\t").drop_duplicates("setting")
        t["par"] = t.setting.map(sid_parse)
        rv["scores"] = t.set_index("setting", drop=False)
        rv["scores_time"] = datetime.datetime.fromtimestamp(os.path.getmtime(ps)).strftime("%Y-%m-%d %H:%M")
        src = t["source"].value_counts().to_dict() if "source" in t else {}
        rv["mode"] = ("every setting re-called from the pipeline's tables (score_verdicts.py --recall_all)"
                      if src and set(src) == {"re-called"} else
                      "the page's own calls for the %s swept settings; %s outside the sweep re-called"
                      % (fi(src.get("sweep", 0)), fi(src.get("re-called", 0))))
    p = "qc/co_review/%s/scores.txt" % S
    if have(p):
        t = open(p).read()
        m = re.search(r"^RE-CALLED\s+(.*)$", t, re.M)
        rv["recalled"] = m.group(1).strip() if m else ""
        m = re.search(r"without the molecules in (.*?);", rv["recalled"])
        rv["recall_excl"] = m.group(1) if m else ""
        m = re.search(r"^table:\s*(\S+)", t, re.M)
        if m and os.path.realpath(m.group(1)) != os.path.realpath(ps):
            WARN.append("%s: %s was written with %s, not with the table this report reads (%s)" % (S, p, m.group(1), ps))
        m = re.search(r"^VERDICTS\s+.*?: (\d+) events in (\d+) (?:reviewed )?cells", t, re.M)
        if m and "events" in rv and (int(m.group(1)), int(m.group(2))) != (rv["events"], rv["cells"]):
            WARN.append("%s: %s scored %s events in %s nuclei, the verdicts file has %s in %s: scores from other "
                        "verdicts?" % (S, p, m.group(1), m.group(2), rv["events"], rv["cells"]))
    d["rv"] = rv


for sp, d in D.items():
    guarded("review %s" % sp, lambda sp=sp, d=d: review(sp, d))


def score_row(d, sid):
    t = d.get("rv", {}).get("scores")
    return t.loc[sid] if t is not None and sid in t.index else None


# =================================================================== numbers out
for sp, d in D.items():
    for k in ("raw", "filtered", "seen", "good_n", "good_per_mb", "clusters", "gap_median_kb", "gap_p95_kb",
              "gap_max_mb", "gap_1mb_share", "good_in_genes", "rdna_bp", "rdna_markers_est", "raw_k",
              "raw_share_cand", "raw_share_dp_low", "raw_share_dp_high", "raw_share_qual_fail", "raw_share_ab_fail",
              "raw_share_pass", "pair_sim_median", "dup_pairs", "dup_cells", "distinct_genotypes", "called_share",
              "ncell_co", "co_total", "co_mean", "co_sd", "co_median", "map_cm", "below_obligate", "below_obligate_exp", "masked_calls",
              "width_median_kb", "U_markers", "U_genes", "rho_markers", "window_mean", "shape", "asm2n", "dead_mb",
              "phi_nuclei", "fw_thr", "excl_dropped", "link_groups"):
        if k in d:
            put(sp, k, d[k])
    if "dup_chance" in d:
        put(sp, "dup_pairs_chance", "%.1f (95%% %.0f–%.0f), p %.3g" % (d["dup_chance"] + (d["dup_p"],)))
    if "land" in d:
        for k in ("U", "U_lo", "U_hi"):
            put(sp, "land.%s" % k, d["land"][k])
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
    for s, v in d.get("per_seq", {}).items():
        put(sp, "cos_per_nucleus.%s" % s, v)
        put(sp, "nuclei_counting.%s" % s, d["nseq_used"][s])
    for c, v in d.get("per_chrom", {}).items():
        put(sp, "cos_per_grain_chromosome.%s" % c, v)
    if "masks" in d:
        for k, v in d["masks"]["class"].value_counts().items():
            put(sp, "dosage_class." + k, v)
        put(sp, "dosage_masked_pairs", len(d.get("masked_pairs", ())))
    if "dead" in d:
        put(sp, "dead_spots", len(d["dead"]))
        put(sp, "dead_spots_chance", "%.1f (95%% %.0f–%.0f)" % d["dead_chance"])
        put(sp, "dead_spots_min_mb", d["dead_thr"])
        for r in d["dead"]:
            put(sp, "dead.%s:%s" % (r[0], r[1]), "%s Mb, inv %s%%, transl %s%%, no homolog %s%%, markers %s, %s"
                % (r[2], r[3], r[4], r[5], r[6], r[8]))
    for r in d.get("busiest", []):
        put(sp, "busiest_5mb.%s:%s" % (r[0], r[1]), "%d crossovers, p %s" % (r[2], r[4]))
    if "intf" in d:
        f = d["intf"]
        put(sp, "interference_pairs", f["pairs"])
        put(sp, "interference_short_ratio", "%.2f (null range %.2f–%.2f)" % (f["ratio"], f["ratio_lo"], f["ratio_hi"]))
    for r in d.get("co_pairs", []):
        put(sp, "cosegregating.%s|%s" % (r["a"], r["b"]), "phi %+.2f at %.0f / %.0f Mb, %d nuclei, chance %.2f, p %.3g; %s"
            % (r["phi"], r["wa"][0], r["wb"][0], r["n"], r["chance"], r["p"], r["why"]))
    for c, x, a, b in d.get("cut_gaps", []):
        put(sp, "cut.%s:%.3f_Mb" % (c, x / 1e6), "no good marker over %.0f kb (%.3f-%.3f Mb)" % ((b - a) / 1e3, a / 1e6,
                                                                                               b / 1e6))
    rv = d.get("rv", {})
    for k in ("chosen", "previous", "page", "events", "cells", "none", "on_masked", "missed", "n_settings"):
        if k in rv:
            put(sp, "review." + k, rv[k])
    for k in ("eye", "auto"):
        for v, n_ in rv.get(k, {}).items():
            put(sp, "review.%s_%s" % (k, v), n_)
    for lab in ("chosen", "previous"):
        r = score_row(d, rv.get(lab))
        if r is not None:
            put(sp, "score_%s" % lab, "%s F1 %.4f (%.4f-%.4f) TP %d FP %d FN %d calls/cell %.2f %s"
                % (r.setting, r.F1, r.F1_lo, r.F1_hi, r.TP, r.FP, r.FN, r.cos_per_cell, r.source))

# ======================================================================== figures


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return p


def cut_txt(d, c):
    return " | ".join(label(d, s).split(" ", 1)[-1] for s in pieces_of(d, c)) if len(pieces_of(d, c)) > 1 else ""


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
            for x in d["cuts"].get(c, []):
                ax.axvline((off + x) / 1e6, color="#2C2C2A", lw=0.7, ls=":")
            ax.text((off + d["L"][c] / 2) / 1e6, -0.08, short(c).replace("chr", ""),
                    transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=7)
            off += d["L"][c]
        ax.set_xlim(0, off / 1e6)
        ax.set_xticks([])
        ax.set_ylabel("good markers / Mb")
        ax.set_title("%s (%s): %s good markers%s" % (NAME[sp], d["S"], fi(d.get("good_n")),
                                                    "; dotted = cuts of the crossover reference" if d["pieces"] else ""),
                     loc="left", fontsize=9)
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
        ax.hist(v, bins=bins, color=COL[sp], alpha=0.85)
        ax.axvline(d["n"] / 2, color="#444441", ls="--", lw=1)
        ax.set_xlabel("crossovers per pollen grain (dashed: n/2 = %s, one crossover per bivalent on average)"
                      % ff(d["n"] / 2), fontsize=8)
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
        if "mprof" in d:
            ax.plot(x, d["mprof"], color="#888780", lw=1.2, label="good markers")
        if "gprof" in d:
            ax.plot(x, d["gprof"], color="#888780", lw=1, ls="--", label="genes")
        ax.axhline(1, color="#D3D1C7", lw=0.8)
        ax.set_xlabel("position along chromosome (% of length, all chromosomes pooled)")
        ax.set_title("%s: %s" % (NAME[sp], HEAD[sp]), loc="left", fontsize=8.5, color=COL[sp])
        ax.legend(frameon=False, fontsize=7, loc="upper center")
    axs[0].set_ylabel("density relative to the mean")
    fig.tight_layout()
    return save(fig, "fig_landscape_relative.png")


def fig_per_chrom():
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [max(len(D[s]["main"]), 7) for s in D]})
    for ax, (sp, d) in zip(axs, D.items()):
        cs = sorted(d["main"], key=lambda c: d["L"][c])
        for i, c in enumerate(cs):
            bottom = 0.0
            pcs = pieces_of(d, c)
            for k, s in enumerate(pcs):
                v = d["per_seq"][s]
                ax.bar(i, v, bottom=bottom, color=COL[sp], alpha=1.0 if k % 2 == 0 else 0.55, edgecolor="white", lw=0.6)
                if len(pcs) > 1 and v > 0.15:
                    ax.text(i, bottom + v / 2, label(d, s).split(" ", 1)[-1], ha="center", va="center", fontsize=6,
                            color="white")
                bottom += v
        ax.axhline(0.5, color="#444441", ls="--", lw=1)
        ax.set_xticks(range(len(cs)))
        ax.set_xticklabels(["%s\n%.0f" % (short(c).replace("chr", ""), d["L"][c] / 1e6) for c in cs], fontsize=6)
        ax.set_xlabel("chromosome (length, Mb), shortest to longest")
        ax.set_title("%s: crossovers per grain per chromosome%s" % (NAME[sp], "\n(its pieces stacked; each "
                     "per nucleus that counts for it)" if d["pieces"] else ""), loc="left", fontsize=9)
    axs[0].set_ylabel("mean crossovers per grain")
    fig.tight_layout()
    return save(fig, "fig_crossovers_per_chromosome.png")


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


PARS = [(1, "marker_num\n(molecules a block needs)"), (2, "block_size, Mb\n(shortest interior block)"),
        (3, "terminal_marker_num\n(molecules a chromosome-end block needs)"), (0, "base_af\n(smoothing threshold)")]


def fig_review():
    fig, axs = plt.subplots(1, 4, figsize=(12.5, 3.6), sharey=True)
    ymin, drawn = 1.0, 0
    for sp, d in D.items():
        rv = d.get("rv", {})
        t, ch = rv.get("scores"), sid_parse(rv.get("chosen"))
        if t is None or ch is None or score_row(d, rv.get("chosen")) is None:
            continue
        prev = score_row(d, rv.get("previous"))
        for ax, (k, lab) in zip(axs, PARS):
            m = [p is not None and all(p[i] == ch[i] for i in range(4) if i != k) for p in t.par]
            sl = t[m].copy()
            if not len(sl):
                continue
            sl["x"] = [p[k] for p in sl.par]
            sl = sl.sort_values("x")
            x, y = sl.x.values, sl.F1.values
            ax.plot(x, y, color=COL[sp], lw=1.2, zorder=2)
            ax.vlines(x, sl.F1_lo.values, sl.F1_hi.values, color=COL[sp], lw=0.8, alpha=0.6, zorder=1)
            un = sl.unreviewed_calls.values if "unreviewed_calls" in sl else np.zeros(len(sl))
            ax.scatter(x[un == 0], y[un == 0], s=18, color=COL[sp], zorder=3)
            ax.scatter(x[un > 0], y[un > 0], s=18, facecolor="white", edgecolor=COL[sp], zorder=3)
            for xx, yy, uu in zip(x, y, un):
                if uu > 0:
                    ax.annotate("+%d" % uu, (xx, yy), xytext=(3, -9), textcoords="offset points", fontsize=6,
                                color=COL[sp])
            ax.scatter([ch[k]], [score_row(d, rv["chosen"]).F1], s=70, facecolor=COL[sp], edgecolor="#2C2C2A",
                       lw=1.2, zorder=4)
            if prev is not None:
                ax.axhline(prev.F1, color=COL[sp], lw=0.9, ls="--", alpha=0.8)
            ymin = min(ymin, float(np.nanmin(sl.F1_lo.values)), float(prev.F1) if prev is not None else 1.0)
            ax.set_xlabel(lab, fontsize=8)
            ax._xs = sorted(set(x) | set(getattr(ax, "_xs", [])))
            drawn += 1
    if not drawn:
        plt.close(fig)
        return None
    for ax in axs:
        ax.set_xticks(getattr(ax, "_xs", []))
        ax.set_xticklabels(["%g" % v for v in getattr(ax, "_xs", [])], fontsize=7)
        ax.set_ylim(max(0.0, ymin - 0.02), 1.0)
        ax.grid(axis="y", color="#F1EFE8", lw=0.8)
    axs[0].set_ylabel("F1 against the review (95% CI)")
    handles = [plt.Line2D([0], [0], color=COL[sp], marker="o", lw=1.2, label=NAME[sp]) for sp in D] + [
        plt.Line2D([0], [0], color="#2C2C2A", marker="o", mfc="#888780", lw=0, ms=8, label="chosen settings"),
        plt.Line2D([0], [0], color="#888780", ls="--", lw=0.9, label="previous settings (all four as before)"),
        plt.Line2D([0], [0], color="#888780", marker="o", mfc="white", lw=0,
                   label="+n: calls that match no reviewed event (not scored)")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=7.5, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("Each setting moved away from the chosen value with the other three held, scored against the review's "
                 "verdicts", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    return save(fig, "fig_settings_from_review.png")


def fig_pairs():
    sps = [sp for sp in D if "pair_obs" in D[sp]]
    if not sps:
        return None
    fig, axs = plt.subplots(1, len(sps), figsize=(6.3 * len(sps), 6.1))
    axs = np.atleast_1d(axs)
    for ax, sp in zip(axs, sps):
        d = D[sp]
        M = d["pair_obs"]
        labs = [label(d, s) for s in d["seqs"]]
        n = len(labs)
        im = ax.imshow(np.ma.masked_invalid(np.where(M > d["fw_thr"], M, np.nan)), cmap=CMAP[sp], vmin=0, vmax=1,
                       interpolation="nearest")
        for i in range(n):
            ax.add_patch(plt.Rectangle((i - 0.5, i - 0.5), 1, 1, color="#B4B2A9", lw=0))
            for j in range(n):
                if i != j and np.isfinite(M[i, j]) and M[i, j] > d["fw_thr"]:
                    ax.text(j, i, "%.2f" % M[i, j], ha="center", va="center", fontsize=6.5 if n > 10 else 8,
                            color="white" if M[i, j] > 0.55 else "#2C2C2A")
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(labs, rotation=90, fontsize=7)
        ax.set_yticklabels(labs, fontsize=7)
        ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
        ax.grid(which="minor", color="white", lw=1.2)
        ax.tick_params(which="minor", length=0)
        for s_ in ax.spines.values():
            s_.set_visible(False)
        ax.set_title("%s: %d %s, %s nuclei\ncoloured: above chance (|φ| > %s, which no pair reaches in 95%% of %d "
                     "shuffles); pale: not" % (NAME[sp], n, "pieces of %d chromosomes" % len(d["main"]) if d["pieces"]
                                               else "chromosomes", fi(d["phi_nuclei"]), ff(d["fw_thr"], 2), NPERM),
                     loc="left", fontsize=8.5)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cb.set_label("strongest |φ| between their %g Mb windows (pairs above chance)" % WIN_MB, fontsize=7)
    fig.tight_layout()
    return save(fig, "fig_cosegregation_pairs.png")


def fig_phi_maps():
    sps = [sp for sp in D if "phi" in D[sp]]
    if not sps:
        return None
    fig, axs = plt.subplots(1, len(sps), figsize=(6.8 * len(sps), 7.0))
    axs = np.atleast_1d(axs)
    for ax, sp in zip(axs, sps):
        d = D[sp]
        P = np.abs(d["phi"]).copy()
        np.fill_diagonal(P, np.nan)
        nw = [d["nwin"][s] for s in d["seqs"]]
        ax.imshow(np.ma.masked_invalid(P), cmap=CMAP[sp], vmin=0, vmax=1, interpolation="nearest")
        edges = np.cumsum([0] + nw) - 0.5
        for k, e in enumerate(edges[1:-1]):
            same = d["src"][d["seqs"][k]][0] == d["src"][d["seqs"][k + 1]][0]
            for f_ in (ax.axhline, ax.axvline):
                f_(e, color="#2C2C2A" if same else "#888780", lw=0.7 if same else 0.4, ls=":" if same else "-")
        mids = (edges[:-1] + edges[1:]) / 2
        ax.set_xticks(mids)
        ax.set_xticklabels([label(d, s) for s in d["seqs"]], rotation=90, fontsize=6.5)
        ax.set_yticks(mids)
        ax.set_yticklabels([label(d, s) for s in d["seqs"]], fontsize=6.5)
        for s_ in ax.spines.values():
            s_.set_visible(False)
        ax.set_title("%s: every two %g Mb windows, %s nuclei%s" % (NAME[sp], WIN_MB, fi(d["phi_nuclei"]),
                     "\n(dotted: a cut between pieces of one chromosome)" if d["pieces"] else ""), loc="left", fontsize=9)
    fig.suptitle("Pollen co-segregation along the genome: |φ| between window genotypes (white 0, independent; dark 1, "
                 "always inherited together; grey: too few nuclei)", x=0.01, ha="left", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save(fig, "fig_cosegregation_windows.png")


def chrom_windows(d, c):
    """The window indices of chromosome c (its pieces in order) and their edges on the chromosome, Mb."""
    W = WIN_MB * 1e6
    idx, edges = [], [0.0]
    for s in pieces_of(d, c):
        o, k = d["src"][s][1], d["nwin"][s]
        for i in range(k):
            idx.append(d["widx"][(s, i)])
            edges.append((o + (i + 1) * W if i < k - 1 else o + d["Ls"][s]) / 1e6)
    return idx, np.array(edges)


def zoom_pairs(d):
    """Chromosome pairs with sequences inherited together above chance (not one chromosome): unexplained first."""
    out = {}
    for r in d.get("co_pairs", []):
        ca, cb = d["src"][r["a"]][0], d["src"][r["b"]][0]
        if ca == cb:
            continue
        k = tuple(sorted((ca, cb), key=natural))
        if k not in out or r["absphi"] > out[k]["absphi"]:
            out[k] = r
    return sorted(out.items(), key=lambda kv: (kv[1]["explained"], -kv[1]["absphi"]))[:2]


def fig_phi_zoom():
    panels = [(sp, k, r) for sp in D if "phi" in D[sp] for k, r in zoom_pairs(D[sp])]
    if not panels:
        return None
    fig, axs = plt.subplots(1, len(panels), figsize=(6.4 * len(panels), 5.8))
    axs = np.atleast_1d(axs)
    for ax, (sp, (c1, c2), r) in zip(axs, panels):
        d = D[sp]
        i1, e1 = chrom_windows(d, c1)
        i2, e2 = chrom_windows(d, c2)
        M = np.abs(d["phi"][np.ix_(i1, i2)])
        ax.pcolormesh(e2, e1, np.ma.masked_invalid(M), cmap=CMAP[sp], vmin=0, vmax=1)
        ax.set_xlim(0, e2[-1])
        ax.set_ylim(e1[-1], 0)
        for x in d["cuts"].get(c2, []):
            ax.axvline(x / 1e6, color="#2C2C2A", lw=0.9, ls=":")
        for y in d["cuts"].get(c1, []):
            ax.axhline(y / 1e6, color="#2C2C2A", lw=0.9, ls=":")
        ax.set_xlabel("%s, Mb%s" % (c2, " (dotted: cut %s)" % cut_txt(d, c2) if cut_txt(d, c2) else ""))
        ax.set_ylabel("%s, Mb%s" % (c1, " (dotted: cut %s)" % cut_txt(d, c1) if cut_txt(d, c1) else ""))
        ax.set_title("%s %s × %s: strongest |φ| %s\n%s" % (NAME[sp], short(c1), short(c2), ff(r["absphi"], 2), r["why"]),
                     loc="left", fontsize=8.5)
    fig.tight_layout()
    return save(fig, "fig_cosegregation_zoom.png")


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


def draw_hic(ax, H, nper, order, agg_bp, d, short_names=True):
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
    names = [(c.replace("_hap", " h") if short_names else c) for c in order]
    ax.set_xticks((edges[:-1] + edges[1:]) / 2)
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_yticks((edges[:-1] + edges[1:]) / 2)
    ax.set_yticklabels(names, fontsize=6)
    for i, c in enumerate(order):        # the cuts of the crossover reference, on the chromosomes it was cut from
        if d["cuts"].get(c) and d.get("L2n", {}).get(c) == d["L"].get(c):
            for x in d["cuts"][c]:
                p_ = edges[i] + x / 1e6
                ax.plot([p_, p_], [0, tot], color="#2C2C2A", lw=0.6, ls=":")
                ax.plot([0, tot], [p_, p_], color="#2C2C2A", lw=0.6, ls=":")
    ax.set_xlim(0, tot)
    ax.set_ylim(tot, 0)
    return im


def hic_rank(H, nper, order):
    """Contacts per bin pair between every two chromosome copies that are not homologs, against the median pair:
    {(a, b): (x the median, rank)}, number of pairs."""
    edges = np.cumsum([0] + list(nper))
    dens = {}
    for i in range(len(order)):
        for j in range(i + 1, len(order)):
            if short(order[i]) == short(order[j]):
                continue
            b = H[edges[i]:edges[i + 1], edges[j]:edges[j + 1]]
            dens[(order[i], order[j])] = float(b.sum()) / max(b.size, 1)
    med = float(np.median(list(dens.values()))) if dens else float("nan")
    rk = sorted(dens, key=lambda k: -dens[k])
    return {k: (v / med if med > 0 else float("nan"), rk.index(k) + 1) for k, v in dens.items()}, len(rk)


def hic_load(sp, name):
    d = D[sp]
    p = os.path.join(A.phd_root, "results", d["species"], "qc/hic_remap", name)
    if not (have(p) and "L2n" in d):
        return None
    z = np.load(p)
    return z["M"], [str(x) for x in z["chroms"]], int(z["bin"]), p


def fig_hic_dual(sp):
    """Hi-C on both haplotypes, all reads with MAPQ >= 1: reads that fit either copy of a sequence put contacts
    between the copies, so homologous sequence shows as a red line off the diagonal. Right: the chromosomes of
    any pair inherited together that the assembly does not explain."""
    d = D[sp]
    got = hic_load(sp, "contacts_all.npz") or hic_load(sp, "contacts_unique.npz")
    if not got:
        return None
    M, chroms, bin_bp, p = got
    order = sorted([c for c in chroms if re.match(r"^chr\d+_hap[12]$", c)],
                   key=lambda c: (int(re.search(r"\d+", c).group()), c[-1]))
    base = sorted({short(d["src"][x][0]) for r in d.get("unexplained", []) for x in (r["a"], r["b"])}, key=natural)
    focus = [c for c in order if short(c) in base]
    H1, n1, b1 = hic_agg(M, chroms, bin_bp, d["L2n"], order, 2e6)
    H2 = None
    if focus:
        H2, n2, b2 = hic_agg(M, chroms, bin_bp, d["L2n"], focus, 1e6)
    del M
    d["hic_rank"], d["hic_npairs"] = hic_rank(H1, n1, order)
    d["hic_src"], d["hic_focus"], d["hic_rank_src"] = p, base, os.path.basename(p).replace(".npz", "")
    pu = os.path.join(A.phd_root, "results", d["species"], "qc/hic_remap", "contacts_unique.npz")
    if os.path.basename(p) != "contacts_unique.npz" and have_opt(pu):   # rank on uniquely mapped reads
        z = np.load(pu)
        Hu, nu, _ = hic_agg(z["M"], [str(x) for x in z["chroms"]], int(z["bin"]), d["L2n"], order, 2e6)
        del z
        d["hic_rank"], d["hic_npairs"] = hic_rank(Hu, nu, order)
        d["hic_rank_src"] = "contacts_unique"
    ncol = 2 if H2 is not None else 1
    fig, axs = plt.subplots(1, ncol, figsize=(7.3 * ncol, 7.0), gridspec_kw={"width_ratios": [1.05, 1][:ncol]})
    axs = np.atleast_1d(axs)
    fig.subplots_adjust(left=0.06, right=0.92, top=0.9, bottom=0.1, wspace=0.18)
    im = draw_hic(axs[0], H1, n1, order, b1, d)
    axs[0].set_title("all %d chromosome copies, hap1 next to hap2 (%g Mb bins)" % (len(order), b1 / 1e6), loc="left",
                     fontsize=9)
    if H2 is not None:
        draw_hic(axs[1], H2, n2, focus, b2, d, short_names=False)
        axs[1].set_title("%s only, both haplotypes (%g Mb bins)" % (" and ".join(base), b2 / 1e6), loc="left",
                         fontsize=9)
    cb = fig.colorbar(im, ax=list(axs), fraction=0.02, pad=0.01)
    cb.set_label("log10 contacts (white = typical between two chromosomes)", fontsize=7)
    fig.suptitle("%s Hi-C on both haplotypes (%s)%s" % (NAME[sp], os.path.basename(p).replace(".npz", ""),
                 "; dotted = cuts of the crossover reference" if any(d["cuts"].values()) else ""),
                 x=0.06, y=0.97, ha="left", fontsize=9)
    return save(fig, "fig_hic_both_haplotypes_%s.png" % sp)


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
            ax.bar((np.arange(nb) + 0.5) * 5, 100.0 * cnt / bin_cells(d, c, nb, 5e6) / 5, width=5, color=COL[sp],
                   alpha=0.85)
            pos = d["good"].pos[d["good"].chrom == c].values
            mk = np.bincount(np.minimum((pos // 5e6).astype(int), nb - 1), minlength=nb)
            ax2 = ax.twinx()
            ax2.plot((np.arange(nb) + 0.5) * 5, mk / 5, color="#888780", lw=0.8)
            ax2.set_yticks([])
            ax2.spines["right"].set_visible(False)
            for x in d["cuts"].get(c, []):
                ax.axvline(x / 1e6, color="#2C2C2A", lw=0.8, ls=":")
            ax.set_title("%s %s%s" % (NAME[sp].split()[1], short(c), " (%s)" % cut_txt(d, c) if cut_txt(d, c) else ""),
                         fontsize=7, loc="left")
            ax.tick_params(labelsize=6)
        for j in range(len(d["main"]), nr * ncol):
            axs[r0 + j // ncol, j % ncol].axis("off")
        r0 += nr
    fig.text(0.0, 0.5, "cM/Mb (bars); good markers per Mb (grey, own scale); dotted = cut of the crossover reference",
             rotation=90, va="center", fontsize=8)
    fig.tight_layout()
    return save(fig, "fig_per_chromosome_landscapes.png")


def fig_structure(sp):
    d = D[sp]
    if "struct" not in d:
        return None
    n = len(d["main"])
    fig = plt.figure(figsize=(11, 2.15 * n + 1.0))
    gs = fig.add_gridspec(4 * n, 1, height_ratios=[3.0, 0.3, 0.42, 0.55] * n, hspace=0.12)
    pcol = ("#ECEAF8", "#D6D2F0")
    big = max(d["L"].values()) / 1e6
    colk = {"collinear": "#D3D1C7", "INVERTED": "#A32D2D", "TRANSLOCATED": "#185FA5"}
    co = d["co"]
    for i, c in enumerate(d["main"]):
        ax, axa, axs_ = fig.add_subplot(gs[4 * i]), fig.add_subplot(gs[4 * i + 1]), fig.add_subplot(gs[4 * i + 2])
        L = d["L"][c] / 1e6
        nb = int(math.ceil(d["L"][c] / 5e6))
        cnt = np.bincount(np.minimum((co.mid[co.chrom == c] // 5e6).astype(int), nb - 1), minlength=nb)
        y = 100.0 * cnt / bin_cells(d, c, nb, 5e6) / 5
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
        for x in d["cuts"].get(c, []):
            for a_ in (ax, axa, axs_):
                a_.axvline(x / 1e6, color="#2C2C2A", lw=0.9, ls=":")
        ax.set_xlim(0, big)
        ax.set_ylabel("cM/Mb", fontsize=6)
        ax.tick_params(labelsize=6, labelbottom=False)
        ax.set_title("%s%s" % (c, " = %s" % " | ".join(label(d, s) for s in pieces_of(d, c)) if cut_txt(d, c) else ""),
                     loc="left", fontsize=8, pad=2)
        for k_, s in enumerate(pieces_of(d, c)):
            a_, b_ = d["src"][s][1] / 1e6, (d["src"][s][1] + d["Ls"][s]) / 1e6
            axa.add_patch(plt.Rectangle((a_, 0), b_ - a_, 1, color=pcol[k_ % 2], lw=0))
            if cut_txt(d, c):
                axa.text((a_ + b_) / 2, 0.5, label(d, s), ha="center", va="center", fontsize=6.5)
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
    handles = [plt.Rectangle((0, 0), 1, 1, color=pcol[0]), plt.Rectangle((0, 0), 1, 1, color=pcol[1])] + \
        [plt.Rectangle((0, 0), 1, 1, color=colk[k]) for k in KINDS] + \
        [plt.Rectangle((0, 0), 1, 1, fc="white", ec="#B4B2A9"), plt.Line2D([0], [0], color="#2C2C2A", lw=3)]
    fig.legend(handles, ["piece of the crossover reference", "next piece", "homolog collinear", "homolog inverted",
                         "homolog from another chromosome", "no homolog alignment",
                         "no crossover over >= %d Mb" % d.get("dead_thr", 20)],
               loc="lower center", ncol=4, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.text(0.5, 0.035, "Mb. Per chromosome: crossovers (bars, cM/Mb; grey line = good markers), the pieces it was cut "
             "into for crossover calling (dotted = cut), and where the other haplotype's copy aligns (hap2 on hap1)",
             ha="center", fontsize=7)
    fig.subplots_adjust(top=0.985, bottom=0.075, left=0.07, right=0.98)
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


FIGS = {}
for name, fn in (("markers", fig_markers), ("raw", fig_raw), ("cells", fig_cells), ("per_grain", fig_per_grain),
                 ("landscape", fig_landscape), ("per_chrom", fig_per_chrom), ("interference", fig_interference),
                 ("review", fig_review), ("pairs", fig_pairs), ("phi_maps", fig_phi_maps), ("phi_zoom", fig_phi_zoom),
                 ("chrom_landscapes", fig_chrom_landscapes)):
    FIGS[name] = guarded("figure " + name, fn)
FIGS["structure"] = guarded("figure structure landscape", lambda: fig_structure("paradoxa"))
FIGS["hic_dual"] = guarded("figure Hi-C both haplotypes", lambda: fig_hic_dual("paradoxa"))
for sp, d in D.items():
    for (a, b), (x, rk) in sorted(d.get("hic_rank", {}).items()):
        if short(a) in d.get("hic_focus", []) and short(b) in d.get("hic_focus", []):
            put(sp, "hic_contacts.%s|%s" % (a, b), "%.2fx the median non-homologous pair, rank %d of %d"
                % (x, rk, d["hic_npairs"]))
EXAMPLES = guarded("example cells", example_cells) or {}
BROWSER = {}
if not A.no_cell_pdfs:
    for sp in D:
        BROWSER[sp] = guarded("cell browser " + sp, lambda sp=sp: cell_browser(sp))

# =========================================================================== html
MOL_GAP, = code_const("workflow/scripts/molecule_tables.py", r"\nGAP = (\d+)", 1)


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


def sel(d, k):
    return d.get("sel", {}).get(k, "–")


def setting_txt(sp):
    c = CC.get(D[sp]["species"]) or {}
    bs = ("%g" % (float(c["block_size"]) / 1e6)) if "block_size" in c else "–"
    return "marker_num %s, block_size %s Mb, terminal_marker_num %s, base_af %s" % (
        c.get("marker_num", "–"), bs, c.get("terminal_marker_num", "–"), c.get("base_af", "–"))


def f1_txt(d, lab, full=False):
    r = score_row(d, d.get("rv", {}).get(lab))
    if r is None:
        return "–"
    t = "F1 %s (95%% %s–%s)" % (ff(r.F1, 3), ff(r.F1_lo, 3), ff(r.F1_hi, 3))
    if full:
        t = "%s: %s; TP %s, FP %s, FN %s; precision %s, recall %s; %s calls per reviewed nucleus" % (
            r.setting, t, fi(r.TP), fi(r.FP), fi(r.FN), ff(r.precision, 3), ff(r.recall, 3), ff(r.cos_per_cell, 2))
    return t


def seq_desc(d):
    if d["pieces"]:
        return "%d pieces of %d chromosomes" % (len(d["seqs"]), len(d["main"]))
    return "%d chromosomes (%s)" % (len(d["main"]), d["row"]["haplotype"])


def excl_txt(d, dropped=True):
    if not d.get("excl"):
        return "none"
    t = ", ".join("%s:%s–%s" % (r["chrom"], fi(r["start"]), fi(r["end"])) for r in d["excl"])
    return t + (" (%s molecules dropped)" % fi(d.get("excl_dropped")) if dropped else "")


def status(sp):
    d, rv = D[sp], D[sp].get("rv", {})
    s = "settings stamped %s from the review of %s nuclei" % (rv.get("stamp_date") or "–", fi(rv.get("cells")))
    opens = []
    if d.get("excl"):
        opens.append("%s left out by hand, to be replaced by an automatic filter" % excl_txt(d, dropped=False))
    for r in d.get("unexplained", []):
        opens.append("%s and %s inherited together, not explained by the assembly" % (label(d, r["a"]), label(d, r["b"])))
    return s + ("; open: " + "; ".join(opens) + " (section 7)" if opens else "")


def pl(n, word, plural=None):
    """'1 pair', '3 pairs'"""
    return "%s %s" % (fi(n), word if ok(n) and float(n) == 1 else (plural or word + "s"))


def dup_txt(d):
    t = "%s among %s" % (pl(d.get("dup_pairs"), "pair"), pl(d.get("dup_cells"), "nucleus", "nuclei"))
    c = d.get("dup_chance")
    return t + ("; %s by chance (95%% %s–%s; p %s)" % (ff(c[0]), fi(c[1]), fi(c[2]), "%.3g" % d["dup_p"]) if c else "")


def dup_above_chance(sp):
    return D[sp].get("dup_p", 1.0) < 0.05


def dup_reading(sp):
    d, c = D[sp], D[sp].get("dup_chance")
    if not c or "dup_pairs" not in d:
        return "<i>%s</i> not tested" % NAME[sp]
    if dup_above_chance(sp):
        ex = d["dup_pairs"] - c[0]
        return ("<i>%s</i>: %s against %s by chance (p %s), an excess of about %s pairs: %s two nuclei of one grain"
                % (NAME[sp], pl(d["dup_pairs"], "pair"), ff(c[0]), "%.3g" % d["dup_p"], fi(ex),
                   "most are likely" if ex > d["dup_pairs"] / 2 else "some are likely"))
    return ("<i>%s</i>: %s, no more than unrelated nuclei give by chance (%s, p %s; the genome segregates as %d "
            "linkage groups here)" % (NAME[sp], pl(d["dup_pairs"], "pair"), ff(c[0]), "%.3g" % d["dup_p"],
                                       d["link_groups"]))


def sel_txt(d):
    s = d.get("sel", {})
    t = "haploidness ≥ %s over ≥ %s windows; ≥ %s molecules on the weakest chromosome" % (
        s.get("min_haploidness", "–"), s.get("min_windows", "–"), s.get("min_chrom_molecules", "–"))
    mm = s.get("min_molecules", "0")
    return t + ("; ≥ %s molecules in total" % mm if mm not in ("0", "0.0", "–") else "; no total-molecule floor")


def masked_classes(d):
    return ", ".join((CC.get(d["species"]) or {}).get("mask_classes", []) or ["none"])


H = []
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
commit = subprocess.run("git rev-parse --short HEAD", shell=True, stdout=subprocess.PIPE,
                        universal_newlines=True).stdout.strip()
GOOD_TXT = ("heterozygous SNPs that segregate cleanly in the pollen: ALT in %s–%s%% of the nuclei that cover them, "
            "and one allele per nucleus" % (ff(100 * GOOD_LO, 0), ff(100 * GOOD_HI, 0)))

H.append("<h1>Meiotic crossover landscapes of <i>Drosera binata</i> and <i>D. paradoxa</i> from single pollen nuclei</h1>")
H.append("<p class='meta'>Version 5, generated %s from CO_smk commit %s by workflow/scripts/supervisor_report.py "
         "(<i>D. binata</i>: %s; <i>D. paradoxa</i>: %s). Every number, table and figure is computed from the "
         "pipeline's outputs (files listed in the appendix).</p>" % (now, commit, B["S"], Pd["S"]))


def run_section(name, fn):
    """Build one section; if it fails, say so in red where it would have been and list it at the top."""
    try:
        fn()
    except Exception as e:
        PROBLEMS.append((name, "%s: %s" % (type(e).__name__, e)))
        traceback.print_exc()
        H.append("<p class='missing'>%s could not be built: %s</p>" % (html.escape(name), html.escape(str(e))))



# ---- 0
def section_0():
    H.append("<h2>0. Summary</h2>")

    def sp_line(sp):
        d = D[sp]
        ref = ("; crossovers called on a reference cut into %d pieces at the %d points where a chromosome's partner in "
               "the other haplotype changes, so no join between translocated pieces is assumed"
               % (len(d["seqs"]), sum(len(v) for v in d["cuts"].values()))) if d["pieces"] else ""
        return ("<b>%s</b> (%s%s): %s nuclei, %s crossovers per grain against %s from one crossover per bivalent "
                "(n/2), a genetic map of %s cM; landscape %s." % (
                    NAME[sp], d["row"]["centromere"], ref, fi(d.get("ncell_co")), ff(d.get("co_mean"), 2), ff(d["n"] / 2),
                    fi(d.get("map_cm")), HEAD[sp]))

    summ = ["Single pollen nuclei of <i>D. binata</i> and <i>D. paradoxa</i> were sequenced (scRNA-seq) and genotyped at "
            "heterozygous SNPs called from HiFi reads; a crossover is called where a nucleus switches between the two "
            "haplotypes.", sp_line("binata"), sp_line("paradoxa")]
    if all(score_row(D[sp], D[sp].get("rv", {}).get("chosen")) is not None for sp in D):
        summ.append("The crossover-calling settings of both were chosen by scoring every swept setting against the "
                    "verdicts of a review of the calls (judged by eye, or accepted automatically where nearly every "
                    "setting agrees) in %s and %s nuclei: %s and %s, against %s and %s for the settings used before "
                    "(section 5)." % (fi(B["rv"].get("cells")), fi(Pd["rv"].get("cells")),
                                                   f1_txt(B, "chosen"), f1_txt(Pd, "chosen"), f1_txt(B, "previous"),
                                                   f1_txt(Pd, "previous")))
    for sp, d in D.items():
        for r in d.get("unexplained", []):
            summ.append("Open in <i>%s</i>: %s and %s are inherited together in the pollen (|φ| up to %s between their "
                        "windows, chance level %s) although they are separate chromosomes in the assembly (section 7)."
                        % (NAME[sp], label(d, r["a"]), label(d, r["b"]), ff(r["absphi"], 2), ff(d["fw_thr"], 2)))
        if d.get("excl"):
            summ.append("<i>%s</i> %s is left out of crossover calling by hand, its alleles following another place in "
                        "the genome; an automatic filter is to replace the hand-picked list (section 7)."
                        % (NAME[sp], excl_txt(d, dropped=False)))
    H.append(P(" ".join(summ)))
    H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
        ["centromere type", B["row"]["centromere"], Pd["row"]["centromere"]],
        ["chromosomes (n)", B["n"], Pd["n"]],
        ["reference sequences crossovers are called on", seq_desc(B), seq_desc(Pd)],
        ["assembly, both haplotypes (2C), chromosomes, Mb", ff(B.get("asm2n_chr_mb"), 0), ff(Pd.get("asm2n_chr_mb"), 0)],
        ["reference used for crossovers (1C), Mb", ff(B["G"] / 1e6, 0), ff(Pd["G"] / 1e6, 0)],
        ["good markers (%s)" % GOOD_TXT, fi(B.get("good_n")), fi(Pd.get("good_n"))],
        ["pollen nuclei used", fi(len(B["cells"])), fi(len(Pd["cells"]))],
        ["informative molecules per nucleus, median", fi(B.get("selected_molecules_median")),
         fi(Pd.get("selected_molecules_median"))],
        ["crossovers per grain, mean (SD)", "%s (%s)" % (ff(B.get("co_mean"), 2), ff(B.get("co_sd"), 2)),
         "%s (%s)" % (ff(Pd.get("co_mean"), 2), ff(Pd.get("co_sd"), 2))],
        ["one crossover per bivalent gives, per grain on average (n / 2)", ff(B["n"] / 2), ff(Pd["n"] / 2)],
        ["genetic map length, cM (crossovers per nucleus summed over sequences, each over the nuclei that count for it)",
         fi(B.get("map_cm")), fi(Pd.get("map_cm"))],
        ["landscape", HEAD["binata"], HEAD["paradoxa"]],
        ["crossover-calling settings", "%s; %s" % (setting_txt("binata"), f1_txt(B, "chosen")),
         "%s; %s" % (setting_txt("paradoxa"), f1_txt(Pd, "chosen"))],
        ["status", status("binata"), status("paradoxa")]]))
    H.append(P("2n assemblies: <i>D. binata</i> <code>%s</code>; <i>D. paradoxa</i> <code>%s</code>."
               % (html.escape(str(B.get("asm2n", "not found"))), html.escape(str(Pd.get("asm2n", "not found"))))))


run_section('section 0', section_0)

# ---- 1
def section_1():
    H.append("<h2>1. Material and data</h2>")

    def ref_desc(d):
        if d["pieces"]:
            return ("%s, cut into %d pieces at the %d points where a chromosome's partner in the other haplotype changes "
                    "(table below)" % (", ".join(d["main"]), len(d["seqs"]), sum(len(v) for v in d["cuts"].values())))
        return "haplotype %s of the phased assembly: %s–%s" % (d["row"]["haplotype"].replace("hap", ""), d["main"][0],
                                                              d["main"][-1])

    star_keys = [("Number of input reads", "reads into STARsolo"), ("Average input read length", "read length, bp"),
                 ("Uniquely mapped reads %", "mapped uniquely"), ("% of reads mapped to multiple loci", "mapped to several loci"),
                 ("% of reads mapped to too many loci", "mapped to too many loci"),
                 ("% of reads unmapped: too short", "unmapped: too short"), ("% of reads unmapped: other", "unmapped: other")]
    rows1 = [["sample in the pipeline", B["S"], Pd["S"]],
             ["reference used", ref_desc(B), ref_desc(Pd)],
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
    for sp, d in D.items():
        if d["pieces"]:
            H.append("<h3><i>%s</i>: the pieces of the crossover reference (%s)</h3>" % (NAME[sp], d["S"]))
            H.append(table(["piece", "said as", "cut from", "Mb on that chromosome", "length, Mb",
                            "how the cut was placed (config/split_plans)"],
                           [[s, label(d, s), d["src"][s][0], "%s–%s" % (ff(d["src"][s][1] / 1e6), ff((d["src"][s][1] +
                             d["Ls"][s]) / 1e6)), ff(d["Ls"][s] / 1e6), d["how"].get(s, "")] for s in d["seqs"]]))
            H.append(P("Each piece is a sequence of its own: crossovers are called inside a piece only, a nucleus is "
                       "checked for carrying each piece once (section 4), and how the pieces travel together is measured "
                       "from the pollen (section 7). Landscapes and the one-crossover-per-chromosome checks put the pieces "
                       "back on the chromosome they were cut from. The cut points come from the hap2-on-hap1 alignment "
                       "(breakpoints_from_paf.py; its notes in the last column)."))


run_section('section 1', section_1)

# ---- 2
def section_2():
    H.append("<h2>2. Pipeline at a glance</h2>")

    def mask_txt(d):
        return "%s nucleus × sequence masked (%s), %s calls dropped" % (fi(len(d.get("masked_pairs", ()))),
                                                                         masked_classes(d), fi(d.get("masked_calls")))

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
        ["6 molecule tables", "good markers only; markers < %s bp apart ride on one read and count once (the deepest); "
                              "stretches listed in config/exclude_regions are left out" % ff(MOL_GAP, 0),
         excl_txt(B), excl_txt(Pd)],
        ["7 call crossovers", "hapCO per nucleus, with the settings chosen in the review (section 5)",
         "%s crossovers kept" % fi(B.get("co_total")), "%s crossovers kept" % fi(Pd.get("co_total"))],
        ["8 dosage masks", "a nucleus that carries a sequence twice, not at all or mixed does not count for it, and its "
                           "calls there are dropped (section 4)", mask_txt(B), mask_txt(Pd)]]))


run_section('section 2', section_2)

# ---- 3
def section_3():
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


run_section('section 3', section_3)

# ---- 4
def section_4():
    H.append("<h2>4. Cells: how many, how deep, how clean</h2>")
    H.append(img(FIGS.get("cells")))
    H.append(P("Doublets: two nuclei in one droplet carry both haplotypes, so their %s-molecule windows sit at 0, 0.5 "
               "and 1 instead of 0 or 1, giving haploidness near 0.5 or below; the haploidness floor removes them, together "
               "with diploid barcodes. This filter replaced the switch-rate filter used for <i>Cuscuta</i> and "
               "<i>Spondias</i>: in these libraries most neighbouring markers sit on one read and agree automatically, "
               "so the switch rate measured read clustering rather than doublets (rule select_cells)." % ff(HT_W, 0)))
    H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
        ["barcodes tested / called", "%s / %s" % (fi(NUM.get("binata.barcodes_tested")), fi(NUM.get("binata.cells_called"))),
         "%s / %s" % (fi(NUM.get("paradoxa.barcodes_tested")), fi(NUM.get("paradoxa.cells_called")))],
        ["called barcodes with allele counts (cellsnp)", sel(B, "barcodes counted by cellsnp"),
         sel(Pd, "barcodes counted by cellsnp")],
        ["scored for haploidness (≥ %s windows of %s molecules)" % (ff(HT_MINWIN, 0), ff(HT_W, 0)),
         sel(B, "scored (>= 5 full windows)"), sel(Pd, "scored (>= 5 full windows)")],
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
        ["near-identical pairs (≥ %s%% of ≥ %d windows called in both agree); expected among unrelated nuclei"
         % (ff(100 * DUP_AGREE, 0), DUP_MIN), dup_txt(B), dup_txt(Pd)],
        ["distinct genotypes if each near-identical group counts once (only where pairs exceed chance)",
         fi(B.get("distinct_genotypes")) if dup_above_chance("binata") else "–",
         fi(Pd.get("distinct_genotypes")) if dup_above_chance("paradoxa") else "–"]]))
    H.append(P("Two nuclei of one pollen grain (vegetative and generative or sperm nuclei) come from one meiotic product and "
               "share its genotype: one meiosis observed twice. Unrelated nuclei also agree by chance when the genome "
               "segregates as few units, so the pairs are compared with shuffles that keep how the genome is linked "
               "(nuclei shuffled independently per linkage group, the sequences inherited together in section 7): %s. "
               "The crossover numbers count every nucleus." % "; ".join(dup_reading(sp) for sp in D)))
    H.append("<h3>Dosage: is each sequence there once in each nucleus?</h3>")
    H.append(P("For every nucleus and reference sequence, from the molecules hapCO reads (co_masks.py): <i>too few</i> = "
               "fewer than %s molecules expected; <i>missing</i> = fewer than %s× the expected molecules; <i>doubled</i> = "
               "neighbouring molecules switch allele in at least %s of steps (two copies interleaved) at ≥ %s× the "
               "expected molecules; <i>mixed</i> = switching that often at normal coverage; otherwise <i>haploid</i>. "
               "Expected = the nucleus's molecules × the sequence's median share across nuclei. A masked nucleus × "
               "sequence (config mask_classes: %s in <i>D. binata</i>, %s in <i>D. paradoxa</i>) does not count for that "
               "sequence and its calls there are dropped." % (ff(EXP_LOW, 0), ff(MISS_COV, 2), ff(SW_MIX, 2), ff(DBL_COV, 1),
                                                              masked_classes(B), masked_classes(Pd))))

    def dosage_rows():
        out = []
        for k in ("haploid", "doubled", "missing", "mixed", "too few"):
            out.append(["nucleus × sequence: %s" % k] + [fi(D[sp]["masks"]["class"].eq(k).sum()) if "masks" in D[sp] else "–"
                                                        for sp in D])
        out.append(["masked"] + [fi(len(D[sp].get("masked_pairs", ()))) if "masks" in D[sp] else "–" for sp in D])
        top = []
        for sp in D:
            d = D[sp]
            if "masks" not in d:
                top.append("–")
                continue
            m = d["masks"][d["masks"].masked == "yes"]
            t = m.groupby("chrom")["class"].agg(lambda x: ", ".join("%s %d" % kv for kv in x.value_counts().items()))
            cnt = m.chrom.value_counts()
            top.append("; ".join("%s: %d (%s)" % (label(d, c), cnt[c], t[c]) for c in cnt.index[:4]) or "none")
        out.append(["sequences with the most masked nuclei"] + top)
        return out

    if any("masks" in D[sp] for sp in D):
        H.append(table(["", NAME["binata"], NAME["paradoxa"]], dosage_rows()))
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


run_section('section 4', section_4)

# ---- 5
def section_5():
    H.append("<h2>5. Crossover calling: settings chosen from a review of the calls</h2>")

    def port_txt(sp):
        pc = D[sp].get("port_check")
        if not pc:
            return "<i>%s</i>: no check log, logs/%s_sweep_check.log" % (NAME[sp], D[sp]["S"])
        return "<i>%s</i> identical in %s of %s nuclei at %s, logs/%s_sweep_check.log" % (
            NAME[sp], fi(pc[0]), fi(pc[1]), pc[2], D[sp]["S"])
    H.append(P("hapCO smooths each nucleus's molecules into haplotype blocks and calls a crossover where the block "
               "genotype switches; the settings in the table below decide which blocks count. For each species a sample "
               "of nuclei was called at every combination of the swept values by co_sweep.py, hapCO re-implemented in "
               "Python (checked against the R's own calls at the settings the pipeline ran with: %s); calls of the same switch "
               "were pooled into candidate crossovers ('events'), and each event was judged on a review page "
               "(review_page.py): accepted automatically when at least %s%% of the settings call it with at least %s "
               "molecules at ≥ %s%% purity on each side (review_page.py's defaults), otherwise judged by eye as real, "
               "false or unsure. Every setting is then scored against the verdicts (score_verdicts.py): TP = real events "
               "it calls, FP = false events it calls, FN = real events it misses (and crossovers no setting called, "
               "marked by hand); precision = TP/(TP+FP), recall = TP/(TP+FN), F1 = 2·precision·recall/(precision+recall); "
               "95%% intervals from resampling the reviewed nuclei."
               % ("; ".join(port_txt(sp) for sp in D), ff(100 * AUTO_FRAC, 0), ff(AUTO_SUP, 0), ff(100 * AUTO_PUR, 0))))

    def rv_cell(d, f):
        try:
            return f(d.get("rv", {}))
        except (KeyError, TypeError, AttributeError, ValueError):
            return "–"

    def grid_txt(rv):
        g = rv["grid"]
        return "%d: base_af %s × marker_num %s × block_size %s Mb × terminal_marker_num %s" % (
            rv["n_settings"], "/".join("%g" % v for v in g["base_af"]), "/".join("%g" % v for v in g["marker_num"]),
            "/".join("%g" % (v / 1e6) for v in g["block_size"]), "/".join("%g" % v for v in g["terminal_marker_num"]))

    def best_share(d):
        r = score_row(d, d.get("rv", {}).get("chosen"))
        return "%s%% of resamples" % ff(100 * r.best_share, 0) if r is not None and "best_share" in r else "–"

    rows5 = [("review page (docs/co_review)", lambda rv: "page %s: %s nuclei with at least one event" % (rv["page"],
                                                                                                      fi(rv["cells"]))),
             ("settings swept", grid_txt),
             ("candidate crossovers (events)", lambda rv: fi(rv["events"])),
             ("accepted automatically as real", lambda rv: fi(rv["auto"].get("real", 0))),
             ("judged by eye: real / false / unsure", lambda rv: "%s / %s / %s" % (
                 fi(rv["eye"].get("real", 0)), fi(rv["eye"].get("false", 0)), fi(rv["eye"].get("unsure", 0)))),
             ("not scored: unsure; no verdict; on a sequence masked in the review or in a nucleus marked bad",
              lambda rv: "%s; %s; %s" % (fi(rv["eye"].get("unsure", 0)), fi(rv["none"]), fi(rv["on_masked"]))),
             ("crossovers no setting called, marked by hand", lambda rv: fi(rv["missed"])),
             ("calls scored", lambda rv: rv["mode"] + ("; the molecules in %s left out" % rv["recall_excl"]
                                                      if rv.get("recall_excl") else "")),
             ("scores written (qc/co_review/<sample>/scores.tsv)", lambda rv: rv["scores_time"])]
    H.append(table(["", NAME["binata"], NAME["paradoxa"]],
                   [[lab, rv_cell(B, f), rv_cell(Pd, f)] for lab, f in rows5] +
                   [["settings used before the review", f1_txt(B, "previous", True), f1_txt(Pd, "previous", True)],
                    ["settings chosen", f1_txt(B, "chosen", True), f1_txt(Pd, "chosen", True)],
                    ["chosen settings tie for the best F1 in", best_share(B), best_share(Pd)]]))
    H.append(img(FIGS.get("review")))
    H.append(P("What to look for: the chosen value (large dot) at or near the top of each panel, and how flat the top is: "
               "a flat top means the result does not hinge on the exact value. The dashed line is the F1 of the settings "
               "used before. Hollow points make calls the review never saw (+n per point): those calls are not scored, "
               "so their F1 is optimistic. Values outside the swept grid were re-called from the pipeline's own tables "
               "and matched to the reviewed events (score_verdicts.py --extra)."))
    cc_rows = [("input_rows", "what a block counts: molecules = informative markers on one read, counted once"),
               ("block_size", "minimum length of an interior haplotype block, bp"),
               ("marker_num", "minimum molecules in a block"),
               ("terminal_marker_num", "molecules needed by a chromosome's first and last block (no length rule there)"),
               ("base_af", "smoothing threshold, pass 1"), ("window_af", "smoothing threshold, pass 2"),
               ("genotype", "ALT share that decides a block's genotype"),
               ("mask_classes", "dosage classes whose nucleus × sequence is masked (section 4)")]

    def cfg_val(sp, k):
        v = (CC.get(D[sp]["species"]) or {}).get(k, "–")
        return ", ".join(v) if isinstance(v, list) else v

    H.append(table(["hapCO setting", "meaning", NAME["binata"], NAME["paradoxa"]],
                   [[k, m, cfg_val("binata", k), cfg_val("paradoxa", k)] for k, m in cc_rows] +
                   [["exclusions", "stretches left out of the molecule tables (config/exclude_regions)",
                     excl_txt(B, dropped=False), excl_txt(Pd, dropped=False)],
                    ["stamped", "config _approved: date of the decision (docs/decisions.md)",
                     B.get("rv", {}).get("stamp_date", "–"), Pd.get("rv", {}).get("stamp_date", "–")]]))
    for sp, d in D.items():
        if d.get("rv", {}).get("notes"):
            H.append(details("%s: verdicts changed after a second look (listed in docs/co_review/%s_verdicts.csv)"
                             % (NAME[sp], d["S"]), text="\n".join(l.lstrip("# ") for l in d["rv"]["notes"])))
    H.append("<h3>Crossovers called</h3>")
    H.append(img(FIGS.get("per_grain")))

    def below_half_txt(d):
        bh = d.get("below_half")
        if bh is None:
            return "–"
        return "%d of %d%s" % (len(bh), len(d["main"]), (": " + ", ".join("%s %s" % (short(c), ff(d["per_chrom"][c], 3))
                                                                             for c in bh)) if bh else "")

    H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
        ["grains", fi(B.get("ncell_co")), fi(Pd.get("ncell_co"))],
        ["crossovers (after the dosage masks)", fi(B.get("co_total")), fi(Pd.get("co_total"))],
        ["calls dropped by the dosage masks", fi(B.get("masked_calls")), fi(Pd.get("masked_calls"))],
        ["per grain, mean (SD); median", "%s (%s); %s" % (ff(B.get("co_mean"), 2), ff(B.get("co_sd"), 2),
                                                          ff(B.get("co_median"), 0)),
         "%s (%s); %s" % (ff(Pd.get("co_mean"), 2), ff(Pd.get("co_sd"), 2), ff(Pd.get("co_median"), 0))],
        ["genetic map length, cM", fi(B.get("map_cm")), fi(Pd.get("map_cm"))],
        ["one crossover per bivalent gives, per grain on average (n / 2)", ff(B["n"] / 2), ff(Pd["n"] / 2)],
        ["grains with fewer than n / 2: observed / expected if every bivalent had exactly one crossover",
         "%s / %s of %s" % (fi(B.get("below_obligate")), ff(B.get("below_obligate_exp"), 1), fi(B.get("ncell_co"))),
         "%s / %s of %s" % (fi(Pd.get("below_obligate")), ff(Pd.get("below_obligate_exp"), 1), fi(Pd.get("ncell_co")))],
        ["chromosomes below 0.5 crossovers per grain (pieces summed)", below_half_txt(B), below_half_txt(Pd)],
        ["crossover interval, median kb (resolution)", fi(B.get("width_median_kb")), fi(Pd.get("width_median_kb"))]]))
    H.append(P("A crossover involves two of the four chromatids, so a bivalent with exactly one crossover passes it to "
               "half of the pollen: with one per bivalent, a grain carries n / 2 on average and many carry fewer by chance "
               "(the expected count above); the test that every bivalent gets one is per chromosome: at least 0.5 per "
               "grain. Crossovers per grain counts every crossover of a nucleus; the genetic map sums, over sequences, the "
               "crossovers per nucleus that counts for each, so a masked sequence does not pull it down. <i>D. paradoxa</i> "
               "nuclei carry %s informative molecules against %s in <i>D. binata</i> (medians), so its crossovers are placed "
               "less precisely (median interval %s kb against %s kb)."
               % (fi(Pd.get("selected_molecules_median")), fi(B.get("selected_molecules_median")),
                  fi(Pd.get("width_median_kb")), fi(B.get("width_median_kb")))))


run_section('section 5', section_5)

# ---- 6
def section_6():
    H.append("<h2>6. Landscapes</h2>")
    H.append(img(FIGS.get("landscape")))
    H.append(table(["", NAME["binata"], NAME["paradoxa"]], [
        ["crossovers: ends vs middle", B["shape"], Pd["shape"]],
        ["good markers: ends vs middle", ff(B.get("U_markers"), 2), ff(Pd.get("U_markers"), 2)],
        ["genes: ends vs middle", ff(B.get("U_genes"), 2), ff(Pd.get("U_genes"), 2)],
        ["crossovers vs good markers per 5 Mb window, Spearman", ff(B.get("rho_markers"), 2), ff(Pd.get("rho_markers"), 2)]]))
    H.append(P("'Ends vs middle' = density in the outer 20%% of each chromosome (10%% at each end) over the inner 60%%; "
               "1 = flat. Positions are along the chromosomes, with the pieces of a cut reference put back in place. "
               "Markers follow genes, so the same ratio for markers and genes shows how much end bias detection alone "
               "could make: the ratio is %s for crossovers against %s for markers in <i>D. binata</i>, and %s against %s "
               "in <i>D. paradoxa</i>; crossovers also track marker density window by window (Spearman %s and %s)."
               % (fg(B.get("land", {}).get("U")), fg(B.get("U_markers")), fg(Pd.get("land", {}).get("U")),
                  fg(Pd.get("U_markers")), ff(B.get("rho_markers"), 2), ff(Pd.get("rho_markers"), 2))))
    H.append("<h3>Busiest 5 Mb windows (Poisson p against the genome-wide mean per window)</h3>")
    for sp, d in D.items():
        if d.get("busiest"):
            nw = sum(int(math.ceil(d["L"][c] / 5e6)) for c in d["main"])
            H.append("<p><b>%s</b> (mean %s crossovers per window; the 6 busiest of %s windows, so p must be below "
                     "%.1e to stand out after testing them all)</p>" % (NAME[sp], ff(d["window_mean"]), fi(nw), 0.05 / nw))
            H.append(table(["chromosome", "Mb", "crossovers", "expected", "p"], d["busiest"]))
    H.append(img(FIGS.get("per_chrom")))
    H.append(P("Each bar: crossovers per grain on one chromosome; for a cut reference the pieces are stacked, each "
               "divided by the nuclei that count for that piece (dosage masks). The dashed line is 0.5: below it, fewer "
               "than one crossover per bivalent (a crossover involves two of the four chromatids, so one per bivalent "
               "reaches half of the pollen)."))
    H.append("<h3>Crossover interference</h3>")
    H.append(img(FIGS.get("interference")))
    H.append(P("For every nucleus and chromosome with two or more crossovers, the distance between neighbouring crossovers "
               "as a share of the chromosome; the line is what the same crossovers would give if placed independently "
               "(drawn from that chromosome's crossover positions, 300 times; band = 95%%). Interference shows as too few "
               "close pairs. Closer than 20%% of the chromosome: <i>D. binata</i> %s× expected (%s pairs), "
               "<i>D. paradoxa</i> %s× (%s pairs). This replaces the coefficient-of-coincidence curve, which needs more "
               "double crossovers than <i>D. paradoxa</i>'s %s nuclei give."
               % (ff(B.get("intf", {}).get("ratio"), 2), B.get("intf", {}).get("pairs", "–"),
                  ff(Pd.get("intf", {}).get("ratio"), 2), Pd.get("intf", {}).get("pairs", "–"), fi(Pd.get("ncell_co")))))
    if FIGS.get("structure"):
        H.append("<h3><i>D. paradoxa</i>: crossovers along each chromosome against the structure of its homolog</h3>")
        H.append(img(FIGS.get("structure")))
        H.append(P("Under each chromosome: the pieces it was cut into for crossover calling, and where the other "
                   "haplotype's copy aligns to it (hap2 on hap1, synteny blocks from translocation_map.py): collinear, "
                   "inverted or from a differently numbered chromosome, and where nothing aligns. Black bars: stretches "
                   "without a crossover in any of the %s nuclei, long enough to be unlikely by chance (table below)."
                   % fi(Pd.get("ncell_co"))))
    if Pd.get("dead") is not None:
        dc = Pd["dead_chance"]
        H.append(P("<b>Stretches without crossovers.</b> With %s crossovers per grain, gaps up to %d Mb arise by chance, "
                   "so only stretches of at least %d Mb are listed: %d, %s Mb in total, against %s expected if the same "
                   "crossovers fell uniformly along each chromosome (95%% range %s–%s). What each overlaps:"
                   % (ff(Pd.get("co_mean"), 2), Pd["dead_thr"] - 5, Pd["dead_thr"], len(Pd["dead"]), ff(Pd["dead_mb"], 0),
                      ff(dc[0]), ff(dc[1], 0), ff(dc[2], 0))))
        H.append(table(["chromosome", "Mb", "length, Mb", "inverted, %", "translocated, %", "no homolog, %",
                        "markers vs chromosome median", "contains a cut", "most likely reason"], Pd["dead"]))


run_section('section 6', section_6)

# ---- 7
def section_7():
    H.append("<h2>7. Is each reference right? What the pollen inherit together</h2>")
    H.append(P("A reference agrees with the pollen when the windows of one chromosome are inherited together (only a "
               "crossover separates them) and windows of different chromosomes independently. For every two %g Mb "
               "windows, φ (phi) is the correlation of their genotypes across nuclei: 0 = inherited independently, ±1 = "
               "always together; the sign only says which allele of one goes with which of the other, and depends on which "
               "haplotype each sequence was taken from. Genotypes as in section 4 (≥ %d molecules, ≥ %s%% of them one "
               "allele), a nucleus × sequence masked for dosage left out, and at least %d nuclei with both windows; %s. "
               "Chance level: the nuclei are shuffled independently for each sequence, which keeps each sequence's own "
               "linkage and breaks any between sequences; in 95%% of %d shuffles the strongest pair of sequences stays "
               "below the level given (so fewer than one pair in 20 reports passes it by chance)."
               % (WIN_MB, MIN_MOL, ff(100 * GT_MAJ, 0), PHI_N, "; ".join(
                   "<i>%s</i>: %s" % (NAME[sp], "one nucleus per group of near-identical nuclei, %s of %s (more pairs "
                                      "than chance, section 4; chance pairs go too, so this loses some real meioses)"
                                      % (fi(D[sp]["phi_nuclei"]), fi(len(D[sp]["cells"])))
                                      if D[sp].get("phi_one_per_group") else "all %s nuclei" % fi(D[sp]["phi_nuclei"]))
                   for sp in D if "phi" in D[sp]), NPERM)))
    H.append(img(FIGS.get("pairs")))
    for sp, d in D.items():
        if "co_pairs" not in d:
            continue
        H.append("<p><b>%s</b>: %d of %d pairs of sequences inherited together above chance (|φ| > %s)%s</p>" % (
            NAME[sp], len(d["co_pairs"]), len(d["seqs"]) * (len(d["seqs"]) - 1) // 2, ff(d["fw_thr"], 2),
            ":" if d["co_pairs"] else "."))
        if d["co_pairs"]:
            H.append(table(["sequences", "strongest φ", "at (Mb on each chromosome)", "nuclei", "chance level for this pair",
                            "shuffles as strong (any pair)", "what the assembly says"],
                           [["%s × %s" % (label(d, r["a"]), label(d, r["b"])), "%+.2f" % r["phi"],
                             "%s–%s / %s–%s" % (ff(r["wa"][0], 0), ff(r["wa"][1], 0), ff(r["wb"][0], 0),
                                                ff(r["wb"][1], 0)), fi(r["n"]), ff(r["chance"], 2),
                             "%d of %d" % (r["as_strong"], NPERM), r["why"]] for r in d["co_pairs"]]))
    exp = []
    for sp, d in D.items():
        if not d.get("expected_pairs"):
            continue
        miss = d.get("expected_missing", [])
        exp.append("<i>%s</i>: %s" % (NAME[sp], "all %d such pairs are above chance" % d["expected_pairs"] if not miss else
                                      "%d of %d are above chance; not: %s" % (
                                          d["expected_pairs"] - len(miss), d["expected_pairs"], "; ".join(
                                              "%s × %s (|φ| %s; %s)" % (label(d, r["a"]), label(d, r["b"]),
                                                                       ff(r["absphi"], 2), r["why"]) for r in miss))))
    H.append(P("How to read it: pieces cut from one chromosome, and a piece with the chromosome its copy in the other "
               "haplotype lies on (a translocation the assembly shows), should be inherited together%s. Pairs linked "
               "only through such a piece follow from the two. A pair the assembly does not explain is a question for "
               "the assembly or the biology; expected linkage that is missing would point at the cut reference."
               % ((": " + "; ".join(exp)) if exp else "")))
    H.append(img(FIGS.get("phi_maps")))
    H.append(P("Every window against every other. What to look for: in a correct reference each chromosome is one dark "
               "block on the diagonal, fading away from the diagonal where crossovers separate distant windows, and "
               "everything off the blocks is pale; a misjoin breaks a block, and two chromosomes inherited together light "
               "up the block where they cross."))
    if FIGS.get("phi_zoom"):
        H.append(img(FIGS.get("phi_zoom")))
        H.append(P("Two chromosomes inherited together, window by window (axes in Mb along each chromosome; pieces put "
                   "back in place). What to look for: two chromosomes joined end to end would be linked most strongly next "
                   "to the joined ends, fading with distance; a translocation the assembly knows of shows the same from "
                   "the exchanged stretch. Linkage between the middles of two chromosomes while all their ends recombine "
                   "freely is what a reciprocal exchange between them in one parent, with breakpoints in the middles, "
                   "would give (the exchanged chromosomes travel together through meiosis); selection against pollen with "
                   "some combinations would also link them. Physical evidence decides: Hi-C (below), HiFi reads across "
                   "the junctions, cytology."))
    if FIGS.get("hic_dual"):
        H.append("<h3><i>D. paradoxa</i>: Hi-C on both haplotypes</h3>")
        H.append(img(FIGS.get("hic_dual"), 90))
        H.append(P("Reads that fit either copy of a sequence put contacts between the two copies, so homologous sequence "
                   "on two chromosomes shows as a red line off the diagonal: for an ordinary pair (chr3 h1 and chr3 h2, "
                   "say) one line through their block. Hi-C reads come from one chromosome molecule at a time, so they "
                   "show how the tissue's chromosomes are built, not how the pollen inherit them."))
        rk = Pd.get("hic_rank", {})
        foc = [(k, v) for k, v in sorted(rk.items(), key=lambda kv: kv[1][1])
               if short(k[0]) in Pd.get("hic_focus", []) and short(k[1]) in Pd.get("hic_focus", [])]
        if foc:
            H.append(P("Contacts between %s copies, per bin pair, against the median of all %d pairs of chromosome copies "
                       "that are not homologs (2 Mb bins, %s): %s. What to look for: if the tissue carried %s material "
                       "joined in a way the assembly does not have, reads would connect them: a red block between the "
                       "copies in the map, and that pair at the top of this ranking. Shared repeats raise contacts too, so "
                       "a high rank alone is not a join."
                       % (" and ".join(Pd["hic_focus"]), Pd["hic_npairs"], Pd.get("hic_rank_src", "?"),
                          "; ".join("%s × %s %s× (rank %d)" % (a.replace("_hap", " h"), b.replace("_hap", " h"), ff(x, 2), r)
                                    for (a, b), (x, r) in foc), " and ".join(Pd["hic_focus"]))))
    H.append("<h3>Assembly flags found while choosing the settings</h3>")
    fl = []
    for sp, d in D.items():
        for r in d.get("excl", []):
            fl.append([NAME[sp], "%s:%s–%s" % (r["chrom"], fi(r["start"]), fi(r["end"])), r.get("reason", ""),
                       "left out of the molecule tables (%s molecules dropped); picked by hand during the review, to be "
                       "replaced by an automatic filter" % fi(d.get("excl_dropped"))])
    if fl:
        H.append(table(["species", "stretch", "evidence (%s)" % ", ".join(d["excl_file"] for d in D.values() if
                                                                           d.get("excl_file")), "what was done"], fl))
    H.append(P("Other places noted in the decision entries below (a few molecules that follow another place in the "
               "genome, a stretch where the assembly's phase switches) are listed there and not masked."))
    dec = "docs/decisions.md"
    if have(dec):
        ents = re.split(r"\n(?=## )", open(dec).read())
        for sp, d in D.items():
            hit = [e for e in ents if e.startswith("## ") and NAME[sp] in e.split("\n", 1)[0] and
                   "crossover calling" in e.split("\n", 1)[0]]
            if hit:
                H.append(details("%s: decision entry for the crossover-calling settings (docs/decisions.md)" % NAME[sp],
                                 text=hit[-1].strip()))
    for sp, d in D.items():
        if d.get("unexplained"):
            p = "qc/linkage/Dparadoxa_std/chr5_chr6_check.txt"
            if sp == "paradoxa" and any({short(d["src"][r["a"]][0]), short(d["src"][r["b"]][0])} == {"chr5", "chr6"}
                                        for r in d["unexplained"]) and os.path.exists(p):
                H.append(details("earlier check of chr5 against chr6 (pollen linkage on Dparadoxa_std, where chr5_hap1 and "
                                 "chr6_hap1 are the same sequences; hap2-on-hap1 alignment)", p))
    prev = []
    for p in [x.strip() for x in A.previous.split(",") if x.strip()]:
        if os.path.exists(p):
            USED.append(p)
            m = re.search(r"Generated ([\d-]+ [\d:]+) from CO_smk commit (\w+)", open(p, errors="replace").read())
            prev.append("%s (generated %s from commit %s)" % (p, m.group(1), m.group(2)) if m else p)
    H.append(P("Not repeated here: earlier versions of this report called <i>D. paradoxa</i> on two joined references, "
               "chr1 and chr2 as the tissue's assembly joins them (A + D, Dparadoxa_std) and as the pollen inherit their "
               "arms (C + B, Dparadoxa_CB), and compared them join by join (arm diagram, crossovers at the joins, rDNA "
               "and Hi-C libraries at the joins), with a pollen-linkage and Hi-C map on each crossover reference (both "
               "species) and the rDNA arrays. The cut reference used here needs no join, so that comparison does not "
               "change the counts; it is kept in %s." % (" and ".join(html.escape(x) for x in prev) if prev else
                                                         "the version-4 report (reports/supervisor_CB/report.html)")))


run_section('section 7', section_7)

# ---- 8, 9
def section_8_9():
    H.append("<h2>8. Limitations</h2>")
    lim = ["Markers come from pollen RNA and sit in expressed genes (%s%% and %s%% of good markers inside annotated genes); "
           "%s%% and %s%% of the chromosome sequence lies in marker gaps over 1 Mb, where crossovers are placed only "
           "roughly." % (ff(100 * B.get("good_in_genes", float("nan")), 0), ff(100 * Pd.get("good_in_genes", float("nan")), 0),
                         ff(100 * B.get("gap_1mb_share", float("nan"))), ff(100 * Pd.get("gap_1mb_share", float("nan")))),
           "<i>D. paradoxa</i> rests on %s nuclei with a median of %s informative molecules, against %s nuclei with %s in "
           "<i>D. binata</i>." % (fi(len(Pd["cells"])), fi(Pd.get("selected_molecules_median")), fi(len(B["cells"])),
                                  fi(B.get("selected_molecules_median"))),
           "Near-identical nuclei are counted separately in the crossover numbers (section 4; <i>D. binata</i> %s, "
           "<i>D. paradoxa</i> %s)." % (dup_txt(B), dup_txt(Pd))]
    for sp, d in D.items():
        if d.get("cut_gaps"):
            lim.append("<i>%s</i>: a crossover between the last good marker before a cut and the first after it cannot be "
                       "called: %s." % (NAME[sp], "; ".join("%s at %s Mb, %s kb" % (short(c), ff(x / 1e6, 1),
                                                                                   fi((b - a) / 1e3))
                                                           for c, x, a, b in d["cut_gaps"])))
        if d.get("excl"):
            lim.append("<i>%s</i>: %s was left out by hand, after the review showed its alleles following another place in "
                       "the genome (reason in %s). One test applied to every sample should replace the hand-picked list "
                       "(section 9)." % (NAME[sp], excl_txt(d, dropped=False), d["excl_file"]))
    lim.append("A nucleus × sequence masked for dosage (section 4) does not count for that sequence: %s in <i>D. binata</i>, "
               "%s in <i>D. paradoxa</i>." % (fi(len(B.get("masked_pairs", ()))), fi(len(Pd.get("masked_pairs", ())))))
    lim.append("hapCO's step that narrows a breakpoint does not check the direction of the switch it narrows to, so an "
               "occasional breakpoint lands on a nearby opposite switch (seen in the <i>D. binata</i> review): it moves "
               "a crossover, it does not add or remove one.")
    H.append("<ul>%s</ul>" % "".join("<li>%s</li>" % x for x in lim))
    H.append("<h2>9. Next steps and asks</h2>")
    nxt = []
    if any(d.get("excl") for d in D.values()):
        nxt.append("Replace the hand-picked exclusion by an automatic filter: per gene cluster, test across nuclei whether "
                   "its alleles follow its neighbours, and leave out the clusters that follow another place instead; run "
                   "it with the next <i>D. paradoxa</i> rerun, so both species are called the same way.")
    for sp, d in D.items():
        for r in d.get("unexplained", []):
            nxt.append("<i>%s</i> %s × %s: HiFi reads that span a junction between them, and Hi-C at finer bins near the "
                       "stretches inherited together; meiotic spreads with paints for both decide."
                       % (NAME[sp], label(d, r["a"]), label(d, r["b"])))
    nxt.append("Cytology: meiotic chromosome spreads with oligo-FISH paints for the translocated arms and 45S rDNA FISH.")
    dup_sp = [sp for sp in D if dup_above_chance(sp)]
    if dup_sp:
        nxt.append("Count near-identical nuclei once in the crossover numbers (%s: more pairs than unrelated nuclei give by "
                   "chance, section 4)." % " and ".join("<i>%s</i>" % NAME[sp] for sp in dup_sp))
    nxt.append("The places flagged in the decision entries (section 7), with HiFi reads or Hi-C, when the filter above is "
               "built.")
    H.append("<ul>%s</ul>" % "".join("<li>%s</li>" % x for x in nxt))


run_section('section 8, 9', section_8_9)

# ---- appendix
def section_appendix():
    H.append("<h2>Appendix</h2><h3>A. Per-chromosome landscapes</h3>")
    H.append(img(FIGS.get("chrom_landscapes")))
    H.append("<h3>B. Per reference sequence</h3>")
    for sp, d in D.items():
        if "per_seq" not in d:
            continue
        mk = d.get("masks")
        mm = mk[mk.masked == "yes"] if mk is not None else None
        rows = []
        for s in d["seqs"]:
            cl = mm[mm.chrom == s]["class"].value_counts() if mm is not None else {}
            rows.append([s, label(d, s), ff(d["Ls"][s] / 1e6, 1), fi(d["nseq_used"][s]),
                         ", ".join("%s %d" % kv for kv in cl.items()) if len(cl) else "", fi(d["co_seq"][s]),
                         ff(d["per_seq"][s], 2), ff(100 * d["per_seq"][s] / (d["Ls"][s] / 1e6), 2)])
        H.append("<p><b>%s</b> (%s)</p>" % (NAME[sp], d["S"]))
        H.append(table(["sequence", "said as", "Mb", "nuclei that count", "masked (class)", "crossovers",
                        "per nucleus", "cM/Mb"], rows))
    H.append("<h3>C. Crossover-calling parameters as run</h3>")
    for sp, d in D.items():
        H.append(details("%s: co_calling_params.txt" % NAME[sp], "results/crossovers/%s/co_calling_params.txt" % d["S"]))
    H.append("<h3>D. Files read</h3><pre>%s</pre>" % html.escape("\n".join(sorted(set(map(str, USED))))))
    if MISSING:
        H.append("<h3>E. Files not found</h3><pre>%s</pre>" % html.escape("\n".join(sorted(set(map(str, MISSING))))))


run_section('section appendix', section_appendix)

TOP = []
if WARN:
    TOP.append("<p class='missing'><b>Checks that failed</b> (the numbers below may not be what the text says): %s</p>"
               % html.escape("; ".join(WARN)))
if PROBLEMS:
    TOP.append("<p class='missing'>Parts that could not be built: %s</p>" % html.escape(
        "; ".join("%s (%s)" % p for p in PROBLEMS)))
H[2:2] = TOP
CSS = """body{font-family:Helvetica,Arial,sans-serif;max-width:1080px;margin:24px auto;padding:0 16px;color:#2C2C2A;
font-size:13px;line-height:1.45}h1{font-size:21px}h2{font-size:16px;border-bottom:1px solid #D3D1C7;margin-top:30px}
h3{font-size:13px;margin-top:18px}table{border-collapse:collapse;margin:8px 0 14px;font-size:12px}td,th{border:1px
solid #D3D1C7;padding:3px 8px;text-align:left;vertical-align:top}th{background:#F1EFE8}table.plain td{border:none}
.meta{color:#5F5E5A;font-size:11px}.missing{color:#A32D2D}.path{color:#888780;font-size:11px}code{font-size:11px}
pre{font-size:10.5px;background:#F7F6F2;padding:8px;overflow-x:auto;white-space:pre-wrap}img{display:block;margin:8px 0;
max-width:100%}details{margin:4px 0}@media print{img{page-break-inside:avoid}table{page-break-inside:avoid}}"""
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
