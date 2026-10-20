#!/usr/bin/env python3
"""rdna_ideogram.py -- all 12 D. paradoxa chromosomes as the soma carries them (the assembly) and as the pollen
inherit them, with every 45S and 5S rDNA locus at its place and each circle's area in proportion to its size.

Nothing typed in. Run from the CO_smk root. Reads:
  results/reference/SAMPLE/MANIFEST.txt   source_fasta = the two-haplotype assembly; its .fai gives the lengths
  qc/reference/SAMPLE/breakpoints.csv     the joins (where a chromosome changes homolog)
  reports/repeat_context/rdna_loci_SAMPLE.csv   every 45S and 5S locus (CARP pieces merged within 20 kb)
  REFS/SPECIES_repeats/rDNA.gff3          the CARP rDNA features: the rDNA bp inside each locus = its size (without
                                          this file the locus span is the size, and the figure says so)
  reports/recon_audit/recon_audit.txt     optional: which way P and Q face, and the chr5/chr6 coupling windows
  reports/cosegregation_2n/karyotype_test_SAMPLE.*   optional: the coupling point
Names: A = chr1_hap1 = L1 + P, B = chr2_hap1 = L2 + P, C = chr1_hap2 = L1 + Q, D = chr2_hap2 = L2 + Q;
X = the arm on chr4_hap1 (X + core) and on chr3_hap2 (chr3 + X).
Soma = the assembly. Pollen: A' = A's L1 + D's Q and D' = D's L2 + A's P, joined at the assembly's joins; B, C and
chr3/chr4 as assembled; chr5/chr6 as assembled with the coupled middles shaded (their structure is unknown).
Writes OUT/rdna_ideogram_SAMPLE.{png,pdf} and OUT/rdna_loci_both_karyotypes_SAMPLE.tsv.
Usage: rdna_ideogram.py [--sample Dparadoxa_split] [--species Drosera_paradoxa]
                        [--refs /netscratch/dep_mercier/grp_marques/Aaryan/refs] [--out reports/ideogram]
                        [--size bp|span] [--label45_kb 5] [--label5_kb 50]
                        [--chr5 chr5_hap1:90:225] [--chr6 chr6_hap1:120:295]
"""
import argparse
import bisect
import csv
import datetime
import glob
import math
import os
import re
import subprocess
import sys
from urllib.parse import unquote

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle

ap = argparse.ArgumentParser()
ap.add_argument("--sample", default="Dparadoxa_split")
ap.add_argument("--species", default="Drosera_paradoxa")
ap.add_argument("--refs", default="/netscratch/dep_mercier/grp_marques/Aaryan/refs")
ap.add_argument("--out", default="reports/ideogram")
ap.add_argument("--size", choices=("bp", "span"), default="bp",
                help="bp: the rDNA bp in each locus (needs rDNA.gff3); span: first to last rDNA base")
ap.add_argument("--label45_kb", type=float, default=5.0, help="45S loci at least this big get a size label")
ap.add_argument("--label5_kb", type=float, default=50.0, help="5S loci at least this big get a size label")
ap.add_argument("--chr5", default="chr5_hap1:90:225", help="coupling window if recon_audit.txt has none")
ap.add_argument("--chr6", default="chr6_hap1:120:295", help="coupling window if recon_audit.txt has none")
A = ap.parse_args()
S = A.sample
MB = 1e6
os.makedirs(A.out, exist_ok=True)
READ, MISSING, NOTES = [], [], []
ARM_COL = {"L1": "#5DCAA5", "L2": "#F0997B", "P": "#AFA9EC", "Q": "#FAC775", "X": "#85B7EB", "core": "#D3D1C7",
           "chr3": "#ECEAE3", "chr5": "#ECEAE3", "chr6": "#ECEAE3"}
RCOL = {"45S": "#A32D2D", "5S": "#E07B39"}
COUP = "#7B6FD0"
INK, MUTED = "#2C2C2A", "#6B6A64"
S_MAX = {"45S": 150.0, "5S": 380.0}      # circle area (pt^2) of the largest locus of each type
S_MIN = 1.5                              # loci whose circle would be smaller than this are drawn as ticks
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.spines.left": False, "pdf.fonttype": 42})


def have(p):
    if p and os.path.exists(p):
        READ.append(p)
        return True
    MISSING.append(p)
    return False


def git_head():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, universal_newlines=True).stdout.strip() or "?"
    except OSError:
        return "?"


def fmt(kb):
    if kb >= 1000:
        return "%.1f Mb" % (kb / 1e3)
    return ("%.0f kb" if kb >= 100 else "%.1f kb" if kb >= 1 else "%.2f kb") % kb


def fmt_group(sizes):
    if len(sizes) > 3:
        return "%s in %d loci" % (fmt(sum(sizes)), len(sizes))
    if all(s_ < 1000 for s_ in sizes):
        return " + ".join(("%.0f" if s_ >= 100 else "%.1f") % s_ for s_ in sizes) + " kb"
    return " + ".join(fmt(s_) for s_ in sizes)


# ------------------------------------------------------------------ lengths and joins
man = "results/reference/%s/MANIFEST.txt" % S
if not have(man):
    sys.exit("STOP: %s is missing" % man)
src = ""
for line in open(man):
    if line.startswith("source_fasta:"):
        src = line.split(":", 1)[1].strip()
if not have(src + ".fai"):
    sys.exit("STOP: %s.fai is missing (source_fasta in %s)" % (src, man))
LEN = {}
for line in open(src + ".fai"):
    f = line.split("\t")
    if len(f) > 1 and re.match(r"^chr\d+_hap\d$", f[0]):
        LEN[f[0]] = int(f[1])
bpf = "qc/reference/%s/breakpoints.csv" % S
if not have(bpf):
    sys.exit("STOP: %s is missing" % bpf)
BP = {}
for r in csv.DictReader(open(bpf)):
    BP.setdefault(r["chrom"], []).append((int(float(r["cut"])), r["left_homolog"].strip(), r["right_homolog"].strip()))
EXPECT = {"chr1_hap1": ("chr1_hap2", "none"), "chr2_hap1": ("chr2_hap2", "none"), "chr1_hap2": ("chr1_hap1", "none"),
          "chr2_hap2": ("chr2_hap1", "none"), "chr3_hap2": ("chr3_hap1", "chr4_hap1"),
          "chr4_hap1": ("chr3_hap2", "chr4_hap2")}
bad = [c for c, (lh, rh) in EXPECT.items() if len(BP.get(c, [])) != 1 or BP[c][0][1:] != (lh, rh)]
bad += sorted(set(BP) - set(EXPECT))
if set(LEN) != {"chr%d_hap%d" % (i, h) for i in range(1, 7) for h in (1, 2)}:
    sys.exit("STOP: %s.fai does not hold exactly chr1_hap1 .. chr6_hap2 (it has: %s)" % (src, ", ".join(sorted(LEN))))
if bad:
    sys.exit("STOP: the joins in %s are not the six this figure draws (A, B, C, D, chr3|X, X|core): check %s"
             % (bpf, ", ".join(bad)))
CUT = {c: BP[c][0][0] for c in EXPECT}
cA, cB, cC, cD = CUT["chr1_hap1"], CUT["chr2_hap1"], CUT["chr1_hap2"], CUT["chr2_hap2"]
lA, lD = LEN["chr1_hap1"], LEN["chr2_hap2"]
c3, c4 = CUT["chr3_hap2"], CUT["chr4_hap1"]
SEGS = {"chr1_hap1": [("L1", 0, cA), ("P", cA, lA)], "chr2_hap1": [("L2", 0, cB), ("P", cB, LEN["chr2_hap1"])],
        "chr1_hap2": [("L1", 0, cC), ("Q", cC, LEN["chr1_hap2"])], "chr2_hap2": [("L2", 0, cD), ("Q", cD, lD)],
        "chr3_hap1": [("chr3", 0, LEN["chr3_hap1"])], "chr3_hap2": [("chr3", 0, c3), ("X", c3, LEN["chr3_hap2"])],
        "chr4_hap1": [("X", 0, c4), ("core", c4, LEN["chr4_hap1"])], "chr4_hap2": [("core", 0, LEN["chr4_hap2"])]}
for c in ("chr5_hap1", "chr5_hap2", "chr6_hap1", "chr6_hap2"):
    SEGS[c] = [(c[:4], 0, LEN[c])]
lAp, lDp = cA + (lD - cD), cD + (lA - cA)
SEGS["Ap"] = [("L1", 0, cA), ("Q", cA, lAp)]
SEGS["Dp"] = [("L2", 0, cD), ("P", cD, lDp)]
OFF = {"chr4_hap2": c4}                   # chr4_hap2 (core only) is drawn under chr4_hap1's core

# ------------------------------------------------------------------ orientation of the moved arms, coupling
ORI = {"hap1": True, "hap2": True}        # True: A's P and D's Q go in as they are
ori_how = "taken as is (no recon_audit.txt)"
WIN, win_how = {}, "--chr5/--chr6 defaults"
aud = "reports/recon_audit/recon_audit.txt"
if have(aud):
    cur, seen = None, []
    for line in open(aud):
        m = re.match(r"^\s+(hap[12]): target ", line)
        if m:
            cur = m.group(1)
        if cur and "-> the join ends" in line:
            ORI[cur] = "do NOT correspond" not in line
            seen.append(cur)
        m = re.search(r"window on (chr[56]_hap[12]): ([\d.]+)-([\d.]+) Mb", line)
        if m:
            WIN[m.group(1)] = (float(m.group(2)) * MB, float(m.group(3)) * MB)
    if len(set(seen)) == 2:
        ori_how = "as P against P and Q against Q say (%s)" % ", ".join(
            "%s %s" % ("A's P" if h == "hap1" else "D's Q", "as is" if ORI[h] else "reversed")
            for h in ("hap1", "hap2"))
    else:
        ori_how = "taken as is (%s has no P/Q alignment)" % aud
    if WIN:
        win_how = aud
if not WIN:
    for spec in (A.chr5, A.chr6):
        c, a, b = spec.split(":")
        WIN[c] = (float(a) * MB, float(b) * MB)
CP = {}
for kt in sorted(glob.glob("reports/cosegregation_2n/karyotype_test_%s.*" % S)):
    if not kt.endswith((".tsv", ".txt", ".csv")):
        continue
    m = re.search(r"(chr5_hap\d) ([\d.]+) Mb x (chr6_hap\d) ([\d.]+) Mb", open(kt, errors="replace").read())
    if m:
        READ.append(kt)
        CP = {m.group(1): float(m.group(2)) * MB, m.group(3): float(m.group(4)) * MB}
        break


def to_pollen(copy, pos):
    """A position on an assembled chromosome -> (row, position) in the pollen karyotype."""
    if copy == "chr1_hap1" and pos > cA:
        return "Dp", (cD + (pos - cA)) if ORI["hap1"] else (cD + (lA - pos) + 1)
    if copy == "chr2_hap2" and pos > cD:
        return "Ap", (cA + (pos - cD)) if ORI["hap2"] else (cA + (lD - pos) + 1)
    if copy == "chr1_hap1":
        return "Ap", pos
    if copy == "chr2_hap2":
        return "Dp", pos
    return copy, pos


# ------------------------------------------------------------------ the loci and their size
rd = "reports/repeat_context/rdna_loci_%s.csv" % S
if not have(rd):
    sys.exit("STOP: %s is missing (repeat_context.py writes it)" % rd)
LOCI = []
for r in csv.DictReader(open(rd)):
    if r["copy"] in LEN and r["type"] in RCOL:
        LOCI.append(dict(type=r["type"], copy=r["copy"], start=int(float(r["start"])), end=int(float(r["end"])),
                         span=float(r["kb"]), pieces=r.get("pieces", "")))
gff = os.path.join(A.refs, "%s_repeats" % A.species, "rDNA.gff3")
IV = {}
if have(gff):
    for line in open(gff):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or "ID=UA_L2_" in f[8]:       # nested members: counted in their container, as in the loci
            continue
        cls = ""
        for kv in f[8].split(";"):
            if kv.startswith("classification="):
                cls = unquote(kv.split("=", 1)[1])
        if cls.startswith("rDNA"):
            IV.setdefault((f[0], "5S" if cls.startswith("rDNA/5S") else "45S"), []).append((int(f[3]) - 1, int(f[4])))
for k in IV:
    IV[k].sort()
STARTS = {k: [a for a, _ in v] for k, v in IV.items()}


def rdna_bp(copy, t, start, end):
    """bp covered by this type's rDNA features inside [start, end] (1-based, inclusive), overlaps counted once."""
    iv, st = IV.get((copy, t), []), STARTS.get((copy, t), [])
    lo, hi = start - 1, end
    i = max(0, bisect.bisect_left(st, lo) - 200)
    tot, cs, ce = 0, None, None
    while i < len(iv) and iv[i][0] < hi:
        a, b = max(iv[i][0], lo), min(iv[i][1], hi)
        i += 1
        if b <= a:
            continue
        if cs is None or a > ce:
            tot += (ce - cs) if cs is not None else 0
            cs, ce = a, b
        else:
            ce = max(ce, b)
    return tot + ((ce - cs) if cs is not None else 0)


for L_ in LOCI:
    L_["bp_kb"] = rdna_bp(L_["copy"], L_["type"], L_["start"], L_["end"]) / 1e3 if IV else None
nzero = sum(1 for L_ in LOCI if L_["bp_kb"] == 0)
USE_BP = A.size == "bp" and bool(IV) and nzero <= 0.02 * max(len(LOCI), 1)
if A.size == "bp" and IV and not USE_BP:
    NOTES.append("%d of %d loci hold no rDNA feature of their type in %s: the table and the GFF disagree, so the "
                 "size shown is the locus span" % (nzero, len(LOCI), gff))
elif USE_BP and nzero:
    NOTES.append("%d loci hold no rDNA feature of their type in %s (drawn as ticks)" % (nzero, gff))
elif A.size == "bp" and not IV:
    NOTES.append("no rDNA features read from %s, so the size shown is the locus span" % gff)
SIZE_WHAT = "the rDNA bp in the locus (%s)" % gff if USE_BP else "the locus span, first to last rDNA base"
for L_ in LOCI:
    L_["size"] = L_["bp_kb"] if USE_BP else L_["span"]
    mid = (L_["start"] + L_["end"]) / 2.0
    L_["soma"] = (L_["copy"], mid)
    L_["pollen"] = to_pollen(L_["copy"], mid)
BIG = {t: max([L_["size"] for L_ in LOCI if L_["type"] == t] or [1.0]) for t in RCOL}


def area(t, size):
    return S_MAX[t] * size / BIG[t] if BIG[t] > 0 else 0.0


SMALLEST = {t: BIG[t] * S_MIN / S_MAX[t] for t in RCOL}      # kb below which a locus is a tick

# ------------------------------------------------------------------ figure
ROWS = {
    "soma": [("chr1_hap1", "A = chr1_hap1", "L1 + P"), ("chr1_hap2", "C = chr1_hap2", "L1 + Q"),
             ("chr2_hap1", "B = chr2_hap1", "L2 + P"), ("chr2_hap2", "D = chr2_hap2", "L2 + Q")],
    "pollen": [("Ap", "A\u2032 = A's L1 + D's Q", "L1 + Q, not in the assembly"),
               ("chr1_hap2", "C = chr1_hap2", "L1 + Q"), ("chr2_hap1", "B = chr2_hap1", "L2 + P"),
               ("Dp", "D\u2032 = D's L2 + A's P", "L2 + P, not in the assembly")]}
for k in ROWS:
    ROWS[k] += [("chr3_hap1", "chr3_hap1", "chr3"), ("chr3_hap2", "chr3_hap2", "chr3 + X"),
                ("chr4_hap1", "chr4_hap1", "X + core"), ("chr4_hap2", "chr4_hap2", "core, drawn under chr4_hap1's"),
                ("chr5_hap1", "chr5_hap1", ""), ("chr5_hap2", "chr5_hap2", ""),
                ("chr6_hap1", "chr6_hap1", ""), ("chr6_hap2", "chr6_hap2", "")]
JOINS = {"soma": {"chr1_hap1": cA, "chr2_hap1": cB, "chr1_hap2": cC, "chr2_hap2": cD, "chr3_hap2": c3,
                  "chr4_hap1": c4},
         "pollen": {"Ap": cA, "Dp": cD, "chr2_hap1": cB, "chr1_hap2": cC, "chr3_hap2": c3, "chr4_hap1": c4}}
NAME = {"Ap": "A'", "Dp": "D'"}         # ASCII names for the table and the printout
ROWH, GAP, BAR = 1.0, 0.5, 0.24
YPOS = []
y = 0.0
for i in range(12):
    if i and i % 2 == 0:
        y -= GAP
    YPOS.append(y)
    y -= ROWH
YLIM = (YPOS[-1] - 0.75, YPOS[0] + 0.75)
XEND = max([LEN[c] + OFF.get(c, 0) for c in LEN] + [lAp, lDp]) / MB
XLIM = (-6, XEND + 40)
fig = plt.figure(figsize=(22, 14))
gs = fig.add_gridspec(2, 2, height_ratios=[11.6, 1.25], hspace=0.07, wspace=0.15, left=0.088, right=0.992,
                      top=0.935, bottom=0.07)
REND = fig.canvas.get_renderer()
LABEL_PT = 6.4


def text_w(txt):
    """width of a size label, points"""
    return REND.get_text_width_height_descent(txt, FontProperties(size=LABEL_PT), ismath=False)[0] * 72.0 / fig.dpi


TITLES = {"soma": "SOMA: the assembly (HiFi contigs + Hi-C)",
          "pollen": "POLLEN: the karyotype the pollen genotypes imply "
                    "(A\u2032 and D\u2032 rebuilt at the assembly's joins)"}
for col, kar in enumerate(("soma", "pollen")):
    ax = fig.add_subplot(gs[0, col])
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    W_in, H_in = fig.get_size_inches()
    p = ax.get_position()
    PPX = p.width * W_in * 72 / (XLIM[1] - XLIM[0])          # points per Mb
    PPY = p.height * H_in * 72 / (YLIM[1] - YLIM[0])         # points per row unit
    labels = []
    for i, (row, name, arms) in enumerate(ROWS[kar]):
        y, off = YPOS[i], OFF.get(row, 0)
        segs = SEGS[row]
        if kar == "pollen" and row in WIN:
            w0, w1 = WIN[row]
            ax.add_patch(Rectangle(((off + w0) / MB, y - 0.46), (w1 - w0) / MB, 0.92, color=COUP, alpha=0.12, lw=0,
                                   zorder=0))
            if row in CP:
                ax.plot([(off + CP[row]) / MB] * 2, [y - 0.46, y + 0.46], color=COUP, lw=1.0, zorder=0.5)
        for arm, a, b in segs:
            ax.add_patch(Rectangle(((off + a) / MB, y - BAR / 2), (b - a) / MB, BAR, color=ARM_COL[arm], lw=0,
                                   zorder=1))
            if (b - a) / MB > 30 and arm not in ("chr3", "chr5", "chr6"):
                ax.text((off + (a + b) / 2) / MB, y, arm, ha="center", va="center", fontsize=6.5, color=INK, zorder=2)
        ax.add_patch(Rectangle((off / MB, y - BAR / 2), segs[-1][2] / MB, BAR, fill=False, ec="#9A988F", lw=0.6,
                               zorder=2))
        if row in JOINS[kar]:
            j = (off + JOINS[kar][row]) / MB
            ax.plot([j, j], [y - BAR / 2 - 0.07, y + BAR / 2 + 0.07], color=INK, lw=1.2, zorder=3,
                    ls=(0, (2, 1.2)) if row in NAME else "-")
        labels.append("%s\n%s" % (name, arms) if arms else name)
        for t, sgn in (("45S", 1), ("5S", -1)):
            here = [L_ for L_ in LOCI if L_["type"] == t and L_[kar][0] == row]
            small = [L_ for L_ in here if area(t, L_["size"]) < S_MIN]
            big = sorted([L_ for L_ in here if area(t, L_["size"]) >= S_MIN], key=lambda L_: -L_["size"])
            if small:
                y0 = y + sgn * BAR / 2
                ax.vlines([(off + L_[kar][1]) / MB for L_ in small], y0, y0 + sgn * 0.09, color=RCOL[t], lw=0.5,
                          alpha=0.55, zorder=2)
            # circles hang below (5S) or stand on (45S) the chromosome; geometry in points, relative to the bar centre
            G = {id(L_): ((off + L_[kar][1]) / MB * PPX, math.sqrt(area(t, L_["size"])) / 2) for L_ in big}
            CY = {k: sgn * (BAR / 2 * PPY + 1.0 + r) for k, (cx, r) in G.items()}
            if big:
                ax.scatter([G[id(L_)][0] / PPX for L_ in big], [y + CY[id(L_)] / PPY for L_ in big],
                           s=[area(t, L_["size"]) for L_ in big], color=RCOL[t], alpha=0.85, ec="white", lw=0.5,
                           zorder=5)
            thr = A.label45_kb if t == "45S" else A.label5_kb
            lab = sorted([L_ for L_ in big if L_["size"] >= thr], key=lambda L_: L_[kar][1])
            groups = []
            for L_ in lab:
                if groups and L_[kar][1] - groups[-1][-1][kar][1] <= 3 * MB:
                    groups[-1].append(L_)
                else:
                    groups.append([L_])
            placed = []                                       # label boxes already drawn: (x0, x1, y)
            for g in groups:
                txt = fmt_group([L_["size"] for L_ in g])
                w_ = text_w(txt)
                right = max(G[id(L_)][0] + G[id(L_)][1] for L_ in g) + 2.5
                left = min(G[id(L_)][0] - G[id(L_)][1] for L_ in g) - 2.5 - w_
                ref = id(max(g, key=lambda L_: L_["size"]))
                y_lab = sgn * max(abs(CY[ref]), BAR / 2 * PPY + 4.5)
                for lev in range(16):          # right of the circles, else left, else one step further out
                    x0, yl = (right, left)[lev % 2], y_lab + sgn * 8.0 * (lev // 2)
                    x1 = x0 + w_
                    if not any(x1 > bx0 - 2 and x0 < bx1 + 2 and abs(yl - by) < 7.5 for bx0, bx1, by in placed) \
                            and not any(x1 > cx - r and x0 < cx + r and abs(yl - CY[k]) < r + 3.5
                                        for k, (cx, r) in G.items()):
                        break
                placed.append((x0, x1, yl))
                ax.text(x0 / PPX, y + yl / PPY, txt, ha="left", va="center", fontsize=LABEL_PT, color=RCOL[t],
                        zorder=6)
    ax.set_yticks(YPOS)
    ax.set_yticklabels(labels, fontsize=7.6)
    ax.tick_params(axis="y", length=0)
    ax.set_xticks(range(0, int(XEND) + 1, 50))
    ax.set_xlabel("Mb from the chromosome's start", fontsize=8)
    ax.set_title(TITLES[kar], loc="left", fontsize=11, color=INK)

# ------------------------------------------------------------------ legend row
lg = fig.add_subplot(gs[1, :])
lg.axis("off")
lg.set_xlim(0, 1)
lg.set_ylim(0, 1)
x = 0.0
for arm in ("L1", "P", "L2", "Q", "X", "core"):
    lg.add_patch(Rectangle((x, 0.66), 0.016, 0.13, color=ARM_COL[arm]))
    lg.text(x + 0.02, 0.725, arm, va="center", fontsize=8)
    x += 0.042
lg.add_patch(Rectangle((x, 0.66), 0.016, 0.13, fc=ARM_COL["chr3"], ec="#9A988F", lw=0.4))
lg.text(x + 0.02, 0.725, "rest of chr3, chr5, chr6", va="center", fontsize=8)
lg.plot([0.0, 0.016], [0.36, 0.36], color=INK, lw=1.2)
lg.text(0.02, 0.36, "join in the assembly", va="center", fontsize=8)
lg.plot([0.118, 0.134], [0.36, 0.36], color=INK, lw=1.2, ls=(0, (2, 1.2)))
lg.text(0.138, 0.36, "join rebuilt for the pollen", va="center", fontsize=8)
lg.add_patch(Rectangle((0.0, 0.0), 0.016, 0.16, color=COUP, alpha=0.18, lw=0))
lg.text(0.02, 0.08, "pollen: chr5 and chr6 middles inherited together (structure unknown, nothing moved); "
        "line = coupling point", va="center", fontsize=8, color=COUP)


def nice(v):
    """the 1-2-5 number nearest to v on a log scale"""
    e = math.floor(math.log10(v))
    return min((m * 10 ** e for m in (1, 2, 5, 10)), key=lambda c: abs(math.log10(c / v)))


x0 = 0.385
for t, where in (("45S", "above"), ("5S", "below")):
    lg.text(x0, 0.93, "%s rDNA, %s each chromosome" % (t, where), fontsize=8.5, color=RCOL[t], va="center")
    vals = sorted({nice(BIG[t] / 20), nice(BIG[t] / 4)} - {BIG[t]}) + [BIG[t]]
    for k, v in enumerate(vals):
        xx = x0 + 0.018 + k * 0.05
        lg.scatter([xx], [0.5], s=[area(t, v)], color=RCOL[t], alpha=0.85, ec="white", lw=0.5)
        lg.text(xx, 0.14, fmt(v) + (" (largest)" if k == len(vals) - 1 else ""), ha="center", fontsize=7,
                color=MUTED)
    lg.text(x0 + 0.018 + len(vals) * 0.05, 0.5, "tick: < %s" % fmt(SMALLEST[t]), va="center", fontsize=7,
            color=MUTED)
    x0 += 0.215
NT = sum(1 for L_ in LOCI if area(L_["type"], L_["size"]) < S_MIN)
for k, s_ in enumerate((
        "Circle area in proportion to %s;" % ("the rDNA bp in the locus" if USE_BP else "the locus span"),
        "45S and 5S on separate scales (45S max %s, 5S max %s)." % (fmt(BIG["45S"]), fmt(BIG["5S"])),
        "%d loci: %d 45S, %d 5S; %d too small to draw are ticks." % (
            len(LOCI), sum(1 for L_ in LOCI if L_["type"] == "45S"), sum(1 for L_ in LOCI if L_["type"] == "5S"), NT),
        "Labelled: 45S >= %g kb, 5S >= %g kb; loci within 3 Mb share a label." % (A.label45_kb, A.label5_kb))):
    lg.text(0.82, 0.93 - k * 0.25, s_, fontsize=7.6, va="center")
fig.text(0.088, 0.972, "D. paradoxa", fontsize=13.5, style="italic", va="baseline")
w_sp = REND.get_text_width_height_descent("D. paradoxa", FontProperties(size=13.5, style="italic"), False)[0]
fig.text(0.088 + w_sp / fig.dpi / fig.get_size_inches()[0], 0.972,
         ": the 45S and 5S rDNA loci on all 12 chromosomes, soma vs pollen", fontsize=13.5, va="baseline")
foot = ("Names: A = chr1_hap1 = L1 + P; B = chr2_hap1 = L2 + P; C = chr1_hap2 = L1 + Q; D = chr2_hap2 = L2 + Q; "
        "X = the arm shared by chr4_hap1 (X + core) and chr3_hap2 (chr3 + X). Pollen: A' = A's L1 + D's Q and "
        "D' = D's L2 + A's P, cut at the assembly's joins (A %.2f Mb, D %.2f Mb); moved arms %s. Loci: CARP 45S/5S "
        "pieces merged within 20 kb (%s); size = %s. Sizes are what is assembled: the largest 45S locus is %s, "
        "far less than a NOR holds, so 45S circles show where 45S sits, not how bright it will be; collapsed "
        "arrays read short. chr5/chr6 windows: %s. CO_smk %s, %s." % (
            cA / MB, cD / MB, ori_how, rd, SIZE_WHAT, fmt(BIG["45S"]), win_how, git_head(),
            datetime.date.today().isoformat()))
fig.text(0.088, 0.008, foot, fontsize=6.8, color=MUTED, wrap=True)
png = os.path.join(A.out, "rdna_ideogram_%s.png" % S)
fig.savefig(png, dpi=200)
fig.savefig(png[:-4] + ".pdf")

# ------------------------------------------------------------------ table and summary
tsv = os.path.join(A.out, "rdna_loci_both_karyotypes_%s.tsv" % S)
with open(tsv, "w") as fo:
    fo.write("type\tcopy\tstart\tend\tspan_kb\trdna_kb\tpieces\tsoma_mb\tpollen_chrom\tpollen_mb\n")
    for L_ in sorted(LOCI, key=lambda L_: (L_["type"], L_["copy"], L_["start"])):
        fo.write("%s\t%s\t%d\t%d\t%.1f\t%s\t%s\t%.3f\t%s\t%.3f\n" % (
            L_["type"], L_["copy"], L_["start"], L_["end"], L_["span"],
            "" if L_["bp_kb"] is None else "%.1f" % L_["bp_kb"], L_["pieces"], L_["soma"][1] / MB,
            NAME.get(L_["pollen"][0], L_["pollen"][0]), L_["pollen"][1] / MB))
print("rDNA IDEOGRAM  %s   size = %s" % (S, SIZE_WHAT))
print("  moved arms: %s; chr5/chr6 windows: %s" % (ori_how, win_how))
for n in NOTES:
    print("  NOTE: " + n)
print("  per chromosome: total size of all loci; [the labelled ones, Mb: size]")
for kar in ("soma", "pollen"):
    print("  %s" % kar.upper())
    for row, name, arms in ROWS[kar]:
        if kar == "pollen" and row not in NAME:
            continue
        parts = []
        for t in ("45S", "5S"):
            h = sorted([L_ for L_ in LOCI if L_["type"] == t and L_[kar][0] == row], key=lambda L_: L_[kar][1])
            thr = A.label45_kb if t == "45S" else A.label5_kb
            top = ", ".join("%.1f: %s" % (L_[kar][1] / MB, fmt(L_["size"])) for L_ in h if L_["size"] >= thr)
            parts.append("%s %s in %d loci%s" % (t, fmt(sum(L_["size"] for L_ in h)), len(h),
                                            (" [%s]" % top) if top else ""))
        print("    %-22s %s" % (name.replace("\u2032", "'"), "   ".join(parts)))
    if kar == "pollen":
        print("    (all other pollen chromosomes as in the soma)")
print("files read:\n   " + "\n   ".join(READ))
if MISSING:
    print("files missing:\n   " + "\n   ".join(MISSING))
print("wrote %s, %s and %s" % (png, png[:-4] + ".pdf", tsv))
