#!/usr/bin/env python3
"""x_arm_dosage.py -- how many copies of the X arm does each pollen nucleus carry, and does that follow from which
chr3 and chr4 it inherited?

In the two-haplotype assembly of D. paradoxa the X arm (~85 Mb) sits on chr4 in hap1 (chr4_hap1 = X + core) and on
chr3 in hap2 (chr3_hap2 = chr3 + X), joined to its chromosome by the same end. A pollen nucleus gets one chr3 and one
chr4, so it should carry X
  once, hap1's X   if it got chr4_hap1 and chr3_hap1   (balanced)
  once, hap2's X   if it got chr4_hap2 and chr3_hap2   (balanced)
  twice, both      if it got chr4_hap1 and chr3_hap2   (unbalanced: X duplicated)
  never            if it got chr4_hap2 and chr3_hap1   (unbalanced: X missing)
Which chr4 and chr3 a nucleus got is read next to the joins: the first NEAR Mb of chr4's core (next to X) and the
last NEAR Mb of chr3 (the end X hangs on in hap2), where a crossover between the window and the join is rare. How
many X copies it carries is read from X itself: molecules over the expected number (co_masks.py's 'coverage': ~1 for
one copy, ~2 for two, ~0 for none) and the alleles (one allele for one copy; both, interleaved, for two).

On the cut reference (Dparadoxa_split) X is the piece chr4_hap1_X: REF = hap1's X (from chr4_hap1), ALT = hap2's X
(from chr3_hap2); chr4_hap1_core: REF = chr4_hap1, ALT = chr4_hap2; chr3_hap1: REF = chr3_hap1, ALT = chr3_hap2.

Writes OUT/x_arm_SAMPLE.tsv (one row per nucleus) and OUT/x_arm_SAMPLE.png; prints the tables. Run from the CO_smk
root; a few seconds.
"""
import argparse
import csv
import os
import statistics

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser()
ap.add_argument("sample", nargs="?", default="Dparadoxa_split")
ap.add_argument("--x", default="chr4_hap1_X", help="the X piece (REF = hap1's X, ALT = hap2's X)")
ap.add_argument("--core", default="chr4_hap1_core", help="the piece X joins in hap1; read at its START")
ap.add_argument("--other", default="chr3_hap1", help="the piece whose hap2 copy carries X; read at its END")
ap.add_argument("--near_n", type=int, default=6, help="molecules nearest the join that say which copy a nucleus got")
ap.add_argument("--near_max_mb", type=float, default=60, help="never look further than this from the join")
ap.add_argument("--win_mb", type=float, default=5, help="windows along X in the figure")
ap.add_argument("--out", default="reports/x_arm")
A = ap.parse_args()
S = A.sample
MB = 1_000_000
NEAR, WIN = int(A.near_max_mb * MB), int(A.win_mb * MB)
MIN_MOL, CALL = 3, 0.2           # a side is called from >= 3 molecules with <= 20% (REF) or >= 80% (ALT) ALT
os.makedirs(A.out, exist_ok=True)
H1, H2, MIX, NONE, GREY = "#378ADD", "#D85A30", "#7F3C8D", "#2C2C2A", "#C9C7BE"
INK, MUTED = "#2C2C2A", "#6B6A64"

ref_dir = "results/reference/%s" % S
L = {}
for r in csv.DictReader(open(os.path.join(ref_dir, "pieces.tsv")), delimiter="\t"):
    L[r["piece"]] = int(r["length"]) if r.get("length") else int(r["end"]) - int(r["start"]) + 1
for p in (A.x, A.core, A.other):
    if p not in L:
        raise SystemExit("%s is not a piece of %s (pieces.tsv): %s" % (p, S, ", ".join(L)))
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
masks = {}
mp = "results/crossovers/%s/masks.tsv" % S
if os.path.exists(mp):
    for r in csv.DictReader(open(mp), delimiter="\t"):
        masks[(r["barcode"], r["chrom"])] = r


def side(mol):
    """Which copy, from the molecules nearest the join: [(distance, alt)] -> (REF|ALT|mixed|few, n, Mb reached).
    Takes molecules nearest-first until --near_n of them (or all within --near_max_mb); <= 20% ALT is REF,
    >= 80% ALT is ALT. A crossover between the join and these molecules would mislabel the copy, so near is better."""
    mol = sorted(mol)[:A.near_n]
    n = len(mol)
    if n < MIN_MOL:
        return "few", n, float("nan")
    s = sum(a for _, a in mol) / n
    return ("ALT" if s >= 1 - CALL else "REF" if s <= CALL else "mixed"), n, mol[-1][0] / MB


rows = []
nwx = -(-L[A.x] // WIN)
for bc in cells:
    p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        continue
    xm, core, oth, tot = [], [], [], 0
    with open(p) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rc, ac = int(float(f[3])), int(float(f[5]))
            if rc == ac:
                continue
            tot += 1
            c, pos, alt = f[0], int(float(f[1])), 1 if ac > rc else 0
            if c == A.x:
                xm.append((pos, alt))
            elif c == A.core and pos <= NEAR:
                core.append((pos, alt))
            elif c == A.other and pos > L[A.other] - NEAR:
                oth.append((L[A.other] - pos, alt))
    xm.sort()
    xa = [a for _, a in xm]
    sw = sum(1 for u, v in zip(xa, xa[1:]) if u != v) / (len(xa) - 1) if len(xa) > 1 else float("nan")
    win_n, win_alt = np.zeros(nwx), np.zeros(nwx)
    for pos, a in xm:
        k = min(pos // WIN, nwx - 1)
        win_n[k] += 1
        win_alt[k] += a
    c4, n4, s4 = side(core)
    c3, n3, s3 = side(oth)
    m = masks.get((bc, A.x), {})
    rows.append(dict(barcode=bc, total=tot, chr4=c4, chr4_n=n4, chr4_mb=s4, chr3=c3, chr3_n=n3, chr3_mb=s3,
                     x_n=len(xa),
                     x_alt=(sum(xa) / len(xa)) if xa else float("nan"), x_switch=sw,
                     x_cov=float(m["coverage"]) if m.get("coverage") not in (None, "", "nan") else float("nan"),
                     x_class=m.get("class", ""), win_n=win_n, win_alt=win_alt))
if not rows:
    raise SystemExit("no per-nucleus tables in results/cell_data_mol/%s" % S)
# coverage when masks.tsv lacks it: X molecules over (the nucleus's molecules x the median share of X)
med = statistics.median([r["x_n"] / r["total"] for r in rows if r["total"]])
for r in rows:
    if r["x_cov"] != r["x_cov"]:
        r["x_cov"] = r["x_n"] / max(r["total"] * med, 1e-9)
    r["expected_x"] = r["total"] * med


def predict(r):
    if r["chr4"] not in ("REF", "ALT") or r["chr3"] not in ("REF", "ALT"):
        return "unresolved"
    return {("REF", "REF"): "1: hap1's X", ("ALT", "ALT"): "1: hap2's X", ("REF", "ALT"): "2: both",
            ("ALT", "REF"): "0: none"}[(r["chr4"], r["chr3"])]


def observe(r):
    """co_masks.py's classes, read here as copies of X: too few | 0 | 1 (with its allele) | 2 | mixed"""
    if r["expected_x"] < 8:
        return "too few"
    if r["x_cov"] < 0.25:
        return "0: none"
    if r["x_switch"] >= 0.25:
        return "2: both" if r["x_cov"] >= 1.2 else "mixed"
    if r["x_alt"] <= CALL:
        return "1: hap1's X"
    if r["x_alt"] >= 1 - CALL:
        return "1: hap2's X"
    return "1: crossover in X"                  # one copy, one allele switch along X


GROUPS = ["1: hap1's X", "1: hap2's X", "2: both", "0: none", "unresolved"]
GCOL = {"1: hap1's X": H1, "1: hap2's X": H2, "2: both": MIX, "0: none": NONE, "unresolved": GREY}
for r in rows:
    r["predicted"], r["observed"] = predict(r), observe(r)
    r["agree"] = ("" if r["predicted"] == "unresolved" or r["observed"] in ("too few",)
                  else "yes" if r["observed"] == r["predicted"]
                  or (r["predicted"].startswith("1") and r["observed"] == "1: crossover in X") else "no")

with open(os.path.join(A.out, "x_arm_%s.tsv" % S), "w") as fo:
    fo.write("barcode\tmolecules\tchr4_next_to_X\tchr4_molecules\tchr4_reach_Mb\tchr3_end\tchr3_molecules\t"
             "chr3_reach_Mb\tpredicted_X\t"
             "X_molecules\tX_expected\tX_coverage\tX_alt_share\tX_switch\tobserved_X\tco_masks_class\tagree\t"
             "X_windows_alt_share_%gMb\n" % A.win_mb)
    for r in rows:
        ws = ",".join("%.2f" % (a / n) if n else "." for a, n in zip(r["win_alt"], r["win_n"]))
        fo.write("%s\t%d\t%s\t%d\t%.1f\t%s\t%d\t%.1f\t%s\t%d\t%.1f\t%.2f\t%.2f\t%.2f\t%s\t%s\t%s\t%s\n" % (
            r["barcode"], r["total"], r["chr4"], r["chr4_n"], r["chr4_mb"], r["chr3"], r["chr3_n"], r["chr3_mb"],
            r["predicted"], r["x_n"],
            r["expected_x"], r["x_cov"], r["x_alt"], r["x_switch"], r["observed"], r["x_class"], r["agree"], ws))

# ------------------------------------------------------------------ tables
n = len(rows)
reach = [r["chr4_mb"] for r in rows if r["chr4"] in ("REF", "ALT")] + \
        [r["chr3_mb"] for r in rows if r["chr3"] in ("REF", "ALT")]
print("%s: %d nuclei; X = %s (%.1f Mb). Which chr4 / chr3: the %d molecules nearest the join on %s (start) and %s "
      "(end), within %g Mb; median reach %.1f Mb" % (S, n, A.x, L[A.x] / MB, A.near_n, A.core, A.other,
                                                     A.near_max_mb, statistics.median(reach) if reach else float("nan")))
print("WHICH chr4 x WHICH chr3 (nuclei; REF = hap1 copy, ALT = hap2 copy; X is on chr4 in hap1, on chr3 in hap2):")
cnt = {}
for r in rows:
    cnt[(r["chr4"], r["chr3"])] = cnt.get((r["chr4"], r["chr3"]), 0) + 1
for a in ("REF", "ALT", "mixed", "few"):
    print("   chr4 %-5s: chr3 REF %3d   chr3 ALT %3d   chr3 mixed %3d   chr3 few %3d" % (
        a, cnt.get((a, "REF"), 0), cnt.get((a, "ALT"), 0), cnt.get((a, "mixed"), 0), cnt.get((a, "few"), 0)))
called = [r for r in rows if r["predicted"] != "unresolved"]
unbal = [r for r in called if r["predicted"][0] in "02"]
print("   %d nuclei with both sides called; %d of them (%.0f%%) got an unbalanced pair (X twice or never)" % (
    len(called), len(unbal), 100 * len(unbal) / max(len(called), 1)))
if called:
    r4 = sum(r["chr4"] == "REF" for r in called) / len(called)
    r3 = sum(r["chr3"] == "REF" for r in called) / len(called)
    print("   if chr3 and chr4 assorted independently: X twice ~%.0f, X never ~%.0f of %d (observed %d and %d);"
          " only balanced pairs = alternate segregation of a chr3-chr4 ring/cross" % (
              len(called) * r4 * (1 - r3), len(called) * (1 - r4) * r3, len(called),
              sum(r["predicted"] == "2: both" for r in called), sum(r["predicted"] == "0: none" for r in called)))
print("PREDICTED (from chr3 and chr4) vs OBSERVED (on X): nuclei, median X coverage, median ALT share on X")
obs_order = ["1: hap1's X", "1: hap2's X", "1: crossover in X", "2: both", "mixed", "0: none", "too few"]
for g in GROUPS:
    rs = [r for r in rows if r["predicted"] == g]
    if not rs:
        continue
    oc = {}
    for r in rs:
        oc[r["observed"]] = oc.get(r["observed"], 0) + 1
    covs = [r["x_cov"] for r in rs if r["x_cov"] == r["x_cov"]]
    alts = [r["x_alt"] for r in rs if r["x_alt"] == r["x_alt"]]
    print("   %-13s %3d nuclei: coverage %.2f, ALT share %.2f; observed %s" % (
        g, len(rs), statistics.median(covs) if covs else float("nan"),
        statistics.median(alts) if alts else float("nan"),
        ", ".join("%s %d" % (o, oc[o]) for o in obs_order if o in oc)))
ag = [r for r in rows if r["agree"] in ("yes", "no")]
print("   the X copies match the prediction in %d of %d nuclei with both read" % (
    sum(r["agree"] == "yes" for r in ag), len(ag)))
if masks:
    mc = {}
    for r in rows:
        mc[r["x_class"]] = mc.get(r["x_class"], 0) + 1
    print("   co_masks.py classes for %s: %s" % (A.x, ", ".join("%s %d" % kv for kv in sorted(mc.items()))))

# ------------------------------------------------------------------ figure
order = sorted(rows, key=lambda r: (GROUPS.index(r["predicted"]), r["x_alt"] if r["x_alt"] == r["x_alt"] else 2,
                                    r["x_cov"]))
H = max(6.0, 0.075 * n + 2.6)
fig = plt.figure(figsize=(15, H))
top, bot = 1 - 1.25 / H, 0.9 / H
axh = fig.add_axes([0.07, bot, 0.40, top - bot])
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("x", [H1, "#B9A9D6", H2])
for i, r in enumerate(order):
    for k in range(nwx):
        nn = r["win_n"][k]
        if nn:
            col = cmap(r["win_alt"][k] / nn)
            axh.add_patch(Rectangle((k, i), 1, 1, color=col, alpha=min(1.0, 0.35 + 0.65 * min(nn, 6) / 6), lw=0))
    for j, (key, cc) in enumerate((("chr4", "chr4"), ("chr3", "chr3"))):
        v = r[key]
        axh.add_patch(Rectangle((nwx + 0.6 + j * 1.2, i), 1, 1, lw=0,
                                color={"REF": H1, "ALT": H2, "mixed": MIX}.get(v, "#EFEDE6")))
    axh.add_patch(Rectangle((nwx + 3.2, i + 0.15), min(r["x_cov"], 3) * 1.5, 0.7, lw=0,
                            color=GCOL.get(r["predicted"], GREY)))
axh.set_xlim(0, nwx + 3.2 + 4.8)
axh.set_ylim(len(order), 0)
axh.set_yticks([])
xt = [k for k in range(0, nwx + 1, max(1, int(round(20 / A.win_mb))))]
axh.set_xticks(xt + [nwx + 1.1, nwx + 2.3, nwx + 3.2, nwx + 3.2 + 1.5, nwx + 3.2 + 3.0])
axh.set_xticklabels(["%d" % (k * A.win_mb) for k in xt] + ["chr4", "chr3", "0", "1", "2"], fontsize=7)
axh.set_xlabel("Mb along X (%s; X joins chr4's core / chr3 at %.0f Mb)            next to the joins     X coverage"
               % (A.x, L[A.x] / MB), fontsize=8)
for x_ in (nwx + 3.2 + 1.5, nwx + 3.2 + 3.0):
    axh.axvline(x_, color=MUTED, lw=0.6, ls=":")
y = 0
for g in GROUPS:
    k = sum(1 for r in order if r["predicted"] == g)
    if k:
        axh.axhline(y, color=INK, lw=0.8)
        axh.text(-0.3, y + k / 2, "%s\n(%d)" % (g, k), ha="right", va="center", fontsize=8, color=GCOL[g])
        y += k
for sp in ("top", "right", "left"):
    axh.spines[sp].set_visible(False)
axh.set_title("a. every nucleus (row), grouped by the X copies its chr4 and chr3 predict", loc="left", fontsize=9.5)

axs = fig.add_axes([0.60, bot + 0.08, 0.37, top - bot - 0.08])
for g in GROUPS:
    rs = [r for r in rows if r["predicted"] == g and r["x_alt"] == r["x_alt"]]
    axs.scatter([r["x_cov"] for r in rs], [r["x_alt"] for r in rs], s=26, color=GCOL[g], alpha=0.85,
                edgecolor="white", lw=0.5, label="%s (%d)" % (g, len(rs)), zorder=3)
for (cx, cy, t) in ((1, 0, "one copy,\nhap1's X"), (1, 1, "one copy,\nhap2's X"), (2, 0.5, "two copies")):
    axs.plot(cx, cy, marker="+", ms=18, color=MUTED, mew=1.2, zorder=2)
    axs.text(cx + 0.06, cy + (0.06 if cy < 0.9 else -0.1), t, fontsize=7.5, color=MUTED)
axs.axvspan(-0.1, 0.25, color="#EFEDE6", lw=0, zorder=1)
axs.text(0.02, 0.5, "no X", rotation=90, fontsize=7.5, color=MUTED, va="center")
axs.set_xlim(-0.1, max(3.0, max(r["x_cov"] for r in rows) * 1.05))
axs.set_ylim(-0.08, 1.08)
axs.set_xlabel("X molecules / expected for one copy (co_masks.py coverage)", fontsize=8.5)
axs.set_ylabel("share of X molecules with hap2's allele", fontsize=8.5)
axs.legend(fontsize=7.5, frameon=False, loc="lower right", title="predicted from chr4 + chr3", title_fontsize=8)
axs.set_title("b. X dosage and allele per nucleus", loc="left", fontsize=9.5)
for sp in ("top", "right"):
    axs.spines[sp].set_visible(False)
fig.suptitle("%s: copies of the X arm per pollen nucleus (%d nuclei), against the chr4 and chr3 each inherited" % (
    S, n), x=0.01, y=1 - 0.15 / H, ha="left", fontsize=11)
fig.text(0.01, 1 - 0.5 / H, "X sits on chr4 in hap1 and on chr3 in hap2. A nucleus with chr4_hap1 + chr3_hap1 or "
         "chr4_hap2 + chr3_hap2 should carry X once (blue: hap1's allele; orange: hap2's); chr4_hap1 + chr3_hap2 twice "
         "(both alleles, interleaved; ~2x molecules); chr4_hap2 + chr3_hap1 never (no X molecules).\nColour in a = share "
         "of hap2's allele per %g Mb window (paler = fewer molecules); the two squares = the copy of chr4 and of chr3 "
         "next to the joins (blue hap1, orange hap2); the bar = X molecules over the number expected for one copy."
         % A.win_mb, fontsize=7.8, color=MUTED, va="top")
p = os.path.join(A.out, "x_arm_%s.png" % S)
fig.savefig(p, dpi=140)
print("wrote %s and %s" % (p, os.path.join(A.out, "x_arm_%s.tsv" % S)))
