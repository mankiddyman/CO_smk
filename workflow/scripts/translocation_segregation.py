#!/usr/bin/env python3
"""translocation_segregation.py -- which cells carry a translocated arm twice or not at all?

A reciprocal translocation heterozygote makes balanced gametes (each arm once)
and, by adjacent segregation, unbalanced ones: one translocated arm twice, its
reciprocal arm not at all. On the mapping reference a duplicated arm reads as
BOTH alleles interleaved, a missing arm as no molecules, and the jump from a
haploid arm into either one looks like a crossover at the breakpoint.

Translocated segments are found from the assembly, not typed in: every 1 Mb
window of a reference chromosome takes the chromosome NUMBER of its homolog on
the other haplotype (all-vs-all hits, 100 kb pieces, >= 50 kb aligned); runs
of >= MIN_MB whose homolog carries a different number are translocated.

Per cell and segment, from the molecule tables (one row per informative molecule):
  coverage   molecules / expected (cell total x the median share across cells)
  switches   allele changes between neighbouring molecules / (molecules - 1):
             ~0 for one haplotype (a crossover adds one), ~0.5 for two interleaved
  class      missing (coverage < 0.25); duplicated (switches >= 0.25 and
             coverage >= 1.2: both copies present); mixed (switches >= 0.25 at
             normal coverage: contamination, not segregation); haploid; or low
             (too few molecules expected to tell)
Crossovers whose midpoint lies within BP_MB of a segment boundary are
breakpoint crossovers, tallied by the cell's class for that segment.

Usage: translocation_segregation.py SAMPLE ALLVSALL_DIR OUTDIR [--min_mb 10] [--bp_mb 5]
"""
import argparse
import collections
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("sample"); ap.add_argument("allvsall"); ap.add_argument("outdir")
ap.add_argument("--min_mb", type=float, default=10)
ap.add_argument("--bp_mb", type=float, default=5)
ap.add_argument("--piece", type=int, default=100000)
A = ap.parse_args()
S, W = A.sample, 1000000
os.makedirs(A.outdir, exist_ok=True)


def cnum(c):
    return int(re.search(r"\d+", c).group())


# --- translocated segments, from the homolog's chromosome number per window
import csv
_rows = list(csv.DictReader(open("config/samples.csv")))
_row = next((r for r in _rows if S in (r.get("sample_id"), r.get("sample"), list(r.values())[0])), None)
if _row is None:
    raise SystemExit("%s not found in config/samples.csv" % S)
fai = _row["assembly_fasta"] + ".fai"
ref = [l.split("\t")[0] for l in open(fai)]
lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
rs = set(ref)
votes = collections.defaultdict(collections.Counter)
for l in open(os.path.join(A.allvsall, "pieces_vs_both.paf")):
    f = l.split("\t")
    if len(f) < 12:
        continue
    q, s0 = f[0].rsplit("__", 1)
    if q not in rs or f[5] in rs or int(f[3]) - int(f[2]) < A.piece // 2:
        continue
    votes[(q, int(s0) // W)][cnum(f[5])] += 1
segments = []
for c in ref:
    nw = lens[c] // W + 1
    lab = [votes[(c, w)].most_common(1)[0][0] if votes.get((c, w)) else None for w in range(nw)]
    known = [(w, k) for w, k in enumerate(lab) if k is not None]
    runs = []                                    # [partner, first window, last window, n windows]
    for w, k in known:
        if runs and runs[-1][0] == k:
            runs[-1][2], runs[-1][3] = w, runs[-1][3] + 1
        else:
            runs.append([k, w, w, 1])
    changed = True
    while changed:                               # absorb short runs (< 5 windows) into a neighbour
        changed = False
        for i, r in enumerate(runs):
            if r[3] < 5 and len(runs) > 1:
                j = i - 1 if i > 0 else i + 1
                runs[j][1], runs[j][2] = min(runs[j][1], r[1]), max(runs[j][2], r[2])
                runs[j][3] += r[3]
                del runs[i]
                merged = []
                for x in runs:
                    if merged and merged[-1][0] == x[0]:
                        merged[-1][2], merged[-1][3] = x[2], merged[-1][3] + x[3]
                    else:
                        merged.append(x)
                runs, changed = merged, True
                break
    for k, a, b, n in runs:
        if k != cnum(c) and (b - a + 1) >= A.min_mb:
            segments.append({"chrom": c, "start": a * W, "end": min((b + 1) * W, lens[c]), "homolog_on": "chr%d" % k})
seg = pd.DataFrame(segments)
if seg.empty:
    raise SystemExit("no translocated segments found on %s" % ", ".join(ref))
seg["name"] = ["%s:%.0f-%.0f" % (r.chrom, r.start / 1e6, r.end / 1e6) for r in seg.itertuples()]
boundaries = []                                  # breakpoints: segment ends that are not chromosome ends
for r in seg.itertuples():
    for x in (r.start, r.end):
        if 0 < x < lens[r.chrom] - W:
            boundaries.append((r.chrom, x, r.name))

# --- per cell, per segment
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
rows = []
for bc in cells:
    p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
    if not os.path.exists(p):
        continue
    m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"])
    m = m[m.rc != m.ac]
    m["alt"] = (m.ac > m.rc).astype(int)
    total = len(m)
    for r in seg.itertuples():
        x = m[(m.chrom == r.chrom) & (m.pos >= r.start) & (m.pos < r.end)].sort_values("pos")
        n = len(x)
        sw = int((np.diff(x.alt.values) != 0).sum()) if n > 1 else 0
        rows.append({"barcode": bc, "segment": r.name, "total": total, "n": n,
                     "alt_frac": x.alt.mean() if n else np.nan, "switch_rate": sw / (n - 1) if n > 1 else np.nan})
t = pd.DataFrame(rows)
share = (t.n / t.total).groupby(t.segment).median()
t["expected"] = t.total * t.segment.map(share)
t["coverage"] = t.n / t.expected.replace(0, np.nan)


def classify(r):
    if r.expected < 8:
        return "low"
    if r.coverage < 0.25:
        return "missing"
    if r.switch_rate >= 0.25:
        return "duplicated" if r.coverage >= 1.2 else "mixed"
    return "haploid"


t["class"] = t.apply(classify, axis=1)
t.to_csv(os.path.join(A.outdir, "segregation_cells.tsv"), sep="\t", index=False, float_format="%.3f")

# --- crossovers at the breakpoints, by the cell's class for that segment
cls = t.set_index(["barcode", "segment"])["class"]
bp_rows, all_cos = [], collections.Counter()
for bc in t.barcode.unique():
    p = "results/crossovers/%s/per_cell/%s_co_pred.txt" % (S, bc)
    if not os.path.exists(p):
        continue
    for l in open(p):
        f = l.split()
        if len(f) < 3 or not f[1].isdigit() or f[0] not in lens:
            continue
        mid = (int(f[1]) + int(f[2])) / 2.0
        all_cos[f[0]] += 1
        for c, x, name in boundaries:
            if f[0] == c and abs(mid - x) <= A.bp_mb * 1e6:
                bp_rows.append((bc, name, c, x, cls.get((bc, name), "n/a")))
bp = pd.DataFrame(bp_rows, columns=["barcode", "segment", "chrom", "breakpoint", "class"])
ncell = t.barcode.nunique()

L = ["TRANSLOCATION SEGREGATION  %s  (%d cells)" % (S, ncell),
     "  translocated segments (homolog on another chromosome number, runs >= %g Mb):" % A.min_mb]
for r in seg.itertuples():
    L.append("    %-22s homolog on %-5s  median share of a cell's molecules %.3f" % (r.name, r.homolog_on, share[r.name]))
L.append("  cells by class, per segment:")
tab = t.groupby(["segment", "class"]).size().unstack(fill_value=0)
for s_, row in tab.iterrows():
    L.append("    %-22s " % s_ + "   ".join("%s %d" % (k, v) for k, v in row.items() if v))
L.append("  combinations across segments (one line per pattern):")
combo = t.pivot(index="barcode", columns="segment", values="class").fillna("n/a")
for pat, n in combo.apply(lambda r: " | ".join("%s=%s" % (k, v) for k, v in r.items()), axis=1).value_counts().items():
    L.append("    %3d cells   %s" % (n, pat))
L.append("  crossovers within %g Mb of a breakpoint, by the cell's class for that segment:" % A.bp_mb)
if bp.empty:
    L.append("    none")
else:
    for (name, x), g in bp.groupby(["segment", "breakpoint"]):
        k = g["class"].value_counts()
        ncls = t[t.segment == name]["class"].value_counts()
        L.append("    %-22s at %5.0f Mb: %3d COs  (%s)" % (name, x / 1e6, len(g), ", ".join(
            "%d in %d %s cells" % (k[c], ncls.get(c, 0), c) for c in k.index)))
L.append("  COs per cell by chromosome: all | without breakpoint COs")
bpc = bp.groupby("chrom").size() if not bp.empty else pd.Series(dtype=int)
for c in ref:
    L.append("    %-11s %.2f | %.2f" % (c, all_cos[c] / ncell, (all_cos[c] - bpc.get(c, 0)) / ncell))
L.append("    total       %.2f | %.2f" % (sum(all_cos.values()) / ncell, (sum(all_cos.values()) - len(bp)) / ncell))
open(os.path.join(A.outdir, "segregation_summary.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))

colors = {"haploid": "#2E7D32", "duplicated": "#D85A30", "missing": "#534AB7", "mixed": "#BA7517", "low": "#B4B2A9"}
fig, ax = plt.subplots(1, len(seg), figsize=(4.2 * len(seg), 3.8), squeeze=False)
for a, r in zip(ax[0], seg.itertuples()):
    x = t[t.segment == r.name]
    for k, g in x.groupby("class"):
        a.scatter(g.coverage, g.switch_rate.fillna(0), s=18, color=colors.get(k, "k"), label="%s (%d)" % (k, len(g)))
    a.axvline(0.25, color="k", lw=.5, ls=":"); a.axvline(1.2, color="k", lw=.5, ls=":"); a.axhline(0.25, color="k", lw=.5, ls=":")
    a.set_title(r.name + "  (homolog on %s)" % r.homolog_on, fontsize=9)
    a.set_xlabel("coverage (observed / expected molecules)"); a.set_ylabel("allele switch rate")
    a.legend(frameon=False, fontsize=7)
fig.suptitle("%s: translocated segments per cell -- haploid, duplicated (both alleles), missing" % S, fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(A.outdir, "segregation.png"), dpi=130)
print("\nwrote %s/{segregation_summary.txt,segregation_cells.tsv,segregation.png}" % A.outdir)
