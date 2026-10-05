#!/usr/bin/env python3
"""rdna_overview.py -- where is the rDNA in D. paradoxa, and does it sit on the crossover
hotspots or the translocation joins?

rDNA: barrnap annotations (results/blacklist/<set>/rrna.gff) for whichever sets exist:
Dparadoxa_std (A, D, chr3-6 hap1), Dparadoxa_hap1 (A, B, chr3-6 hap1), hap2_rest
(C = chr1_hap2, chr3-6 hap2). Features < 50 kb apart form one array; 45S = the array
holds 18S/5.8S/28S, else 5S.
Crossovers: results/crossovers/Dparadoxa_std/co_intervals.bed (the crossover reference
A + D + chr3-6 hap1); interval midpoints counted in 5 Mb windows; hotspot = a window at
least 3 SD above the genome-wide mean.
Joins: the translocation joins, and the chr5 copy of the repeat at A's join.
Context: the share of the annotated genome within 5 Mb of any array, i.e. how often a
random point would land that close by chance.

Run from the CO_smk root. Usage: rdna_overview.py DUAL.fa.fai OUTDIR
Writes OUTDIR/rdna_overview.txt and OUTDIR/rdna_overview.png
"""
import collections
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FAI, OUT = sys.argv[1:3]
os.makedirs(OUT, exist_ok=True)
SETS = ["Dparadoxa_std", "Dparadoxa_hap1", "hap2_rest"]
CO_BED = "results/crossovers/Dparadoxa_std/co_intervals.bed"
ORDER = ["chr1_hap1", "chr1_hap2", "chr2_hap1", "chr2_hap2", "chr3_hap1", "chr3_hap2",
         "chr4_hap1", "chr4_hap2", "chr5_hap1", "chr5_hap2", "chr6_hap1", "chr6_hap2"]
NAME = {"chr1_hap1": "A  L1\u00b7P", "chr1_hap2": "C  L1\u00b7Q", "chr2_hap1": "B  L2\u00b7P", "chr2_hap2": "D  L2\u00b7Q"}
JOINS = [("chr1_hap1", 262.93, "A join L1|P"), ("chr2_hap2", 215.0, "D join L2|Q"),
         ("chr2_hap1", 211.0, "B join L2|P (approx.)"), ("chr1_hap2", 262.0, "C join L1|Q (approx.)"),
         ("chr4_hap1", 85.53, "chr4 piece|rest (hap1)"), ("chr3_hap2", 290.0, "chr3|chr4 piece (hap2, approx.)"),
         ("chr5_hap1", 61.1, "chr5 copy of the A-join repeat")]
WIN, NEAR, GAP = 5e6, 5.0, 50000

lens = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(FAI)}
annotated, arrays = set(), collections.defaultdict(list)       # chrom -> [(start, end, Counter)]
used = []
for S in SETS:
    p = "results/blacklist/%s/rrna.gff" % S
    if not os.path.exists(p):
        continue
    used.append(S)
    feats = collections.defaultdict(list)
    seqs = set()
    for l in open(p):
        if l.startswith("##sequence-region"):
            seqs.add(l.split()[1])
        if l.startswith("#") or not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        name = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv).get("Name", "?").replace("_rRNA", "")
        feats[f[0]].append((int(f[3]), int(f[4]), name))
    for c in set(feats) | seqs:
        if c in annotated or c not in lens:
            continue
        annotated.add(c)
        for s, e, n in sorted(feats.get(c, [])):
            if arrays[c] and s - arrays[c][-1][1] <= GAP:
                a = arrays[c][-1]
                arrays[c][-1] = (a[0], max(a[1], e), a[2] + collections.Counter({n: 1}))
            else:
                arrays[c].append((s, e, collections.Counter({n: 1})))
# a set's annotation covers every chromosome of its reference, with or without hits
for S, cs in (("Dparadoxa_std", ["chr1_hap1", "chr2_hap2", "chr3_hap1", "chr4_hap1", "chr5_hap1", "chr6_hap1"]),
              ("Dparadoxa_hap1", ["chr1_hap1", "chr2_hap1", "chr3_hap1", "chr4_hap1", "chr5_hap1", "chr6_hap1"]),
              ("hap2_rest", ["chr1_hap2", "chr3_hap2", "chr4_hap2", "chr5_hap2", "chr6_hap2"])):
    if S in used:
        annotated |= set(cs)


def kind(cnt):
    return "45S" if any(t in cnt for t in ("18S", "28S", "5_8S", "5.8S")) else ("5S" if "5S" in cnt else "other")


def nearest(c, mb):
    best = None
    for s, e, cnt in arrays.get(c, []):
        d = 0.0 if s / 1e6 <= mb <= e / 1e6 else min(abs(mb - s / 1e6), abs(mb - e / 1e6))
        if best is None or d < best[0]:
            best = (d, s, e, cnt)
    return best


out = ["rDNA IN D. PARADOXA: arrays, crossover hotspots, translocation joins",
       "  rDNA annotations used: %s; chromosomes annotated: %d of %d" % (", ".join(used) or "none", len(annotated), len(ORDER)), ""]
out.append("1. rDNA ARRAYS (features < 50 kb apart = one array)")
tot = collections.Counter()
for c in ORDER:
    if c not in annotated:
        out.append("  %-10s not annotated yet" % c)
        continue
    if not arrays.get(c):
        out.append("  %-10s no rDNA" % c)
    for s, e, cnt in arrays.get(c, []):
        tot[kind(cnt)] += sum(cnt.values())
        out.append("  %-10s %8.2f-%8.2f Mb %7.0f kb  %-4s %s" % (c, s / 1e6, e / 1e6, (e - s) / 1e3, kind(cnt),
                                                             ", ".join("%s x%d" % kv for kv in sorted(cnt.items()))))
out += ["  features in total: %s" % ", ".join("%s %d" % kv for kv in sorted(tot.items())), ""]

# crossover hotspots on the crossover reference
hot, dens = [], {}
if os.path.exists(CO_BED):
    mids = collections.defaultdict(list)
    for l in open(CO_BED):
        f = l.split("\t")
        if len(f) >= 3 and f[0] in lens:
            mids[f[0]].append((int(f[1]) + int(f[2])) / 2)
    allw = []
    for c, v in mids.items():
        h = np.bincount((np.array(v) // WIN).astype(int), minlength=int(lens[c] // WIN) + 1)
        dens[c] = h
        allw += [(n, c, i) for i, n in enumerate(h)]
    counts = np.array([n for n, _, _ in allw])
    mu, sd = counts.mean(), counts.std()
    out.append("2. CROSSOVER HOTSPOTS on the crossover reference (%d crossovers; 5 Mb windows; mean %.1f, SD %.1f per window)"
               % (sum(len(v) for v in mids.values()), mu, sd))
    for n, c, i in sorted(allw, reverse=True)[:10]:
        z = (n - mu) / sd if sd else 0
        nb = nearest(c, (i + 0.5) * WIN / 1e6)
        near = ("nearest rDNA: %s array %.1f Mb away" % (kind(nb[3]), nb[0])) if nb else (
            "no rDNA on this chromosome" if c in annotated else "chromosome not annotated")
        tag = "HOTSPOT" if z >= 3 else ""
        if z >= 3:
            hot.append((c, i))
        out.append("  %-10s %5.0f-%5.0f Mb  %3d crossovers  z %5.1f  %-7s  %s" % (c, i * WIN / 1e6, (i + 1) * WIN / 1e6, n, z,
                                                                            tag, near))
    out.append("")
else:
    out += ["2. CROSSOVER HOTSPOTS: %s not found" % CO_BED, ""]

out.append("3. TRANSLOCATION JOINS and the nearest rDNA")
for c, mb, lab in JOINS:
    if c not in annotated:
        out.append("  %-34s %-10s %7.1f Mb   chromosome not annotated yet" % (lab, c, mb))
        continue
    nb = nearest(c, mb)
    out.append("  %-34s %-10s %7.1f Mb   %s" % (lab, c, mb, ("nearest rDNA: %s, %s, %.1f Mb away" % (
        kind(nb[3]), ", ".join("%s x%d" % kv for kv in sorted(nb[3].items())), nb[0])) if nb else "no rDNA on this chromosome"))
cover = sum(sum(min(lens[c] / 1e6, e / 1e6 + NEAR) - max(0, s / 1e6 - NEAR) for s, e, _ in arrays.get(c, []))
            for c in annotated)
glen = sum(lens[c] / 1e6 for c in annotated)
nj = [nearest(c, mb) for c, mb, _ in JOINS if c in annotated]
nh = [nearest(c, (i + 0.5) * WIN / 1e6) for c, i in hot]
out += ["", "4. BY CHANCE?  share of the annotated genome within %g Mb of an rDNA array: %.1f%% (upper bound, overlaps "
        "counted twice)" % (NEAR, 100 * min(cover / glen, 1) if glen else float("nan")),
        "  joins within %g Mb of an array: %d of %d" % (NEAR, sum(1 for x in nj if x and x[0] <= NEAR), len(nj)),
        "  hotspots within %g Mb of an array: %d of %d" % (NEAR, sum(1 for x in nh if x and x[0] <= NEAR), len(nh))]
open(os.path.join(OUT, "rdna_overview.txt"), "w").write("\n".join(out) + "\n")
print("\n".join(out))

# figure: one bar per scaffold, rDNA ticks, crossover density on the crossover reference, joins
fig, ax = plt.subplots(figsize=(12, 7.2))
ymax = len(ORDER)
for k, c in enumerate(ORDER):
    y = ymax - k
    L = lens.get(c, 0) / 1e6
    ax.add_patch(plt.Rectangle((0, y - 0.12), L, 0.24, color="#D3D1C7" if c in annotated else "#F1EFE8",
                               ec="#888780", lw=0.5, hatch=None if c in annotated else "///"))
    ax.text(-8, y, "%s  %s" % (c, NAME.get(c, "")), ha="right", va="center", fontsize=8)
    for s, e, cnt in arrays.get(c, []):
        h = 0.14 + 0.06 * np.log10(sum(cnt.values()))
        ax.add_patch(plt.Rectangle((s / 1e6, y - h), max((e - s) / 1e6, 1.5), 2 * h,
                                   color="#534AB7" if kind(cnt) == "5S" else "#D85A30", lw=0))
    if c in dens:
        h = dens[c]
        x = (np.arange(len(h)) + 0.5) * WIN / 1e6
        ax.plot(x, y + 0.18 + 0.32 * h / max(1, max(v.max() for v in dens.values())), color="#0F6E56", lw=1)
        for cc, i in hot:
            if cc == c:
                ax.add_patch(plt.Rectangle((i * WIN / 1e6, y + 0.16), WIN / 1e6, 0.36, color="#EF9F27", alpha=0.35, lw=0))
    for cc, mb, lab in JOINS:
        if cc == c:
            ax.plot([mb, mb], [y - 0.42, y - 0.2], color="#2C2C2A", lw=1)
            ax.plot([mb], [y - 0.42], marker="^", color="#2C2C2A", ms=5)
ax.set_xlim(-140, max(lens.get(c, 0) for c in ORDER) / 1e6 + 10)
ax.set_ylim(0.3, ymax + 0.7)
ax.set_yticks([])
ax.set_xticks(range(0, int(max(lens.get(c, 0) for c in ORDER) / 1e6) + 1, 100))
ax.set_xlabel("Mb")
for s in ("left", "right", "top"):
    ax.spines[s].set_visible(False)
ax.set_title("D. paradoxa rDNA (purple 5S, orange 45S), crossover density on the crossover reference (green; "
             "hotspots shaded), joins (black triangles)", fontsize=9)
fig.text(0.01, 0.01, "Hatched: rDNA not annotated yet. Crossovers exist only for the crossover reference "
         "(chr1_hap1, chr2_hap2, chr3-6 hap1).", fontsize=8, color="#5F5E5A")
fig.tight_layout(rect=(0, 0.03, 1, 1))
fig.savefig(os.path.join(OUT, "rdna_overview.png"), dpi=150)
print("wrote %s/rdna_overview.txt and .png" % os.path.abspath(OUT))
