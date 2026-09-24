#!/usr/bin/env python3
"""Does the switch-rate floor come from the markers themselves?

Same cells, same reads: switch rate computed on assembly-confirmed markers
(a hap1-vs-hap2 difference at that position with those alleles) versus
unconfirmed ones. If the unconfirmed markers carry the noise, most "dirty"
cells are real nuclei damaged by bad markers, not ambient droplets.

Per-cell tables (cellsnp_to_per_cell): chrom pos ref ref_count alt alt_count,
no header.
"""
import collections
import math
import os
import random
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SAMPLE = sys.argv[1] if len(sys.argv) > 1 else "Dbinata_hap1"
ASM = sys.argv[2] if len(sys.argv) > 2 else "/tmp/binata_asm.vcf"
OUT = sys.argv[3] if len(sys.argv) > 3 else "qc/cell_qc/%s" % SAMPLE
N_PER = int(sys.argv[4]) if len(sys.argv) > 4 else 150
os.makedirs(OUT, exist_ok=True)

# ---- assembly SNVs: hap1 base, hap2 base
asm = collections.defaultdict(dict)
n_asm = 0
with open(ASM) as f:
    for l in f:
        if l.startswith("#"):
            continue
        t = l.split("\t", 6)
        if len(t) > 4 and len(t[3]) == 1 and len(t[4]) == 1:
            asm[t[0]][int(t[1])] = (t[3].upper(), t[4].upper())
            n_asm += 1
print("assembly SNVs: %s" % format(n_asm, ","))

# ---- cells, by their pipeline switch rate
sw = {}
with open("results/cell_qc/%s/switches.tsv" % SAMPLE) as f:
    h = f.readline().rstrip("\n").split("\t")
    ib, im, ir = h.index("barcode"), h.index("total_markers"), h.index("switch_rate")
    for l in f:
        t = l.rstrip("\n").split("\t")
        try:
            sw[t[ib]] = (int(t[im]), float(t[ir]))
        except (ValueError, IndexError):
            pass
groups = {
    "clean (<=0.10)": [b for b, (m, r) in sw.items() if r <= 0.10 and m >= 2000],
    "middle": [b for b, (m, r) in sw.items() if 0.10 < r < 0.13 and m >= 2000],
    "dirty (>=0.13)": [b for b, (m, r) in sw.items() if r >= 0.13 and m >= 2000],
}
rng = random.Random(1)
for g in groups:
    groups[g] = rng.sample(groups[g], min(N_PER, len(groups[g])))
    print("  %-15s %d cells sampled (>=2000 markers)" % (g, len(groups[g])))


def load(bc):
    rows = []
    with open("results/cell_data/%s/%s.tsv" % (SAMPLE, bc)) as f:
        for l in f:
            t = l.rstrip("\n").split("\t")
            if len(t) < 6:
                continue
            try:
                c, p, rb, rc, ab, ac = t[0], int(t[1]), t[2].upper(), int(t[3]), t[4].upper(), int(t[5])
            except ValueError:
                continue
            dp = rc + ac
            if dp == 0:
                continue
            fr = ac / float(dp)
            g = 0 if fr <= 0.2 else 1 if fr >= 0.8 else None
            if g is None:
                continue
            rows.append((c, p, g, asm.get(c, {}).get(p) == (rb, ab), dp))
    return rows


def rate(rows, keep):
    by = collections.defaultdict(list)
    for c, p, g, conf, dp in rows:
        if keep(conf, dp):
            by[c].append((p, g))
    s = n = 0
    for v in by.values():
        v.sort()
        for i in range(1, len(v)):
            n += 1
            s += v[i][1] != v[i - 1][1]
    return (s / float(n) if n >= 50 else float("nan")), n


TESTS = [("all", lambda c, d: True),
         ("confirmed", lambda c, d: c),
         ("unconfirmed", lambda c, d: not c),
         ("all, DP>=2", lambda c, d: d >= 2),
         ("confirmed, DP>=2", lambda c, d: c and d >= 2)]
res = []
for g, cells in groups.items():
    for bc in cells:
        rows = load(bc)
        if not rows:
            continue
        fc = sum(1 for r in rows if r[3]) / float(len(rows))
        rec = {"group": g, "barcode": bc, "called": len(rows), "frac_conf": fc}
        for name, keep in TESTS:
            rec[name], rec[name + " n"] = rate(rows, keep)
        res.append(rec)


def med(v):
    v = sorted(x for x in v if x == x)
    return v[len(v) // 2] if v else float("nan")


print("\nmedian switch rate per cell  (markers used, median)")
print("  %-15s %7s" % ("group", "confirmed%") + "".join("%18s" % t[0] for t in TESTS))
for g in groups:
    sub = [r for r in res if r["group"] == g]
    if not sub:
        continue
    line = "  %-15s %9.1f%%" % (g, 100 * med(r["frac_conf"] for r in sub))
    for name, _ in TESTS:
        line += "   %.3f (%6d)" % (med(r[name] for r in sub), med(r[name + " n"] for r in sub))
    print(line)

with open(os.path.join(OUT, "marker_switch_test.tsv"), "w") as f:
    cols = ["group", "barcode", "called", "frac_conf"] + [x for t in TESTS for x in (t[0], t[0] + " n")]
    f.write("\t".join(cols) + "\n")
    for r in res:
        f.write("\t".join(str(r[c]) for c in cols) + "\n")

# ---- plots
C = {"clean (<=0.10)": "#2E7D32", "middle": "#888780", "dirty (>=0.13)": "#D85A30"}
fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))

a = ax[0, 0]
bins = [i * 0.005 for i in range(0, 61)]
for name, col in (("all", "#D85A30"), ("confirmed", "#2E7D32"), ("unconfirmed", "#534AB7")):
    a.hist([r[name] for r in res if r[name] == r[name]], bins=bins, histtype="step", lw=2,
           color=col, label=name)
a.axvline(0.10, color="black", ls="--", lw=0.8)
a.set_xlabel("switch rate"); a.set_ylabel("cells"); a.legend(frameon=False)
a.set_title("same cells, three marker sets")

a = ax[0, 1]
for g, col in C.items():
    sub = [r for r in res if r["group"] == g]
    a.scatter([r["all"] for r in sub], [r["confirmed"] for r in sub], s=10, color=col, label=g, alpha=.7)
a.plot([0, 0.3], [0, 0.3], color="grey", lw=0.8, ls=":")
a.axhline(0.10, color="black", ls="--", lw=0.8)
a.set_xlabel("switch rate, all markers"); a.set_ylabel("switch rate, confirmed markers only")
a.set_title("does confirmation rescue the dirty cells?"); a.legend(frameon=False, fontsize=8)

a = ax[1, 0]
names = [t[0] for t in TESTS]
pos = 0
for g, col in C.items():
    sub = [r for r in res if r["group"] == g]
    data = [[r[n] for r in sub if r[n] == r[n]] for n in names]
    bp = a.boxplot(data, positions=[pos + i for i in range(len(names))], widths=0.6,
                   patch_artist=True, showfliers=False)
    for b in bp["boxes"]:
        b.set_facecolor(col); b.set_alpha(.6)
    pos += len(names) + 1
a.set_xticks([i for i in range(pos) if i % (len(names) + 1) != len(names)])
a.set_xticklabels(names * len(C), rotation=60, ha="right", fontsize=7)
a.axhline(0.10, color="black", ls="--", lw=0.8)
a.set_ylabel("switch rate"); a.set_title("by group: clean | middle | dirty")

a = ax[1, 1]
for g, col in C.items():
    sub = [r for r in res if r["group"] == g]
    a.scatter([100 * r["frac_conf"] for r in sub], [r["all"] for r in sub], s=10, color=col, alpha=.7)
a.set_xlabel("% of a cell's called markers that are assembly-confirmed")
a.set_ylabel("switch rate, all markers"); a.set_title("more confirmed markers, cleaner cell?")

fig.tight_layout()
pdf = os.path.join(OUT, "marker_switch_test.pdf")
fig.savefig(pdf); fig.savefig(pdf.replace(".pdf", ".png"), dpi=120)
print("\nwrote %s (+ .png, .tsv)" % pdf)
