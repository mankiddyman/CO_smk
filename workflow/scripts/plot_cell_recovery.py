#!/usr/bin/env python3
"""plot_cell_recovery.py -- how many usable cells each library gives, side by side.

For each sample: the recovery funnel (called -> counted -> kept at each
switch-rate threshold -> enough markers), the switch-rate distribution with
its valley, cells kept as a function of the threshold, markers per kept cell,
and per-chromosome marker coverage of kept cells -- the numbers cell_qc's
min_total_markers and min_per_chrom_markers are chosen from.

Usage: plot_cell_recovery.py OUTDIR SAMPLE [SAMPLE ...]
"""
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = sys.argv[1]
SAMPLES = sys.argv[2:] or ["Dbinata_hap1", "Dparadoxa_hap1"]
COL = ["#D85A30", "#534AB7", "#2E7D32", "#888780"]
THR = [0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.15]
PER_CHROM = [20, 30, 40, 60, 100]
KEEP = 0.10
os.makedirs(OUT, exist_ok=True)


def load(s):
    rows, chroms = [], []
    with open(os.path.join("results/cell_qc", s, "switches.tsv")) as f:
        h = f.readline().rstrip("\n").split("\t")
        mk = [i for i, c in enumerate(h) if c.startswith("markers_")]
        chroms = [h[i][len("markers_"):] for i in mk]
        im, ir = h.index("total_markers"), h.index("switch_rate")
        for l in f:
            t = l.rstrip("\n").split("\t")
            try:
                m, r = int(t[im]), float(t[ir])
            except (ValueError, IndexError):
                continue
            if m > 0:
                rows.append((m, r, [int(float(t[i])) for i in mk]))
    called = None
    try:
        with open(os.path.join("results/cells", s, "qc_summary.tsv")) as f:
            h = f.readline().rstrip("\n").split("\t")
            called = int(float(f.readline().rstrip("\n").split("\t")[h.index("n_cells_default")]))
    except (OSError, ValueError, IndexError):
        pass
    return rows, chroms, called


def valley(rates, lo=0.05, hi=0.13, step=0.005):
    nb = int(0.5 / step) + 1
    h = [0] * nb
    for r in rates:
        h[min(nb - 1, int(r / step))] += 1
    sm = [(h[max(0, i - 1)] + h[i] + h[min(nb - 1, i + 1)]) / 3.0 for i in range(nb)]
    cand = [i for i in range(nb) if lo <= (i + 0.5) * step <= hi]
    return (min(cand, key=lambda i: sm[i]) + 0.5) * step if cand else float("nan")


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else float("nan")


data = {s: load(s) for s in SAMPLES}

# ---------------------------------------------------------------- table
lines = []
say = lambda s="": (lines.append(s), print(s))
say("CELL RECOVERY" + "".join("%18s" % s.replace("_hap1", "") for s in SAMPLES))
def row(label, f):
    say("  %-30s" % label + "".join("%18s" % f(s) for s in SAMPLES))
fmt = lambda x: format(x, ",") if isinstance(x, int) else ("%.3f" % x if x == x else "n/a")
row("called by emptyDrops", lambda s: fmt(data[s][2]) if data[s][2] is not None else "n/a")
row("counted by cellsnp", lambda s: fmt(len(data[s][0])))
row("median markers / cell", lambda s: fmt(int(med([r[0] for r in data[s][0]]))))
row("median switch rate", lambda s: fmt(med([r[1] for r in data[s][0]])))
row("histogram valley (0.05-0.13)", lambda s: fmt(valley([r[1] for r in data[s][0]])))
say()
say("  cells kept at switch rate <=")
for t in THR:
    row("    %.2f" % t, lambda s, t=t: fmt(sum(1 for r in data[s][0] if r[1] <= t)))
say()
say("  at <= %.2f, cells also passing:" % KEEP)
for lo in (800, 1500, 2500):
    row("    >= %d markers total" % lo,
        lambda s, lo=lo: fmt(sum(1 for r in data[s][0] if r[1] <= KEEP and r[0] >= lo)))
for pc in PER_CHROM:
    row("    >= %d markers on EVERY chrom" % pc,
        lambda s, pc=pc: fmt(sum(1 for r in data[s][0] if r[1] <= KEEP and min(r[2]) >= pc)))
row("  median markers, kept cells", lambda s: fmt(int(med([r[0] for r in data[s][0] if r[1] <= KEEP]) or 0)))
with open(os.path.join(OUT, "cell_recovery.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")

# ---------------------------------------------------------------- plots
fig, ax = plt.subplots(2, 3, figsize=(15, 9))

a = ax[0, 0]
bins = [i * 0.005 for i in range(0, 61)]
for k, s in enumerate(SAMPLES):
    rs = [r[1] for r in data[s][0]]
    a.hist(rs, bins=bins, histtype="step", lw=2, color=COL[k], density=True,
           label="%s (n=%s)" % (s.replace("_hap1", ""), format(len(rs), ",")))
    v = valley(rs)
    if v == v:
        a.axvline(v, color=COL[k], ls=":", lw=1)
a.axvline(KEEP, color="black", ls="--", lw=0.8)
a.set_xlabel("switch rate"); a.set_ylabel("density"); a.legend(frameon=False, fontsize=8)
a.set_title("switch-rate distributions (dotted: each valley)")

a = ax[0, 1]
xs = [0.02 + 0.0025 * i for i in range(73)]
for k, s in enumerate(SAMPLES):
    rs = sorted(r[1] for r in data[s][0])
    ys, j = [], 0
    for x in xs:
        while j < len(rs) and rs[j] <= x:
            j += 1
        ys.append(j)
    a.plot(xs, ys, color=COL[k], lw=2, label=s.replace("_hap1", ""))
    kept = sum(1 for r in rs if r <= KEEP)
    a.annotate(format(kept, ","), (KEEP, kept), textcoords="offset points", xytext=(6, -4),
               color=COL[k], fontsize=9)
a.axvline(KEEP, color="black", ls="--", lw=0.8)
a.set_xlabel("max_switch_rate"); a.set_ylabel("cells kept"); a.set_yscale("log")
a.legend(frameon=False, fontsize=8); a.set_title("cells kept as the threshold moves")

a = ax[0, 2]
stages = ["called", "counted", "<=0.12", "<=0.10", "<=0.10\n>=800", "<=0.10\n>=2500"]
w = 0.8 / len(SAMPLES)
for k, s in enumerate(SAMPLES):
    rows = data[s][0]
    vals = [data[s][2] or 0, len(rows),
            sum(1 for r in rows if r[1] <= 0.12), sum(1 for r in rows if r[1] <= 0.10),
            sum(1 for r in rows if r[1] <= 0.10 and r[0] >= 800),
            sum(1 for r in rows if r[1] <= 0.10 and r[0] >= 2500)]
    a.bar([i + k * w for i in range(len(stages))], vals, width=w, color=COL[k],
          label=s.replace("_hap1", ""))
    for i, v in enumerate(vals):
        a.text(i + k * w, v, format(v, ","), ha="center", va="bottom", fontsize=6, rotation=90)
a.set_xticks([i + w * (len(SAMPLES) - 1) / 2 for i in range(len(stages))])
a.set_xticklabels(stages, fontsize=8); a.set_yscale("log")
a.set_ylabel("cells"); a.legend(frameon=False, fontsize=8); a.set_title("recovery funnel")

a = ax[1, 0]
lb = [2 + 0.05 * i for i in range(61)]
for k, s in enumerate(SAMPLES):
    v = [math.log10(r[0]) for r in data[s][0] if r[1] <= KEEP]
    if v:
        a.hist(v, bins=lb, histtype="step", lw=2, color=COL[k], label=s.replace("_hap1", ""))
for x in (800, 2500):
    a.axvline(math.log10(x), color="grey", ls=":", lw=1)
a.set_xlabel("log10(markers), kept cells"); a.set_ylabel("cells")
a.legend(frameon=False, fontsize=8); a.set_title("markers per kept cell (dotted: 800, 2500)")

for k, s in enumerate(SAMPLES[:2]):
    a = ax[1, 1 + k]
    rows, chroms = data[s][0], data[s][1]
    kept = [r for r in rows if r[1] <= KEEP]
    if not kept:
        a.set_axis_off(); continue
    meds = [med([r[2][j] for r in kept]) for j in range(len(chroms))]
    lows = [sorted(r[2][j] for r in kept)[len(kept) // 10] for j in range(len(chroms))]
    xs = range(len(chroms))
    a.bar(xs, meds, color=COL[k], alpha=.8, label="median")
    a.plot(xs, lows, "k_", markersize=12, mew=2, label="10th percentile")
    for pc in (30, 100):
        a.axhline(pc, color="grey", ls=":", lw=1)
    a.set_xticks(list(xs)); a.set_xticklabels([c.replace("_hap1", "").replace("chr", "")
                                               for c in chroms], fontsize=8)
    a.set_yscale("log"); a.set_xlabel("chromosome"); a.set_ylabel("markers in a kept cell")
    a.legend(frameon=False, fontsize=8)
    a.set_title("%s: per-chromosome markers, %s kept cells" % (s.replace("_hap1", ""),
                                                               format(len(kept), ",")))

fig.tight_layout()
png = os.path.join(OUT, "cell_recovery.png")
fig.savefig(png, dpi=120); fig.savefig(png.replace(".png", ".pdf"))
print("\nwrote %s (+ .pdf, cell_recovery.txt)" % png)
