#!/usr/bin/env python3
"""plot_prereview.py -- plots for the cell-calling pre-review.

Draws, from the per-barcode table prereview_cells.py writes (and, where given,
a reference library's UMI totals and the pre-review report):
  01 barcode-rank curve, log-log, against the reference
  02 log10(UMI) histogram, both libraries, with the ambient hump marked
  03 genes detected vs UMIs
  04 rDNA-gene share and intronic share by UMI bin
  05 ambient share and expected wrong-allele rate vs nucleus UMIs
  06 where reads sit along gene bodies (read out of the report)
PNGs plus one combined PDF. The caveat line goes on every panel, so a plot
can never be read out of context later.
"""
import argparse
import gzip
import math
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

THRESH = (100, 300, 500, 1000)
C1, C2, GREY = "#D85A30", "#534AB7", "#888780"


def read_table(path):
    gf, ge, ng, rd = [], [], [], []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        head = f.readline().rstrip("\n").split("\t")
        col = {k: i for i, k in enumerate(head)}
        for l in f:
            t = l.rstrip("\n").split("\t")
            gf.append(float(t[col["genefull_umi"]]))
            ge.append(float(t[col["gene_umi"]]))
            ng.append(float(t[col["genes_detected"]]))
            rd.append(float(t[col["rdna_umi"]]))
    return gf, ge, ng, rd


def read_totals(path):
    with open(path) as f:
        return [float(l.split()[0]) for l in f if l.strip()]


def hump(vals, step=0.1, lo=1.5, hi=3.5):
    nb = int(6 / step)
    h = [0] * nb
    for v in vals:
        if v > 0:
            h[max(0, min(nb - 1, int(math.log10(v) / step)))] += 1
    sm = [(h[max(0, i - 1)] + h[i] + h[min(nb - 1, i + 1)]) / 3.0 for i in range(nb)]
    a, b = int(lo / step), int(hi / step)
    k = max(range(a, b), key=lambda i: sm[i])
    return 10 ** ((k + 0.5) * step)


def profile_from_report(path):
    """The 'where reads sit along genes' block of the report: labels and percents."""
    labs, vals = [], []
    try:
        for l in open(path):
            m = re.match(r"^\s{5}(5' flank 1 kb|gene\s+\d+-\s*\d+%|3' flank 1 kb)\s+([\d.]+)%", l)
            if m:
                labs.append(re.sub(r"\s+", " ", m.group(1)).replace(" 1 kb", "").replace("gene ", ""))
                vals.append(float(m.group(2)))
    except OSError:
        pass
    return (labs, vals) if len(vals) >= 3 else (None, None)


def finish(fig, caveat):
    if caveat:
        fig.text(0.5, 0.005, caveat, ha="center", va="bottom", fontsize=7, color=C1, wrap=True)
    fig.tight_layout(rect=(0, 0.045, 1, 1))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--table", required=True)
    ap.add_argument("--label", default="sample")
    ap.add_argument("--ref-totals")
    ap.add_argument("--ref-label", default="reference")
    ap.add_argument("--report")
    ap.add_argument("--caveat", default="")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    gf, ge, ng, rd = read_table(args.table)
    ref = read_totals(args.ref_totals) if args.ref_totals and os.path.exists(args.ref_totals) else None
    A = hump(gf)
    srt = sorted(gf, reverse=True)
    figs = []

    # ---- 01 rank curve
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.loglog(range(1, len(srt) + 1), srt, color=C1, lw=1.6, label=args.label)
    if ref:
        rs = sorted(ref, reverse=True)
        ax.loglog(range(1, len(rs) + 1), rs, color=C2, lw=1.6, label=args.ref_label)
    for t in THRESH:
        ax.axhline(t, color=GREY, lw=0.6, ls=":")
        ax.text(1.2, t * 1.05, "%d UMIs" % t, fontsize=7, color=GREY)
    ax.axhline(A, color=C1, lw=0.9, ls="--")
    ax.text(len(srt), A * 1.15, "ambient hump ~%d " % A, fontsize=7, color=C1, ha="right")
    ax.set_xlabel("barcode rank")
    ax.set_ylabel("UMIs (GeneFull)")
    ax.set_title("Barcode-rank curve")
    ax.legend(frameon=False)
    ax.grid(alpha=0.15, which="both")
    finish(fig, args.caveat)
    figs.append(("01_rank_curve", fig))

    # ---- 02 histogram
    n = 2 if ref else 1
    fig, axes = plt.subplots(n, 1, figsize=(7.5, 3 + 2 * n), sharex=True)
    axes = axes if n > 1 else [axes]
    bins = [i * 0.05 for i in range(int(6 / 0.05) + 1)]
    axes[0].hist([math.log10(v) for v in gf if v > 0], bins=bins, color=C1)
    axes[0].set_title("%s: barcodes per 0.05 log10(UMI) bin (table starts at 30 UMIs)" % args.label)
    axes[0].axvline(math.log10(A), color="black", lw=1, ls="--")
    axes[0].text(math.log10(A), axes[0].get_ylim()[1] * 0.9, " ambient hump ~%d" % A, fontsize=8)
    if ref:
        axes[1].hist([math.log10(v) for v in ref if v > 0], bins=bins, color=C2)
        axes[1].set_title("%s: barcodes per 0.05 log10(UMI) bin" % args.ref_label)
        axes[1].axvline(math.log10(hump(ref)), color="black", lw=1, ls="--")
    for a in axes:
        for t in THRESH:
            a.axvline(math.log10(t), color=GREY, lw=0.6, ls=":")
        a.set_yscale("log")
        a.grid(alpha=0.15)
    axes[-1].set_xlabel("log10(UMIs per barcode)   [dotted: 100, 300, 500, 1000]")
    axes[-1].set_xlim(0.5, max(4.0, math.log10(max(srt)) + 0.2))
    finish(fig, args.caveat)
    figs.append(("02_umi_histogram", fig))

    # ---- 03 genes vs UMIs
    fig, ax = plt.subplots(figsize=(7.5, 5))
    x = [math.log10(v) for v in gf if v > 0]
    y = [math.log10(max(1.0, g)) for v, g in zip(gf, ng) if v > 0]
    hb = ax.hexbin(x, y, gridsize=70, bins="log", cmap="magma_r", mincnt=1)
    fig.colorbar(hb, ax=ax, label="barcodes (log)")
    ax.set_xlabel("log10(UMIs)")
    ax.set_ylabel("log10(genes detected)")
    ax.set_title("Genes detected vs UMIs")
    for t in THRESH:
        ax.axvline(math.log10(t), color=GREY, lw=0.6, ls=":")
    finish(fig, args.caveat)
    figs.append(("03_genes_vs_umis", fig))

    # ---- 04 rDNA and intronic share by UMI bin
    edges = [30, 100, 300, 500, 1000, 2000, 5000, 10 ** 12]
    labs, rdna, intr, cnt = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = [(a, b, c) for a, b, c in zip(gf, ge, rd) if lo <= a < hi]
        if not sel:
            continue
        tf = sum(a for a, _, _ in sel)
        labs.append("%d-%d" % (lo, hi) if hi < 10 ** 12 else "%d+" % lo)
        rdna.append(100.0 * sum(c for _, _, c in sel) / tf)
        intr.append(100.0 * (1 - sum(b for _, b, _ in sel) / tf))
        cnt.append(len(sel))
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.bar(range(len(labs)), rdna, color=C1, label="rDNA-gene share of UMIs")
    ax.plot(range(len(labs)), intr, color=C2, marker="o", lw=1.6, label="intronic share")
    for i, c in enumerate(cnt):
        ax.text(i, max(rdna + intr) * 1.02, format(c, ","), ha="center", fontsize=7, color=GREY)
    ax.set_xticks(range(len(labs)))
    ax.set_xticklabels(labs, rotation=30, ha="right")
    ax.set_ylabel("% of UMIs in the bin")
    ax.set_xlabel("UMIs per barcode")
    ax.set_title("What the counted UMIs are, by barcode size")
    ax.legend(frameon=False)
    ax.grid(alpha=0.15, axis="y")
    finish(fig, args.caveat)
    figs.append(("04_rdna_intronic", fig))

    # ---- 05 contamination
    fig, ax = plt.subplots(figsize=(7.5, 5))
    xs = [10 ** (2 + 0.01 * k) for k in range(301)]
    ax.semilogx(xs, [100 * min(1.0, A / t) for t in xs], color=C1, lw=1.8, label="ambient share")
    ax.semilogx(xs, [100 * min(1.0, A / t) / 2 for t in xs], color=C2, lw=1.4, ls="--",
                label="expected wrong-allele reads")
    ax.axhline(20, color=GREY, lw=0.8, ls=":")
    ax.axvline(5 * A, color="black", lw=0.9, ls="--")
    ax.text(5 * A * 1.05, 60, "  %d UMIs: 20%% ambient\n  %s barcodes reach it"
            % (5 * A, format(sum(1 for v in gf if v >= 5 * A), ",")), fontsize=8)
    ax.set_xlabel("UMIs in a nucleus")
    ax.set_ylabel("%")
    ax.set_ylim(0, 100)
    ax.set_title("Ambient contamination at ~%d ambient UMIs per droplet" % A)
    ax.legend(frameon=False)
    ax.grid(alpha=0.15, which="both")
    finish(fig, args.caveat)
    figs.append(("05_contamination", fig))

    # ---- 06 read position along genes
    if args.report:
        labs, vals = profile_from_report(args.report)
        if labs:
            fig, ax = plt.subplots(figsize=(7.5, 5))
            ax.bar(range(len(vals)), vals, color=[C2] + [C1] * (len(vals) - 2) + [C2])
            ax.set_xticks(range(len(labs)))
            ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=8)
            ax.set_ylabel("% of reads")
            ax.set_xlabel("position along the gene, 5' -> 3' (flanks in blue)")
            ax.set_title("Where reads sit along gene bodies")
            ax.grid(alpha=0.15, axis="y")
            finish(fig, args.caveat)
            figs.append(("06_read_position", fig))

    with PdfPages(os.path.join(args.outdir, "prereview_%s.pdf" % args.label)) as pdf:
        for name, fig in figs:
            fig.savefig(os.path.join(args.outdir, "%s_%s.png" % (name, args.label)), dpi=130)
            pdf.savefig(fig)
            plt.close(fig)
    for name, _ in figs:
        print("  %s" % os.path.join(args.outdir, "%s_%s.png" % (name, args.label)))
    print("  %s" % os.path.join(args.outdir, "prereview_%s.pdf" % args.label))


if __name__ == "__main__":
    main()
