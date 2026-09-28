#!/usr/bin/env python3
"""allvsall_dotplot.py -- every chromosome of both haplotypes against every other, in one dot plot.

The hap2-on-hap1 map cannot show one hap1 chromosome against another, and it
hid anything with two equally good homologs (MAPQ 0). Here, pieces (100 kb,
one every 250 kb) from BOTH haplotypes are mapped against BOTH haplotypes
together (minimap2 asm20, approximate mapping, secondary hits kept), each
piece's hit to its own position is dropped, and everything else is plotted:
homologous pairs, translocations, inversions, and sequence shared between
chromosomes -- one square, laid out like a Hi-C map.

Two layouts: homologs interleaved (chr1_hap1, chr1_hap2, chr2_hap1, ...) and
by haplotype (all hap1, then all hap2). Plus a matrix: for each chromosome,
the share of its pieces with a hit on each other chromosome.

Usage: allvsall_dotplot.py SAMPLE HAP1_FASTA OUTDIR [--threads 24] [--piece 100000] [--step 250000] [--remap]
Reads: results/hap_align/SAMPLE/hap2.fa. Re-running replots from the saved hits unless --remap.
"""
import argparse
import collections
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("sample"); ap.add_argument("hap1_fasta"); ap.add_argument("outdir")
ap.add_argument("--threads", type=int, default=24)
ap.add_argument("--piece", type=int, default=100000)
ap.add_argument("--step", type=int, default=250000)
ap.add_argument("--hap2")
ap.add_argument("--remap", action="store_true")
A = ap.parse_args()
S, H1, OUT, PIECE, STEP = A.sample, A.hap1_fasta, A.outdir, A.piece, A.step
H2 = A.hap2 or "results/hap_align/%s/hap2.fa" % S
os.makedirs(OUT, exist_ok=True)
TGT, PCS = os.path.join(OUT, "both_haplotypes.fa"), os.path.join(OUT, "pieces.fa")
PAF, LENS = os.path.join(OUT, "pieces_vs_both.paf"), os.path.join(OUT, "chrom_lengths.tsv")


def is_chrom(n):
    return re.match(r"^chr\d+_hap[12]$", n) is not None


def cnum(c):
    return int(re.search(r"\d+", c).group())


if A.remap or not (os.path.exists(PAF) and os.path.getsize(PAF) > 0 and os.path.exists(LENS)):
    mm = sorted(glob.glob(".snakemake/conda/*/bin/minimap2"))
    MM = os.path.abspath(mm[0]) if mm else shutil.which("minimap2")
    if not MM:
        sys.exit("minimap2 not found")
    lens = {}
    npc = 0
    with open(TGT, "w") as ft, open(PCS, "w") as fp:
        def done(name, seq):
            global npc
            lens[name] = len(seq)
            for s0 in range(0, len(seq) - PIECE + 1, STEP):
                p = seq[s0:s0 + PIECE]
                if p.count("N") + p.count("n") > PIECE // 2:
                    continue
                fp.write(">%s__%d\n%s\n" % (name, s0, p))
                npc += 1
        for fa in (H1, H2):
            name, keep, parts = None, False, []
            with open(fa) as fh:
                for l in fh:
                    if l.startswith(">"):
                        if keep:
                            done(name, "".join(parts))
                        name = l[1:].split()[0]
                        keep, parts = is_chrom(name), []
                        if keep:
                            ft.write(">%s\n" % name)
                        continue
                    if keep:
                        ft.write(l)
                        parts.append(l.strip())
            if keep:
                done(name, "".join(parts))
    pd.Series(lens).to_csv(LENS, sep="\t", header=False)
    print("%d chromosomes (%.0f Mb), %d pieces of %d kb every %d kb; mapping with %d threads ..."
          % (len(lens), sum(lens.values()) / 1e6, npc, PIECE // 1000, STEP // 1000, A.threads), flush=True)
    cmd = [MM, "-x", "asm20", "-t", str(A.threads), "--secondary=yes", "-N", "20", "-p", "0.3", "-I", "16G",
           TGT, PCS]
    with open(PAF, "w") as fo, open(os.path.join(OUT, "minimap2.log"), "w") as fl:
        r = subprocess.run(cmd, stdout=fo, stderr=fl)
    if r.returncode:
        sys.exit("minimap2 failed -- see %s/minimap2.log" % OUT)
    os.remove(TGT)
lens = pd.read_csv(LENS, sep="\t", header=None, index_col=0)[1].to_dict()

# ---- hits, minus each piece's own position
MINHIT = PIECE // 5
hits, pieces_of = [], collections.Counter()
for l in open(PCS):
    if l.startswith(">"):
        pieces_of[l[1:].rsplit("__", 1)[0]] += 1
for l in open(PAF):
    f = l.split("\t")
    if len(f) < 12:
        continue
    qn, q0 = f[0].rsplit("__", 1)
    q0 = int(q0)
    t, ts, te = f[5], int(f[7]), int(f[8])
    if t == qn and abs(ts - (q0 + int(f[2]))) < PIECE:
        continue
    if int(f[3]) - int(f[2]) < MINHIT or t not in lens:
        continue
    hits.append((qn, q0 + (int(f[2]) + int(f[3])) // 2, f[4], t, (ts + te) // 2, q0))
h = pd.DataFrame(hits, columns=["qchrom", "qpos", "strand", "tchrom", "tpos", "piece"])
print("hits kept: %s (>= %d kb of a piece, self removed)" % (format(len(h), ","), MINHIT // 1000))

order_i = sorted(lens, key=lambda c: (cnum(c), c))
order_h = sorted(lens, key=lambda c: (c.endswith("hap2"), cnum(c)))

# ---- matrix: share of each chromosome's pieces with a hit on each other chromosome
pairs = h.drop_duplicates(["qchrom", "piece", "tchrom"]).groupby(["qchrom", "tchrom"]).size()
lines = ["ALL-VS-ALL  %s  (%d kb pieces every %d kb from both haplotypes, mapped to both; self removed)"
         % (S, PIECE // 1000, STEP // 1000),
         "  share of the ROW chromosome's pieces with a hit on the COLUMN chromosome (%)",
         "  %-10s " % "" + "".join("%9s" % c.replace("_hap", "h") for c in order_i)]
for q in order_i:
    lines.append("  %-10s " % q.replace("_hap", "h") + "".join(
        "%9s" % ("." if q == t else "%.0f" % (100.0 * pairs.get((q, t), 0) / max(pieces_of[q], 1)))
        for t in order_i))
open(os.path.join(OUT, "%s_allvsall.txt" % S), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))


def plot(order, fname, title):
    off, acc = {}, 0
    for c in order:
        off[c] = acc
        acc += lens[c]
    x = (h.qchrom.map(off) + h.qpos) / 1e6
    y = (h.tchrom.map(off) + h.tpos) / 1e6
    col = np.where(h.strand == "+", "#185FA5", "#D85A30")
    fig, ax = plt.subplots(figsize=(14, 14))
    ax.scatter(x, y, s=1.2, c=col, linewidths=0, alpha=.8, rasterized=True)
    for c in order:
        ax.axvline(off[c] / 1e6, color="#B4B2A9", lw=.6)
        ax.axhline(off[c] / 1e6, color="#B4B2A9", lw=.6)
        mid = (off[c] + lens[c] / 2.0) / 1e6
        ax.text(mid, -acc / 1e6 * 0.006, c.replace("_hap", " h"), ha="center", va="bottom", fontsize=8, rotation=90)
        ax.text(-acc / 1e6 * 0.006, mid, c.replace("_hap", " h"), ha="right", va="center", fontsize=8)
    ax.set_xlim(0, acc / 1e6); ax.set_ylim(acc / 1e6, 0)
    ax.set_xlabel("piece position (Mb)", labelpad=8); ax.set_ylabel("hit position (Mb)", labelpad=70)
    ax.set_title(title, pad=48)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, fname), dpi=150)


plot(order_i, "%s_allvsall_interleaved.png" % S,
     "%s: both haplotypes against both, homologs adjacent; blue + strand, orange - strand" % S)
plot(order_h, "%s_allvsall_by_haplotype.png" % S,
     "%s: both haplotypes against both, hap1 then hap2; blue + strand, orange - strand" % S)
print("wrote %s/%s_allvsall_{interleaved,by_haplotype}.png and %s_allvsall.txt" % (OUT, S, S))
