#!/usr/bin/env python3
"""co_sweep.py -- hapCO's calling over a grid of settings for every selected cell, in seconds.

Uses hapco_port.py (hapCO_identification.R transcribed, no plots). A cell is smoothed once per baseAF and cut
into genotype blocks once; every block setting (marker_num x block_size x terminal_marker_num) then only
filters those blocks, which is cheap.

  --check   first, the port at the settings the pipeline ran with (co_calling_params.txt) against the
            pipeline's own per-cell files (results/crossovers/SAMPLE/per_cell/*_co_pred.txt and
            *_co_block_pred.txt), cell by cell. Any difference -> exit 1: the sweep is not to be trusted.
Outputs (OUTDIR, default qc/co_sweep/SAMPLE):
  settings.tsv     id, base_af, marker_num, block_size, terminal_marker_num, label
  cells.tsv        barcode, molecules (informative rows), called (1 = at least cell_markers, so hapCO calls it)
  calls.tsv.gz     setting, barcode, chrom, start, end, pre, post, cstart, cend: one row per crossover.
                   start-end: the two molecules hapCO puts the crossover between (its fine breakpoint);
                   cstart-cend: the gap between the two blocks that survived the filter (where the crossover
                   could be); pre/post: genotype before and after (1 = REF, 0 = ALT)
  blocks.tsv.gz    the final blocks of the current setting and the --show settings (setting, barcode, chrom,
                   start, end, genotype)
  summary.tsv      per setting: cells, crossovers per cell, per chromosome (pieces summed into their source
                   chromosome), chromosomes below 0.5, close doubles, clean vs less clean cells

Usage: co_sweep.py SAMPLE [--outdir D] [--threads 8] [--check] [--base_af 0.2,0.3] [--marker_num 3,4,5,6,8]
                   [--block_mb 1,2,3,5] [--terminal 4,8,12] [--show b0.3_mn4_bs2_t8,...]
Run from the CO_smk directory.
"""
import argparse
import csv
import gzip
import itertools
import os
import re
import sys
import time
from multiprocessing import get_context

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hapco_port as H  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("--outdir")
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--check", action="store_true")
ap.add_argument("--base_af", default="0.2,0.3")
ap.add_argument("--marker_num", default="3,4,5,6,8")
ap.add_argument("--block_mb", default="1,2,3,5")
ap.add_argument("--terminal", default="4,8,12")
ap.add_argument("--show", default="")
A = ap.parse_args()
S = A.sample
OUT = A.outdir or "qc/co_sweep/%s" % S
os.makedirs(OUT, exist_ok=True)
T0 = time.time()


def as_run():
    """The settings the pipeline's calls were made with, and the per-cell input it read."""
    p = "results/crossovers/%s/co_calling_params.txt" % S
    d = {}
    for l in open(p):
        m = re.match(r"\s*([A-Za-z_ ]+?)(?:\s*\(-\w\))?:\s*(.*)", l)
        if m:
            d[m.group(1).strip()] = m.group(2).strip()
    prm = dict(base_af=float(d["baseAF"]), window_af=float(d["windowAF"]), genotype=float(d["genotype"]),
               block_size=float(d["block_size"]), marker_num=int(float(d["marker_num"])),
               terminal_marker_num=int(float(d.get("terminal_marker_num", 0))), cell_markers=int(float(d["cell_markers"])))
    indir = d.get("input tables") or "results/cell_data/%s" % S
    return prm, indir


def sid(b, mn, bs, t):
    return "b%g_mn%d_bs%g_t%d" % (b, mn, bs / 1e6, t)


RUN, INDIR = as_run()
CMAP = H.read_chrom_map("results/cell_data/%s/chrom_map.tsv" % S)
CELLS = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
CURRENT = sid(RUN["base_af"], RUN["marker_num"], RUN["block_size"], RUN["terminal_marker_num"])
MISSING = [bc for bc in CELLS if not os.path.exists("%s/%s.tsv" % (INDIR, bc))]
if MISSING:
    print("WARNING  %d of %d selected cells have no table in %s (first: %s); they are left out"
          % (len(MISSING), len(CELLS), INDIR, MISSING[0]))
    CELLS = [bc for bc in CELLS if bc not in set(MISSING)]


def load_pipeline(bc):
    out = []
    for suf in ("_co_pred.txt", "_co_block_pred.txt"):
        p = "results/crossovers/%s/per_cell/%s%s" % (S, bc, suf)
        rows = []
        if os.path.exists(p):
            for l in open(p):
                f = l.split()
                if f:
                    rows.append(tuple([f[0]] + [float(x) for x in f[1:]]))
        out.append(rows)
    return out


def check_cell(bc):
    res = H.call_cell(H.read_cell("%s/%s.tsv" % (INDIR, bc)), CMAP, RUN)
    co, bk = [], []
    if res:
        for c, (bps, final) in res.items():
            co += [tuple([c] + [float(x) for x in b]) for b in bps]
            bk += [tuple([c] + [float(x) for x in b]) for b in final]
    pc, pb = load_pipeline(bc)
    return bc, (co == pc and bk == pb), len(pc), (co, pc)


if A.check:
    n_files = len([f for f in os.listdir("results/crossovers/%s/per_cell" % S) if f.endswith("_co_block_pred.txt")]) \
        if os.path.isdir("results/crossovers/%s/per_cell" % S) else 0
    if n_files == 0:
        print("CHECK  no per-cell hapCO output in results/crossovers/%s/per_cell: nothing to check against" % S)
        sys.exit(1)
    with get_context("fork").Pool(A.threads) as pool:   # fork: no re-import of this script
        r = pool.map(check_cell, CELLS)
    same = [x for x in r if x[1]]
    print("CHECK  port vs the pipeline's hapCO output at the settings it ran with (%s; input %s)" % (CURRENT, INDIR))
    print("  identical in %d of %d cells (%d crossovers in the pipeline's calls)" % (len(same), len(r), sum(x[2] for x in r)))
    for bc, ok, n, (co, pc) in r:
        if not ok:
            print("  DIFFERS %s\n    port:     %s\n    pipeline: %s" % (bc, co[:4], pc[:4]))
            break
    if len(same) != len(r):
        sys.exit(1)

BAFS = [float(x) for x in A.base_af.split(",")]
MNS = [int(x) for x in A.marker_num.split(",")]
BSS = [float(x) * 1e6 for x in A.block_mb.split(",")]
TS = [int(x) for x in A.terminal.split(",")]
GRID = [(b, mn, bs, t) for b, mn, bs, t in itertools.product(BAFS, MNS, BSS, TS)]
if (RUN["base_af"], RUN["marker_num"], RUN["block_size"], RUN["terminal_marker_num"]) not in GRID:
    GRID.append((RUN["base_af"], RUN["marker_num"], RUN["block_size"], RUN["terminal_marker_num"]))
IDS = [sid(*g) for g in GRID]
SHOW = [CURRENT] + [s for s in A.show.split(",") if s and s != CURRENT]
for s in SHOW:
    if s not in IDS:
        print("WARNING  --show %s is not in the grid; no blocks written for it" % s)


def sweep_cell(bc):
    """hapCO's blocks_to_calls for every setting, with the smoothing and genotype blocks shared per baseAF."""
    rows = H.read_cell("%s/%s.tsv" % (INDIR, bc))
    prep = H.prepare(rows, CMAP, RUN["cell_markers"])
    n_inf = sum(1 for r in rows if r[2] > 0 or r[3] > 0)
    calls, blocks = [], []
    if prep is None:
        return bc, n_inf, calls, blocks
    for b in sorted(set(g[0] for g in GRID)):
        sm = H.smooth_cell(prep, b, RUN["window_af"])
        gb = {c: H.genotype_blocks(pos, smt, RUN["genotype"]) for c, (pos, smt) in sm.items()}
        for (bb, mn, bs, t) in GRID:
            if bb != b:
                continue
            s = sid(b, mn, bs, t)
            for c, size in CMAP:
                if c not in sm:
                    continue
                pos, smt = sm[c]
                filtered = H.filter_blocks(gb[c], bs, mn, t)
                if not filtered:
                    continue
                bps = H.get_breakpoints(filtered)
                bpu = H.fine_breakpoints(bps, pos, smt) if bps else []
                for coarse, fine in zip(bps, bpu):
                    calls.append((s, bc, c, int(fine[0]), int(fine[1]), int(fine[2]), int(fine[3]),
                                  int(coarse[0]), int(coarse[1])))
                if s in SHOW:
                    final = H.make_final_blocks(filtered, size)
                    if bps:
                        final = H.fine_blocks(bpu, final)
                    for st, en, g, _ in final:
                        blocks.append((s, bc, c, int(st), int(en), int(g)))
    return bc, n_inf, calls, blocks


with get_context("fork").Pool(A.threads) as pool:   # fork: no re-import of this script
    out = pool.map(sweep_cell, CELLS)
with open(os.path.join(OUT, "settings.tsv"), "w") as fo:
    fo.write("id\tbase_af\tmarker_num\tblock_size\tterminal_marker_num\tlabel\n")
    for b, mn, bs, t in GRID:
        s = sid(b, mn, bs, t)
        fo.write("%s\t%g\t%d\t%d\t%d\t%s\n" % (s, b, mn, bs, t, "current" if s == CURRENT else ""))
with open(os.path.join(OUT, "cells.tsv"), "w") as fo:
    fo.write("barcode\tmolecules\tcalled\n")
    for bc, n_inf, _, _ in out:
        fo.write("%s\t%d\t%d\n" % (bc, n_inf, int(n_inf >= RUN["cell_markers"])))
with gzip.open(os.path.join(OUT, "calls.tsv.gz"), "wt") as fo:
    fo.write("setting\tbarcode\tchrom\tstart\tend\tpre\tpost\tcstart\tcend\n")
    for _, _, calls, _ in out:
        for r in calls:
            fo.write("\t".join(map(str, r)) + "\n")
with gzip.open(os.path.join(OUT, "blocks.tsv.gz"), "wt") as fo:
    fo.write("setting\tbarcode\tchrom\tstart\tend\tgenotype\n")
    for _, _, _, blocks in out:
        for r in blocks:
            fo.write("\t".join(map(str, r)) + "\n")

# ---- summary per setting (over the cells hapCO calls)
src = {c: c for c, _ in CMAP}
pf = "results/reference/%s/pieces.tsv" % S
if os.path.exists(pf):
    src.update({r["piece"]: r["source"] for r in csv.DictReader(open(pf), delimiter="\t")})
L = dict(CMAP)
hap = {}
ht = "qc/haplotypes/%s/haplotype_tracks.tsv.gz" % S
if os.path.exists(ht):
    with gzip.open(ht, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                hap[r["barcode"]] = float(r["haploidness"])
            except (ValueError, KeyError, TypeError):
                pass
called = [bc for bc, n_inf, _, _ in out if n_inf >= RUN["cell_markers"]]
clean = {bc for bc in called if hap.get(bc, 0) >= 0.9}
n = max(len(called), 1)
per_set = {s: [] for s in IDS}
for _, _, calls, _ in out:
    for r in calls:
        per_set[r[0]].append(r)
units = sorted(set(src.values()), key=lambda c: [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", c)])
with open(os.path.join(OUT, "summary.tsv"), "w") as fo:
    fo.write("setting\tcells\tcos_per_cell\tmin_unit\tunits_below_0.5\tclose_double_share\tclean_cos_per_cell\t"
             "other_cos_per_cell\t" + "\t".join(units) + "\n")
    for s in IDS:
        cs = per_set[s]
        per_unit = {u: 0 for u in units}
        for r in cs:
            per_unit[src.get(r[2], r[2])] += 1
        pu = {u: v / n for u, v in per_unit.items()}
        by = {}
        for r in cs:
            by.setdefault((r[1], r[2]), []).append((r[3] + r[4]) / 2)
        pairs = close = 0
        for (bc, c), mids in by.items():
            mids.sort()
            for a, b2 in zip(mids, mids[1:]):
                pairs += 1
                close += (b2 - a) < 0.1 * L[c]
        cc = sum(1 for r in cs if r[1] in clean) / max(len(clean), 1)
        oc = sum(1 for r in cs if r[1] not in clean) / max(n - len(clean), 1)
        fo.write("%s\t%d\t%.3f\t%.3f\t%d\t%.3f\t%.3f\t%.3f\t%s\n" % (
            s, len(called), len(cs) / n, min(pu.values()) if pu else 0, sum(v < 0.5 for v in pu.values()),
            close / pairs if pairs else 0, cc, oc, "\t".join("%.3f" % pu[u] for u in units)))
print("SWEEP  %s: %d cells called (of %d selected) x %d settings -> %s in %.0f s" % (
    S, len(called), len(CELLS), len(GRID), OUT, time.time() - T0))
print("  crossovers per cell: current setting %s %.2f; range over the grid %.2f to %.2f" % (
    CURRENT, len(per_set[CURRENT]) / n, min(len(v) for v in per_set.values()) / n,
    max(len(v) for v in per_set.values()) / n))
