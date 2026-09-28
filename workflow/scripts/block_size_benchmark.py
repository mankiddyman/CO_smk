#!/usr/bin/env python3
"""block_size_benchmark.py -- what does hapCO's interior block-length minimum cost, and what does it buy?

hapCO drops any INTERIOR block shorter than block_size, so two crossovers
closer together than that -- an island -- cannot be called. With molecule input
(one row per molecule, informative markers only), marker_num already sets a
minimum that scales with marker density: 4 molecules span a few hundred kb in
gene-rich ends and several Mb in deserts. This measures on real cells whether
the fixed length rule is still needed on top.

1. PLANT: in N real cells, flip the alleles of every molecule inside islands of
   known length (0.25-3 Mb by default), each placed inside an existing called
   block, >= 1 Mb from that block's edges and from the chromosome's ends. Real
   noise and real coverage, known truth.
2. CALL: hapCO on the planted AND the untouched inputs at each block_size, every
   other setting exactly as the pipeline ran (read from co_calling_params.txt).
   The untouched run at 2 Mb must reproduce the pipeline's calls: checked.
3. SCORE:
   recall of planted islands, by length and by the molecules inside them
   untouched cells: islands that appear only below 2 Mb, by orientation and
   whether they sit at a recurrent structural hotspot. Real close doubles go
   both ways (REF-in-ALT ~ ALT-in-REF) and are scattered; leftover artefacts
   are one-sided or pile up at the hotspots.

Usage: block_size_benchmark.py SAMPLE [--cells 100] [--seed 1] [--threads 16]
         [--sizes 2000000,1000000,500000,250000,1] [--lengths 0.25,0.5,1,1.5,2,3]
         [--hotspots CHROM:START-END,...]
Output: qc/benchmark/SAMPLE/{block_size_summary.txt, planted.tsv, recall.tsv, new_islands.tsv}
"""
import argparse
import collections
import itertools
import os
import random
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("--cells", type=int, default=100)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--threads", type=int, default=16)
ap.add_argument("--sizes", default="2000000,1000000,500000,250000,1")
ap.add_argument("--lengths", default="0.25,0.5,1,1.5,2,3")
ap.add_argument("--per_cell", type=int, default=4)
ap.add_argument("--hotspots", default="")
ap.add_argument("--rscript", default=".snakemake/conda/738e17a23d1043909c1a53d6c039531c_/bin/Rscript")
ap.add_argument("--hapco", default="workflow/scripts/hapCO_identification.R")
A = ap.parse_args()
S = A.sample
RS, HAP = os.path.abspath(A.rscript), A.hapco
MOL = "results/cell_data_mol/%s" % S
PIPE = "results/crossovers/%s/per_cell" % S
CMAP = "results/cell_data/%s/chrom_map.tsv" % S
OUT = "qc/benchmark/%s" % S
SIZES = [int(x) for x in A.sizes.split(",")]
LENS = [float(x) for x in A.lengths.split(",")]
REF_BS = 2000000
HOT = []
for h in filter(None, A.hotspots.split(",")):
    c, r = h.split(":")
    HOT.append((c, int(r.split("-")[0]), int(r.split("-")[1])))
FLIP = {"0": "1", "1": "0"}
os.makedirs(OUT, exist_ok=True)

# ---- the pipeline's own settings
P = {}
for l in open("results/crossovers/%s/co_calling_params.txt" % S):
    if ":" in l:
        k, v = l.split(":", 1)
        if k.strip():
            P[k.strip().split()[0]] = v.strip()
need = ["cell_markers", "marker_num", "baseAF", "windowAF", "genotype", "terminal_marker_num"]
if [k for k in need if k not in P]:
    sys.exit("co_calling_params.txt lacks %s" % [k for k in need if k not in P])
if "cell_data_mol" not in P.get("input", ""):
    sys.exit("the pipeline run did not use molecule input (%s) -- nothing to benchmark" % P.get("input"))


def read_tab(p):
    return pd.read_csv(p, sep="\t", header=None, names=["chrom", "pos", "ref", "rc", "alt", "ac"])


def read_blocks(p):
    rows = []
    if os.path.exists(p):
        for l in open(p):
            f = l.split()
            if len(f) >= 4 and f[1].isdigit() and f[2].isdigit():
                rows.append((f[0], int(f[1]), int(f[2]), f[3]))
    return rows


def n_cos(p):
    return sum(1 for l in open(p) if len(l.split()) >= 3 and l.split()[1].isdigit()) if os.path.exists(p) else None


def islands(blocks):
    by = collections.defaultdict(list)
    for r in blocks:
        by[r[0]].append(r)
    out = []
    for bl in by.values():
        bl.sort(key=lambda r: r[1])
        for i in range(1, len(bl) - 1):
            if bl[i - 1][3] == bl[i + 1][3] != bl[i][3]:
                out.append(bl[i])
    return out


good_cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
rng = random.Random(A.seed)
cells = rng.sample(good_cells, min(A.cells, len(good_cells)))

# ---- 1. plant
pdir = os.path.join(OUT, "inputs", "planted")
shutil.rmtree(pdir, ignore_errors=True)
os.makedirs(pdir)
planted, mols = [], {}
cyc = itertools.cycle(LENS)
for bc in cells:
    d = read_tab(os.path.join(MOL, bc + ".tsv"))
    dp = (d.rc + d.ac).clip(lower=1)
    fr = d.ac / dp
    mols[bc] = d.assign(call=np.where(fr <= 0.2, 0, np.where(fr >= 0.8, 1, -1)))[["chrom", "pos", "call"]]
    by = collections.defaultdict(list)
    for r in read_blocks(os.path.join(PIPE, bc + "_co_block_pred.txt")):
        by[r[0]].append(r)
    free = list(by)
    rng.shuffle(free)
    for _ in range(A.per_cell):
        L = int(next(cyc) * 1e6)
        for ch in list(free):
            bl = sorted(by[ch], key=lambda r: r[1])
            c0, c1 = bl[0][1], bl[-1][2]
            cands = []
            for _, s0, e0, gt in bl:
                lo = max(s0 + 1000000, c0 + 1000000)
                hi = min(e0 - 1000000, c1 - 1000000) - L
                if hi > lo:
                    cands.append((lo, hi, gt))
            if not cands:
                continue
            lo, hi, gt = rng.choices(cands, weights=[h - l for l, h, _ in cands])[0]
            a = rng.randint(lo, hi)
            b = a + L
            sel = ((d.chrom == ch) & (d.pos >= a) & (d.pos < b)).values
            d.loc[sel, ["rc", "ac"]] = d.loc[sel, ["ac", "rc"]].values
            planted.append((bc, ch, a, b, L, int(sel.sum()), gt))
            free.remove(ch)
            break
    d.to_csv(os.path.join(pdir, bc + ".tsv"), sep="\t", header=False, index=False)
pt = pd.DataFrame(planted, columns=["barcode", "chrom", "start", "end", "length", "molecules", "block_gt"])
pt.to_csv(os.path.join(OUT, "planted.tsv"), sep="\t", index=False)
print("planted %d islands in %d cells (%s)" % (len(pt), len(cells),
      ", ".join("%.2g Mb x%d" % (k / 1e6, v) for k, v in sorted(pt.length.value_counts().items()))), flush=True)


# ---- 2. call
def run(job):
    inp, bc, out, bs = job
    cmd = [RS, HAP, "--input", inp, "--prefix", bc, "--chrom_map", CMAP, "--outpath", out,
           "--cell_markers", P["cell_markers"], "--block_size", str(bs), "--marker_num", P["marker_num"],
           "--baseAF", P["baseAF"], "--windowAF", P["windowAF"], "--genotype", P["genotype"],
           "--terminal_marker_num", P["terminal_marker_num"]]
    r = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, universal_newlines=True)
    return bc, out, r.returncode, r.stderr[-400:]


jobs = []
for bs in SIZES:
    for kind, indir in (("planted", pdir), ("base", MOL)):
        o = os.path.join(OUT, "runs", "%s_bs%d" % (kind, bs))
        shutil.rmtree(o, ignore_errors=True)
        os.makedirs(o)
        jobs += [(os.path.join(indir, bc + ".tsv"), bc, o, bs) for bc in cells]
print("calling: %d hapCO runs on %d threads ..." % (len(jobs), A.threads), flush=True)
res = []
with ThreadPoolExecutor(A.threads) as ex:
    for i, r in enumerate(ex.map(run, jobs), 1):
        res.append(r)
        if i % 100 == 0 or i == len(jobs):
            print("  %d / %d" % (i, len(jobs)), flush=True)
bad = [r for r in res if r[2] != 0]
print("  done; %d failed" % len(bad))
for r in bad[:3]:
    print("  FAILED %s in %s:\n%s" % (r[0], r[1], r[3]))

lines = ["BLOCK SIZE BENCHMARK  %s  (%d cells, seed %d; marker_num %s molecules, terminal %s, baseAF %s, "
         "windowAF %s, genotype %s)" % (S, len(cells), A.seed, P["marker_num"], P["terminal_marker_num"],
                                         P["baseAF"], P["windowAF"], P["genotype"])]
if REF_BS in SIZES:
    same = sum(1 for bc in cells if n_cos(os.path.join(PIPE, bc + "_co_pred.txt")) ==
               n_cos(os.path.join(OUT, "runs", "base_bs%d" % REF_BS, bc + "_co_pred.txt")))
    lines.append("  check: untouched cells at 2 Mb reproduce the pipeline's CO count in %d of %d cells" % (same, len(cells)))

# which hapCO genotype code carries ALT? learned from the pipeline's big blocks
fits = collections.defaultdict(list)
for bc in cells:
    m = mols[bc][mols[bc].call >= 0]
    for ch, s0, e0, gt in read_blocks(os.path.join(PIPE, bc + "_co_block_pred.txt")):
        x = m.call[(m.chrom == ch) & (m.pos >= s0) & (m.pos <= e0)]
        if len(x) >= 10:
            fits[gt].append(x.mean())
alt_code = max(fits, key=lambda k: np.mean(fits[k]))

# ---- 3a. recall of planted islands
rec = []
for bs in SIZES:
    cache = {}
    for bc, ch, a, b, L, nm, gt in planted:
        if bc not in cache:
            cache[bc] = read_blocks(os.path.join(OUT, "runs", "planted_bs%d" % bs, bc + "_co_block_pred.txt"))
        hit, err = False, np.nan
        for c2, s0, e0, g2 in cache[bc]:
            if c2 == ch and g2 == FLIP.get(gt, gt) and min(e0, b) - max(s0, a) >= 0.5 * L and e0 - s0 <= 2 * L + 1e6:
                hit, err = True, max(abs(s0 - a), abs(e0 - b))
                break
        rec.append((bs, bc, ch, L, nm, hit, err))
rt = pd.DataFrame(rec, columns=["block_size", "barcode", "chrom", "length", "molecules", "recovered", "boundary_err"])
rt.to_csv(os.path.join(OUT, "recall.tsv"), sep="\t", index=False)
rt["mol_bin"] = pd.cut(rt.molecules, [-1, 3, 7, 15, 31, 1e9], labels=["0-3", "4-7", "8-15", "16-31", "32+"])

lens = sorted(rt.length.unique())
lines.append("  RECALL of planted islands, by island length (planted per length: %s)"
             % ", ".join("%d" % (pt.length == L).sum() for L in lens))
lines.append("    %-11s" % "block_size" + "".join("%9s" % ("%.2g Mb" % (L / 1e6)) for L in lens)
             + "   boundary err (median)")
for bs in SIZES:
    x = rt[rt.block_size == bs]
    lines.append("    %-11s" % ("%.2g Mb" % (bs / 1e6) if bs > 1 else "off")
                 + "".join("%8.0f%%" % (100 * x[x.length == L].recovered.mean()) for L in lens)
                 + "   %6.0f kb" % (x[x.recovered].boundary_err.median() / 1e3 if x.recovered.any() else np.nan))
bins = ["0-3", "4-7", "8-15", "16-31", "32+"]
lines.append("  RECALL by informative molecules inside the island (islands per bin: %s)"
             % ", ".join("%s %d" % (k, (rt[rt.block_size == SIZES[0]].mol_bin == k).sum()) for k in bins))
for bs in SIZES:
    x = rt[rt.block_size == bs]
    lines.append("    %-11s" % ("%.2g Mb" % (bs / 1e6) if bs > 1 else "off")
                 + "".join("%8.0f%%" % (100 * x[x.mol_bin == k].recovered.mean()) if (x.mol_bin == k).any()
                           else "%9s" % "-" for k in bins))

# ---- 3b. untouched cells: islands that only appear below 2 Mb
ref = {bc: islands(read_blocks(os.path.join(OUT, "runs", "base_bs%d" % REF_BS, bc + "_co_block_pred.txt")))
       for bc in cells}
new_rows, per = [], []
for bs in SIZES:
    tot_cos = tot_isl = 0
    for bc in cells:
        blk = read_blocks(os.path.join(OUT, "runs", "base_bs%d" % bs, bc + "_co_block_pred.txt"))
        tot_cos += n_cos(os.path.join(OUT, "runs", "base_bs%d" % bs, bc + "_co_pred.txt")) or 0
        isl = islands(blk)
        tot_isl += len(isl)
        if bs == REF_BS:
            continue
        m = mols[bc][mols[bc].call >= 0]
        for ch, s0, e0, gt in isl:
            if any(c == ch and min(e, e0) - max(s, s0) >= 0.5 * (e0 - s0) for c, s, e, _ in ref[bc]):
                continue
            x = m.call[(m.chrom == ch) & (m.pos >= s0) & (m.pos <= e0)]
            allele = "ALT" if gt == alt_code else "REF"
            hot = any(c == ch and min(e, e0) > max(s, s0) for c, s, e in HOT)
            new_rows.append((bs, bc, ch, s0, e0, e0 - s0, len(x), allele, hot))
    per.append((bs, tot_cos / float(len(cells)), tot_isl / float(len(cells))))
nt = pd.DataFrame(new_rows, columns=["block_size", "barcode", "chrom", "start", "end", "length", "molecules",
                                     "allele", "hotspot"])
nt.to_csv(os.path.join(OUT, "new_islands.tsv"), sep="\t", index=False)
lines.append("  UNTOUCHED cells: all calls, and islands that appear only below 2 Mb (per cell)")
lines.append("    %-11s %8s %9s %9s %11s %11s %9s %13s" % ("block_size", "COs", "islands", "new", "REF-in-ALT",
                                                            "ALT-in-REF", "hotspot", "median length"))
N = float(len(cells))
for bs, cos, isl in per:
    x = nt[nt.block_size == bs]
    lines.append("    %-11s %8.2f %9.2f %9.2f %11.2f %11.2f %9.2f %10s"
                 % ("%.2g Mb" % (bs / 1e6) if bs > 1 else "off", cos, isl, len(x) / N,
                    (x.allele == "REF").sum() / N, (x.allele == "ALT").sum() / N, x.hotspot.sum() / N,
                    ("%.0f kb" % (x.length.median() / 1e3)) if len(x) else "-"))
lines += ["  READ: pick the smallest block_size whose new islands are (i) balanced between REF-in-ALT and",
          "  ALT-in-REF and (ii) not piled on the hotspots, while short planted islands are recovered."]
open(os.path.join(OUT, "block_size_summary.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
