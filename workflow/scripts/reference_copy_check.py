#!/usr/bin/env python3
"""reference_copy_check.py -- is every piece of the genome in the mapping reference exactly once?

The crossover chain maps all reads to one haploid reference. Pieces (100 kb,
one every 250 kb) of every chromosome of BOTH haplotypes of the dual assembly
are mapped to it (minimap2 asm20, approximate; secondary hits scoring >= 80% of
the best kept). A piece's copies = its distinct hits covering >= half of it.
Per 1 Mb window, the median over its pieces:
  0 copies   the reference lacks this sequence: its reads have nowhere to go
  1 copy     as required
  2+ copies  both homologs are in the reference: reads split between them and
             every site looks homozygous
Either failure leaves the region without heterozygous markers. Runs >= 5 Mb
are listed. FAIL when more than --max_bad_mb Mb are at 0 or 2+ copies; exit 1
on FAIL unless --report_only. Writes OUTDIR/copy_check.txt, a per-window
table, and OUTDIR/copy_check.pass on PASS.

Usage: reference_copy_check.py DUAL_FASTA REFERENCE_FASTA OUTDIR [--threads 24] [--max_bad_mb 20]
                               [--report_only] [--label NAME]
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

ap = argparse.ArgumentParser()
ap.add_argument("dual"); ap.add_argument("reference"); ap.add_argument("outdir")
ap.add_argument("--threads", type=int, default=24)
ap.add_argument("--max_bad_mb", type=float, default=20)
ap.add_argument("--report_only", action="store_true")
ap.add_argument("--label", default="")
A = ap.parse_args()
PIECE, STEP, W, RUN = 100000, 250000, 1000000, 5
os.makedirs(A.outdir, exist_ok=True)
for f in ("copy_check.pass",):
    if os.path.exists(os.path.join(A.outdir, f)):
        os.remove(os.path.join(A.outdir, f))
mm = sorted(glob.glob(".snakemake/conda/*/bin/minimap2"))
MM = os.path.abspath(mm[0]) if mm else shutil.which("minimap2")
if not MM:
    sys.exit("minimap2 not found")


def is_chrom(n):
    return re.match(r"^chr\d+_hap[12]$", n) is not None


ref_names = [l[1:].split()[0] for l in open(A.reference) if l.startswith(">")]
pcs = os.path.join(A.outdir, "pieces.fa")
lens = {}
with open(pcs, "w") as fp:
    name, parts = None, []

    def done(n, seq):
        lens[n] = len(seq)
        for s0 in range(0, len(seq) - PIECE + 1, STEP):
            p = seq[s0:s0 + PIECE]
            if p.count("N") + p.count("n") <= PIECE // 2:
                fp.write(">%s__%d\n%s\n" % (n, s0, p))

    with open(A.dual) as fh:
        for l in fh:
            if l.startswith(">"):
                if name and is_chrom(name):
                    done(name, "".join(parts))
                name, parts = l[1:].split()[0], []
            elif name and is_chrom(name):
                parts.append(l.strip())
    if name and is_chrom(name):
        done(name, "".join(parts))
print("%s: %d dual-assembly chromosomes, reference %s (%s)" % (A.label or A.reference, len(lens), A.reference,
                                                               ", ".join(ref_names)), flush=True)
paf = os.path.join(A.outdir, "pieces_vs_reference.paf")
with open(paf, "w") as fo, open(os.path.join(A.outdir, "minimap2.log"), "w") as fl:
    r = subprocess.run([MM, "-x", "asm20", "-t", str(A.threads), "--secondary=yes", "-N", "5", "-p", "0.8",
                        "-I", "16G", A.reference, pcs], stdout=fo, stderr=fl)
if r.returncode:
    sys.exit("minimap2 failed -- see %s/minimap2.log" % A.outdir)

hits = collections.defaultdict(set)
pieces = []
for l in open(pcs):
    if l.startswith(">"):
        q, s0 = l[1:].strip().rsplit("__", 1)
        pieces.append((q, int(s0)))
for l in open(paf):
    f = l.split("\t")
    if len(f) < 12 or int(f[3]) - int(f[2]) < PIECE // 2:
        continue
    q, s0 = f[0].rsplit("__", 1)
    hits[(q, int(s0))].add((f[5], int(f[7]) // PIECE))
rows = []
for q, s0 in pieces:
    rows.append((q, s0 // W, len(hits.get((q, s0), ()))))
t = pd.DataFrame(rows, columns=["chrom", "window", "copies"])
wt = t.groupby(["chrom", "window"]).copies.median().round().astype(int).reset_index()
wt["start_mb"], wt["end_mb"] = wt.window, wt.window + 1
wt["class"] = np.where(wt.copies == 0, "0 copies", np.where(wt.copies == 1, "1 copy", "2+ copies"))
wt.to_csv(os.path.join(A.outdir, "copy_check_windows.tsv"), sep="\t", index=False)

bad_mb = int((wt["class"] != "1 copy").sum())
ok = bad_mb <= A.max_bad_mb
lines = ["REFERENCE COPY CHECK  %s  (%s)" % (A.label or os.path.basename(os.path.dirname(A.reference)), A.reference),
         "  reference chromosomes: %s" % ", ".join(ref_names),
         "  every 1 Mb window of both haplotypes, copies present in the reference (median of its 100 kb pieces):",
         "    %-12s %10s %10s %10s" % ("chromosome", "0 copies", "1 copy", "2+ copies")]
for c in sorted(wt.chrom.unique(), key=lambda c: (c.endswith("hap2"), int(re.search(r"\d+", c).group()))):
    x = wt[wt.chrom == c]["class"].value_counts()
    lines.append("    %-12s %8d Mb %8d Mb %8d Mb" % (c, x.get("0 copies", 0), x.get("1 copy", 0), x.get("2+ copies", 0)))
lines.append("  runs >= %d Mb not at exactly 1 copy:" % RUN)
nrun = 0
for c, x in wt.sort_values(["chrom", "window"]).groupby("chrom"):
    x = x.reset_index(drop=True)
    i = 0
    while i < len(x):
        if x["class"][i] != "1 copy":
            j = i
            while j + 1 < len(x) and x["class"][j + 1] == x["class"][i] and x.window[j + 1] == x.window[j] + 1:
                j += 1
            if j - i + 1 >= RUN:
                lines.append("    %-12s %6d-%6d Mb  %s" % (c, x.window[i], x.window[j] + 1, x["class"][i]))
                nrun += 1
            i = j + 1
        else:
            i += 1
if nrun == 0:
    lines.append("    none")
lines.append("  %s: %d Mb of the genome not at exactly 1 copy (limit %.0f Mb)" % ("PASS" if ok else "FAIL", bad_mb,
                                                                                 A.max_bad_mb))
open(os.path.join(A.outdir, "copy_check.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
if ok:
    open(os.path.join(A.outdir, "copy_check.pass"), "w").write("PASS %d Mb\n" % bad_mb)
os.remove(pcs)
sys.exit(0 if ok or A.report_only else 1)
