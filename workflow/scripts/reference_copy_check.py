#!/usr/bin/env python3
"""reference_copy_check.py -- is every piece of the genome in the mapping reference exactly once?

The crossover chain maps all reads to one haploid reference, chosen from the
chromosomes of a dual (hap1 + hap2) assembly. Pieces (100 kb, one every
250 kb) of every chromosome of BOTH haplotypes are mapped against the WHOLE
dual assembly (minimap2 asm20, approximate, secondary hits kept); hits on the
piece's own chromosome are ignored, and a hit must cover >= half the piece.
For each piece, its copies in the reference = (1 if its own chromosome is in
the reference) + the number of OTHER reference chromosomes it hits:
  ok          exactly 1
  duplicated  2 or more: both homologs are in the reference, reads split
              between them and every site looks homozygous
  absent      0, but the piece does hit a non-reference chromosome: this
              sequence exists in the assembly, twice, and not in the reference
              -- its reads have nowhere to go
  unmatched   0, and no counterpart anywhere: divergent or hemizygous
              sequence. The same for every choice of reference, so it is
              reported but never fails the check
Each 1 Mb window takes the class of at least half its pieces (else ok). The
failure this guards against is large and contiguous (a whole arm), so the
check counts windows in runs of >= --run duplicated or absent windows (gaps of
one window allowed); FAIL when those runs exceed --max_bad_mb. Exit 1 on FAIL
unless --report_only. Writes OUTDIR/copy_check.txt, the window table and, on
PASS, OUTDIR/copy_check.pass. The piece-vs-assembly mapping depends only on
the dual assembly: --paf reuses one (pieces named NAME__START).

Usage: reference_copy_check.py DUAL_FASTA REFERENCE_FASTA OUTDIR [--threads 24] [--max_bad_mb 20]
                               [--report_only] [--label NAME] [--paf PIECES_VS_DUAL.paf]
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
ap.add_argument("--paf")
ap.add_argument("--piece", type=int, default=100000)
ap.add_argument("--step", type=int, default=250000)
ap.add_argument("--window", type=int, default=1000000)
ap.add_argument("--run", type=int, default=5)
A = ap.parse_args()
PIECE, STEP, W, RUN = A.piece, A.step, A.window, A.run
MB = W / 1e6
os.makedirs(A.outdir, exist_ok=True)
if os.path.exists(os.path.join(A.outdir, "copy_check.pass")):
    os.remove(os.path.join(A.outdir, "copy_check.pass"))


def is_chrom(n):
    return re.match(r"^chr\d+_hap[12]$", n) is not None


def cnum(c):
    return int(re.search(r"\d+", c).group())


ref = set(l[1:].split()[0] for l in open(A.reference) if l.startswith(">"))
lens = {}
if A.paf:
    paf, name = A.paf, None
    for l in open(A.dual):
        if l.startswith(">"):
            name = l[1:].split()[0]
            lens[name] = 0
        elif name is not None:
            lens[name] += len(l.strip())
    lens = {k: v for k, v in lens.items() if is_chrom(k)}
else:
    mm = sorted(glob.glob(".snakemake/conda/*/bin/minimap2"))
    MM = os.path.abspath(mm[0]) if mm else shutil.which("minimap2")
    if not MM:
        sys.exit("minimap2 not found")
    tgt, pcs = os.path.join(A.outdir, "dual_chromosomes.fa"), os.path.join(A.outdir, "pieces.fa")
    with open(tgt, "w") as ft, open(pcs, "w") as fp:
        name, keep, parts = None, False, []

        def done(n, seq):
            lens[n] = len(seq)
            for s0 in range(0, len(seq) - PIECE + 1, STEP):
                p = seq[s0:s0 + PIECE]
                if p.count("N") + p.count("n") <= PIECE // 2:
                    fp.write(">%s__%d\n%s\n" % (n, s0, p))

        with open(A.dual) as fh:
            for l in fh:
                if l.startswith(">"):
                    if keep:
                        done(name, "".join(parts))
                    name = l[1:].split()[0]
                    keep, parts = is_chrom(name), []
                    if keep:
                        ft.write(">%s\n" % name)
                elif keep:
                    ft.write(l)
                    parts.append(l.strip())
        if keep:
            done(name, "".join(parts))
    paf = os.path.join(A.outdir, "pieces_vs_dual.paf")
    print("mapping pieces of %d chromosomes against the whole dual assembly ..." % len(lens), flush=True)
    with open(paf, "w") as fo, open(os.path.join(A.outdir, "minimap2.log"), "w") as fl:
        r = subprocess.run([MM, "-x", "asm20", "-t", str(A.threads), "--secondary=yes", "-N", "20", "-p", "0.3",
                            "-I", "16G", tgt, pcs], stdout=fo, stderr=fl)
    os.remove(tgt)
    if r.returncode:
        sys.exit("minimap2 failed -- see %s/minimap2.log" % A.outdir)

missing = sorted(ref - set(lens))
if missing:
    sys.exit("reference chromosomes not in the dual assembly: %s" % missing)

# every piece, and the OTHER chromosomes it hits
hits = collections.defaultdict(set)
for l in open(paf):
    f = l.split("\t")
    if len(f) < 12:
        continue
    q, s0 = f[0].rsplit("__", 1)
    if f[5] == q or f[5] not in lens or int(f[3]) - int(f[2]) < PIECE // 2:
        continue
    hits[(q, int(s0))].add(f[5])
pieces = [(q, s0) for q, n in lens.items() for s0 in range(0, n - PIECE + 1, STEP)]
if not A.paf:
    os.remove(os.path.join(A.outdir, "pieces.fa"))

rows = []
for q, s0 in pieces:
    h = hits.get((q, s0), set())
    copies = (1 if q in ref else 0) + sum(1 for c in h if c in ref)
    if copies == 1:
        cls = "ok"
    elif copies >= 2:
        cls = "duplicated"
    elif any(c not in ref for c in h):
        cls = "absent"
    else:
        cls = "unmatched"
    rows.append((q, s0 // W, cls))
t = pd.DataFrame(rows, columns=["chrom", "window", "cls"])
frac = t.groupby(["chrom", "window"]).cls.value_counts(normalize=True).unstack(fill_value=0)
for c in ("ok", "duplicated", "absent", "unmatched"):
    if c not in frac.columns:
        frac[c] = 0.0
wclass = np.where(frac["duplicated"] >= 0.5, "duplicated",
                  np.where(frac["absent"] >= 0.5, "absent", np.where(frac["unmatched"] >= 0.5, "unmatched", "ok")))
wt = frac.reset_index()[["chrom", "window"]].assign(**{"class": wclass})
wt.to_csv(os.path.join(A.outdir, "copy_check_windows.tsv"), sep="\t", index=False)

# runs of duplicated / absent windows (one-window gaps allowed)
runs = []
for c, x in wt.sort_values(["chrom", "window"]).groupby("chrom"):
    cls = dict(zip(x.window, x["class"]))
    for kind in ("duplicated", "absent"):
        ws = sorted(w for w, k in cls.items() if k == kind)
        i = 0
        while i < len(ws):
            j = i
            while j + 1 < len(ws) and ws[j + 1] - ws[j] <= 2:
                j += 1
            n = sum(1 for w in range(ws[i], ws[j] + 1) if cls.get(w) == kind)
            if n >= RUN:
                runs.append((c, ws[i], ws[j] + 1, kind, n))
            i = j + 1
bad_run_mb = sum(r[4] for r in runs) * MB
ok = bad_run_mb <= A.max_bad_mb
lines = ["REFERENCE COPY CHECK  %s  (%s)" % (A.label or os.path.basename(os.path.dirname(A.reference)), A.reference),
         "  reference chromosomes: %s" % ", ".join(sorted(ref, key=lambda c: (cnum(c), c))),
         "  every %g Mb window of both haplotypes, by how often its sequence is in the reference:" % MB,
         "    %-12s %9s %12s %9s %11s" % ("chromosome", "ok", "duplicated", "absent", "unmatched")]
for c in sorted(lens, key=lambda c: (c.endswith("hap2"), cnum(c))):
    x = wt[wt.chrom == c]["class"].value_counts()
    lines.append("    %-12s %6.0f Mb %9.0f Mb %6.0f Mb %8.0f Mb   %s"
                 % (c, x.get("ok", 0) * MB, x.get("duplicated", 0) * MB, x.get("absent", 0) * MB,
                    x.get("unmatched", 0) * MB, "in reference" if c in ref else ""))
lines.append("  runs of >= %d duplicated or absent windows (the failure this check exists for):" % RUN)
for c, a, b, kind, n in sorted(runs, key=lambda r: (r[0].endswith("hap2"), cnum(r[0]), r[1])):
    lines.append("    %-12s %6.1f-%6.1f Mb  %-10s (%d windows)" % (c, a * MB, b * MB, kind, n))
if not runs:
    lines.append("    none")
lines.append("  unmatched (no counterpart anywhere in the assembly; the same for any reference): %.0f Mb"
             % ((wt["class"] == "unmatched").sum() * MB))
lines.append("  %s: %.0f Mb in duplicated/absent runs (limit %.0f Mb)" % ("PASS" if ok else "FAIL", bad_run_mb,
                                                                          A.max_bad_mb))
open(os.path.join(A.outdir, "copy_check.txt"), "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
if ok:
    open(os.path.join(A.outdir, "copy_check.pass"), "w").write("PASS %.0f Mb\n" % bad_run_mb)
sys.exit(0 if ok or A.report_only else 1)
