#!/usr/bin/env python3
"""co_masks.py -- per cell and chromosome (or piece): there once, twice, not at all, or mixed? And masked?

A nucleus that carries a sequence twice reads both alleles interleaved; hapCO calls noise there, not
crossovers. The classes are calling_overview.py's, from the per-cell tables hapCO reads (molecules with
both alleles equally often are left out; a molecule is ALT when ALT reads outnumber REF reads):
  expected   the cell's molecules x the median share of this sequence across the selected cells
  switch     allele changes between neighbouring molecules / (molecules - 1): ~0 for one haplotype (a
             crossover adds one), ~0.5 when two copies are interleaved
  class      too few (expected < 8) | missing (< 0.25 x expected) | doubled (switch >= 0.25 and
             >= 1.2 x expected) | mixed (switch >= 0.25 at normal coverage) | haploid
A pair is masked when its class is in --classes (config co_calling.<species>.mask_classes); co_aggregate.R
then drops its calls and the cell does not count for that sequence.

Usage: co_masks.py SAMPLE GOOD_CELLS CHROM_MAP TABLES_DIR OUT [--classes doubled,missing,mixed]
Writes OUT (barcode, chrom, molecules, total, switch, expected, coverage, class, masked) and prints a summary;
if qc/calling_overview/SAMPLE_dosage.tsv exists, says how many classes agree with it.
"""
import argparse
import csv
import os
import statistics

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("good_cells")
ap.add_argument("chrom_map")
ap.add_argument("tables")
ap.add_argument("out")
ap.add_argument("--classes", default="")
A = ap.parse_args()
CLASSES = set(x for x in A.classes.replace(" ", "").split(",") if x)
CHROMS = [r["name"] for r in csv.DictReader(open(A.chrom_map), delimiter="\t")]
CELLS = [l.strip() for l in open(A.good_cells) if l.strip()]

rows, missing = [], 0
for bc in CELLS:
    p = os.path.join(A.tables, bc + ".tsv")
    if not os.path.exists(p):
        missing += 1
        continue
    by = {c: [] for c in CHROMS}
    tot = 0
    with open(p) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[0] not in by:
                continue
            rc, ac = int(float(f[3])), int(float(f[5]))
            if rc != ac:
                by[f[0]].append((int(float(f[1])), 1 if ac > rc else 0))
                tot += 1
    for c in CHROMS:
        a = [g for _, g in sorted(by[c])]
        n = len(a)
        sw = sum(1 for x, y in zip(a, a[1:]) if x != y) / (n - 1) if n > 1 else float("nan")
        rows.append([bc, c, n, tot, sw])
share = {c: statistics.median([r[2] / max(r[3], 1) for r in rows if r[1] == c]) for c in CHROMS} if rows else {}
out, counts = [], {}
for bc, c, n, tot, sw in rows:
    exp = tot * share[c]
    cov = n / exp if exp else float("nan")
    if exp < 8:
        cl = "too few"
    elif cov < 0.25:
        cl = "missing"
    elif sw >= 0.25:
        cl = "doubled" if cov >= 1.2 else "mixed"
    else:
        cl = "haploid"
    m = "yes" if cl in CLASSES else "no"
    out.append((bc, c, n, tot, sw, exp, cov, cl, m))
    counts[(c, cl)] = counts.get((c, cl), 0) + 1
os.makedirs(os.path.dirname(A.out) or ".", exist_ok=True)
with open(A.out, "w") as fo:
    fo.write("barcode\tchrom\tmolecules\ttotal\tswitch\texpected\tcoverage\tclass\tmasked\n")
    for bc, c, n, tot, sw, exp, cov, cl, m in out:
        fo.write("%s\t%s\t%d\t%d\t%s\t%.2f\t%s\t%s\t%s\n" % (bc, c, n, tot, "NA" if sw != sw else "%.4f" % sw, exp,
                                                         "NA" if cov != cov else "%.4f" % cov, cl, m))
nm = sum(1 for r in out if r[8] == "yes")
print("co_masks %s: %d cells (%d without a table), %d sequences; masked classes: %s -> %d cell x sequence masked"
      % (A.sample, len(CELLS) - missing, missing, len(CHROMS), ", ".join(sorted(CLASSES)) or "none", nm))
order = ("haploid", "doubled", "missing", "mixed", "too few")
print("  %-18s %s" % ("sequence", "".join("%9s" % k for k in order)))
for c in CHROMS:
    print("  %-18s %s" % (c, "".join("%9d" % counts.get((c, k), 0) for k in order)))
ov = "qc/calling_overview/%s_dosage.tsv" % A.sample
if os.path.exists(ov):
    mine = {(r[0], r[1]): r[7] for r in out}
    same = tot = 0
    for r in csv.DictReader(open(ov), delimiter="\t"):
        k = (r["barcode"], r["chrom"])
        if k in mine:
            tot += 1
            same += mine[k] == r["class"]
    print("  classes identical to %s for %d of %d cell x sequence" % (ov, same, tot))
