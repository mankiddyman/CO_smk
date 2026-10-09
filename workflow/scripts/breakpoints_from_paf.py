#!/usr/bin/env python3
"""breakpoints_from_paf.py -- where each chromosome changes homolog, read from the hap2-on-hap1 alignment,
and where to cut the crossover mapping reference so that no translocation join is left in it.

From first principles; nothing typed in (no expected positions, no tables from other scripts):
 1. ALIGNMENTS. Every hap2 chromosome was aligned to all of hap1 (rule hap_align: minimap2 asm20,
    results/hap_align/SAMPLE/hap2_on_hap1.paf). One PAF line says: this stretch of a hap2 chromosome is
    the same sequence as that stretch of a hap1 chromosome. Used: primary lines with MAPQ >= --min_mapq,
    >= --min_kb long on both sides.
 2. HOMOLOGY = COLLINEAR CHAINS. Two homologous stretches align as a long run of alignments in the same
    order on both sides; repeats align as scattered hits. Alignments between one hap2 and one hap1
    chromosome on one strand are chained while each next one is <= --chain_gap_mb further on BOTH sides;
    only chains with >= --min_chain_mb aligned count as homology.
 3. PAINT. Every chromosome of both haplotypes, --bin_mb at a time, takes the chromosome of the other
    haplotype whose chains cover most of that stretch ('none' when they cover < --min_cover of it).
    A chromosome untouched by a translocation is painted by its homolog end to end; a translocated one by
    two homologs (one per arm), or by one and then 'none' (an arm whose other copy is not in the other
    haplotype).
 4. CLEAN-UP. 'none' between two runs of the same homolog is an alignment gap, not a breakpoint. Runs of
    another homolog shorter than --min_run_mb are islands (small duplications): listed, not breakpoints.
    'none' shorter than --min_arm_mb at a chromosome end is an unaligned tip.
 5. BREAKPOINT = where the homolog changes. Exactly: between the last aligned base of the left homolog's
    chains and the first aligned base of the right one's; the cut goes in the middle of the stretch
    between them (nothing aligns there), or right after the left homolog's last base when the right arm
    has no homolog in the other haplotype.
 6. PLAN. The mapping reference (BASE_FAI) is cut at every breakpoint on its chromosomes and nowhere else.
    Pieces are named with --names CHROM=LEFT,RIGHT[,...]. A breakpoint on a reference chromosome without
    names, names that do not match the breakpoints found, or a reference chromosome whose length differs
    from the alignment's, stops the plan: look before cutting.
Breakpoints on chromosomes outside the reference are reported too: they are the same joins seen in the
other haplotype, and what the Hi-C map shows as well. --names may also label the arms of chromosomes
outside the reference (labels and colours only). A name is a claim: an arm and its homolog in the other
haplotype get the same name (passed on to unnamed homologs), and a name given to two different pairs of
homologs, or to the wrong number of arms, stops the plan. Arms with no homolog in the other haplotype
(here P and Q) cannot be matched by this alignment; giving two of them one name is not checked.
Figures colour every arm by its name, the same wherever it sits; hatched = no homolog in the other haplotype.

Usage: breakpoints_from_paf.py PAF BASE_FAI OUTDIR PLAN_CSV [--names CHROM=A,B ...] [--label NAME]
         [--min_mapq 5] [--min_kb 20] [--chain_gap_mb 5] [--min_chain_mb 1] [--bin_mb 1]
         [--min_cover 0.1] [--min_run_mb 5] [--min_arm_mb 10] [--zoom_mb 15]
Writes OUTDIR/{breakpoints.txt, breakpoints.csv, homolog_paint.png, breakpoints_zoom.png} and, unless the
plan is rejected, PLAN_CSV (one '#' line, then piece, source, start, end, length, how; 1-based, inclusive;
together exactly the reference's sequence). A plan identical to the one already there is left untouched.
"""
import argparse
import collections
import datetime
import io
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("paf"); ap.add_argument("base_fai"); ap.add_argument("outdir"); ap.add_argument("plan")
ap.add_argument("--names", action="append", default=[])
ap.add_argument("--label", default="")
ap.add_argument("--min_mapq", type=int, default=5)
ap.add_argument("--min_kb", type=float, default=20)
ap.add_argument("--chain_gap_mb", type=float, default=5)
ap.add_argument("--min_chain_mb", type=float, default=1)
ap.add_argument("--bin_mb", type=float, default=1)
ap.add_argument("--min_cover", type=float, default=0.1)
ap.add_argument("--min_run_mb", type=float, default=5)
ap.add_argument("--min_arm_mb", type=float, default=10)
ap.add_argument("--zoom_mb", type=float, default=15)
A = ap.parse_args()
MB = 1000000
BIN, GAP, MINLEN = int(A.bin_mb * MB), int(A.chain_gap_mb * MB), int(A.min_kb * 1000)
MIN_RUN, MIN_ARM = max(1, int(round(A.min_run_mb * MB / BIN))), max(1, int(round(A.min_arm_mb * MB / BIN)))
GAP_BINS = max(1, int(round(GAP / BIN)))
os.makedirs(A.outdir, exist_ok=True)
REPORT = []


def say(s=""):
    REPORT.append(s)
    print(s, flush=True)


def mb(x):
    return "%.3f" % (x / MB)


# ---------------------------------------------------------------- 1. alignments
L, side = {}, {}                 # sequence -> length; sequence -> 'hap1' (PAF target) or 'hap2' (PAF query)
alns, n_all, n_sec, n_weak = [], 0, 0, 0
with open(A.paf) as fh:
    for line in fh:
        f = line.split("\t", 12)
        if len(f) < 12:
            continue
        n_all += 1
        L[f[0]], L[f[5]] = int(f[1]), int(f[6])
        side[f[0]], side[f[5]] = "hap2", "hap1"
        if len(f) > 12 and "tp:A:S" in f[12]:
            n_sec += 1
            continue
        qs, qe, ts, te, mq = int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[11])
        if mq < A.min_mapq or qe - qs < MINLEN or te - ts < MINLEN:
            n_weak += 1
            continue
        alns.append((f[0], qs, qe, f[4], f[5], ts, te, mq))
if not alns:
    sys.exit("no usable alignments in %s" % A.paf)

# ---------------------------------------------------------------- 2. collinear chains
groups = collections.defaultdict(list)
for a in alns:
    groups[(a[0], a[4], a[3])].append(a)
chains = []
for lst in groups.values():
    lst.sort(key=lambda a: (a[1], a[2]))
    open_ = []
    for a in lst:
        q, qs, qe, st, t, ts, te, mq = a
        best, best_d = None, None
        for ch in open_:
            p = ch[-1]
            dq = qs - p[2]
            dt = (ts - p[6]) if st == "+" else (p[5] - te)
            if -GAP <= dq <= GAP and -GAP <= dt <= GAP:
                d = abs(dq) + abs(dt)
                if best is None or d < best_d:
                    best, best_d = ch, d
        if best is None:
            open_.append([a])
        else:
            best.append(a)
        keep = []
        for ch in open_:
            (chains if qs - ch[-1][2] > GAP else keep).append(ch)
        open_ = keep
    chains.extend(open_)
homol = [ch for ch in chains if sum(a[6] - a[5] for a in ch) >= A.min_chain_mb * MB]
used = [a for ch in homol for a in ch]
say("WHERE EACH CHROMOSOME CHANGES HOMOLOG  %s  (%s)" % (A.label, datetime.datetime.now().isoformat(timespec="seconds")))
say("  alignment: %s (each hap2 chromosome aligned to all of hap1)" % A.paf)
say("  PAF lines %s: secondary %s, weak (MAPQ < %d or < %.0f kb) %s, confident %s (%.0f Mb on hap1)"
    % (format(n_all, ","), format(n_sec, ","), A.min_mapq, A.min_kb, format(n_weak, ","), format(len(alns), ","),
       sum(a[6] - a[5] for a in alns) / MB))
say("  homology = collinear chains with >= %.1f Mb aligned (next alignment <= %.0f Mb further on both sides):"
    % (A.min_chain_mb, A.chain_gap_mb))
say("    %d chains, %s alignments, %.0f Mb on hap1; left out as repeats or short duplications: %s alignments, %.0f Mb"
    % (len(homol), format(len(used), ","), sum(a[6] - a[5] for a in used) / MB, format(len(alns) - len(used), ","),
       (sum(a[6] - a[5] for a in alns) - sum(a[6] - a[5] for a in used)) / MB))

# per sequence: its chain alignments in its own coordinates
HOM = collections.defaultdict(list)   # seq -> [(start, end, partner, p_start, p_end, strand, mapq)]
for q, qs, qe, st, t, ts, te, mq in used:
    HOM[t].append((ts, te, q, qs, qe, st, mq))
    HOM[q].append((qs, qe, t, ts, te, st, mq))


# ---------------------------------------------------------------- 3. paint
def order(c):
    m = re.match(r"^chr(\d+)_hap([12])$", c)
    return (0, int(m.group(1)), int(m.group(2))) if m else (1, 0, c)


ROWS = sorted([c for c in L if re.match(r"^chr\d+_hap[12]$", c) or L[c] >= 10 * MB], key=order)


def nbins(c):
    return (L[c] + BIN - 1) // BIN


def paint(c):
    n = nbins(c)
    cov = collections.defaultdict(lambda: np.zeros(n))
    for s, e, p, _, _, _, _ in HOM.get(c, []):
        for b in range(s // BIN, (e - 1) // BIN + 1):
            cov[p][b] += min(e, (b + 1) * BIN) - max(s, b * BIN)
    size = np.minimum(BIN, L[c] - np.arange(n) * BIN).astype(float)
    if not cov:
        return [None] * n, {}
    parts = list(cov)
    M = np.vstack([cov[p] for p in parts])
    best, val = M.argmax(0), M.max(0)
    return [parts[best[b]] if val[b] >= A.min_cover * size[b] else None for b in range(n)], cov


def rle(lab):
    runs = []
    for b, x in enumerate(lab):
        if runs and runs[-1][0] == x:
            runs[-1][2] = b
        else:
            runs.append([x, b, b])
    return runs


def clean(lab):
    """Steps 3-4: fill alignment gaps, drop islands, absorb unaligned tips. Returns labels, islands, gaps."""
    lab, islands, gaps = list(lab), [], []
    while True:
        runs = rle(lab)
        changed = False
        for i in range(1, len(runs) - 1):          # 'none' between two runs of one homolog: a gap
            x, b0, b1 = runs[i]
            if x is None and runs[i - 1][0] is not None and runs[i - 1][0] == runs[i + 1][0]:
                if b1 - b0 + 1 >= 2:
                    gaps.append((runs[i - 1][0], b0, b1))
                lab[b0:b1 + 1] = [runs[i - 1][0]] * (b1 - b0 + 1)
                changed = True
        if changed:
            continue
        homologs = [r for r in runs if r[0] is not None]
        for x, b0, b1 in homologs:                 # short runs of a homolog: islands
            if len(homologs) > 1 and b1 - b0 + 1 < MIN_RUN:
                islands.append((x, b0, b1))
                lab[b0:b1 + 1] = [None] * (b1 - b0 + 1)
                changed = True
        if changed:
            continue
        for i in (0, len(runs) - 1):                # short unaligned chromosome tips
            x, b0, b1 = runs[i]
            if x is None and len(runs) > 1 and b1 - b0 + 1 < MIN_ARM:
                nb = runs[1][0] if i == 0 else runs[-2][0]
                lab[b0:b1 + 1] = [nb] * (b1 - b0 + 1)
                changed = True
        if not changed:
            return lab, islands, gaps


def arms(lab):
    """Runs that are arms: homolog runs, and 'none' runs at an end or >= MIN_ARM long. Shorter interior
    'none' runs (between two different homologs) are the unaligned junction of a breakpoint."""
    runs, out = rle(lab), []
    for i, (x, b0, b1) in enumerate(runs):
        if x is None and 0 < i < len(runs) - 1 and b1 - b0 + 1 < MIN_ARM:
            continue
        out.append((x, b0, b1))
    return out


PAINT, COV, ISL, GAPS, ARMS = {}, {}, {}, {}, {}
for c in ROWS:
    raw, COV[c] = paint(c)
    PAINT[c], ISL[c], GAPS[c] = clean(raw)
    ARMS[c] = arms(PAINT[c])


# ---------------------------------------------------------------- 5. breakpoints
def partner_pos(h, at_end):
    s, e, p, ps, pe, st, mq = h
    return (pe if st == "+" else ps) if at_end else (ps if st == "+" else pe)


BPS = []
for c in ROWS:
    ar = ARMS[c]
    for k_arm, ((x, a0, a1), (y, c0, c1)) in enumerate(zip(ar, ar[1:])):
        hits = HOM.get(c, [])
        left = [h for h in hits if x is not None and h[2] == x and a0 * BIN <= h[0] < (c0 + 1) * BIN]
        right = [h for h in hits if y is not None and h[2] == y and a1 * BIN < h[1] <= (c1 + 1) * BIN]
        lh = max(left, key=lambda h: h[1]) if left else None
        rh = min(right, key=lambda h: h[0]) if right else None
        lend, rstart = (lh[1] if lh else None), (rh[0] if rh else None)
        if lend is not None and rstart is not None:
            cut = (lend + rstart) // 2
        else:
            cut = lend if lend is not None else rstart
        if cut is None:                              # no chain hit near the boundary: fall back to the bins
            cut = c0 * BIN
        BPS.append({"chrom": c, "k": k_arm, "left": x, "right": y, "lend": lend, "rstart": rstart, "cut": cut,
                    "lh": lh, "rh": rh, "left_hits": sorted(left, key=lambda h: h[1])[-3:],
                    "right_hits": sorted(right, key=lambda h: h[0])[:3]})


def homolog_name(x, c):
    return x if x is not None else "none (no copy in %s)" % ("hap2" if side.get(c) == "hap1" else "hap1")


base = []
for l in open(A.base_fai):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 2:
        base.append((f[0], int(f[1])))
NAMES, bad = {}, []
for spec in A.names:
    if "=" not in spec:
        sys.exit("--names %r: want CHROM=LEFT,RIGHT[,...]" % spec)
    c, v = spec.split("=", 1)
    NAMES[c] = [x for x in v.split(",") if x]
inref = {c for c, _ in base}

# arm names: given ones, passed on to the homolog in the other haplotype, else chrN when both share a number
ARMNAME = {}
for c, nms in NAMES.items():
    if c not in ARMS:
        bad.append("--names %s: %s is not in the alignment" % (c, c))
    elif len(nms) != len(ARMS[c]):
        bad.append("--names %s=%s: %d names, but %s has %d arm(s) here" % (c, ",".join(nms), len(nms), c, len(ARMS[c])))
    else:
        for k, nm in enumerate(nms):
            ARMNAME[(c, k)] = nm


def homolog_arms(c, k):
    p = ARMS[c][k][0]
    return [(p, j) for j, r in enumerate(ARMS.get(p, [])) if r[0] == c] if p is not None else []


for (c, k), nm in sorted(ARMNAME.items()):
    for pj in homolog_arms(c, k):
        if pj in ARMNAME and ARMNAME[pj] != nm:
            bad.append("%s arm %s pairs with %s arm %s: one arm, two names" % (c, nm, pj[0], ARMNAME[pj]))
        ARMNAME.setdefault(pj, nm)
for c in ROWS:
    for k, (p, b0, b1) in enumerate(ARMS[c]):
        if (c, k) not in ARMNAME:
            mc, mp = re.match(r"^(chr\d+)_", c), re.match(r"^(chr\d+)_", p or "")
            ARMNAME[(c, k)] = ("" if p is None else mc.group(1) if (mc and mp and mc.group(1) == mp.group(1))
                               else "%s~%s" % (c, p))
pairs_of = collections.defaultdict(set)
for (c, k), nm in ARMNAME.items():
    if nm:
        p = ARMS[c][k][0]
        pairs_of[nm].add(frozenset((c, p)) if p is not None else ("no homolog",))
for nm, ps in sorted(pairs_of.items()):
    paired = [x for x in ps if x != ("no homolog",)]
    if len(paired) > 1 or (paired and ("no homolog",) in ps):
        bad.append("arm name %s is given to different stretches: %s" % (nm, "; ".join(
            " with ".join(sorted(x)) if x != ("no homolog",) else "one with no homolog" for x in ps)))


def edges(c):
    return [0] + sorted(b["cut"] for b in BPS if b["chrom"] == c) + [L[c]]


say("")
say("PAINT  (each run = one homolog in the other haplotype, [arm]; Mb)")
for c in ROWS:
    say("  %-11s %7.1f Mb  %s" % (c, L[c] / MB, "  |  ".join(
        "%s-%s %s%s" % ("%.1f" % (b0 * BIN / MB), "%.1f" % (min((b1 + 1) * BIN, L[c]) / MB), x or "none",
                        " [%s]" % ARMNAME[(c, k)] if ARMNAME.get((c, k)) else "")
        for k, (x, b0, b1) in enumerate(ARMS[c]))))
say("")
say("ARMS  (where each arm sits, between the cuts; an arm and its homolog share a name)")
seen_arm = collections.OrderedDict()
for c in ROWS:
    ed = edges(c)
    for k, (p, b0, b1) in enumerate(ARMS[c]):
        nm = ARMNAME.get((c, k)) or "(unnamed)"
        seen_arm.setdefault(nm, []).append("%s %.1f-%.1f%s" % (c, ed[k] / MB, ed[k + 1] / MB,
                                                               "" if p is not None else " (no homolog in %s)" %
                                                               ("hap2" if side.get(c) == "hap1" else "hap1")))
for nm, where in seen_arm.items():
    say("  %-10s %s" % (nm, ";  ".join(where)))
say("  (two copies with no homolog in the other haplotype, like P on two hap1 chromosomes, are not matched to each "
    "other by this alignment)")

# ---------------------------------------------------------------- 6. plan
for c, n in base:
    if c not in L:
        bad.append("%s (reference) is not in the alignment" % c)
    elif L[c] != n:
        bad.append("%s is %d bp in the reference but %d bp in the alignment: not the same assembly" % (c, n, L[c]))

say("")
say("BREAKPOINTS  (cut = middle of the unaligned stretch between the two homologs; 1 Mb / 2 Mb bin as on a Hi-C map)")
bp_of = collections.defaultdict(list)
for b in BPS:
    c = b["chrom"]
    bp_of[c].append(b)
    where = ("in the reference: CUT HERE" if c in inref else "not in the reference: nothing to cut")
    gap = (b["rstart"] - b["lend"]) if (b["lend"] is not None and b["rstart"] is not None) else None
    say("  %-11s %s Mb   (1 Mb bin %d, 2 Mb bin %d)   %s" % (c, mb(b["cut"]), b["cut"] // MB, b["cut"] // (2 * MB), where))
    say("     left : %s, last aligned base %s%s" % (
        homolog_name(b["left"], c), mb(b["lend"]) if b["lend"] is not None else "-",
        "  (= %s:%s)" % (b["left"], mb(partner_pos(b["lh"], True))) if b["lh"] else ""))
    say("     right: %s, first aligned base %s%s" % (
        homolog_name(b["right"], c), mb(b["rstart"]) if b["rstart"] is not None else "-",
        "  (= %s:%s)" % (b["right"], mb(partner_pos(b["rh"], False))) if b["rh"] else ""))
    if gap is not None:
        say("     between them: %s of %.0f kb" % ("unaligned stretch" if gap >= 0 else "overlap", abs(gap) / 1000.0))
    for tag, hs in (("last left", b["left_hits"]), ("first right", b["right_hits"])):
        for h in hs:
            say("       %-11s %s-%s  %-10s %s-%s %s  MAPQ %d" % (tag, mb(h[0]), mb(h[1]), h[2], mb(h[3]), mb(h[4]), h[5], h[6]))
if not BPS:
    say("  none: every chromosome is painted by one homolog")

say("")
say("ISLANDS  (runs < %.0f Mb painted by a chromosome other than their surroundings; not breakpoints)" % A.min_run_mb)
n_isl = 0
for c in ROWS:
    for x, b0, b1 in ISL[c]:
        n_isl += 1
        say("  %-11s %7.1f-%7.1f Mb  %s" % (c, b0 * BIN / MB, min((b1 + 1) * BIN, L[c]) / MB, x))
if not n_isl:
    say("  none")
say("GAPS  (>= 3 Mb with no homology chain inside one homolog's run: satellites, rDNA, assembly gaps; not breakpoints)")
n_gap, small = 0, []
for c in ROWS:
    for x, b0, b1 in GAPS[c]:
        if (b1 - b0 + 1) * BIN >= 3 * MB:
            n_gap += 1
            say("  %-11s %7.1f-%7.1f Mb  inside %s" % (c, b0 * BIN / MB, min((b1 + 1) * BIN, L[c]) / MB, x))
        else:
            small.append((b1 - b0 + 1) * BIN)
if not n_gap:
    say("  none")
if small:
    say("  plus %d shorter gaps (2-3 Mb), %.0f Mb in all" % (len(small), sum(small) / MB))

rows = []
for c, n in base:
    bps = sorted(bp_of.get(c, []), key=lambda b: b["cut"])
    names = NAMES.get(c)
    if bps and not names:
        bad.append("%s (reference) changes homolog at %s Mb but has no --names: look at it before cutting"
                   % (c, ", ".join(mb(b["cut"]) for b in bps)))
        continue
    if names and len(names) != len(bps) + 1:
        bad.append("%s: --names gives %d pieces (%s), the alignment shows %d breakpoint(s)%s"
                   % (c, len(names), ",".join(names), len(bps),
                      (" at " + ", ".join(mb(b["cut"]) + " Mb" for b in bps)) if bps else ""))
        continue
    if not bps:
        rows.append((c, c, 1, n, n, "whole: one homolog end to end"))
        continue
    starts = [1] + [b["cut"] + 1 for b in bps]
    ends = [b["cut"] for b in bps] + [n]
    for k, (s, e) in enumerate(zip(starts, ends)):
        b_left = bps[k - 1] if k > 0 else None
        b_right = bps[k] if k < len(bps) else None
        homolog = b_right["left"] if b_right else b_left["right"]
        how = "homolog %s" % homolog_name(homolog, c)
        if b_right:
            how += "; cut %s Mb" % mb(b_right["cut"])
        rows.append(("%s_%s" % (c, names[k]), c, s, e, e - s + 1, how))
say("")
say("PLAN for %s  (%d sequences)" % (A.base_fai, len(base)))
for r in rows:
    say("  %-16s %-11s %12d %12d  %7.1f Mb  %s" % (r[0], r[1], r[2], r[3], r[4] / MB, r[5]))
if rows:
    say("  %d pieces, %.1f Mb (reference %.1f Mb)" % (len(rows), sum(r[4] for r in rows) / MB, sum(n for _, n in base) / MB))
say("")
say("SUMMARY")
for b in sorted(BPS, key=lambda b: order(b["chrom"])):
    say("  %-4s %-11s %8s Mb   %-5s | %-5s  (homolog %s | %s)" % (
        "CUT" if b["chrom"] in inref else "join", b["chrom"], mb(b["cut"]), ARMNAME.get((b["chrom"], b["k"]), ""),
        ARMNAME.get((b["chrom"], b["k"] + 1), ""), b["left"] or "none", b["right"] or "none"))
if bad:
    say("")
    say("PLAN REJECTED:")
    for b in bad:
        say("  " + b)


# ---------------------------------------------------------------- figures
NUMCOL = {1: "#5DCAA5", 2: "#F0997B", 3: "#6C9BD2", 4: "#D97FB0", 5: "#C9A227", 6: "#8E8E3C"}
NONE_COL = "#E3E0D8"


KNOWN = {"L1": "#5DCAA5", "L2": "#F0997B", "P": "#AFA9EC", "Q": "#FAC775"}   # the arm colours of the Hi-C figures
EXTRA = ["#E4572E", "#D97FB0", "#2A9D8F", "#8D6E63", "#7E57C2", "#5B8E7D", "#B07AA1"]
ARMCOL = {}


def pcol(p):
    m = re.match(r"^chr(\d+)_", p or "")
    return NUMCOL.get(int(m.group(1)), "#999999") if m else ("#999999" if p else NONE_COL)


def acol(nm):
    """One colour per arm name, the same wherever the arm sits."""
    if not nm:
        return NONE_COL
    if nm not in ARMCOL:
        m = re.match(r"^chr(\d+)$", nm)
        if nm in KNOWN:
            ARMCOL[nm] = KNOWN[nm]
        elif m and NUMCOL.get(int(m.group(1))) and NUMCOL[int(m.group(1))] not in ARMCOL.values():
            ARMCOL[nm] = NUMCOL[int(m.group(1))]
        else:
            ARMCOL[nm] = next((x for x in EXTRA if x not in ARMCOL.values()), "#999999")
    return ARMCOL[nm]


fig, ax = plt.subplots(figsize=(15, 0.62 * len(ROWS) + 1.8))
for i, c in enumerate(ROWS):
    y = len(ROWS) - 1 - i
    ed = edges(c)
    for k, (x, b0, b1) in enumerate(ARMS[c]):
        s, e, nm = ed[k] / MB, ed[k + 1] / MB, ARMNAME.get((c, k), "")
        ax.add_patch(plt.Rectangle((s, y - 0.3), e - s, 0.6, fc=acol(nm), lw=0,
                                   ec="#8A857A" if x is None else acol(nm), hatch="///" if x is None else None))
        if nm:
            ax.text((s + e) / 2, y, nm, ha="center", va="center", fontsize=8, fontweight="bold")
    ax.add_patch(plt.Rectangle((0, y - 0.3), L[c] / MB, 0.6, fill=False, ec="#444444", lw=0.6))
    for x, b0, b1 in ISL[c]:
        ax.plot([b0 * BIN / MB, (b1 + 1) * BIN / MB], [y + 0.38] * 2, color=pcol(x), lw=3, solid_capstyle="butt")
    for b in bp_of.get(c, []):
        ref = c in inref
        ax.plot([b["cut"] / MB] * 2, [y - 0.42, y + 0.42], color="black" if ref else "#555555",
                lw=1.6 if ref else 1.0, ls="-" if ref else "--")
        ax.text(b["cut"] / MB + 1.5, y + 0.33, "%s%.2f" % ("CUT " if ref else "", b["cut"] / MB), ha="left", va="bottom",
                fontsize=7.5, fontweight="bold" if ref else "normal")
    ax.text(-4, y, c + ("  (ref)" if c in inref else ""), ha="right", va="center", fontsize=8.5,
            fontweight="bold" if c in inref else "normal")
names_used = [nm for nm in ARMCOL]
handles = [plt.Rectangle((0, 0), 1, 1, fc=ARMCOL[nm], ec=ARMCOL[nm]) for nm in names_used] + \
          [plt.Rectangle((0, 0), 1, 1, fc="white", hatch="///", ec="#8A857A")]
ax.legend(handles, names_used + ["hatched: no homolog in the other haplotype"],
          ncol=min(len(handles), 11), loc="upper center", bbox_to_anchor=(0.5, -0.09), frameon=False, fontsize=8)
ax.set_xlim(-2, max(L[c] for c in ROWS) / MB * 1.02); ax.set_ylim(-0.8, len(ROWS) - 0.2)
ax.set_yticks([]); ax.set_xlabel("position (Mb)")
for s in ("left", "right", "top"):
    ax.spines[s].set_visible(False)
ax.set_title("Each chromosome split into arms where its homolog in the other haplotype changes (hap2-on-hap1 alignment, "
             "collinear chains >= %.0f Mb)\none colour per arm wherever it sits; black = cut in the mapping reference; "
             "dashed = the same joins on chromosomes outside it; strips above a bar = islands < %.0f Mb"
             % (A.min_chain_mb, A.min_run_mb), fontsize=10, loc="left")
fig.tight_layout()
fig.savefig(os.path.join(A.outdir, "homolog_paint.png"), dpi=130)
plt.close(fig)

if BPS:
    nc = min(3, len(BPS))
    nr = (len(BPS) + nc - 1) // nc
    fig, axs = plt.subplots(nr, nc, figsize=(5.6 * nc, 3.6 * nr), squeeze=False)
    usedset = set(used)
    other = collections.defaultdict(list)            # confident alignments outside homology chains (repeats)
    for a in alns:
        if a not in usedset:
            other[a[4]].append((a[5], a[6]))
            other[a[0]].append((a[1], a[2]))
    Z = A.zoom_mb * MB
    for k, b in enumerate(BPS):
        ax = axs[k // nc][k % nc]
        c, lo, hi = b["chrom"], b["cut"] - Z, b["cut"] + Z
        hits = [h for h in HOM.get(c, []) if h[1] > lo and h[0] < hi]
        partners = []
        for h in sorted(hits, key=lambda h: h[0]):
            if h[2] not in partners:
                partners.append(h[2])
        nb = len(partners)
        for j, p in enumerate(partners):             # band j from the top: the homolog's own coordinates
            y0 = nb - j
            vis = []                                  # the part inside the window, partner side interpolated
            for s_, e_, _, ps, pe, st, _ in (h for h in hits if h[2] == p):
                ya, yb = (ps, pe) if st == "+" else (pe, ps)
                xa, xb = max(s_, lo), min(e_, hi)
                fa, fb = (xa - s_) / float(e_ - s_), (xb - s_) / float(e_ - s_)
                vis.append((xa, xb, ya + fa * (yb - ya), ya + fb * (yb - ya)))
            pmin, pmax = min(min(v[2], v[3]) for v in vis), max(max(v[2], v[3]) for v in vis)
            span = max(pmax - pmin, 1)
            arm_nm = (ARMNAME.get((c, b["k"]), "") if p == b["left"] else
                      ARMNAME.get((c, b["k"] + 1), "") if p == b["right"] else "")
            for xa, xb, ya, yb in vis:
                ax.plot([xa / MB, xb / MB], [y0 - 0.85 + 0.7 * (ya - pmin) / span, y0 - 0.85 + 0.7 * (yb - pmin) / span],
                        color=acol(arm_nm) if arm_nm else pcol(p), lw=2.4, solid_capstyle="butt")
            ax.text(lo / MB + 0.2, y0 - 0.08, "%s%s  %.1f-%.1f Mb" % ("%s on " % arm_nm if arm_nm else "", p,
                                                                     pmin / MB, pmax / MB),
                    fontsize=7.5, va="top", color="#333333")
            ax.axhline(y0 - 1, color="#DDDDDD", lw=0.6)
        for s_, e_ in other.get(c, []):              # repeats: x only
            if e_ > lo and s_ < hi:
                ax.plot([s_ / MB, e_ / MB], [-0.45, -0.45], color="#9A9A9A", lw=5, solid_capstyle="butt")
        ax.text(lo / MB + 0.2, -0.12, "other confident alignments (repeats), position only", fontsize=7, va="top",
                color="#777777")
        for side_, x_ in (("left", b["left"]), ("right", b["right"])):
            if x_ is None:
                ax.text((lo + Z * 0.5) / MB if side_ == "left" else (hi - Z * 0.5) / MB, nb / 2.0 + 0.1,
                        "no homolog in %s" % ("hap2" if side.get(c) == "hap1" else "hap1"), ha="center", fontsize=8.5,
                        color="#777777", style="italic")
        if b["lend"] is not None and b["rstart"] is not None:
            ax.axvspan(min(b["lend"], b["rstart"]) / MB, max(b["lend"], b["rstart"]) / MB, color="#FFD966", alpha=0.6, lw=0)
        ref = c in inref
        ax.axvline(b["cut"] / MB, color="black" if ref else "#555555", lw=1.4, ls="-" if ref else "--")
        ax.set_xlim(lo / MB, hi / MB); ax.set_ylim(-0.75, max(nb, 1) + 0.05)
        ax.set_yticks([]); ax.tick_params(labelsize=8)
        ax.set_xlabel("%s (Mb)" % c, fontsize=9)
        ax.set_title("%s  %s %s Mb%s\n%s | %s   (homolog %s | %s)" % (
            c, "CUT at" if ref else "join at", mb(b["cut"]), "" if ref else "  (not in the reference)",
            ARMNAME.get((c, b["k"]), "") or "?", ARMNAME.get((c, b["k"] + 1), "") or "?",
            b["left"] or "none", b["right"] or "none"), fontsize=9.5, loc="left")
    for k in range(len(BPS), nr * nc):
        axs[k // nc][k % nc].axis("off")
    fig.suptitle("Every breakpoint up close (+-%.0f Mb). Each band = one homolog, its alignments drawn in its own "
                 "coordinates (range at the left); line = cut; yellow = unaligned stretch between the two homologs"
                 % A.zoom_mb, fontsize=10, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(A.outdir, "breakpoints_zoom.png"), dpi=120)
    plt.close(fig)

with open(os.path.join(A.outdir, "breakpoints.csv"), "w") as fo:
    fo.write("chrom,haplotype,cut,left_homolog,left_last_base,right_homolog,right_first_base,unaligned_kb,in_reference\n")
    for b in BPS:
        fo.write("%s,%s,%d,%s,%s,%s,%s,%s,%s\n" % (
            b["chrom"], side.get(b["chrom"], "?"), b["cut"], b["left"] or "none",
            b["lend"] if b["lend"] is not None else "", b["right"] or "none",
            b["rstart"] if b["rstart"] is not None else "",
            "%.0f" % ((b["rstart"] - b["lend"]) / 1000.0) if (b["lend"] is not None and b["rstart"] is not None) else "",
            "yes" if b["chrom"] in inref else "no"))
open(os.path.join(A.outdir, "breakpoints.txt"), "w").write("\n".join(REPORT) + "\n")
if bad:
    sys.exit(1)
buf = io.StringIO()
buf.write("# split plan %s, made by workflow/scripts/breakpoints_from_paf.py from %s and %s; report: %s\n"
          % (A.label, A.paf, A.base_fai, os.path.join(A.outdir, "breakpoints.txt")))
buf.write("piece,source,start,end,length,how\n")
for r in rows:
    buf.write("%s,%s,%d,%d,%d,%s\n" % (r[0], r[1], r[2], r[3], r[4], r[5].replace(",", ";")))
if os.path.exists(A.plan) and open(A.plan).read() == buf.getvalue():
    print("plan unchanged: %s left as it is" % A.plan)
else:
    os.makedirs(os.path.dirname(os.path.abspath(A.plan)), exist_ok=True)
    open(A.plan, "w").write(buf.getvalue())
    print("wrote %s" % A.plan)
