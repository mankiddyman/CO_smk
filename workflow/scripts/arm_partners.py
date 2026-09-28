#!/usr/bin/env python3
"""arm_partners.py -- which sequence does an unaligned chromosome tail pair with?

The genome-wide hap2-on-hap1 map calls a region 'no homolog' when nothing
aligns there with MAPQ >= 5. That hides two different things: sequence with
no counterpart, and sequence with TWO equally good counterparts (MAPQ 0 by
definition -- e.g. a right arm shared by two chromosomes). Genome-wide repeat
masking of minimizers can also starve a repetitive arm of anchors.

This takes every unaligned terminal tail (>= 10 Mb) of both haplotypes and
aligns the tails directly to one another, pair by pair, against a small
target (so nothing is lost to multi-mapping or genome-wide masking):
  hap2 tail vs hap1 tail   which hap1 tail is each hap2 tail's homolog?
  hap1 tail vs hap1 tail   do two hap1 chromosomes share a tail?
  hap2 tail vs hap2 tail   the same within hap2
Reported per pair: share of the query covered, share of the target covered,
identity, strand, collinearity (Spearman of positions). Queries are cut into
5 Mb pieces for minimap2 and put back together.
Also the raw HiFi calls in each hap1 tail against its chromosome's body: het
SNPs per Mb, depth, ALT ratio -- no calls means reads that could not be
placed; calls at ~2x depth or with skewed ratios mean collapsed copies.

--sample_mb N aligns an N Mb piece from the middle of each query tail instead
of the whole tail (minutes instead of hours; identity is what it measures).
The typical hap1-hap2 identity of ordinary homologs, from the genome-wide PAF,
is printed alongside: two tails that are the TWO HAPLOTYPES of one arm align
at that identity; a recent duplication aligns near 100%.
--bam HIFI.bam (reads on hap1): mean depth in each hap1 tail against its
chromosome's body. A tail holding one haplotype's copy of an arm gets only
that haplotype's reads: about half the body's depth.

Usage: arm_partners.py SAMPLE HAP1_FASTA OUTDIR [--threads 4] [--min_tail_mb 10] [--sample_mb 0] [--bam HIFI.bam]
Reads: qc/translocations/SAMPLE/SAMPLE_{synteny_blocks,hap1_windows}.tsv,
       results/hap_align/SAMPLE/hap2.fa, results/markers/SAMPLE/hifi_raw.vcf.gz
"""
import argparse
import glob
import itertools
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("sample"); ap.add_argument("hap1_fasta"); ap.add_argument("outdir")
ap.add_argument("--threads", type=int, default=4)
ap.add_argument("--min_tail_mb", type=float, default=10)
ap.add_argument("--hap2")
ap.add_argument("--sample_mb", type=float, default=0)
ap.add_argument("--bam")
A = ap.parse_args()
S, H1, OUT = A.sample, A.hap1_fasta, A.outdir
T = "qc/translocations/%s" % S
H2 = A.hap2 or "results/hap_align/%s/hap2.fa" % S
VCF = "results/markers/%s/hifi_raw.vcf.gz" % S
CHUNK = 5000000
os.makedirs(OUT, exist_ok=True)


def tool(name):
    c = sorted(glob.glob(".snakemake/conda/*/bin/%s" % name))
    return os.path.abspath(c[0]) if c else shutil.which(name)


MM, BCF, SAMT = tool("minimap2"), tool("bcftools"), tool("samtools")
if not MM:
    sys.exit("minimap2 not found (conda envs or PATH)")


def base(c):
    return re.sub(r"_hap[12]$", "", c)


def cnum(c):
    m = re.search(r"(\d+)", base(c))
    return int(m.group(1)) if m else 10 ** 6


def lengths(fa):
    out, name, n = {}, None, 0
    if os.path.exists(fa + ".fai"):
        for l in open(fa + ".fai"):
            f = l.split("\t")
            out[f[0]] = int(f[1])
        return out
    for l in open(fa):
        if l.startswith(">"):
            if name:
                out[name] = n
            name, n = l[1:].split()[0], 0
        else:
            n += len(l.strip())
    if name:
        out[name] = n
    return out


L1, L2 = lengths(H1), lengths(H2)
bt = pd.read_csv(os.path.join(T, "%s_synteny_blocks.tsv" % S), sep="\t")
wt = pd.read_csv(os.path.join(T, "%s_hap1_windows.tsv" % S), sep="\t")
MIN = A.min_tail_mb * 1e6

# ---- unaligned terminal tails
tails = []   # (hap, chrom, start, end, side)
for q in sorted([c for c in L2 if cnum(c) < 10 ** 6], key=cnum):
    x = bt[bt.hap2 == q]
    if len(x):
        lo, hi = x.q_start.min(), x.q_end.max()
        if lo >= MIN:
            tails.append(("hap2", q, 0, int(lo), "left"))
        if L2[q] - hi >= MIN:
            tails.append(("hap2", q, int(hi), L2[q], "right"))
for c in sorted([c for c in L1 if cnum(c) < 10 ** 6], key=cnum):
    x = wt[wt.chrom == c].sort_values("start")
    cov = x[x.aligned >= 0.2]
    if len(cov):
        lo, hi = cov.start.min(), cov.end.max()
        if lo >= MIN:
            tails.append(("hap1", c, 0, int(lo), "left"))
        if L1[c] - hi >= MIN:
            tails.append(("hap1", c, int(hi), L1[c], "right"))
if not tails:
    sys.exit("no unaligned tails >= %.0f Mb" % A.min_tail_mb)
lines = ["ARM PARTNERS  %s  (unaligned terminal tails >= %.0f Mb, aligned to each other directly; minimap2 asm20)"
         % (S, A.min_tail_mb), "  tails:"]
for h, c, s0, e0, side in tails:
    lines.append("    %s  %-11s %-5s %7.1f-%7.1f Mb  (%.1f Mb)" % (h, c, side, s0 / 1e6, e0 / 1e6, (e0 - s0) / 1e6))
print("\n".join(lines), flush=True)


# ---- extract tails (streaming, one pass per assembly)
def extract(fa, wanted, outdir):
    want = {}
    for h, c, s0, e0, side in wanted:
        want.setdefault(c, []).append((s0, e0, "%s_%s" % (c, side)))
    handles = {}
    for c, lst in want.items():
        for s0, e0, nm in lst:
            handles[nm] = open(os.path.join(outdir, nm + ".fa"), "w")
            handles[nm].write(">%s\n" % nm)
    cur, pos = None, 0
    with open(fa) as fh:
        for l in fh:
            if l.startswith(">"):
                cur, pos = l[1:].split()[0], 0
                continue
            if cur not in want:
                continue
            s = l.strip()
            n = len(s)
            for s0, e0, nm in want[cur]:
                a, b = max(s0, pos), min(e0, pos + n)
                if a < b:
                    handles[nm].write(s[a - pos:b - pos] + "\n")
            pos += n
    for h in handles.values():
        h.close()


seqdir = os.path.join(OUT, "tails")
os.makedirs(seqdir, exist_ok=True)
extract(H1, [t for t in tails if t[0] == "hap1"], seqdir)
extract(H2, [t for t in tails if t[0] == "hap2"], seqdir)


def chunk(name):
    src, dst = os.path.join(seqdir, name + ".fa"), os.path.join(seqdir, name + ".chunks.fa")
    seq = "".join(l.strip() for l in open(src) if not l.startswith(">"))
    if A.sample_mb > 0 and len(seq) > A.sample_mb * 1e6:
        half = int(A.sample_mb * 1e6 / 2)
        mid = len(seq) // 2
        seq = seq[mid - half:mid + half]
    with open(dst, "w") as f:
        for off in range(0, len(seq), CHUNK):
            f.write(">%s__%d\n%s\n" % (name, off, seq[off:off + CHUNK]))
    return dst, len(seq)


def align(pair):
    qn, tn, kind = pair
    qfa, qlen = chunk(qn)
    tfa = os.path.join(seqdir, tn + ".fa")
    paf = os.path.join(OUT, "%s__vs__%s.paf" % (qn, tn))
    r = subprocess.run([MM, "-x", "asm20", "-c", "-t", str(A.threads), tfa, qfa], stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, universal_newlines=True)
    rows = []
    for l in r.stdout.splitlines():
        f = l.split("\t")
        if len(f) < 12:
            continue
        nm, off = f[0].rsplit("__", 1)
        off = int(off)
        f[0], f[1], f[2], f[3] = nm, str(qlen), str(int(f[2]) + off), str(int(f[3]) + off)
        rows.append("\t".join(f))
    open(paf, "w").write("\n".join(rows) + ("\n" if rows else ""))
    return pair, paf, r.returncode, r.stderr[-300:]


def union(iv):
    iv.sort()
    tot, cs, ce = 0, None, None
    for s0, e0 in iv:
        if cs is None or s0 > ce:
            if cs is not None:
                tot += ce - cs
            cs, ce = s0, e0
        else:
            ce = max(ce, e0)
    return tot + (ce - cs if cs is not None else 0)


def summarize(paf):
    rows = []
    for l in open(paf):
        f = l.rstrip("\n").split("\t")
        if len(f) >= 12 and "tp:A:P" in l and int(f[3]) - int(f[2]) >= 5000:
            rows.append((int(f[1]), int(f[2]), int(f[3]), f[4], int(f[6]), int(f[7]), int(f[8]), int(f[9]),
                         int(f[10])))
    if not rows:
        return 0.0, 0.0, float("nan"), float("nan"), float("nan"), 0
    ql, tl = rows[0][0], rows[0][4]
    qcov = union([(r[1], r[2]) for r in rows]) / float(ql)
    tcov = union([(r[5], r[6]) for r in rows]) / float(tl)
    ident = sum(r[7] for r in rows) / float(sum(r[8] for r in rows))
    plus = sum(r[2] - r[1] for r in rows if r[3] == "+") / float(sum(r[2] - r[1] for r in rows))
    big = [r for r in rows if r[2] - r[1] >= 50000]
    rho = (pd.Series([(r[1] + r[2]) / 2 for r in big]).corr(pd.Series([(r[5] + r[6]) / 2 for r in big]),
                                                             method="spearman") if len(big) >= 5 else float("nan"))
    return qcov, tcov, ident, plus, rho, len(rows)


names = {t: "%s_%s" % (t[1], t[4]) for t in tails}
printed = len(lines)

# ---- the genome-wide PAF, with NO MAPQ filter: where did the hap2 tails' pieces go?
import collections
PAF = "results/hap_align/%s/hap2_on_hap1.paf" % S
agg = collections.defaultdict(float)
t2 = [(t[1], t[2], t[3], names[t]) for t in tails if t[0] == "hap2"]
if os.path.exists(PAF):
    with open(PAF) as fh:
        for l in fh:
            f = l.split("\t", 12)
            if len(f) < 12:
                continue
            for c, s0, e0, nm in t2:
                if f[0] == c and int(f[2]) >= s0 and int(f[3]) <= e0 + 1:
                    if int(f[8]) - int(f[7]) >= 20000:
                        mq = int(f[11])
                        agg[(nm, f[5], "MAPQ 0" if mq == 0 else "MAPQ 1-4" if mq < 5 else "MAPQ 5+",
                             "secondary" if "tp:A:S" in l else "primary")] += int(f[3]) - int(f[2])
                    break
    lines.append("  genome-wide PAF, NO MAPQ filter, alignments >= 20 kb: where the hap2 tails' pieces landed")
    for k in sorted(agg, key=lambda k: (k[0], -agg[k])):
        lines.append("    %-20s -> %-11s %-9s %-10s %7.1f Mb" % (k[0], k[1], k[2], k[3], agg[k] / 1e6))
    if not agg:
        lines.append("    nothing: no alignment of any quality from these tails")
    print("\n".join(lines[printed:]), flush=True)
    printed = len(lines)
h1 = [names[t] for t in tails if t[0] == "hap1"]
h2 = [names[t] for t in tails if t[0] == "hap2"]
pairs = [(q, t, "hap2 -> hap1") for q in h2 for t in h1]
pairs += [(a, b, "hap1 -> hap1") for a, b in itertools.combinations(h1, 2)]
pairs += [(a, b, "hap2 -> hap2") for a, b in itertools.combinations(h2, 2)]
print("aligning %d tail pairs, %d threads each ..." % (len(pairs), A.threads), flush=True)
with ThreadPoolExecutor(min(len(pairs), 6)) as ex:
    res = list(ex.map(align, pairs))

lines.append("  pairs (query -> target): query covered, target covered, identity, + strand share, collinearity rho")
rows = []
for (qn, tn, kind), paf, rc, err in res:
    if rc != 0:
        lines.append("    %-13s %-20s -> %-20s  minimap2 FAILED: %s" % (kind, qn, tn, err.strip()[-120:]))
        continue
    qc, tc, idn, plus, rho, n = summarize(paf)
    rows.append((kind, qn, tn, qc, tc, idn, plus, rho, n))
    lines.append("    %-13s %-20s -> %-20s  %5.0f%%  %5.0f%%  %s  %s  %s  (%d alignments)"
                 % (kind, qn, tn, 100 * qc, 100 * tc, "id %.3f" % idn if n else "id   -  ",
                    "+%3.0f%%" % (100 * plus) if n else "   -  ", "rho %+.2f" % rho if rho == rho else "rho   - ", n))
pd.DataFrame(rows, columns=["kind", "query", "target", "query_cov", "target_cov", "identity", "plus_share",
                            "rho", "n_aln"]).to_csv(os.path.join(OUT, "%s_arm_pairs.tsv" % S), sep="\t", index=False)

# ---- raw HiFi calls: hap1 tails against their chromosome's body
if BCF and os.path.exists(VCF) and (os.path.exists(VCF + ".tbi") or os.path.exists(VCF + ".csi")):
    lines.append("  raw HiFi het SNPs (biallelic), hap1 tail vs the rest of its chromosome:")
    lines.append("    %-20s %12s %12s %14s" % ("region", "het SNP/Mb", "median DP", "median ALT ratio"))
    for t in [t for t in tails if t[0] == "hap1"]:
        h, c, s0, e0, side = t
        for label, regs in (("%s tail" % names[t], ["%s:%d-%d" % (c, s0 + 1, e0)]),
                            ("%s body" % c, ["%s:%d-%d" % (c, e0 + 1, L1[c])] if side == "left"
                             else ["%s:1-%d" % (c, s0)])):
            r = subprocess.run([BCF, "query", "-r", ",".join(regs), "-i", 'TYPE="snp" && GT="het"',
                                "-f", "[%DP]\t[%AD]\n", VCF], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               universal_newlines=True)
            dp, ar = [], []
            for l in r.stdout.splitlines():
                f = l.split("\t")
                try:
                    d = int(f[0]); ad = [int(v) for v in f[1].split(",")]
                except (ValueError, IndexError):
                    continue
                if d > 0 and len(ad) == 2:
                    dp.append(d); ar.append(ad[1] / float(d))
            span = sum(int(x.split(":")[1].split("-")[1]) - int(x.split(":")[1].split("-")[0]) + 1 for x in regs)
            lines.append("    %-20s %12.0f %12s %14s" % (label, len(dp) / (span / 1e6),
                         "%.0f" % np.median(dp) if dp else "-", "%.2f" % np.median(ar) if ar else "-"))
else:
    lines.append("  raw HiFi calls: skipped (bcftools, %s or its index not found)" % VCF)
# ---- baseline: identity of ordinary homologs in the genome-wide PAF (primary, MAPQ >= 5, >= 50 kb)
if os.path.exists(PAF):
    mt = bl = 0
    with open(PAF) as fh:
        for l in fh:
            f = l.split("\t", 12)
            if len(f) >= 12 and "tp:A:S" not in l and int(f[11]) >= 5 and int(f[3]) - int(f[2]) >= 50000 \
                    and base(f[0]) == base(f[5]):
                mt += int(f[9]); bl += int(f[10])
    lines.append("  baseline: ordinary hap1-hap2 homologs align at identity %.3f (genome-wide PAF, primary, "
                 "MAPQ >= 5, >= 50 kb)" % (mt / float(max(bl, 1))))

# ---- HiFi depth: hap1 tails against their chromosome's body
if A.bam and SAMT and os.path.exists(A.bam):
    def depth(reg):
        r = subprocess.run([SAMT, "coverage", "-r", reg, A.bam], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           universal_newlines=True)
        rows = [l.split("\t") for l in r.stdout.splitlines() if l and not l.startswith("#")]
        return float(rows[0][6]) if rows else float("nan")
    lines.append("  HiFi mean depth (%s), 10 Mb from the middle of each hap1 tail vs 10 Mb from its body:" % A.bam)
    for t in [t for t in tails if t[0] == "hap1"]:
        h, c, s0, e0, side = t
        tm = (s0 + e0) // 2
        bs0, be0 = (e0, L1[c]) if side == "left" else (0, s0)
        bm = (bs0 + be0) // 2
        dt = depth("%s:%d-%d" % (c, tm - 5000000, tm + 5000000))
        db = depth("%s:%d-%d" % (c, bm - 5000000, bm + 5000000))
        lines.append("    %-20s tail %6.1fx   body %6.1fx   ratio %.2f" % (names[t], dt, db, dt / db if db else float("nan")))
elif A.bam:
    lines.append("  HiFi depth: skipped (samtools or %s not found)" % A.bam)
open(os.path.join(OUT, "%s_arm_partners.txt" % S), "w").write("\n".join(lines) + "\n")
print("\n".join(lines[printed:]))
print("wrote %s/%s_arm_partners.txt, %s_arm_pairs.tsv and one PAF per pair" % (OUT, S, S))
