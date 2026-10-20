#!/usr/bin/env python3
"""recon_audit.py -- the exact coordinates a pollen-consistent D. paradoxa karyotype would rest on, read from files.

Nothing typed in except the chr5/chr6 coupling windows, which are printed next to what the karyotype test wrote so
they can be checked. Run from the CO_smk root.
  1. ASSEMBLY: the two-haplotype assembly the crossover reference was cut from (source_fasta in the reference's
     MANIFEST.txt) and the length of every chromosome.
  2. JOINS: every place a chromosome changes homolog (qc/reference/SAMPLE/breakpoints.csv: the last base aligned
     to the left homolog, the first base aligned to the right one) and the cut the crossover reference used (its
     split plan).
  3. P AND Q: the arms with no homolog in the other haplotype (P twice in hap1, Q twice in hap2) aligned to each
     other (minimap2 asm20, chains only): where each starts, which strand, how collinear. That closes the junction
     stretch on the P/Q side and says which way Q must face on A's L1, and P on D's L2.
  4. GAPS: every run of >= --min_gap N in the assembly (gaps.tsv) and the ones at each join. A join at a gap is a
     scaffolding join (two contigs put side by side by the Hi-C); a join with no gap runs through assembled (HiFi)
     sequence, i.e. one molecule.
  5. RECIPE: the two pollen chromosomes the assembly lacks, A' = A's L1 + D's Q and D' = D's L2 + A's P, as exact
     coordinate ranges, or why they cannot be written yet.
  6. chr5/chr6: scaffold gaps inside the pollen coupling windows on every copy (hap2 windows translated through
     hap2_on_hap1.paf). In the haplotype that carries a swap, its chr5 and its chr6 must both be broken, i.e. carry
     a scaffold gap inside the window (a contig running through is one molecule).
  7. rDNA arrays near the joins and inside the windows (reports/repeat_context/rdna_loci_SAMPLE.csv).
Names (definitions, not data): A = chr1_hap1 = L1 + P, B = chr2_hap1 = L2 + P, C = chr1_hap2 = L1 + Q,
D = chr2_hap2 = L2 + Q; X = the arm on chr4_hap1 (X + core) and on chr3_hap2 (chr3 + X).
Writes OUT/{recon_audit.txt, gaps.tsv, pq_hap1_*.paf, pq_hap2_*.paf} (the PAFs are reused when the arms are unchanged).
Usage: recon_audit.py [--sample Dparadoxa_split] [--out reports/recon_audit] [--threads 16] [--min_gap 10]
                      [--near_mb 2] [--min_aln_kb 20] [--min_mapq 5] [--chr5 chr5_hap1:90:225] [--chr6 chr6_hap1:120:295]
                      [--no_align]
"""
import argparse
import csv
import datetime
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--sample", default="Dparadoxa_split")
ap.add_argument("--out", default="reports/recon_audit")
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--min_gap", type=int, default=10)
ap.add_argument("--near_mb", type=float, default=2.0)
ap.add_argument("--min_aln_kb", type=float, default=20.0)
ap.add_argument("--min_mapq", type=int, default=5)
ap.add_argument("--chr5", default="chr5_hap1:90:225")
ap.add_argument("--chr6", default="chr6_hap1:120:295")
ap.add_argument("--no_align", action="store_true")
A = ap.parse_args()
MB = 1000000
S = A.sample
NEAR = int(A.near_mb * MB)
MINALN = int(A.min_aln_kb * 1000)
os.makedirs(os.path.join(A.out, "arms"), exist_ok=True)
OUT, READ, MISSING = [], [], []
LET = {"chr1_hap1": "A (L1|P)", "chr2_hap1": "B (L2|P)", "chr1_hap2": "C (L1|Q)", "chr2_hap2": "D (L2|Q)",
       "chr4_hap1": "X|core", "chr3_hap2": "chr3|X"}


def say(s=""):
    OUT.append(s)
    print(s, flush=True)


def mb(x):
    return "%.3f" % (x / MB)


def num(c):
    m = re.search(r"chr(\d+)", c)
    return int(m.group(1)) if m else 99


def hapof(c):
    m = re.search(r"_(hap\d)$", c)
    return m.group(1) if m else "?"


def natural(cs):
    return sorted(cs, key=lambda c: (num(c), hapof(c), c))


def ival(s):
    s = (s or "").strip()
    try:
        return int(float(s))
    except ValueError:
        return None


def have(p):
    if os.path.exists(p):
        st = os.stat(p)
        READ.append("%s  (%s, %.1f MB)" % (p, datetime.datetime.fromtimestamp(st.st_mtime).strftime("%Y-%m-%d %H:%M"),
                                          st.st_size / 1e6))
        return True
    MISSING.append(p)
    return False


def git_head():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL, universal_newlines=True).stdout.strip() or "?"
    except OSError:
        return "?"


# ------------------------------------------------------------------ 1. assembly
man = "results/reference/%s/MANIFEST.txt" % S
if not have(man):
    sys.exit("STOP: %s is missing" % man)
meta = {}
for line in open(man):
    if ":" in line:
        k, v = line.split(":", 1)
        meta[k.strip()] = v.strip()
SRC = meta.get("source_fasta", "")
if not SRC or not have(SRC):
    sys.exit("STOP: %s names source_fasta %r, which does not exist" % (man, SRC))
if SRC.endswith((".gz", ".bgz", ".bgzf")):
    sys.exit("STOP: %s is compressed; this script reads a plain FASTA through its .fai" % SRC)
FAI = SRC + ".fai"
if not have(FAI):
    sys.exit("STOP: %s is missing" % FAI)
FX = {}
for line in open(FAI):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 5:
        FX[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))
LEN = {c: v[0] for c, v in FX.items()}
FA = open(SRC, "rb")


def fetch(c, start=1, end=None):
    """1-based inclusive slice of a sequence, read through the .fai (no samtools needed)."""
    n, off, lb, lw = FX[c]
    end = n if end is None else min(end, n)
    s0, e0 = start - 1, end - 1
    b0 = off + (s0 // lb) * lw + s0 % lb
    b1 = off + (e0 // lb) * lw + e0 % lb + 1
    FA.seek(b0)
    return FA.read(b1 - b0).replace(b"\n", b"").replace(b"\r", b"")


say("RECONSTRUCTION AUDIT  %s  (%s, CO_smk commit %s)" % (S, datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
                                                        git_head()))
say("Names: A = chr1_hap1 = L1 + P; B = chr2_hap1 = L2 + P; C = chr1_hap2 = L1 + Q; D = chr2_hap2 = L2 + Q.")
say("       Pollen: A' = A's L1 + D's Q and D' = D's L2 + A's P (not in the assembly); B and C as assembled.")
say("       X = the arm on chr4_hap1 (X + core) and on chr3_hap2 (chr3 + X). Coordinates 1-based, inclusive.")
say("")
say("1. ASSEMBLY  %s" % SRC)
CHR = natural([c for c in LEN if re.match(r"^chr\d+_hap\d$", c)])
if not CHR:
    sys.exit("STOP: no chrN_hapM sequences in %s" % FAI)
say("   %-11s %13s %10s" % ("chromosome", "length (bp)", "Mb"))
tot = {}
for c in CHR:
    say("   %-11s %13d %10s" % (c, LEN[c], mb(LEN[c])))
    tot[hapof(c)] = tot.get(hapof(c), 0) + LEN[c]
say("   " + "; ".join("%s %s Mb" % (h, mb(v)) for h, v in sorted(tot.items())))
rest = [c for c in LEN if c not in CHR]
if rest:
    say("   plus %d other sequence(s), %s Mb (unplaced; not scanned)" % (len(rest), mb(sum(LEN[c] for c in rest))))

# ------------------------------------------------------------------ gaps (scanned first, reported in section 4)
print("scanning %s for runs of >= %d N ..." % (SRC, A.min_gap), file=sys.stderr, flush=True)
PAT = re.compile(rb"[Nn]{%d,}" % A.min_gap)
GAPS = {}
for c in CHR:
    seq = fetch(c)
    GAPS[c] = [(m.start() + 1, m.end()) for m in PAT.finditer(seq)]
    del seq
with open(os.path.join(A.out, "gaps.tsv"), "w") as fo:
    fo.write("chrom\tstart\tend\tlength\n")
    for c in natural(GAPS):
        for s_, e_ in GAPS[c]:
            fo.write("%s\t%d\t%d\t%d\n" % (c, s_, e_, e_ - s_ + 1))


def gaps_in(c, a, b):
    return [(s_, e_) for s_, e_ in GAPS.get(c, []) if e_ >= a and s_ <= b]


def gap_txt(g, ref=None):
    s_, e_ = g
    t = "%s-%s (%d N, at %s Mb)" % (s_, e_, e_ - s_ + 1, mb(s_))
    if ref is not None:
        d = 0 if s_ <= ref <= e_ else (s_ - ref if s_ > ref else ref - e_)
        t += ", %.1f kb %s the cut" % (d / 1000.0, "at" if d == 0 else ("right of" if s_ > ref else "left of"))
    return t


# ------------------------------------------------------------------ 2. joins and the plan
bpf = "qc/reference/%s/breakpoints.csv" % S
if not have(bpf):
    sys.exit("STOP: %s is missing" % bpf)
BPS = []
for r in csv.DictReader(open(bpf)):
    BPS.append(dict(chrom=r["chrom"], hap=r.get("haplotype") or hapof(r["chrom"]), cut=ival(r["cut"]),
                    left=(r.get("left_homolog") or "none").strip(), lend=ival(r.get("left_last_base")),
                    right=(r.get("right_homolog") or "none").strip(), rstart=ival(r.get("right_first_base")),
                    inref=(r.get("in_reference") or "").strip()))
BPS.sort(key=lambda b: (num(b["chrom"]), hapof(b["chrom"]), b["cut"] or 0))
say("")
say("2. JOINS: where each chromosome changes homolog  (%s)" % bpf)
say("   %-11s %-9s %12s   %-22s %-22s %s" % ("chromosome", "join", "cut", "left homolog: last base", "right homolog: first base", "in crossover ref"))
for b in BPS:
    say("   %-11s %-9s %12s   %-22s %-22s %s" % (
        b["chrom"], LET.get(b["chrom"], ""), b["cut"],
        "%s: %s" % (b["left"], b["lend"] if b["lend"] is not None else "-"),
        "%s: %s" % (b["right"], b["rstart"] if b["rstart"] is not None else "-"), b["inref"]))
plan = "config/split_plans/%s.csv" % S
if have(plan):
    rows = [l for l in open(plan) if l.strip() and not l.startswith("#")]
    say("   crossover reference pieces (%s):" % plan)
    for r in csv.DictReader(rows):
        say("     %-18s %-11s %12s %12s  %s Mb" % (r["piece"], r["source"], r["start"], r["end"], mb(int(r["length"]))))

# ------------------------------------------------------------------ 3. P and Q against their second copies
ARMS = []
for b in BPS:
    if b["cut"] is None or b["chrom"] not in LEN:
        continue
    if b["right"] == "none":
        ARMS.append(dict(chrom=b["chrom"], start=b["cut"] + 1, end=LEN[b["chrom"]], join_low=True, bp=b))
    if b["left"] == "none":
        ARMS.append(dict(chrom=b["chrom"], start=1, end=b["cut"], join_low=False, bp=b))
MM = shutil.which("minimap2")
PQ = {}       # hap -> summary
say("")
say("3. ARMS WITH NO HOMOLOG IN THE OTHER HAPLOTYPE, aligned to each other (minimap2 asm20, MAPQ >= %d, >= %g kb)"
    % (A.min_mapq, A.min_aln_kb))


def write_arm(a):
    p = os.path.join(A.out, "arms", "%s_%d_%d.fa" % (a["chrom"], a["start"], a["end"]))
    if not (os.path.exists(p) and os.path.getsize(p) > a["end"] - a["start"]):
        seq = fetch(a["chrom"], a["start"], a["end"])
        with open(p + ".tmp", "wb") as fo:
            fo.write((">%s:%d-%d\n" % (a["chrom"], a["start"], a["end"])).encode())
            for i in range(0, len(seq), 80):
                fo.write(seq[i:i + 80] + b"\n")
        os.replace(p + ".tmp", p)
    return p


def merged(iv):
    out = []
    for s_, e_ in sorted(iv):
        if out and s_ <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e_)
        else:
            out.append([s_, e_])
    return sum(e_ - s_ for s_, e_ in out)


def spearman(x, y):
    if len(x) < 3:
        return float("nan")
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


for h in ("hap1", "hap2"):
    arms = sorted([a for a in ARMS if hapof(a["chrom"]) == h], key=lambda a: (num(a["chrom"]), a["start"]))
    if len(arms) != 2:
        say("   %s: %d such arm(s) (%s); expected 2, not aligned" % (h, len(arms), ", ".join(
            "%s:%d-%d" % (a["chrom"], a["start"], a["end"]) for a in arms) or "none"))
        continue
    t, q = arms
    say("   %s: target %s:%d-%d (%s Mb) <- query %s:%d-%d (%s Mb)" % (
        h, t["chrom"], t["start"], t["end"], mb(t["end"] - t["start"] + 1), q["chrom"], q["start"], q["end"],
        mb(q["end"] - q["start"] + 1)))
    paf = os.path.join(A.out, "pq_%s_%s_%d_%d__%s_%d_%d.paf" % (h, t["chrom"], t["start"], t["end"],
                                                              q["chrom"], q["start"], q["end"]))
    if not os.path.exists(paf) or os.path.getsize(paf) == 0:
        if A.no_align or not MM:
            say("     not aligned (%s)" % ("--no_align" if A.no_align else "minimap2 not on PATH"))
            continue
        tf, qf = write_arm(t), write_arm(q)
        print("minimap2 %s: %s <- %s ..." % (h, tf, qf), file=sys.stderr, flush=True)
        with open(paf + ".tmp", "w") as fo, open(paf + ".log", "w") as fl:
            rc = subprocess.run([MM, "-x", "asm20", "-t", str(A.threads), "--secondary=no", tf, qf],
                                stdout=fo, stderr=fl).returncode
        if rc != 0:
            say("     minimap2 failed (exit %d, see %s.log)" % (rc, paf))
            continue
        os.replace(paf + ".tmp", paf)
        for f_ in (tf, qf):
            os.remove(f_)
    al = []
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12 or int(f[11]) < A.min_mapq or int(f[10]) < MINALN or "tp:A:S" in f[12:]:
            continue
        al.append(dict(qs=int(f[2]), qe=int(f[3]), st=f[4], ts=int(f[7]), te=int(f[8]), al=int(f[10])))
    if not al:
        say("     no alignment passes the filters: the two arms are not homologous here, or too diverged for asm20")
        continue
    byst = {s_: sum(x["qe"] - x["qs"] for x in al if x["st"] == s_) for s_ in "+-"}
    st = max(byst, key=byst.get)
    share = byst[st] / float(sum(byst.values()))
    main = [x for x in al if x["st"] == st]
    rho = spearman([(x["qs"] + x["qe"]) / 2 for x in main], [(x["ts"] + x["te"]) / 2 for x in main])
    qcov = merged([(x["qs"], x["qe"]) for x in al]) / float(q["end"] - q["start"] + 1)
    tcov = merged([(x["ts"], x["te"]) for x in al]) / float(t["end"] - t["start"] + 1)
    first = {t["chrom"]: t["start"] + min(x["ts"] for x in al), q["chrom"]: q["start"] + min(x["qs"] for x in al)}
    last = {t["chrom"]: t["start"] + max(x["te"] for x in al) - 1, q["chrom"]: q["start"] + max(x["qe"] for x in al) - 1}
    corr = (st == "+") == (q["join_low"] == t["join_low"])
    PQ[h] = dict(t=t, q=q, strand=st, share=share, rho=rho, first=first, last=last, corr=corr)
    say("     %d alignments; strand %s carries %.0f%% of the aligned bp; collinearity (Spearman) %s" % (
        len(al), st, 100 * share, "n/a (fewer than 3 alignments)" if rho != rho else "%.2f" % rho))
    say("     covered: %.0f%% of %s, %.0f%% of %s" % (100 * tcov, t["chrom"], 100 * qcov, q["chrom"]))
    for a in (t, q):
        say("     %s: first aligned base %d (%s Mb), last %d (%s Mb); its join end is at %s" % (
            a["chrom"], first[a["chrom"]], mb(first[a["chrom"]]), last[a["chrom"]], mb(last[a["chrom"]]),
            "%d (the arm's low end)" % a["start"] if a["join_low"] else "%d (the arm's high end)" % a["end"]))
    say("     -> the join ends %s" % ("correspond: the two copies face the same way relative to their joins"
                                     if corr else "do NOT correspond: one copy is turned round relative to its join"))

# ------------------------------------------------------------------ 4. gaps at each join
say("")
say("4. SCAFFOLD GAPS (runs of >= %d N; all in %s)" % (A.min_gap, os.path.join(A.out, "gaps.tsv")))
ng = sum(len(v) for v in GAPS.values())
gsize = sorted(e_ - s_ + 1 for v in GAPS.values() for s_, e_ in v)
say("   %d gaps in %s Mb of chromosomes%s" % (ng, mb(sum(LEN[c] for c in CHR)), "; lengths: %s" % ", ".join(
    "%d N x %d" % (L_, gsize.count(L_)) for L_ in sorted(set(gsize))[:8]) if gsize else ""))
for c in CHR:
    v = GAPS[c]
    say("   %-11s %4d gaps%s" % (c, len(v), ("; median spacing %s Mb" % mb(float(np.median(np.diff(
        [s_ for s_, _ in v]))))) if len(v) > 1 else ""))
say("")
say("   AT EACH JOIN. Junction stretch = from the last base aligned to the left homolog to the first base aligned")
say("   to the right one (for P and Q: the first base aligned to the arm's other copy, section 3).")
JUNC = {}
for b in BPS:
    c, cut = b["chrom"], b["cut"]
    if c not in LEN or cut is None:
        continue
    lo = (b["lend"] + 1) if b["lend"] is not None else None
    hi = (b["rstart"] - 1) if b["rstart"] is not None else None
    for h, s in PQ.items():
        if hi is None and b["right"] == "none" and c in s["first"]:
            hi = s["first"][c] - 1
        if lo is None and b["left"] == "none" and c in s["last"]:
            lo = s["last"][c] + 1
    say("   %s %s, cut %d (%s Mb)" % (c, LET.get(c, ""), cut, mb(cut)))
    if lo is not None and hi is not None:
        if hi >= lo:
            say("     junction stretch %d-%d (%.1f kb not aligned to either homolog)" % (lo, hi, (hi - lo + 1) / 1000.0))
        else:
            say("     the two homologs overlap by %.1f kb (%d-%d)" % ((lo - hi - 1) / 1000.0, hi + 1, lo - 1))
        inside = gaps_in(c, min(lo, hi), max(lo, hi))
    else:
        say("     junction stretch open on one side (%s-%s): the other copy of that arm was not aligned" % (lo, hi))
        inside = []
    near = gaps_in(c, cut - NEAR, cut + NEAR)
    if len(inside) == 1:
        verdict = "a scaffold gap: %s" % gap_txt(inside[0], cut)
    elif len(inside) > 1:
        verdict = "%d gaps in the stretch: %s" % (len(inside), "; ".join(gap_txt(g, cut) for g in inside))
    elif near:
        verdict = "no gap in the stretch; within %g Mb: %s" % (A.near_mb, "; ".join(gap_txt(g, cut) for g in near))
    else:
        verdict = "NO gap within %g Mb: the homolog switch runs through contiguous assembled sequence" % A.near_mb
    say("     " + verdict)
    JUNC[c] = dict(lo=lo, hi=hi, inside=inside, near=near)

# ------------------------------------------------------------------ 5. recipe for A' and D'
say("")
say("5. RECIPE: the two pollen chromosomes the assembly lacks")


def recipe(name, base, moved, h, what):
    """base keeps its left arm (1..gap start - 1); moved contributes its arm right of its own gap."""
    jb, jm = JUNC.get(base), JUNC.get(moved)
    why = []
    bb = [b for b in BPS if b["chrom"] == base]
    if len(bb) != 1 or bb[0]["left"] == "none" or bb[0]["right"] != "none":
        why.append("%s: expected one join with its homologous arm on the left and the P/Q arm on the right" % base)
    for c, j in ((base, jb), (moved, jm)):
        if j is None:
            why.append("%s: no join found" % c)
        elif len(j["inside"]) != 1:
            why.append("%s: %d scaffold gaps in its junction stretch (need exactly 1)" % (c, len(j["inside"])))
    if h not in PQ:
        why.append("%s arms not aligned: orientation unknown" % h)
    if why:
        say("   %s = %s: cannot be written exactly yet: %s" % (name, what, "; ".join(why)))
        return None
    gb, gm = jb["inside"][0], jm["inside"][0]
    left = (base, 1, gb[0] - 1)
    s = PQ[h]
    arm = s["q"] if s["q"]["chrom"] == moved else s["t"]
    fwd = (s["corr"] == arm["join_low"])
    right = (moved, gm[1] + 1, LEN[moved]) if arm["join_low"] else (moved, 1, gm[0] - 1)
    n = gb[1] - gb[0] + 1
    total = (left[2] - left[1] + 1) + n + (right[2] - right[1] + 1)
    say("   %s = %s:" % (name, what))
    say("        %s:%d-%d   (%s Mb, forward)" % (left[0], left[1], left[2], mb(left[2] - left[1] + 1)))
    say("      + %d N         (%s's own join gap; the true junction sequence is not assembled)" % (n, base))
    say("      + %s:%d-%d   (%s Mb, %s)" % (right[0], right[1], right[2], mb(right[2] - right[1] + 1),
                                          "forward" if fwd else "REVERSE-COMPLEMENT"))
    say("      = %d bp (%s Mb)" % (total, mb(total)))
    return total


ta = recipe("A'", "chr1_hap1", "chr2_hap2", "hap2", "A's L1 + D's Q")
td = recipe("D'", "chr2_hap2", "chr1_hap1", "hap1", "D's L2 + A's P")
if ta and td:
    ad = LEN["chr1_hap1"] + LEN["chr2_hap2"]
    say("   check: A' + D' = %d bp; A + D = %d bp: %s" % (ta + td, ad, "equal, nothing added or lost" if ta + td == ad
                                                       else "DIFFER by %d bp" % (ta + td - ad)))
say("   B (chr2_hap1), C (chr1_hap2), chr3/chr4 and chr5/chr6: as assembled.")

# ------------------------------------------------------------------ 6. chr5/chr6
say("")
say("6. chr5/chr6: SCAFFOLD GAPS INSIDE THE POLLEN COUPLING WINDOWS")
pafs = sorted(glob.glob("results/hap_align/*/hap2_on_hap1.paf"))
pafs = [p for p in pafs if "paradoxa" in p.lower()] or pafs
H2H1 = pafs[0] if pafs else ""
WIN = {}
for spec in (A.chr5, A.chr6):
    c, a, b = spec.split(":")
    WIN[c] = (int(float(a) * MB), int(float(b) * MB))
    say("   window on %s: %s-%s Mb (from the karyotype test's 95%% interval; check it below)" % (c, a, b))
if H2H1 and have(H2H1):
    want = {c: c.replace("_hap1", "_hap2") for c in list(WIN) if c.endswith("_hap1")}
    hits = {c: [] for c in want}
    for line in open(H2H1):
        f = line.split("\t", 12)
        if len(f) < 12 or f[5] not in want or f[0] != want[f[5]]:
            continue
        if int(f[11]) < A.min_mapq or int(f[10]) < MINALN:
            continue
        qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
        a, b = WIN[f[5]]
        lo, hi = max(ts, a), min(te, b)
        if lo >= hi or te <= ts:
            continue
        r0, r1 = (lo - ts) / float(te - ts), (hi - ts) / float(te - ts)     # where the window cuts this alignment
        if f[4] == "+":
            hits[f[5]].append((qs + r0 * (qe - qs), qs + r1 * (qe - qs)))
        else:
            hits[f[5]].append((qe - r1 * (qe - qs), qe - r0 * (qe - qs)))
    for c1, c2 in want.items():
        if hits[c1]:
            WIN[c2] = (int(min(s_ for s_, _ in hits[c1])), int(max(e_ for _, e_ in hits[c1])))
            say("   window on %s: %s-%s Mb (the part aligned to %s's window, %s)" % (
                c2, mb(WIN[c2][0]), mb(WIN[c2][1]), c1, H2H1))
        else:
            say("   window on %s: not translated (no alignment to %s's window in %s)" % (c2, c1, H2H1))
else:
    say("   hap2 windows not translated: no results/hap_align/*/hap2_on_hap1.paf")
CNT = {}
for c in natural(WIN):
    a, b = WIN[c]
    g = gaps_in(c, a, b)
    CNT[c] = len(g)
    say("   %-11s %s-%s Mb: %d gap(s)%s" % (c, mb(a), mb(b), len(g), (": " + ", ".join(
        "%s (%d N)" % (mb(s_), e_ - s_ + 1) for s_, e_ in g[:40]) + (" ..." if len(g) > 40 else "")) if g else ""))
for h in ("hap1", "hap2"):
    c5, c6 = "chr5_" + h, "chr6_" + h
    if c5 in CNT and c6 in CNT:
        say("   if %s carries the swap: %d x %d gap pairs x 2 orientations = %d candidate joins%s" % (
            h, CNT[c5], CNT[c6], 2 * CNT[c5] * CNT[c6],
            "" if CNT[c5] and CNT[c6] else "  (none: that haplotype cannot carry it unless a contig is chimeric)"))
coup = []
for p in sorted(glob.glob("reports/cosegregation_2n/*")):
    if p.endswith((".txt", ".tsv", ".csv")):
        for line in open(p, errors="replace"):
            if "coupl" in line.lower():
                coup.append("%s: %s" % (os.path.basename(p), line.rstrip()[:200]))
say("   what the karyotype test wrote about the coupling point (%d line(s)):" % len(coup))
for line in coup[:15]:
    say("     " + line)
bpt = "qc/reference/%s/breakpoints.txt" % S
if os.path.exists(bpt):
    isl = [l.rstrip() for l in open(bpt) if "island" in l.lower() and re.search(r"chr[56]_", l)]
    say("   islands on chr5/chr6 in %s: %s" % (bpt, "none" if not isl else ""))
    for line in isl[:20]:
        say("     " + line.strip()[:200])

# ------------------------------------------------------------------ 7. rDNA
say("")
rd = "reports/repeat_context/rdna_loci_%s.csv" % S
say("7. rDNA ARRAYS  (%s)" % rd)
if have(rd):
    R = list(csv.DictReader(open(rd)))
    if R and all(k in R[0] for k in ("type", "copy", "start", "end")):
        joins = [(b["chrom"], b["cut"]) for b in BPS if b["cut"] is not None]
        say("   45S (all):")
        for r in [r for r in R if r["type"].startswith("45")]:
            say("     %-11s %12s-%-12s %6s kb  %s" % (r["copy"], r["start"], r["end"], r.get("kb", ""), r.get("at_join", "")))
        say("   5S within %g Mb of a join:" % (2.5 * A.near_mb))
        for r in R:
            if not r["type"].startswith("5S"):
                continue
            s_, e_ = int(float(r["start"])), int(float(r["end"]))
            d = [abs((s_ + e_) / 2 - p) for c, p in joins if c == r["copy"]]
            if d and min(d) <= 2.5 * NEAR:
                say("     %-11s %12d-%-12d %6s kb  %.2f Mb from the join" % (r["copy"], s_, e_, r.get("kb", ""), min(d) / MB))
        say("   any rDNA inside the chr5/chr6 windows:")
        for r in R:
            c = r["copy"]
            if c in WIN:
                s_, e_ = int(float(r["start"])), int(float(r["end"]))
                if e_ >= WIN[c][0] and s_ <= WIN[c][1]:
                    say("     %-4s %-11s %12d-%-12d %6s kb" % (r["type"], c, s_, e_, r.get("kb", "")))
        say("   (%d arrays in the table: %s)" % (len(R), ", ".join("%s %d" % (t_, sum(1 for r in R if r["type"] == t_))
                                                                  for t_ in sorted(set(r["type"] for r in R)))))
    else:
        say("   unexpected columns: %s" % (", ".join(R[0].keys()) if R else "empty table"))

say("")
say("files read:")
for p in READ:
    say("   " + p)
if MISSING:
    say("files missing:")
    for p in MISSING:
        say("   " + p)
open(os.path.join(A.out, "recon_audit.txt"), "w").write("\n".join(OUT) + "\n")
print("wrote %s and gaps.tsv" % os.path.join(A.out, "recon_audit.txt"), file=sys.stderr)
