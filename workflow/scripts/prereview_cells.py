#!/usr/bin/env python3
"""prereview_cells.py -- read-only look at a STARsolo library before cell calling.

Checks, in the order they matter:
  1. Chemistry vs config. Are the cDNA reads sense or antisense to their genes
     (that decides --soloStrand), and do they pile up at 3' or 5' ends? Reads
     just past annotated gene ends mean 3' UTRs are missing from the
     annotation: STARsolo cannot count those.
  2. Barcodes. log10(UMI) histogram (a valley or a smear), counts above
     100/300/500/1000, UMI share of the top barcodes, side by side with a
     reference library.
  3. Ambient. UMIs in an empty droplet (the hump), its top genes, and whether
     the 300-500 UMI band looks like ambient or like the top barcodes.
  4. Contamination. Ambient share, and the expected wrong-allele rate, at each
     UMI total: the practical floor for genotyping a nucleus.
  5. Intronic fraction (GeneFull vs Gene) and rDNA-gene share by UMI bin.

Read-only. Python standard library only; samtools reads the BAM (check 1).
"""
import argparse
import array
import bisect
import gzip
import math
import os
import random
import re
import shutil
import subprocess
import sys
from collections import defaultdict

RE_CIG = re.compile(r"(\d+)([MIDNSHP=X])")
LINES = []


def say(s=""):
    LINES.append(s)
    print(s, flush=True)


def note(s):
    print("  .. " + s, file=sys.stderr, flush=True)


def pct(x):
    return "%5.1f%%" % (100.0 * x) if x == x else "   n/a"


def bar(x, width):
    return "#" * int(round(x * width))


# ------------------------------------------------------------------ inputs
def star_settings(solo):
    """--soloStrand, whitelist and features as STAR actually ran them (Log.out)."""
    log = os.path.join(os.path.dirname(os.path.abspath(solo.rstrip("/"))), "Log.out")
    got = {}
    try:
        with open(log, errors="replace") as f:
            for i, l in enumerate(f):
                for key in ("soloStrand", "soloCBwhitelist", "soloFeatures"):
                    m = re.search(r"--%s\s+(.+?)(?:\s+--|$)" % key, l)
                    if m and key not in got:
                        got[key] = m.group(1).strip()
                if i > 400 or len(got) == 3:
                    break
    except OSError:
        pass
    return got


def read_genes(gff):
    """(chrom, start, end, strand, id) of every gene; mRNA/transcript if no genes."""
    genes, tx = [], []
    with open(gff, errors="replace") as f:
        for l in f:
            if l.startswith("#"):
                continue
            t = l.rstrip("\n").split("\t")
            if len(t) < 9 or t[6] not in "+-":
                continue
            m = re.search(r"(?:^|;)ID=([^;]+)", t[8])
            if not m:
                continue
            rec = (t[0], int(t[3]), int(t[4]), t[6], m.group(1))
            if t[2] == "gene":
                genes.append(rec)
            elif t[2] in ("mRNA", "transcript"):
                tx.append(rec)
    return genes or tx


def read_features(solo, feature):
    with open(os.path.join(solo, feature, "raw", "features.tsv")) as f:
        return [l.split("\t")[0].strip() for l in f]


def read_barcodes(solo, feature):
    with open(os.path.join(solo, feature, "raw", "barcodes.tsv")) as f:
        return [l.strip() for l in f]


def read_mtx(path, keep=False):
    """Column totals and genes detected from a MatrixMarket file (STARsolo raw).
    With keep, also the (row, col, value) triples for a second look."""
    note("reading %s" % path)
    with open(path) as f:
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        nr, nc, _ = (int(x) for x in line.split()[:3])
        tot = array.array("d", bytes(8 * (nc + 1)))
        ngen = array.array("i", bytes(4 * (nc + 1)))
        R, C, V = (array.array("i"), array.array("i"), array.array("f")) if keep else (None, None, None)
        for l in f:
            a, b, c = l.split()[:3]
            b = int(b)
            v = int(c) if c.isdigit() else float(c)
            tot[b] += v
            ngen[b] += 1
            if keep:
                R.append(int(a))
                C.append(b)
                V.append(v)
    return nr, nc, tot, ngen, (R, C, V)


def rdna_rows(genes, bed, features):
    """Matrix rows of genes overlapping the rDNA blacklist BED."""
    if not bed or not os.path.exists(bed):
        return set(), 0, 0
    iv = defaultdict(list)
    with open(bed) as f:
        for l in f:
            t = l.split()
            if len(t) >= 3 and not l.startswith(("#", "track")):
                iv[t[0]].append((int(t[1]), int(t[2])))
    hit = {g[4] for g in genes if any(g[1] <= e and g[2] > s for s, e in iv.get(g[0], ()))}
    row = {fid: i + 1 for i, fid in enumerate(features)}
    strip = {fid.split(":")[-1]: i + 1 for i, fid in enumerate(features)}
    rows = set()
    for gid in hit:
        r = row.get(gid) or strip.get(gid.split(":")[-1])
        if r:
            rows.add(r)
    return rows, len(hit), len(rows)


# ------------------------------------------------------------------ 1. strand and 3'/5'
def isolated(genes, gap):
    """Genes with no other gene (either strand) within `gap` bp."""
    by = defaultdict(list)
    for g in genes:
        by[g[0]].append(g)
    keep = []
    for gs in by.values():
        gs.sort(key=lambda g: g[1])
        maxend = -10 ** 15
        for i, g in enumerate(gs):
            nxt = gs[i + 1][1] if i + 1 < len(gs) else 10 ** 15
            if maxend < g[1] - gap and nxt > g[2] + gap:
                keep.append(g)
            maxend = max(maxend, g[2])
    return keep


def strand_check(args, genes, configured):
    say("1. CHEMISTRY vs CONFIG")
    say("   STAR ran with: --soloStrand %s | whitelist %s | features %s"
        % (configured.get("soloStrand", "?"), os.path.basename(configured.get("soloCBwhitelist", "?")),
           configured.get("soloFeatures", "?")))
    st = args.samtools or shutil.which("samtools")
    if not st or not os.access(st, os.X_OK):
        say("   samtools not found -- strand and 3'/5' checks skipped (pass --samtools)")
        say()
        return None
    flank = args.flank
    pool = [g for g in isolated(genes, 2 * flank) if 2000 <= g[2] - g[1] + 1 <= 50000]
    rng = random.Random(args.seed)
    pick = rng.sample(pool, min(args.genes, len(pool)))
    if not pick:
        say("   no isolated 2-50 kb genes to sample -- strand check skipped")
        return None
    idx = defaultdict(list)
    for g in pick:
        idx[g[0]].append((max(1, g[1] - flank), g[2] + flank, g))
    for c in idx:
        idx[c].sort()
    starts = {c: [x[0] for x in v] for c, v in idx.items()}
    regions = ["%s:%d-%d" % (g[0], max(1, g[1] - flank), g[2] + flank) for g in sorted(pick)]
    note("strand check: %d isolated genes, up to %d reads" % (len(pick), args.max_reads))
    cmd = [st, "view", "-F", "0x904", "-q", "255", args.bam] + regions
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    nb = 12                                   # 5' flank | gene in 10 bins | 3' flank
    prof = {True: [0] * nb, False: [0] * nb}
    per_gene = defaultdict(lambda: [0, 0])
    n = 0
    for line in p.stdout:
        t = line.split("\t", 6)
        if len(t) < 6:
            continue
        flag, chrom, pos = int(t[1]), t[2], int(t[3])
        span = sum(int(k) for k, op in RE_CIG.findall(t[5]) if op in "MDN=X")
        mid = pos + (span - 1) / 2.0
        v = idx.get(chrom)
        if not v:
            continue
        i = bisect.bisect_right(starts[chrom], mid) - 1
        if i < 0 or mid > v[i][1]:
            continue
        g = v[i][2]
        sense = ("-" if flag & 16 else "+") == g[3]
        length = g[2] - g[1]
        rel = (mid - g[1]) / length if g[3] == "+" else (g[2] - mid) / length
        b = 0 if rel < 0 else nb - 1 if rel > 1 else 1 + min(9, int(rel * 10))
        prof[sense][b] += 1
        per_gene[g[4]][0 if sense else 1] += 1
        n += 1
        if n >= args.max_reads:
            break
    p.terminate()
    err = p.stderr.read().strip() if p.poll() not in (None, 0, -15) else ""
    if n == 0:
        say("   no reads returned by samtools%s" % ((": " + err[:200]) if err else ""))
        return None
    ns, na = sum(prof[True]), sum(prof[False])
    frac = [s / float(s + a) for s, a in per_gene.values() if s + a >= 20]
    frac.sort()
    med = frac[len(frac) // 2] if frac else float("nan")
    say("   sampled %d isolated genes (2-50 kb, no neighbour within %d kb); %s unique reads"
        % (len(pick), 2 * flank // 1000, format(n, ",")))
    say("   sense to their gene   %s   (per-gene median %s over %d genes with >=20 reads)"
        % (pct(ns / float(n)), pct(med), len(frac)))
    say("   antisense             %s" % pct(na / float(n)))
    if frac:
        say("   genes >=80%% sense: %d | >=80%% antisense: %d | mixed: %d"
            % (sum(1 for f in frac if f >= 0.8), sum(1 for f in frac if f <= 0.2),
               sum(1 for f in frac if 0.2 < f < 0.8)))
    share = ns / float(n)
    if share >= 0.8:
        verdict, right = "SENSE", "Forward"
    elif share <= 0.2:
        verdict, right = "ANTISENSE", "Reverse"
    else:
        verdict, right = "MIXED", "Unstranded?"
    conf = configured.get("soloStrand", "?")
    ok = right == conf
    say("   -> reads are %s to their genes: --soloStrand %s counts them.  %s"
        % (verdict, right, "[matches the config]" if ok else "[CONFIG HAS %s -- MISMATCH]" % conf))
    major = share >= 0.5
    pr = prof[major]
    tot = float(sum(pr)) or 1.0
    labels = ["5' flank 1 kb"] + ["gene %3d-%3d%%" % (10 * k, 10 * k + 10) for k in range(10)] + ["3' flank 1 kb"]
    say("   where %s reads sit along genes, 5' -> 3':" % ("sense" if major else "antisense"))
    top = max(pr) / tot or 1.0
    for lab, c in zip(labels, pr):
        say("     %-14s %s %s" % (lab, pct(c / tot), bar(c / tot / top, 40)))
    five, three, past = (pr[0] + pr[1]) / tot, (pr[10] + pr[11]) / tot, pr[11] / tot
    end = "3'" if three > 2 * five else "5'" if five > 2 * three else "neither end"
    say("   -> piles up at the %s end: %s in the 3' tenth + 3' flank, %s in the 5' tenth + 5' flank"
        % (end, pct(three), pct(five)))
    if past >= 0.15:
        say("   -> %s of reads lie past the annotated 3' end: 3' UTRs are missing from the annotation,"
            % pct(past).strip())
        say("      and STARsolo counts none of those reads")
    say()
    countable = share if conf == "Forward" else (1 - share) if conf == "Reverse" else 1.0
    return dict(ok=ok, right=right, conf=conf, verdict=verdict, end=end, past=past, countable=countable)


# ------------------------------------------------------------------ 2-5. barcodes
def hist(tot, lo=0.0, hi=6.0, step=0.1):
    nb = int(round((hi - lo) / step))
    h = [0] * nb
    for v in tot:
        if v > 0:
            k = int((math.log10(v) - lo) / step)
            h[max(0, min(nb - 1, k))] += 1
    return h


def hump(h, step=0.1, lo=1.5, hi=3.5):
    """Centre (UMIs) of the tallest smoothed bin between 10^lo and 10^hi."""
    sm = [(h[max(0, i - 1)] + h[i] + h[min(len(h) - 1, i + 1)]) / 3.0 for i in range(len(h))]
    a, b = int(lo / step), int(hi / step)
    k = max(range(a, b), key=lambda i: sm[i])
    edge = k in (a, b - 1)
    return 10 ** ((k + 0.5) * step), edge


def stats(tot):
    v = sorted((x for x in tot if x > 0), reverse=True)
    s = sum(v) or 1.0
    above = {t: sum(1 for x in v if x >= t) for t in (100, 300, 500, 1000, 2000, 5000)}
    share = {k: sum(v[:k]) / s for k in (1000, 5000, 10000, 20000)}
    return len(v), above, share, v


def affinity(p, q):
    """Hellinger affinity of two count profiles: 1 = identical, 0 = disjoint."""
    sp, sq = float(sum(p.values())) or 1.0, float(sum(q.values())) or 1.0
    return sum(math.sqrt(p[k] / sp * q.get(k, 0) / sq) for k in p)


def barcodes(args, genes):
    solo = args.solo
    feats = read_features(solo, "GeneFull")
    rrna, n_hit, n_rows = rdna_rows(genes, args.rrna_bed, feats)
    nr, nc, gf, ngen, (R, C, V) = read_mtx(os.path.join(solo, "GeneFull", "raw", "matrix.mtx"), keep=True)
    try:
        _, _, ge, _, _ = read_mtx(os.path.join(solo, "Gene", "raw", "matrix.mtx"))
    except OSError:
        ge = None
    ref = None
    if args.ref_solo:
        try:
            ref = read_mtx(os.path.join(args.ref_solo, "GeneFull", "raw", "matrix.mtx"))[2]
        except OSError as e:
            note("reference skipped: %s" % e)
    name, rname = args.name, args.ref_name
    n1, ab1, sh1, sorted1 = stats(gf)
    say("2. BARCODES (GeneFull UMIs per barcode)")
    cols = [(name, n1, ab1, sh1)]
    if ref is not None:
        n2, ab2, sh2, _ = stats(ref)
        cols.append((rname, n2, ab2, sh2))
    say("   %-34s" % "" + "".join("%16s" % c[0] for c in cols))
    say("   %-34s" % "barcodes with >=1 UMI" + "".join("%16s" % format(c[1], ",") for c in cols))
    for t in (100, 300, 500, 1000, 2000, 5000):
        say("   %-34s" % (">= %s UMIs" % format(t, ",")) + "".join("%16s" % format(c[2][t], ",") for c in cols))
    for k in (1000, 5000, 10000, 20000):
        say("   %-34s" % ("UMI share of the top %s" % format(k, ",")) + "".join("%16s" % pct(c[3][k]) for c in cols))
    h1 = hist(gf)
    A, edge = hump(h1)
    hs = [h1] + ([hist(ref)] if ref is not None else [])
    say("   log10(UMI) histogram: barcodes per 0.1 bin, each library scaled to its own tallest bin")
    say("     %8s  %-30s %s" % ("UMIs", name, rname if ref is not None else ""))
    mx = [max(h[5:]) or 1 for h in hs]
    last = max(i for i in range(len(h1)) if any(h[i] for h in hs))
    for i in range(5, last + 1):
        umi = 10 ** (i * 0.1)
        mark = " <- ambient hump" if abs(math.log10(A) - (i + 0.5) * 0.1) < 0.05 else ""
        say("     %8s  %-30s %s%s" % (format(int(round(umi)), ","),
                                        bar(hs[0][i] / float(mx[0]), 30),
                                        bar(hs[1][i] / float(mx[1]), 30) if ref is not None else "", mark))
    if ref is not None:
        Ar, _ = hump(hs[1])
        say("   ambient hump peaks at ~%s UMIs (%s), ~%s (%s)%s" % (format(int(A), ","), name, format(int(Ar), ","),
                                                                 rname, "  [no clear hump: edge of range]" if edge else ""))
    else:
        say("   ambient hump peaks at ~%s UMIs%s" % (format(int(A), ","), "  [no clear hump]" if edge else ""))
    say()

    # ---- 3. ambient profile vs the 300-500 band vs the top barcodes
    lo = max(30.0, A / 2.0)
    topk = min(2000, len(sorted1))
    cut_top = max(1000.0, sorted1[topk - 1] if topk else float("inf"))
    grp = {}
    for c in range(1, nc + 1):
        t = gf[c]
        if t >= cut_top:
            grp[c] = "top"
        elif 300 <= t <= 500:
            grp[c] = "band"
        elif lo <= t <= A:
            grp[c] = "amb"
    prof = {"amb": defaultdict(float), "band": defaultdict(float), "top": defaultdict(float)}
    rr = array.array("d", bytes(8 * (nc + 1)))
    for r, c, v in zip(R, C, V):
        if r in rrna:
            rr[c] += v
            continue
        g = grp.get(c)
        if g:
            prof[g][r] += v
    n_g = {k: sum(1 for x in grp.values() if x == k) for k in prof}
    say("3. AMBIENT (rDNA gene models excluded from the profiles)")
    say("   groups: ambient = %s-%s UMIs, outside the band (n=%s) | band = 300-500 UMIs (n=%s)"
        % (format(int(lo), ","), format(int(A), ","), format(n_g["amb"], ","), format(n_g["band"], ",")))
    say("           top = the %s barcodes with >= %s UMIs" % (format(n_g["top"], ","), format(int(cut_top), ",")))
    if min(n_g.values()) < 20:
        say("   a group has fewer than 20 barcodes -- profile comparison skipped")
    else:
        a_bt, a_ba, a_at = affinity(prof["band"], prof["top"]), affinity(prof["band"], prof["amb"]), \
            affinity(prof["amb"], prof["top"])
        say("   profile similarity (Hellinger affinity, 1 = identical):")
        say("     300-500 band vs ambient       %.3f" % a_ba)
        say("     300-500 band vs top barcodes  %.3f" % a_bt)
        say("     ambient vs top barcodes       %.3f" % a_at)
        if a_ba > a_bt + 0.02:
            say("   -> the 300-500 band looks like AMBIENT, not like the top barcodes")
        elif a_bt > a_ba + 0.02:
            say("   -> the 300-500 band looks like the TOP barcodes: real low-RNA nuclei are likely in it")
        else:
            say("   -> the 300-500 band sits between ambient and the top barcodes: a mixture")
    for k, lab in (("amb", "ambient"), ("top", "top barcodes")):
        s = float(sum(prof[k].values())) or 1.0
        best = sorted(prof[k].items(), key=lambda kv: -kv[1])[:6]
        say("   top genes, %-13s %s" % (lab + ":", "  ".join("%s %.1f%%" % (feats[r - 1], 100 * v / s) for r, v in best)))
    say()

    # ---- 4. contamination
    say("4. CONTAMINATION (~%s ambient UMIs in every droplet, ~50/50 hap1/hap2)" % format(int(A), ","))
    say("   nucleus UMIs   ambient share   expected wrong-allele reads")
    for T in (500, 1000, 2000, 4000, 8000, 16000):
        s = min(1.0, A / T)
        say("   %12s   %13s   %s" % (format(T, ","), pct(s), pct(s / 2)))
    floor = 5 * A
    say("   -> <=20%% ambient needs >= %s UMIs; %s barcodes reach it (%s)"
        % (format(int(floor), ","), format(sum(1 for x in sorted1 if x >= floor), ","), name))
    say()

    # ---- 5. by UMI bin
    say("5. BY UMI BIN (GeneFull)        barcodes   intronic share      rDNA-gene share")
    edges = [30, 100, 300, 500, 1000, 2000, 5000, 10 ** 12]
    for lo_, hi_ in zip(edges[:-1], edges[1:]):
        cs = [c for c in range(1, nc + 1) if lo_ <= gf[c] < hi_]
        if not cs:
            continue
        tf = sum(gf[c] for c in cs)
        intr = (1 - sum(ge[c] for c in cs) / tf) if ge is not None and tf else float("nan")
        rd = sum(rr[c] for c in cs) / tf if tf else float("nan")
        lab = "%s-%s" % (format(lo_, ","), format(hi_, ",")) if hi_ < 10 ** 12 else "%s+" % format(lo_, ",")
        say("   %-26s %12s   %s              %s" % (lab, format(len(cs), ","), pct(intr), pct(rd)))
    say("   (intronic share = 1 - Gene/GeneFull UMIs; rDNA genes: %d models overlap the blacklist, %d found in the matrix)"
        % (n_hit, n_rows))
    say()

    def write_table(path):
        bcs = read_barcodes(solo, "GeneFull")
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with gzip.open(path, "wt") as f:
            f.write("barcode\tgenefull_umi\tgene_umi\tgenes_detected\trdna_umi\n")
            for c in range(1, nc + 1):
                if gf[c] >= 30:
                    f.write("%s\t%d\t%d\t%d\t%d\n" % (bcs[c - 1], gf[c], ge[c] if ge is not None else -1,
                                                      ngen[c], rr[c]))
    return A, floor, edge, write_table


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--solo", required=True, help="Solo.out of the library")
    ap.add_argument("--name", default="sample")
    ap.add_argument("--ref-solo", help="Solo.out of a reference library, for comparison")
    ap.add_argument("--ref-name", default="reference")
    ap.add_argument("--bam", help="STAR's coordinate-sorted, indexed BAM (for check 1)")
    ap.add_argument("--gff", required=True, help="annotation the STAR index was built from (GFF3)")
    ap.add_argument("--rrna-bed", help="rDNA blacklist BED")
    ap.add_argument("--samtools")
    ap.add_argument("--genes", type=int, default=400)
    ap.add_argument("--max-reads", type=int, default=2000000)
    ap.add_argument("--flank", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--report", help="write the report here as well")
    ap.add_argument("--table", help="per-barcode table (tsv.gz), barcodes with >=30 UMIs")
    args = ap.parse_args()

    genes = read_genes(args.gff)
    note("%d gene models" % len(genes))
    configured = star_settings(args.solo)
    say("PRE-REVIEW  %s   (%s)" % (args.name, os.path.abspath(args.solo)))
    say()
    chem = strand_check(args, genes, configured) if args.bam else None
    A, floor, edge, write_table = barcodes(args, genes)

    say("DECISION INPUTS")
    if chem and not chem["ok"]:
        say("   !! STRAND: the matrices were counted with --soloStrand %s, but the reads are %s to their genes."
            % (chem["conf"], chem["verdict"].lower()))
        say("   !! The counts above come from the %s of gene reads on the other strand. Recount with"
            % pct(chem["countable"]).strip())
        say("   !! --soloStrand %s before cell calling; genotypes (allele counts) are unaffected." % chem["right"])
    elif chem:
        say("   strand: config and reads agree (--soloStrand %s); %s library" % (chem["conf"], chem["end"]))
    else:
        say("   !! strand: NOT CHECKED -- settle it before cell calling")
    say("   default_lower: ~%s UMIs, the top of the ambient hump%s" % (format(int(A), ","),
                                                                  "  [hump unclear: look at the histogram]" if edge else ""))
    say("   contamination floor for genotyping: >= %s UMIs (<=20%% ambient) -> record it for cell_qc"
        % format(int(floor), ","))
    if args.report:
        os.makedirs(os.path.dirname(os.path.abspath(args.report)), exist_ok=True)
        with open(args.report, "w") as f:
            f.write("\n".join(LINES) + "\n")
        note("report written to %s" % args.report)
    if args.table:
        write_table(args.table)
        note("per-barcode table written to %s" % args.table)


if __name__ == "__main__":
    main()
