#!/usr/bin/env python3
"""repeat_recombination.py -- do the repeats explain where the pollen crossovers fall, and where the chromosomes
were rearranged?

Q1 CROSSOVERS. Crossover density along every chromosome of the crossover reference (centimorgans per Mb in
   WIN-Mb windows; each call spread evenly over its interval; nuclei masked for dosage on a sequence do not count
   for it), against the CARP repeat composition of the same chromosome copy, next to the other known shapers of
   a crossover landscape: distance to the nearest chromosome end, genes, informative markers (a crossover is only
   seen where markers flank it) and the structure against the other haplotype (collinear, inverted,
   translocated, no homolog; the synteny blocks of translocation_map.py).
   - Spearman rho per variable, raw and partial (after distance to the end and markers). p from rotating the
     variable along its own chromosome (NPERM rotations): a rotation keeps the variable's patchiness and its
     chromosome, and breaks only its alignment with the crossovers. Bonferroni over every variable tested.
   - How much of the ups and downs each set of variables explains: deviance explained by a Poisson model of the
     window counts; the repeats' share against the share that rotated repeat tracks reach.
   - Crossover density within ARRAY Mb of the big tandem arrays (45S and 5S rDNA >= 10 kb, satellite arrays
     >= 50 kb), over the chromosome's own density, against the same arrays rotated along their chromosome.
   - Dead spots (runs of 5 Mb windows without a crossover midpoint, longer than chance; as in
     supervisor_report.py) and what they are made of.
Q2 REARRANGEMENTS. Every join of the two-haplotype assembly where one arm meets another (with the pollen's verdict,
   cosegregation_2n.py) and every inversion or translocation breakpoint between the two haplotypes (edges of
   the synteny blocks): the repeats within NEAR Mb, against random positions on the same chromosome copy.

Reads (run from the CO_smk root): config/samples.csv; the crossover reference's pieces.tsv and
config/split_plans/SAMPLE.csv; the CARP annotation as reproducible_phd publishes it (REFS/SPECIES_repeats/);
results/crossovers/SAMPLE/{co_intervals.bed, co_per_cell.tsv, cells_per_chrom.tsv}; qc/markers/SAMPLE/
marker_classes.tsv.gz; the reference's gene annotation (samples.csv); qc/translocations/HAP1/HAP1_synteny_blocks.tsv;
reports/cosegregation_2n/karyotype_test_SAMPLE.tsv and logs/cosegregation_2n_SAMPLE.log;
qc/reference/SAMPLE/breakpoints.csv. Only the annotation and the crossovers are required; a missing input drops
its part, with a note.
Writes OUT/: recombination_vs_repeats_SAMPLE.png (the tracks), recombination_tests_SAMPLE.png,
breakpoints_vs_repeats_SAMPLE.png, windows_SAMPLE.csv, assoc_SAMPLE.csv, tests_SAMPLE.csv,
breakpoints_SAMPLE.csv, summary_SAMPLE.txt. The annotation, binned, is cached in
results/repeat_recombination/carp_bins_SPECIES.npz and rebuilt when the annotation changes.
Usage: python workflow/scripts/repeat_recombination.py SAMPLE [--win_mb 5]. 1-3 min; 0.5 GB of memory.
"""
import argparse
import array
import csv
import math
import os
import re
import textwrap
import time
import traceback
from urllib.parse import unquote

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

ap = argparse.ArgumentParser(description="Do the repeats explain the crossover landscape and the rearrangements?")
ap.add_argument("sample")
ap.add_argument("--out", default="reports/repeat_recombination")
ap.add_argument("--refs", default="/netscratch/dep_mercier/grp_marques/Aaryan/refs")
ap.add_argument("--coseg", default="reports/cosegregation_2n")
ap.add_argument("--win_mb", type=float, default=0.0,
                help="window of the crossover tests; 0 = the median chromosome / 50, to 0.5 Mb, within 1-5 Mb")
ap.add_argument("--near_mb", type=float, default=1.0, help="Q2: repeats within this distance of a breakpoint")
ap.add_argument("--array_mb", type=float, default=2.5, help="Q1: crossover density within this distance of an array")
ap.add_argument("--nperm", type=int, default=5000, help="rotations / random positions per test")
ap.add_argument("--nboot", type=int, default=300, help="resamples of the nuclei for the density bands")
ap.add_argument("--min_lineage", type=float, default=0.005,
                help="LTR lineages tested: those covering at least this share of the crossover reference")
ap.add_argument("--seed", type=int, default=1)
A = ap.parse_args()
S = A.sample
MB = 1_000_000
FINE = 10_000               # every track is kept in 10 kb bins
MERGE = 20_000              # pieces of one family (or of rDNA) closer than this are one array
RDNA_MIN, SAT_TICK, SAT_BIG = 10_000, 20_000, 50_000
RNG = np.random.default_rng(A.seed)
os.makedirs(A.out, exist_ok=True)
USED, NOTES, FAILED = [], [], []
T0 = time.time()

# repeat classes as in repeat_context.py (same colours), bottom to top
CLASSES = ["Ty3/gypsy", "Ty1/copia", "other LTR", "LINE", "DNA transposon", "simple / low complexity",
           "other repeat", "satellite", "5S rDNA", "45S rDNA"]
CCOL = ["#185FA5", "#0F6E56", "#7FA7D4", "#8DA33B", "#B07A10", "#C9C7BE", "#E3E1DA", "#534AB7", "#E07B39", "#A32D2D"]
I45, I5, ISAT = CLASSES.index("45S rDNA"), CLASSES.index("5S rDNA"), CLASSES.index("satellite")
NK = len(CLASSES)
INK, MUTED, RED = "#2C2C2A", "#6B6A64", "#A32D2D"
KIND_COL = {"collinear": "#B4B2A9", "inverted": "#E07B39", "translocated": "#1D9E75"}


def coarse(cls):
    if cls.startswith("rDNA"):
        return I5 if cls.startswith("rDNA/5S") else I45
    if cls.startswith("Satellite"):
        return ISAT
    if cls.startswith("Class_I/LTR/Ty3_gypsy"):
        return 0
    if cls.startswith("Class_I/LTR/Ty1_copia"):
        return 1
    if cls.startswith("Class_I/LTR"):
        return 2
    if cls.startswith("Class_I/LINE"):
        return 3
    if cls.startswith("Class_II"):
        return 4
    if cls.startswith(("Simple_repeat", "Low_complexity")):
        return 5
    return 6


def lineage(cls):
    """An LTR element's lineage as DANTE names it: the last level under the superfamily (Ty3/gypsy Tekay)."""
    p = [x for x in cls.split("/") if x]
    if len(p) < 3:
        return "LTR, no superfamily"
    sup = {"Ty3_gypsy": "Ty3/gypsy", "Ty1_copia": "Ty1/copia"}.get(p[2], p[2])
    return sup + " " + p[-1] if len(p) > 3 else sup + ", no lineage"


def num(c):
    m = re.match(r"chr(\d+)", c)
    return int(m.group(1)) if m else 0


def hap(c):
    return c[-1]


def fai_lens(p):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(p) if l.strip()}


def note(msg):
    NOTES.append(msg)
    print("note: " + msg, flush=True)


def section(name, fn):
    t = time.time()
    try:
        r = fn()
        print("[%s: %.0f s]" % (name, time.time() - t), flush=True)
        return r
    except Exception as e:
        FAILED.append("%s: %s: %s" % (name, type(e).__name__, e))
        traceback.print_exc()
        return None


def rankavg(x):
    """Ranks 1..n, ties averaged."""
    x = np.asarray(x, float)
    o = np.argsort(x, kind="mergesort")
    xs = x[o]
    b = np.flatnonzero(np.r_[True, xs[1:] != xs[:-1], True])
    r = np.empty(len(x))
    r[o] = np.repeat((b[:-1] + b[1:] - 1) / 2.0 + 1, np.diff(b))
    return r


def corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = math.sqrt(float((a * a).sum() * (b * b).sum()))
    return float((a * b).sum() / d) if d > 0 else float("nan")


def cover(s, e, nb, B, w=None):
    """Sum over intervals [s, e) of (bp in each B-bp bin) x w: a coverage track (w = 1) or a spread track."""
    cov = np.zeros(nb)
    if len(s) == 0:
        return cov
    s = np.asarray(s, np.int64)
    e = np.asarray(e, np.int64)
    w = np.ones(len(s)) if w is None else np.asarray(w, float)
    b0 = np.minimum(s // B, nb - 1)
    b1 = np.minimum((e - 1) // B, nb - 1)
    one = b0 == b1
    np.add.at(cov, b0[one], ((e - s) * w)[one])
    m = ~one
    if m.any():
        np.add.at(cov, b0[m], ((b0[m] + 1) * B - s[m]) * w[m])
        np.add.at(cov, b1[m], (e[m] - b1[m] * B) * w[m])
        d = np.zeros(nb + 1)
        np.add.at(d, b0[m] + 1, B * w[m])
        np.add.at(d, b1[m], -B * w[m])
        cov += np.cumsum(d)[:nb]
    return cov


def merge_iv(s, e, gap):
    """Merge intervals closer than gap: starts, ends, pieces per merged interval."""
    s = np.asarray(s, np.int64)
    e = np.asarray(e, np.int64)
    if not len(s):
        return s, e, np.zeros(0, int)
    o = np.argsort(s, kind="mergesort")
    s, e = s[o], e[o]
    run = np.maximum.accumulate(e)
    new = np.r_[True, s[1:] - run[:-1] > gap]
    st = np.flatnonzero(new)
    return s[new], np.maximum.reduceat(e, st), np.diff(np.r_[st, len(s)])


def ff(x, d=2):
    return "nan" if x is None or not np.isfinite(x) else ("%%.%df" % d) % x


def ptxt(p, n):
    """A rotation / draw p: never printed as 0; '<1/(n+1)' when no draw reached the observed."""
    return ("< %.0e" % (1.0 / (n + 1))) if p <= 1.0 / (n + 1) + 1e-12 else ("%.2g" % p)


# ============================================================ the sample, its assemblies, the crossover reference
rows = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}
if S not in rows:
    raise SystemExit("%s is not in config/samples.csv" % S)
row = rows[S]
SPECIES = row["species"]
SP = SPECIES.replace("_", " ")
refdir = os.path.dirname(row["assembly_fasta"])
L2n = None
for man in (os.path.join(refdir, "MANIFEST.txt"), os.path.join(A.refs, "%s_hap1" % SPECIES, "MANIFEST.txt")):
    if os.path.exists(man):
        m = re.search(r"^source_fasta:\s*(\S+)", open(man).read(), re.M)
        if m and os.path.exists(m.group(1) + ".fai"):
            USED.append(man)
            L2n = fai_lens(m.group(1) + ".fai")
            break
if L2n is None:
    raise SystemExit("no two-haplotype assembly index found via MANIFEST.txt for %s" % S)
COPIES = sorted([c for c in L2n if re.match(r"^chr\d+_hap[12]$", c)], key=lambda c: (num(c), c))
CIDX = {c: i for i, c in enumerate(COPIES)}

PIECES = {}                                     # reference sequence -> (copy, start0, end)
pp = os.path.join(refdir, "pieces.tsv")
if os.path.exists(pp):
    USED.append(pp)
    for r in csv.DictReader(open(pp), delimiter="\t"):
        PIECES[r["piece"]] = (r["source"], int(r["start"]) - 1, int(r["end"]))
else:
    fai = row["assembly_fasta"] + ".fai"
    if not os.path.exists(fai):
        raise SystemExit("no index of the crossover reference: %s" % fai)
    USED.append(fai)
    for c, n in fai_lens(fai).items():
        if c in CIDX:
            if n != L2n[c]:
                raise SystemExit("%s is %d bp in the crossover reference but %d bp in the two-haplotype assembly"
                                 % (c, n, L2n[c]))
            PIECES[c] = (c, 0, n)
bad = sorted({v[0] for v in PIECES.values() if v[0] not in CIDX})
if bad or not PIECES:
    raise SystemExit("reference sequences not on the two-haplotype assembly's chromosomes: %s" % (", ".join(bad) or
                                                                                                    "none found"))
for p_, (c, a, b) in PIECES.items():
    if b > L2n[c]:
        raise SystemExit("piece %s ends at %d, beyond %s (%d bp)" % (p_, b, c, L2n[c]))
REFC = sorted({v[0] for v in PIECES.values()}, key=lambda c: (num(c), c))
CUTS = {c: sorted(a for (src, a, b) in PIECES.values() if src == c and a > 0) for c in REFC}
PSPAN = {c: sorted((a, b, p_) for p_, (src, a, b) in PIECES.items() if src == c) for c in REFC}


def piece_at(c, x):
    for a, b, p_ in PSPAN[c]:
        if a <= x < b:
            return p_
    return PSPAN[c][-1][2]


HOW = {}
plan = "config/split_plans/%s.csv" % S
if os.path.exists(plan):
    USED.append(plan)
    HOW = {r["piece"]: r.get("how", "") or "" for r in csv.DictReader(l for l in open(plan) if not l.startswith("#"))}
NOHOM = sorted(p_ for p_, h in HOW.items() if "homolog none" in h and p_ in PIECES)

# ============================================================ the CARP annotation, in 10 kb bins (cached)
REP = os.path.join(A.refs, "%s_repeats" % SPECIES)
GFF = os.path.join(REP, "Repeat_Annotation_Unified.gff3")
if not os.path.exists(GFF):
    raise SystemExit("no CARP annotation published for %s yet (%s)" % (SPECIES, GFF))
USED.append(GFF)
seqf = os.path.join(REP, "sequences.tsv")
if os.path.exists(seqf):
    LA = {r["seqid"]: int(r["length"]) for r in csv.DictReader(open(seqf), delimiter="\t")}
    bad = [c for c in COPIES if LA.get(c) != L2n[c]]
    if bad:
        raise SystemExit("the CARP annotation was made on a different assembly (lengths differ for %s)" %
                         ", ".join(bad[:4]))
    USED.append(seqf)
carp_version = ""
if os.path.exists(os.path.join(REP, "MANIFEST.txt")):
    m = re.search(r"^carp_version:\s*(\S+)", open(os.path.join(REP, "MANIFEST.txt")).read(), re.M)
    carp_version = m.group(1) if m else ""
MONO = {}
for tab, pref in (("trc_table_default.tsv", "TRC_"), ("trc_table_short.tsv", "TRC_S_")):
    p_ = os.path.join(REP, tab)
    if os.path.exists(p_):
        for r in csv.DictReader(open(p_), delimiter="\t"):
            tid = r.get("TRC_ID", "")
            name = pref + tid.split("TRC_", 1)[1] if tid.startswith("TRC_") else tid
            for col in ("monomer_tarean", "monomer_kite", "prevalent_founder"):
                v = (r.get(col) or "").strip()
                if v and v.replace(".", "", 1).isdigit() and float(v) > 0:
                    MONO[name] = int(float(v))
                    break

NF = {c: -(-L2n[c] // FINE) for c in COPIES}
OFF, o_ = {}, 0
for c in COPIES:
    OFF[c] = o_
    o_ += NF[c]
NFT = o_
BINLEN = np.full(NFT, float(FINE))
for c in COPIES:
    BINLEN[OFF[c] + NF[c] - 1] = L2n[c] - (NF[c] - 1) * FINE
os.makedirs("results/repeat_recombination", exist_ok=True)        # results/ is not tracked by git
CACHE = os.path.join("results/repeat_recombination", "carp_bins_%s.npz" % SPECIES)
KEY = "v1|%s|%d|%d|%d|%s" % (os.path.realpath(GFF), os.path.getsize(GFF), int(os.path.getmtime(GFF)), FINE,
                             ",".join("%s:%d" % (c, L2n[c]) for c in COPIES))


def read_carp():
    cs, ss, es = array.array("h"), array.array("q"), array.array("q")
    ks, ls, fs = array.array("b"), array.array("h"), array.array("i")
    lin_idx, fam_idx = {}, {}
    with open(GFF) as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 9:
                continue
            ci = CIDX.get(f[0])
            if ci is None:
                continue
            at = f[8]
            if "ID=UA_L2_" in at:               # members nested in a tandem LTR container: counted in the container
                continue
            k = at.find("classification=")
            cls = at[k + 15:].split(";", 1)[0].rstrip() if k >= 0 else "Unknown"
            if "%" in cls:
                cls = unquote(cls)
            kc = coarse(cls)
            li = lin_idx.setdefault(lineage(cls), len(lin_idx)) if kc <= 2 else -1
            fi = -1
            if kc in (ISAT, I45, I5):
                k2 = at.find("Name=")
                nm = at[k2 + 5:].split(";", 1)[0].rstrip() if k2 >= 0 else ""
                if nm.startswith("TRC_"):
                    fi = fam_idx.setdefault(nm, len(fam_idx))
            cs.append(ci)
            ss.append(int(f[3]) - 1)
            es.append(int(f[4]))
            ks.append(kc)
            ls.append(li)
            fs.append(fi)
    C = np.frombuffer(cs, np.int16).astype(np.int64)
    Sx = np.frombuffer(ss, np.int64).copy()
    Ex = np.frombuffer(es, np.int64).copy()
    K = np.frombuffer(ks, np.int8)
    LI = np.frombuffer(ls, np.int16)
    FI = np.frombuffer(fs, np.int32)
    Lc = np.array([L2n[c] for c in COPIES], np.int64)[C]
    Ex = np.minimum(Ex, Lc)
    Sx = np.maximum(Sx, 0)
    ok = Ex > Sx
    C, Sx, Ex, K, LI, FI = C[ok], Sx[ok], Ex[ok], K[ok], LI[ok], FI[ok]
    g0 = np.array([OFF[c] for c in COPIES], np.int64)[C] * FINE
    covk = np.zeros((NFT, NK), np.float32)
    for k in range(NK):
        m = K == k
        covk[:, k] = cover(Sx[m] + g0[m], Ex[m] + g0[m], NFT, FINE)
    lbp = np.bincount(LI[LI >= 0], weights=(Ex - Sx)[LI >= 0], minlength=len(lin_idx))
    tot = float(sum(L2n[c] for c in COPIES))
    names = {i: n for n, i in lin_idx.items()}
    keep = [i for i in np.argsort(-lbp) if lbp[i] >= 0.001 * tot]
    covl = np.zeros((NFT, len(keep)), np.float32)
    for j, i in enumerate(keep):
        m = LI == i
        covl[:, j] = cover(Sx[m] + g0[m], Ex[m] + g0[m], NFT, FINE)
    # tandem-repeat arrays: each family's pieces on a copy merged within MERGE; kind = the family's main class
    fams = {i: n for n, i in fam_idx.items()}
    arr = []
    o = np.argsort(FI, kind="mergesort")
    FIo = FI[o]
    starts = np.searchsorted(FIo, np.arange(len(fam_idx)))
    ends = np.searchsorted(FIo, np.arange(len(fam_idx)), side="right")
    for fi in range(len(fam_idx)):
        idx = o[starts[fi]:ends[fi]]
        kk = np.bincount(K[idx], minlength=NK)
        kind = max((ISAT, I5, I45), key=lambda k: kk[k])
        for ci in np.unique(C[idx]):
            j = idx[C[idx] == ci]
            for a, b, n in zip(*merge_iv(Sx[j], Ex[j], MERGE)):
                arr.append((ci, a, b, fi, kind, n))
    # rDNA loci: every 45S (5S) piece, TideCluster arrays and RepeatMasker subunits, merged within MERGE
    loci = []
    for k in (I45, I5):
        for ci in range(len(COPIES)):
            m = (K == k) & (C == ci)
            for a, b, n in zip(*merge_iv(Sx[m], Ex[m], MERGE)):
                loci.append((ci, a, b, k, n))
    np.savez_compressed(CACHE, key=np.array(KEY), covk=covk, covl=covl,
                        lin=np.array([names[i] for i in keep] or ["-"]), nlin=np.array(len(keep)),
                        fam=np.array([fams[i] for i in range(len(fams))] or ["-"]), nfam=np.array(len(fams)),
                        arr=np.array(arr, np.int64).reshape(-1, 6), loci=np.array(loci, np.int64).reshape(-1, 5),
                        nfeat=np.array(len(C)))
    return CACHE


if not (os.path.exists(CACHE) and str(np.load(CACHE)["key"]) == KEY):
    print("reading the CARP annotation (%s, %.0f MB); cached in %s" % (GFF, os.path.getsize(GFF) / 1e6, CACHE),
          flush=True)
    t_ = time.time()
    read_carp()
    print("  read in %.0f s" % (time.time() - t_), flush=True)
else:
    print("CARP annotation: from the cache %s" % CACHE, flush=True)
Z_ = np.load(CACHE)
COVK = Z_["covk"].astype(np.float64)
NLIN = int(Z_["nlin"])
COVL = Z_["covl"].astype(np.float64)[:, :NLIN]
LIN = [str(x) for x in Z_["lin"]][:NLIN]
FAMS = [str(x) for x in Z_["fam"]][:int(Z_["nfam"])]
ARR = pd.DataFrame(Z_["arr"], columns=["ci", "s", "e", "fi", "kind", "pieces"])
ARR["chrom"] = [COPIES[i] for i in ARR.ci]
ARR["family"] = [FAMS[i] for i in ARR.fi]
LOCI = pd.DataFrame(Z_["loci"], columns=["ci", "s", "e", "kind", "pieces"])
LOCI["chrom"] = [COPIES[i] for i in LOCI.ci]
LOCI["type"] = ["45S" if k == I45 else "5S" for k in LOCI.kind]
NFEAT = int(Z_["nfeat"])
del Z_
CSK = np.vstack([np.zeros((1, NK)), np.cumsum(COVK, axis=0)])        # prefix sums for any stretch
CSL = np.vstack([np.zeros((1, NLIN)), np.cumsum(COVL, axis=0)]) if NLIN else np.zeros((NFT + 1, 0))
CSB = np.r_[0.0, np.cumsum(BINLEN)]


def famlabel(f):
    return "%s%s" % (f, " (%d bp)" % MONO[f] if f in MONO else "")


def comp(c, a, b):
    """Share of [a, b) of copy c in each repeat class (10 kb resolution)."""
    i0 = OFF[c] + max(0, int(a // FINE))
    i1 = OFF[c] + min(NF[c], int(-(-b // FINE)))
    d = CSB[i1] - CSB[i0]
    return (CSK[i1] - CSK[i0]) / d if d > 0 else np.full(NK, np.nan)


print("%s (%s): %s features of the CARP annotation%s on %d chromosome copies; crossover reference: %s" % (
    S, SP, format(NFEAT, ","), " (CARP %s)" % carp_version if carp_version else "", len(COPIES),
    ", ".join("%s (%s)" % (c, " + ".join(p_.replace(c + "_", "") for _, _, p_ in PSPAN[c])) if len(PSPAN[c]) > 1
              else c for c in REFC)), flush=True)

# ============================================================ crossovers, markers, genes, structure
cof = "results/crossovers/%s/co_intervals.bed" % S
perf = "results/crossovers/%s/co_per_cell.tsv" % S
for p_ in (cof, perf):
    if not os.path.exists(p_):
        raise SystemExit("no crossover calls for %s: %s" % (S, p_))
USED += [cof, perf]
co = pd.read_csv(cof, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["seq", "start", "end", "barcode"],
                 dtype={0: str, 3: str})
per = pd.read_csv(perf, sep="\t", dtype={"barcode": str})
NUC = list(per.barcode)
NN = len(NUC)
nuc_i = {b: i for i, b in enumerate(NUC)}
n0 = len(co)
co = co[co.seq.isin(PIECES) & co.barcode.isin(nuc_i)].copy()
if len(co) < n0:
    note("%d of %d crossover calls are on sequences outside the reference pieces or in nuclei not in co_per_cell.tsv;"
         " left out" % (n0 - len(co), n0))
co["c"] = [PIECES[p_][0] for p_ in co.seq]
off = np.array([PIECES[p_][1] for p_ in co.seq], np.int64)
pend = np.array([PIECES[p_][2] for p_ in co.seq], np.int64)
co["s"] = np.minimum(np.maximum(co.start.values - 1, 0) + off, pend - 1)
co["e"] = np.minimum(co.end.values + off, pend)
co.loc[co.e <= co.s, "e"] = co.s + 1
co["mid"] = (co.s + co.e) / 2.0
co["nuc"] = co.barcode.map(nuc_i).astype(int)
co["g0"] = co.c.map(OFF).values * FINE
NCELL = {p_: NN for p_ in PIECES}
cpc = "results/crossovers/%s/cells_per_chrom.tsv" % S
if os.path.exists(cpc):
    USED.append(cpc)
    for r in csv.DictReader(open(cpc), delimiter="\t"):
        if r["chrom"] in NCELL:
            NCELL[r["chrom"]] = int(r["n_cells"])
WIDTH_MB = float(np.median(co.e - co.s)) / MB if len(co) else float("nan")
print("crossovers: %d calls in %d nuclei (%.2f per nucleus); median interval %.2f Mb" % (
    len(co), NN, len(co) / max(NN, 1), WIDTH_MB), flush=True)


def co_track(weights=None):
    """Crossover mass per 10 kb bin: each call spread evenly over its interval."""
    w = 1.0 / (co.e.values - co.s.values)
    if weights is not None:
        w = w * weights[co.nuc.values]
    return cover(co.s.values + co.g0.values, co.e.values + co.g0.values, NFT, FINE, w)


COF = co_track()

MK = np.zeros(NFT)
mkf = "qc/markers/%s/marker_classes.tsv.gz" % S
HAVE_MK = os.path.exists(mkf)
if HAVE_MK:
    USED.append(mkf)
    mk = pd.read_csv(mkf, sep="\t", usecols=["chrom", "pos", "class"], dtype={"chrom": str})
    mk = mk[(mk["class"] == "good") & mk.chrom.isin(PIECES)]
    g = np.array([OFF[PIECES[p_][0]] * FINE + PIECES[p_][1] for p_ in mk.chrom], np.int64) + mk.pos.values - 1
    MK = np.bincount(g // FINE, minlength=NFT).astype(float)[:NFT]
    del mk
else:
    note("no marker classes (%s): no detectability control" % mkf)

GN = np.zeros(NFT)
gff = row.get("annotation_gff", "")
HAVE_GN = bool(gff) and os.path.exists(gff)
if HAVE_GN:
    USED.append(gff)
    gg = pd.read_csv(gff, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4], names=["seq", "type", "s", "e"],
                     dtype={0: str}, low_memory=False)
    gg = gg[(gg["type"] == "gene") & gg.seq.isin(PIECES)]
    g = np.array([OFF[PIECES[p_][0]] * FINE + PIECES[p_][1] for p_ in gg.seq], np.int64) + \
        ((gg.s.values + gg.e.values) // 2).astype(np.int64)
    GN = np.bincount(g // FINE, minlength=NFT).astype(float)[:NFT]
    HAVE_GN = len(gg) > 0
    del gg
else:
    note("no gene annotation for %s (%s): genes left out" % (S, gff or "none in samples.csv"))

# structure: where the other haplotype's copy aligns (hap2-on-hap1 synteny blocks), per 10 kb bin
STRUCT = {k: np.full(NFT, np.nan) for k in ("collinear", "inverted", "translocated", "aligned")}
BLOCKS = {}
HAVE_ST = False
sid = next((k for k, r in rows.items() if r["species"] == SPECIES and r["haplotype"] == "hap1"), None)
stf = "qc/translocations/%s/%s_synteny_blocks.tsv" % (sid, sid) if sid else ""
if stf and os.path.exists(stf):
    USED.append(stf)
    HAVE_ST = True
    bt = pd.read_csv(stf, sep="\t")
    bt["kind"] = bt["kind"].fillna("").astype(str).replace("", "collinear").str.lower()
    for c in COPIES:
        x = bt[bt.hap1 == c] if hap(c) == "1" else bt[bt.hap2 == c]
        s0, e0 = (x.t_start.values, x.t_end.values) if hap(c) == "1" else (x.q_start.values, x.q_end.values)
        BLOCKS[c] = list(zip(s0.astype(np.int64), e0.astype(np.int64), x["kind"].values,
                             (x.hap2 if hap(c) == "1" else x.hap1).values))
        sl = slice(OFF[c], OFF[c] + NF[c])
        for k in ("collinear", "inverted", "translocated", "aligned"):
            sel = [(a, b) for a, b, kk, _ in BLOCKS[c] if kk == k or k == "aligned"]
            if sel:
                ms, me, _ = merge_iv([a for a, _ in sel], [b for _, b in sel], 0)
                STRUCT[k][sl] = cover(ms, me, NF[c], FINE) / BINLEN[sl]
            else:
                STRUCT[k][sl] = 0.0
    for p_ in NOHOM:                       # the homolog is in the same haplotype assembly: structure unknown here
        c, a, b = PIECES[p_]
        for k in STRUCT:
            STRUCT[k][OFF[c] + a // FINE:OFF[c] + -(-b // FINE)] = np.nan
else:
    note("no synteny blocks between the haplotypes (%s): structure left out" % (stf or "no hap1 sample"))

# ============================================================ joins and the coupling (cosegregation_2n.py)
JOINS, COUPLES = [], []
kt = os.path.join(A.coseg, "karyotype_test_%s.tsv" % S)
if os.path.exists(kt):
    USED.append(kt)
    for r in csv.DictReader(open(kt), delimiter="\t"):
        if r["kind"] == "join":
            m = re.match(r"^(chr\d+_hap[12])(?: \(([A-D])\))?: (\S+) \| (\S+) at ([\d.]+) Mb", r["what"])
            if m and m.group(1) in CIDX:
                JOINS.append(dict(copy=m.group(1), letter=m.group(2) or "", left=m.group(3), right=m.group(4),
                                  pos=float(m.group(5)) * MB, phi=float(r["phi"]),
                                  lo=float(r["ci_lo"] or "nan"), hi=float(r["ci_hi"] or "nan")))
        elif r["kind"] == "couple":
            m = re.match(r"^(\S+) ([\d.]+) Mb x (\S+) ([\d.]+) Mb", r["what"])
            if m:
                COUPLES.append(dict(p=m.group(1), pp=float(m.group(2)) * MB, q=m.group(3), qp=float(m.group(4)) * MB))
bpf = "qc/reference/%s/breakpoints.csv" % S
if JOINS and os.path.exists(bpf):
    USED.append(bpf)
    for r in csv.DictReader(open(bpf)):
        for j in JOINS:
            if r["chrom"] == j["copy"] and abs(int(r["cut"]) - j["pos"]) < 2 * MB:
                j["pos"] = float(int(r["cut"]))
for j in JOINS:
    j["name"] = "%s%s %s|%s" % (j["letter"] + " " if j["letter"] else "", j["copy"], j["left"], j["right"])
    j["joined"] = j["phi"] >= 0.5
lg = "logs/cosegregation_2n_%s.log" % S
SPANS = {}                                       # reference copy -> (a, b, point) of the coupling interval
if COUPLES:
    ci_ = {}
    if os.path.exists(lg):
        USED.append(lg)
        for l in open(lg):
            m = re.search(r"COUPLING (\S+) x (\S+): point at ([\d.]+) Mb \(95% ([\d.]+)-([\d.]+)\) x ([\d.]+) Mb "
                          r"\(95% ([\d.]+)-([\d.]+)\)", l)
            if m:
                ci_[(m.group(1), m.group(2))] = ((float(m.group(4)) * MB, float(m.group(5)) * MB),
                                                 (float(m.group(7)) * MB, float(m.group(8)) * MB))
    for cp in COUPLES:
        cis = ci_.get((cp["p"], cp["q"]))
        for k, ix in (("p", 0), ("q", 1)):
            src, o2, _ = PIECES.get(cp[k], (cp[k], 0, 0))
            if src in REFC:
                pt = o2 + cp[k + "p"]
                a, b = (o2 + cis[ix][0], o2 + cis[ix][1]) if cis else (pt - 10 * MB, pt + 10 * MB)
                SPANS[src] = (a, b, pt)
        if not cis:
            note("no 95%% interval for the coupling %s x %s in %s: +-10 Mb around the point" % (cp["p"], cp["q"], lg))

# ============================================================ windows
if A.win_mb > 0:
    W = int(A.win_mb * MB)
else:
    med = float(np.median([L2n[c] for c in REFC])) / MB
    W = int(min(5.0, max(1.0, round(med / 50 * 2) / 2)) * MB)
rws = []
for c in REFC:
    L = L2n[c]
    n = max(1, int(round(L / W)))
    eb = np.unique(np.round(np.linspace(0, NF[c], n + 1)).astype(int))
    for k in range(len(eb) - 1):
        a, b = int(eb[k]) * FINE, min(int(eb[k + 1]) * FINE, L)
        mid = (a + b) / 2.0
        p_ = piece_at(c, mid)
        rws.append(dict(chrom=c, k=k, start=a, end=b, g0=OFF[c] + int(eb[k]), g1=OFF[c] + int(eb[k + 1]),
                        len_mb=(b - a) / MB, piece=p_, n_cells=NCELL[p_],
                        at_cut=any(a - W / 2 < x < b + W / 2 for x in CUTS[c]),
                        end_dist_mb=min(mid, L - mid) / MB, end_rel=min(mid, L - mid) / L,
                        coupling=bool(c in SPANS and min(b, SPANS[c][1]) - max(a, SPANS[c][0]) >= 0.5 * (b - a))))
WT = pd.DataFrame(rws)
G0, G1 = WT.g0.values, WT.g1.values
NW = len(WT)
CS_CO = np.r_[0.0, np.cumsum(COF)]
WT["co_mass"] = CS_CO[G1] - CS_CO[G0]
WT["co_mid"] = 0
for i, c in enumerate(REFC):
    w_ = WT.index[WT.chrom == c]
    edges = np.r_[WT.start.values[w_], WT.end.values[w_][-1]]
    mids = co.mid.values[co.c.values == c]
    WT.loc[w_, "co_mid"] = np.histogram(mids, bins=edges)[0]
WT["cM_per_Mb"] = 100.0 * WT.co_mass / WT.n_cells.clip(lower=1) / WT.len_mb
boot = np.zeros((A.nboot, NW))
for bi in range(A.nboot):
    wts = np.bincount(RNG.integers(0, NN, NN), minlength=NN).astype(float)
    cb = np.r_[0.0, np.cumsum(co_track(wts))]
    boot[bi] = 100.0 * (cb[G1] - cb[G0]) / WT.n_cells.clip(lower=1).values / WT.len_mb.values
WT["cM_lo"] = np.percentile(boot, 2.5, axis=0)
WT["cM_hi"] = np.percentile(boot, 97.5, axis=0)
del boot
BPW = CSB[G1] - CSB[G0]
for k in range(NK):
    WT[CLASSES[k]] = (CSK[G1, k] - CSK[G0, k]) / BPW
WT["all repeats"] = WT[CLASSES].sum(axis=1)
for j, ln in enumerate(LIN):
    WT[ln] = (CSL[G1, j] - CSL[G0, j]) / BPW
CS_MK, CS_GN = np.r_[0.0, np.cumsum(MK)], np.r_[0.0, np.cumsum(GN)]
WT["markers_per_mb"] = (CS_MK[G1] - CS_MK[G0]) / WT.len_mb if HAVE_MK else np.nan
WT["genes_per_mb"] = (CS_GN[G1] - CS_GN[G0]) / WT.len_mb if HAVE_GN else np.nan
for k in ("collinear", "inverted", "translocated", "aligned"):
    v = STRUCT[k]
    known = np.r_[0.0, np.cumsum(np.isfinite(v) * BINLEN)]
    vs = np.r_[0.0, np.cumsum(np.nan_to_num(v) * BINLEN)]
    kb = known[G1] - known[G0]
    WT["struct_" + k] = np.where(kb >= 0.5 * BPW, (vs[G1] - vs[G0]) / np.maximum(kb, 1), np.nan)
WT["struct_none"] = 1 - WT.struct_aligned
AN = (~WT.at_cut.values) & (WT.n_cells.values >= 10)
REFBP = float(sum(L2n[c] for c in REFC))
CPCOPY = {c: np.flatnonzero(WT.chrom.values == c) for c in REFC}
print("windows: %d of %.1f Mb on %d chromosomes; %d tested (%d within %.1f Mb of a cut of the reference left out)" % (
    NW, W / MB, len(REFC), int(AN.sum()), int(WT.at_cut.sum()), W / 2 / MB), flush=True)

# rotations: per chromosome, every window moves by the same random number of windows (circular)
PERM = np.empty((A.nperm, NW), np.int64)
for c, ix in CPCOPY.items():
    n = len(ix)
    sh = RNG.integers(1, n, A.nperm) if n > 1 else np.zeros(A.nperm, np.int64)
    PERM[:, ix] = ix[0] + (np.arange(n)[None, :] + sh[:, None]) % n

# ============================================================ Q1a: one variable at a time
VARS = [("distance to the nearest end (Mb)", "position", WT.end_dist_mb.values)]
if HAVE_MK:
    VARS.append(("informative markers per Mb", "markers", WT.markers_per_mb.values))
if HAVE_GN:
    VARS.append(("genes per Mb", "genes", WT.genes_per_mb.values))
if HAVE_ST:
    for k, lab in (("collinear", "collinear"), ("inverted", "inverted"), ("translocated", "translocated"),
                   ("none", "not aligned")):
        VARS.append(("other haplotype: %s" % lab, "structure", WT["struct_" + k].values))
COUP = WT.coupling.values.astype(bool)
COUPLING_VAR = "in the coupling interval (pollen)"
if COUP.any():
    VARS.append((COUPLING_VAR, "structure", COUP.astype(float)))
for k in range(NK):
    if np.nanstd(WT[CLASSES[k]].values[AN]) > 0:
        VARS.append((CLASSES[k], "repeat class", WT[CLASSES[k]].values))
VARS.append(("all repeats", "repeat class", WT["all repeats"].values))
LIN_T = [ln for j, ln in enumerate(LIN) if CSL.shape[1] and
         sum(CSL[OFF[c] + NF[c], j] - CSL[OFF[c], j] for c in REFC) >= A.min_lineage * REFBP]
for ln in LIN_T:
    VARS.append((ln, "LTR lineage", WT[ln].values))
MT = len(VARS)
ALPHA = 0.05 / MT
y_all = WT.cM_per_Mb.values


def ctrl_cols(mask):
    cols = [np.ones(mask.sum()), rankavg(WT.end_dist_mb.values[mask])]
    if HAVE_MK:
        cols.append(rankavg(WT.markers_per_mb.values[mask]))
    Q, _ = np.linalg.qr(np.column_stack(cols))
    return Q


def assoc_core(x, base, is_ctrl):
    """Spearman rho (and partial rho) of x with crossover density over the windows in base, and their rotation
    nulls: x rotated along its chromosome NPERM times."""
    m0 = base & np.isfinite(x)
    if m0.sum() < 10 or np.nanstd(x[m0]) == 0:
        return None
    ry = rankavg(y_all[m0])
    rx = rankavg(x[m0])
    o = dict(n=int(m0.sum()), rho=corr(rx, ry))
    Q = ctrl_cols(m0)
    o["rho_partial"] = np.nan if is_ctrl else corr(rx - Q @ (Q.T @ rx), ry - Q @ (Q.T @ ry))
    XR = x[PERM]
    if np.isfinite(x).all():
        R = np.vstack([rankavg(v) for v in XR[:, m0]])
        Rc = R - R.mean(1, keepdims=True)
        ryc = ry - ry.mean()
        nr = np.sqrt((Rc * Rc).sum(1)) * math.sqrt((ryc * ryc).sum())
        rn = (Rc @ ryc) / np.where(nr > 0, nr, np.nan)
        pn = np.full(A.nperm, np.nan)
        if not is_ctrl:
            Rr = R - (R @ Q) @ Q.T
            yr = ry - Q @ (Q.T @ ry)
            Rr -= Rr.mean(1, keepdims=True)
            yr = yr - yr.mean()
            nr2 = np.sqrt((Rr * Rr).sum(1)) * math.sqrt((yr * yr).sum())
            pn = (Rr @ yr) / np.where(nr2 > 0, nr2, np.nan)
    else:
        rn, pn = np.full(A.nperm, np.nan), np.full(A.nperm, np.nan)
        for i in range(A.nperm):
            xi = XR[i]
            mi = base & np.isfinite(xi)
            if mi.sum() < 10:
                continue
            a_, b_ = rankavg(xi[mi]), rankavg(y_all[mi])
            rn[i] = corr(a_, b_)
            if not is_ctrl:
                Qi = ctrl_cols(mi)
                pn[i] = corr(a_ - Qi @ (Qi.T @ a_), b_ - Qi @ (Qi.T @ b_))
    rn = rn[np.isfinite(rn)]
    o["p_rot"] = (1 + np.sum(np.abs(rn) >= abs(o["rho"]) - 1e-12)) / (1 + len(rn))
    o["null_lo"], o["null_hi"] = np.percentile(rn, [2.5, 97.5])
    o["p_partial"], o["pnull_lo"], o["pnull_hi"] = np.nan, np.nan, np.nan
    pn = pn[np.isfinite(pn)]
    if not is_ctrl and len(pn):
        o["p_partial"] = (1 + np.sum(np.abs(pn) >= abs(o["rho_partial"]) - 1e-12)) / (1 + len(pn))
        o["pnull_lo"], o["pnull_hi"] = np.percentile(pn, [2.5, 97.5])
    return o


def assoc_one(name, group, x):
    is_ctrl = group in ("position", "markers")
    base = AN & np.isfinite(y_all)
    out = dict(variable=name, group=group)
    o = assoc_core(x, base, is_ctrl)
    if o is None:
        out.update(n_windows=int((base & np.isfinite(x)).sum()), rho=np.nan, p_rot=np.nan, rho_partial=np.nan,
                   p_partial=np.nan)
        return out
    out["n_windows"] = o.pop("n")
    out.update(o)
    if COUP.any() and name != COUPLING_VAR:
        o2 = assoc_core(x, base & ~COUP, is_ctrl)        # the same, with the coupled middles left out
        if o2 is not None:
            out["outside_coupling_rho"] = o2["rho"] if is_ctrl else o2["rho_partial"]
            out["outside_coupling_p"] = o2["p_rot"] if is_ctrl else o2["p_partial"]
    m0 = base & np.isfinite(x)
    pc, agree, nchr = [], 0, 0
    for c, ix in CPCOPY.items():
        mm = m0[ix]
        if mm.sum() >= 8 and np.std(x[ix][mm]) > 0:
            r_ = corr(rankavg(x[ix][mm]), rankavg(y_all[ix][mm]))
            pc.append((c, r_))
            nchr += 1
            agree += int(np.sign(r_) == np.sign(out["rho"]))
    out["per_chromosome"] = "; ".join("%s %+.2f" % (c, r_) for c, r_ in pc)
    out["same_sign"] = "%d of %d" % (agree, nchr)
    out["_pc"] = pc
    return out


def run_assoc():
    res = [assoc_one(n_, g_, x_) for n_, g_, x_ in VARS]
    T = pd.DataFrame(res)
    T["p_used"] = np.where(T.group.isin(["position", "markers"]), T.p_rot, T.p_partial)
    T["beyond_chance"] = T.p_used < ALPHA
    return T


TA = section("crossovers x one variable at a time", run_assoc)


# ============================================================ Q1b: how much each set explains (Poisson deviance)
def irls(X, y, offs):
    b = np.zeros(X.shape[1])
    b[0] = math.log(max(y.sum(), 1e-9) / np.exp(offs).sum())
    for _ in range(100):
        eta = np.clip(X @ b + offs, -30, 30)
        mu = np.exp(eta)
        z = eta - offs + (y - mu) / mu
        XtW = X.T * mu
        H = XtW @ X
        H[np.diag_indices_from(H)] += 1e-8
        bn = np.linalg.solve(H, XtW @ z)
        done = np.max(np.abs(bn - b)) < 1e-7
        b = bn
        if done:
            break
    mu = np.exp(np.clip(X @ b + offs, -30, 30))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(y > 0, y * np.log(np.where(y > 0, y, 1) / mu), 0.0)
    return float(2 * np.sum(t - (y - mu)))


def zcols(cols):
    out = []
    for v in cols:
        v = np.asarray(v, float)
        sd = v.std()
        if sd > 0 and np.isfinite(sd):
            out.append((v - v.mean()) / sd)
    return np.column_stack(out) if out else np.zeros((len(cols[0]) if cols else 0, 0))


def run_glm():
    m = AN & np.isfinite(WT.co_mass.values)
    y = WT.co_mass.values[m]
    offs = np.log(WT.n_cells.values[m] * WT.len_mb.values[m])
    u = WT.end_rel.values[m]
    blocks = {"position": zcols([u, u ** 2, np.log1p(WT.end_dist_mb.values[m])])}
    if HAVE_MK:
        blocks["markers"] = zcols([rankavg(WT.markers_per_mb.values[m])])
    if HAVE_GN:
        blocks["genes"] = zcols([rankavg(WT.genes_per_mb.values[m])])
    cols = []
    if HAVE_ST:
        for k in ("collinear", "inverted", "translocated"):
            v = WT["struct_" + k].values[m].copy()
            v[~np.isfinite(v)] = np.nanmean(v) if np.isfinite(v).any() else 0
            cols.append(v)
        cols.append((~np.isfinite(WT.struct_collinear.values[m])).astype(float))
    if COUP.any():
        cols.append(COUP[m].astype(float))
    if cols:
        blocks["structure"] = zcols(cols)
    rep_cols = [k for k in range(NK) if (WT[CLASSES[k]].values[m] > 0).sum() >= 3]
    REPX = WT[[CLASSES[k] for k in rep_cols]].values

    def repblock(perm=None):
        X_ = REPX if perm is None else REPX[perm]
        return zcols([rankavg(X_[m, j]) for j in range(X_.shape[1])])

    one = np.ones((m.sum(), 1))
    d0 = irls(one, y, offs)

    def r2(*bl):
        X = np.hstack([one] + [b_ for b_ in bl if b_.shape[1]])
        return 1 - irls(X, y, offs) / d0

    R = {}
    for k, b_ in blocks.items():
        R[k] = r2(b_)
    rb = repblock()
    R["repeats"] = r2(rb)
    pm = [blocks["position"]] + ([blocks["markers"]] if "markers" in blocks else [])
    rest = pm + [blocks[k] for k in ("genes", "structure") if k in blocks]
    R["position + markers"] = r2(*pm)
    R["position + markers + repeats"] = r2(*(pm + [rb]))
    R["everything but the repeats"] = r2(*rest)
    R["everything"] = r2(*(rest + [rb]))
    nn = min(A.nperm, 1000)
    nul = np.zeros((nn, 3))
    for i in range(nn):
        rr = repblock(PERM[i])
        nul[i] = (r2(rr), r2(*(pm + [rr])) - R["position + markers"], r2(*(rest + [rr])) - R["everything but the repeats"])
    obs = np.array([R["repeats"], R["position + markers + repeats"] - R["position + markers"],
                    R["everything"] - R["everything but the repeats"]])
    p = (1 + (nul >= obs[None, :] - 1e-12).sum(0)) / (1 + nn)
    return dict(R=R, n=int(m.sum()), d0=d0, classes=[CLASSES[k] for k in rep_cols], nul=nul, obs=obs, p=p, nn=nn,
                labels=["repeats alone", "repeats after position + markers", "repeats after everything else"])


GLM = section("crossovers x all variables together (Poisson)", run_glm)

# ============================================================ Q1c: crossover density at the big tandem arrays
REFSET = set(REFC)
MASKF = np.zeros(NFT, bool)                     # bins near a cut of the crossover reference: calls not possible there
for c in REFC:
    for x in CUTS[c]:
        MASKF[OFF[c] + max(0, int((x - W / 2) // FINE)):OFF[c] + min(NF[c], int((x + W / 2) // FINE) + 1)] = True


def array_sets():
    sets = []
    lo = LOCI[LOCI.chrom.isin(REFSET)]
    for t in ("45S", "5S"):
        x = lo[(lo.type == t) & (lo.e - lo.s >= RDNA_MIN)]
        sets.append(("%s rDNA arrays >= %d kb" % (t, RDNA_MIN // 1000), list(zip(x.chrom, x.s, x.e))))
    sa = ARR[ARR.chrom.isin(REFSET) & (ARR.kind == ISAT)]
    big = sa[sa.e - sa.s >= SAT_BIG]
    sets.append(("satellite arrays >= %d kb" % (SAT_BIG // 1000), list(zip(big.chrom, big.s, big.e))))
    mid = sa[sa.e - sa.s >= SAT_TICK]
    tops = mid.assign(l=mid.e - mid.s).groupby("family").agg(n=("l", "size"), bp=("l", "sum"))
    tops = tops[tops.n >= 3].sort_values("bp", ascending=False).head(4)
    for f, r in tops.iterrows():
        x = mid[mid.family == f]
        sets.append(("%s arrays >= %d kb" % (famlabel(f), SAT_TICK // 1000), list(zip(x.chrom, x.s, x.e))))
    return sets


def run_arrays():
    D = int(A.array_mb * MB // FINE)
    nn = min(A.nperm, 2000)
    out = []
    massm = np.where(MASKF, 0.0, COF)
    for name, arrs in array_sets():
        if len(arrs) < 2:
            out.append(dict(set=name, n=len(arrs), note="fewer than 2 arrays"))
            continue
        O = np.zeros(nn + 1)
        E = np.zeros(nn + 1)
        zone_mb = 0.0
        for c in REFC:
            a_c = [(s, e) for cc, s, e in arrs if cc == c]
            if not a_c:
                continue
            n = NF[c]
            sl = slice(OFF[c], OFF[c] + n)
            mass, unm = massm[sl], (~MASKF[sl]).astype(float)
            dens = mass.sum() / max(unm.sum(), 1)
            iv = sorted((max(0, s // FINE - D), min(n, -(-e // FINE) + D)) for s, e in a_c)
            mi = []
            for a, b in iv:
                if mi and a <= mi[-1][1]:
                    mi[-1][1] = max(mi[-1][1], b)
                else:
                    mi.append([a, b])
            a_ = np.array([v[0] for v in mi])
            b_ = np.array([v[1] for v in mi])
            zone_mb += float((b_ - a_).sum()) * FINE / MB
            Pm = np.r_[0.0, np.cumsum(np.r_[mass, mass])]
            Pu = np.r_[0.0, np.cumsum(np.r_[unm, unm])]
            ks = np.r_[0, RNG.integers(1, n, nn)]
            ia, ib = a_[None, :] + ks[:, None], b_[None, :] + ks[:, None]
            O += (Pm[ib] - Pm[ia]).sum(1)
            E += dens * (Pu[ib] - Pu[ia]).sum(1)
        ratio = np.where(E > 0, O / np.where(E > 0, E, 1), np.nan)
        r0, rn = ratio[0], ratio[1:][np.isfinite(ratio[1:])]
        lr0, lrn = np.log(max(r0, 1e-9)), np.log(np.maximum(rn, 1e-9))
        p = (1 + np.sum(np.abs(lrn - np.median(lrn)) >= abs(lr0 - np.median(lrn)) - 1e-12)) / (1 + len(lrn))
        out.append(dict(set=name, n=len(arrs), zone_mb=zone_mb, obs_cM=100 * O[0] / NN, exp_cM=100 * E[0] / NN,
                        ratio=r0, null_lo=np.percentile(rn, 2.5), null_hi=np.percentile(rn, 97.5),
                        null_median=np.median(rn), p=p, nn=len(rn)))
    return pd.DataFrame(out)


TARR = section("crossovers at the tandem arrays", run_arrays)


# ============================================================ Q1d: dead spots (as supervisor_report.py)
def run_dead(w=5 * MB):
    rng = np.random.default_rng(A.seed + 7)

    def runs(cnt, L):
        out, i = [], 0
        while i < len(cnt):
            if cnt[i] == 0:
                j = i
                while j + 1 < len(cnt) and cnt[j + 1] == 0:
                    j += 1
                out.append((i * w, min((j + 1) * w, L)))
                i = j + 1
            else:
                i += 1
        return out
    sim = []
    for _ in range(200):
        lens = []
        for c in REFC:
            L = L2n[c]
            nb = int(math.ceil(L / w))
            pos = rng.uniform(0, L, int((co.c == c).sum()))
            lens += [b - a for a, b in runs(np.bincount(np.minimum((pos // w).astype(int), nb - 1), minlength=nb), L)]
        sim.append(np.array(lens))
    thr = next((t for t in range(20, 205, 5) if np.mean([(x >= t * MB).sum() for x in sim]) < 1.0), 200)
    chance = np.array([(x >= thr * MB).sum() for x in sim])
    rows_ = []
    for c in REFC:
        L = L2n[c]
        nb = int(math.ceil(L / w))
        cnt = np.bincount(np.minimum((co.mid.values[co.c.values == c] // w).astype(int), nb - 1), minlength=nb)
        whole = comp(c, 0, L)
        for a, b in runs(cnt, L):
            if b - a < thr * MB:
                continue
            cc = comp(c, a, b)
            sat = ARR[(ARR.chrom == c) & (ARR.kind == ISAT) & (ARR.e > a) & (ARR.s < b)]
            rl = LOCI[(LOCI.chrom == c) & (LOCI.e > a) & (LOCI.s < b) & (LOCI.e - LOCI.s >= RDNA_MIN)]
            gi0, gi1 = OFF[c] + int(a // FINE), OFF[c] + int(-(-b // FINE))
            mk_ = MK[gi0:gi1].sum() / ((b - a) / MB) / max(MK[OFF[c]:OFF[c] + NF[c]].sum() / (L / MB), 1e-9) \
                if HAVE_MK else np.nan
            st = {k: np.nanmean(STRUCT[k][gi0:gi1]) if HAVE_ST and np.isfinite(STRUCT[k][gi0:gi1]).any() else np.nan
                  for k in ("collinear", "inverted", "translocated", "aligned")}
            cpl = (min(b, SPANS[c][1]) - max(a, SPANS[c][0])) / (b - a) if c in SPANS else 0.0
            rows_.append(dict(chrom=c, start_mb=a / MB, end_mb=b / MB, mb=(b - a) / MB,
                              gypsy=cc[0], gypsy_chrom=whole[0], copia=cc[1], copia_chrom=whole[1],
                              repeats=np.nansum(cc), repeats_chrom=np.nansum(whole),
                              satellite_kb=float((np.minimum(sat.e, b) - np.maximum(sat.s, a)).sum()) / 1e3,
                              rdna_arrays=len(rl), markers_vs_chrom=mk_, collinear=st["collinear"],
                              inverted=st["inverted"], translocated=st["translocated"],
                              not_aligned=1 - st["aligned"] if np.isfinite(st["aligned"]) else np.nan,
                              coupling_share=max(0.0, cpl),
                              cut_inside=any(a < x < b for x in CUTS[c])))
    return dict(thr=thr, chance=(float(chance.mean()), float(np.percentile(chance, 2.5)),
                                 float(np.percentile(chance, 97.5))), rows=pd.DataFrame(rows_))


DEAD = section("dead spots", run_dead)
WT["dead"] = False
if DEAD is not None and len(DEAD["rows"]):
    for _, r in DEAD["rows"].iterrows():
        WT.loc[(WT.chrom == r.chrom) & (WT.end > r.start_mb * MB) & (WT.start < r.end_mb * MB), "dead"] = True

# ============================================================ Q2: breakpoints against the repeats
R_NEAR = int(A.near_mb * MB)


def sv_breakpoints():
    out = []
    if not HAVE_ST:
        return out
    for c, bl in BLOCKS.items():
        for a, b, kk, other in bl:
            if kk not in ("inverted", "translocated"):
                continue
            for x in (a, b):
                if min(x, L2n[c] - x) >= 2 * MB and not any(j["copy"] == c and abs(j["pos"] - x) < 2 * MB
                                                            for j in JOINS):
                    out.append(dict(copy=c, pos=float(x), kind=kk, other=other))
    out.sort(key=lambda r: (CIDX[r["copy"]], r["pos"]))
    merged = []
    for r in out:
        if merged and merged[-1]["copy"] == r["copy"] and r["pos"] - merged[-1]["pos"] < 2 * MB:
            if r["kind"] not in merged[-1]["kind"]:
                merged[-1]["kind"] += "+" + r["kind"]
            continue
        merged.append(dict(r))
    for r in merged:
        r["name"] = "%s %.1f Mb (%s%s)" % (r["copy"], r["pos"] / MB, r["kind"],
                                           ", other haplotype: " + r["other"] if r["kind"] == "translocated" else "")
    return merged


def near_tables(R):
    """Per copy: bins within R of a 45S / 5S rDNA array (>= 10 kb) or a satellite array (>= 20 kb)."""
    t = {}
    Rb = R // FINE
    for lab, df in (("45S", LOCI[(LOCI.type == "45S") & (LOCI.e - LOCI.s >= RDNA_MIN)]),
                    ("5S", LOCI[(LOCI.type == "5S") & (LOCI.e - LOCI.s >= RDNA_MIN)]),
                    ("satellite", ARR[(ARR.kind == ISAT) & (ARR.e - ARR.s >= SAT_TICK)])):
        z = np.zeros(NFT, bool)
        for c, s, e in zip(df.chrom, df.s, df.e):
            z[OFF[c] + max(0, s // FINE - Rb):OFF[c] + min(NF[c], -(-e // FINE) + Rb)] = True
        t[lab] = z
    t["rDNA"] = t["45S"] | t["5S"]
    return t


def nearest(c, x, df):
    if not len(df):
        return np.nan, ""
    d = np.maximum(0, np.maximum(df.s.values - x, x - df.e.values))
    i = int(np.argmin(d))
    return float(d[i]), df.iloc[i]


def run_breakpoints():
    SV = sv_breakpoints()
    BPS = [dict(set="join", chrom=j["copy"], pos=j["pos"], name=j["name"],
                verdict="the pollen inherit it" if j["joined"] else "the pollen do NOT inherit it",
                joined=j["joined"], phi=j["phi"]) for j in JOINS]
    BPS += [dict(set="haplotype breakpoint", chrom=r["copy"], pos=r["pos"], name=r["name"], verdict="", joined=None,
                 phi=np.nan) for r in SV]
    if not BPS:
        return None
    nt = {R: near_tables(R) for R in sorted({R_NEAR, int(0.5 * MB), int(2.5 * MB)})}
    l45 = LOCI[(LOCI.type == "45S") & (LOCI.e - LOCI.s >= RDNA_MIN)]
    l5 = LOCI[(LOCI.type == "5S") & (LOCI.e - LOCI.s >= RDNA_MIN)]
    sat = ARR[(ARR.kind == ISAT) & (ARR.e - ARR.s >= SAT_TICK)]
    rows_ = []
    for b in BPS:
        c, x = b["chrom"], b["pos"]
        cc = comp(c, x - R_NEAR, x + R_NEAR)
        r = dict(b)
        for k in range(NK):
            r[CLASSES[k]] = cc[k]
        for lab in ("45S", "5S"):                      # rDNA of any size (single units, fragments) as well
            d, hit = nearest(c, x, LOCI[(LOCI.type == lab) & (LOCI.chrom == c)])
            r["nearest_%s_any_mb" % lab] = d / MB if np.isfinite(d) else np.nan
            r["nearest_%s_any_kb" % lab] = (hit.e - hit.s) / 1e3 if np.isfinite(d) else np.nan
        for lab, df in (("45S", l45), ("5S", l5), ("satellite", sat)):
            dfc = df[df.chrom == c]
            d, hit = nearest(c, x, dfc)
            r["nearest_%s_mb" % lab] = d / MB if np.isfinite(d) else np.nan
            r["nearest_%s_side" % lab] = ("" if not np.isfinite(d) else "spans it" if d == 0 else
                                          "left" if hit.e <= x else "right")
            if lab == "satellite" and np.isfinite(d):
                r["nearest_satellite_family"] = famlabel(hit.family)
        inn = sat[(sat.chrom == c) & (sat.e > x - R_NEAR) & (sat.s < x + R_NEAR)]
        r["satellite_families_within"] = "; ".join("%s %.0f kb" % (famlabel(f), v / 1e3) for f, v in
                                                   (np.minimum(inn.e, x + R_NEAR) - np.maximum(inn.s, x - R_NEAR))
                                                   .groupby(inn.family).sum().sort_values(ascending=False).items())
        rows_.append(r)
    BT = pd.DataFrame(rows_)
    # chance: each breakpoint moved to random places on its own copy (away from the ends)
    nn = A.nperm
    sets = [("all joins", BT.set == "join"), ("joins the pollen inherit", (BT.set == "join") & (BT.joined == True)),
            ("joins the pollen do not inherit", (BT.set == "join") & (BT.joined == False)),
            ("haplotype breakpoints", BT.set == "haplotype breakpoint")]
    RP = {}
    for i, b in BT.iterrows():
        c = b.chrom
        lo_, hi_ = R_NEAR, max(R_NEAR + FINE, L2n[c] - R_NEAR)
        RP[i] = RNG.integers(lo_ // FINE, hi_ // FINE, nn)
    tests = []
    for sname, mask in sets:
        idx = list(BT.index[mask.values])
        if not idx:
            continue
        Rb = R_NEAR // FINE
        obs = np.nanmean(BT.loc[idx, CLASSES].values, axis=0)
        nul = np.zeros((nn, NK))
        for i in idx:
            c = BT.chrom[i]
            p0 = OFF[c] + RP[i]
            i0, i1 = np.maximum(p0 - Rb, OFF[c]), np.minimum(p0 + Rb, OFF[c] + NF[c])
            nul += (CSK[i1] - CSK[i0]) / (CSB[i1] - CSB[i0])[:, None]
        nul /= len(idx)
        for k in range(NK):
            mu = nul[:, k].mean()
            p = (1 + np.sum(np.abs(nul[:, k] - mu) >= abs(obs[k] - mu) - 1e-12)) / (1 + nn)
            tests.append(dict(test="breakpoints: share within %g Mb" % A.near_mb, set=sname, what=CLASSES[k],
                              n=len(idx), obs=obs[k], null_mean=mu, null_lo=np.percentile(nul[:, k], 2.5),
                              null_hi=np.percentile(nul[:, k], 97.5), ratio=obs[k] / mu if mu > 0 else np.nan, p=p))
        for R, tabs in nt.items():
            for lab, z in tabs.items():
                o_ = sum(bool(z[OFF[BT.chrom[i]] + int(BT.pos[i] // FINE)]) for i in idx)
                cnt = np.zeros(nn)
                for i in idx:
                    cnt += z[OFF[BT.chrom[i]] + RP[i]]
                p = (1 + np.sum(cnt >= o_)) / (1 + nn)
                tests.append(dict(test="breakpoints: with %s within %g Mb" % (lab, R / MB), set=sname, what=lab,
                                  n=len(idx), obs=o_, null_mean=cnt.mean(), null_lo=np.percentile(cnt, 2.5),
                                  null_hi=np.percentile(cnt, 97.5), ratio=o_ / cnt.mean() if cnt.mean() > 0 else np.nan,
                                  p=p, R_mb=R / MB))
    return dict(BT=BT, tests=pd.DataFrame(tests), nn=nn)


BPR = section("breakpoints x repeats", run_breakpoints)

# ============================================================ outputs: tables
WT.drop(columns=["g0", "g1"]).to_csv(os.path.join(A.out, "windows_%s.csv" % S), index=False, float_format="%.5g")
OUTS = [os.path.join(A.out, "windows_%s.csv" % S)]
if TA is not None:
    TA.drop(columns=["_pc"], errors="ignore").to_csv(os.path.join(A.out, "assoc_%s.csv" % S), index=False, float_format="%.4g")
    OUTS.append(os.path.join(A.out, "assoc_%s.csv" % S))
trows = []
if GLM is not None:
    for k, v in GLM["R"].items():
        trows.append(dict(test="deviance explained (Poisson, %d windows)" % GLM["n"], set="", what=k, obs=v))
    for lab, o_, nl, p in zip(GLM["labels"], GLM["obs"], GLM["nul"].T, GLM["p"]):
        trows.append(dict(test="deviance explained: rotation null", set="", what=lab, obs=o_,
                          null_mean=nl.mean(), null_lo=np.percentile(nl, 2.5), null_hi=np.percentile(nl, 97.5), p=p,
                          n=GLM["nn"]))
if TARR is not None:
    for _, r in TARR.iterrows():
        trows.append(dict(test="crossover density within %g Mb of the arrays, observed / expected" % A.array_mb,
                          set="", what=r["set"], n=r.get("n"), obs=r.get("ratio"), null_mean=r.get("null_median"),
                          null_lo=r.get("null_lo"), null_hi=r.get("null_hi"), p=r.get("p")))
if BPR is not None:
    trows += BPR["tests"].to_dict("records")
    BPR["BT"].to_csv(os.path.join(A.out, "breakpoints_%s.csv" % S), index=False, float_format="%.4g")
    OUTS.append(os.path.join(A.out, "breakpoints_%s.csv" % S))
pd.DataFrame(trows).to_csv(os.path.join(A.out, "tests_%s.csv" % S), index=False, float_format="%.4g")
OUTS.append(os.path.join(A.out, "tests_%s.csv" % S))

# ============================================================ figure 1: the tracks
SHORT = {}
for j in JOINS:
    SHORT[(j["copy"], round(j["pos"]))] = j


def arms_of(c):
    """Arm names along a reference copy, from its pieces (chr1_hap2_L1 -> L1) or its join."""
    if len(PSPAN[c]) > 1:
        return [(a, b, p_.replace(c + "_", "")) for a, b, p_ in PSPAN[c]]
    js = [j for j in JOINS if j["copy"] == c]
    if js:
        j = js[0]
        return [(0, j["pos"], j["left"]), (j["pos"], L2n[c], j["right"])]
    return []


def top_repeat_vars(n=2):
    if TA is None:
        return []
    t = TA[TA.group.isin(["repeat class", "LTR lineage"]) & (TA.variable != "all repeats")].copy()
    t["key"] = np.where(np.isfinite(t.rho_partial), t.rho_partial.abs(), t.rho.abs())
    return list(t.sort_values("key", ascending=False).variable[:n])


def zcolors(names):
    """Colours of the repeat lines in the tracks: each in its class colour, a later one in coral if taken."""
    out = []
    for n in names:
        c = CCOL[CLASSES.index(n)] if n in CLASSES else ("#185FA5" if n.startswith("Ty3") else
                                                         "#0F6E56" if n.startswith("Ty1") else "#7B2D8E")
        out.append("#D85A30" if c in out else c)
    return out


def fig_tracks():
    ncol = 1 if len(REFC) <= 8 else 2
    nrow = int(math.ceil(len(REFC) / ncol))
    PERH = 2.25
    H = nrow * PERH + 2.1
    fig = plt.figure(figsize=(16 if ncol == 1 else 20, H))
    tv = top_repeat_vars(2)
    zvars = [(v, c_, "-", 1.4) for v, c_ in zip(tv, zcolors(tv))]
    if HAVE_GN:
        zvars.append(("genes per Mb", "#8E6E53", "-", 1.0))
    if HAVE_MK:
        zvars.append(("informative markers per Mb", "#9A9890", "--", 1.0))
    col_of = {"genes per Mb": "genes_per_mb", "informative markers per Mb": "markers_per_mb"}
    Zs = {}
    for v, _, _, _ in zvars:
        x = WT[col_of.get(v, v)].values.astype(float)
        mu, sd = np.nanmean(x), np.nanstd(x)
        Zs[v] = np.clip((x - mu) / sd, -3, 3) if sd > 0 else np.zeros_like(x)
    xmax = max(L2n[c] for c in REFC) / MB
    ymax = np.nanpercentile(WT.cM_hi.values[AN], 99) if AN.any() else 1
    ymax = max(ymax * 1.15, 0.1)
    gs = GridSpec(nrow * 3, ncol, figure=fig, height_ratios=[1.0, 0.72, 0.16] * nrow,
                  left=0.075 if ncol == 1 else 0.06, right=0.99, top=1 - 1.55 / H, bottom=0.95 / H, hspace=0.28,
                  wspace=0.12)
    for i, c in enumerate(REFC):
        rr, cc = (i % nrow) if ncol > 1 else i, (i // nrow) if ncol > 1 else 0
        ax = fig.add_subplot(gs[rr * 3, cc])
        az = fig.add_subplot(gs[rr * 3 + 1, cc], sharex=ax)
        ast = fig.add_subplot(gs[rr * 3 + 2, cc], sharex=ax)
        w_ = WT[WT.chrom == c]
        xs = np.r_[w_.start.values, w_.end.values[-1]] / MB
        ax.fill_between(xs, np.r_[w_.cM_lo.values, w_.cM_lo.values[-1]], np.r_[w_.cM_hi.values, w_.cM_hi.values[-1]],
                        step="post", color="#B4B2A9", alpha=0.55, lw=0)
        ax.step(xs, np.r_[w_.cM_per_Mb.values, w_.cM_per_Mb.values[-1]], where="post", color=INK, lw=1.2)
        for _, r in w_[w_.at_cut].iterrows():
            ax.axvspan(r.start / MB, r.end / MB, facecolor="none", edgecolor="#9A9890", hatch="////", lw=0)
        if DEAD is not None and len(DEAD["rows"]):
            for _, r in DEAD["rows"][DEAD["rows"].chrom == c].iterrows():
                ax.axvspan(r.start_mb, r.end_mb, color=RED, alpha=0.10, lw=0)
        if c in SPANS:
            for a_ in (ax, az):
                a_.axvspan(SPANS[c][0] / MB, SPANS[c][1] / MB, color="#534AB7", alpha=0.08, lw=0)
            ax.axvline(SPANS[c][2] / MB, color="#534AB7", lw=1.0)
        for x in CUTS[c]:
            for a_ in (ax, az):
                a_.axvline(x / MB, color=INK, lw=0.9, ls="--")
        lo = LOCI[(LOCI.chrom == c) & (LOCI.e - LOCI.s >= RDNA_MIN)]
        for _, r in lo.iterrows():
            ax.plot((r.s + r.e) / 2 / MB, ymax * 0.97, marker="v", ms=6, color=CCOL[I45 if r.type == "45S" else I5],
                    mec="white", mew=0.4, zorder=5, clip_on=False)
        sa = ARR[(ARR.chrom == c) & (ARR.kind == ISAT) & (ARR.e - ARR.s >= SAT_BIG)]
        for _, r in sa.iterrows():
            ax.add_patch(Rectangle((r.s / MB, 0), max((r.e - r.s) / MB, xmax / 500), ymax * 0.05, color=CCOL[ISAT],
                                   lw=0, zorder=4))
        ax.set_ylim(0, ymax)
        ax.set_xlim(0, xmax)
        ax.axvline(L2n[c] / MB, color=MUTED, lw=0.6)
        lab = c
        js = [j for j in JOINS if j["copy"] == c and j["letter"]]
        if js:
            lab += " (%s)" % js[0]["letter"]
        arms = arms_of(c)
        ax.set_ylabel("%s\n%s\ncM/Mb" % (lab, " | ".join(a_[2] for a_ in arms)) if arms else "%s\ncM/Mb" % lab,
                      fontsize=7.5)
        for v, col, ls, lw in zvars:
            z = Zs[v][WT.chrom.values == c]
            az.step(xs, np.r_[z, z[-1]], where="post", color=col, ls=ls, lw=lw)
        az.axhline(0, color="#D3D1C7", lw=0.6)
        az.set_ylim(-3.2, 3.2)
        az.set_yticks([-2, 0, 2])
        az.set_ylabel("z", fontsize=7)
        # structure strip
        ast.set_ylim(0, 1)
        ast.set_yticks([])
        for a_, b_, kk, oth in BLOCKS.get(c, []):
            ast.add_patch(Rectangle((a_ / MB, 0), (b_ - a_) / MB, 1, color=KIND_COL.get(kk, "#B4B2A9"), lw=0))
        for p_ in NOHOM:
            if PIECES[p_][0] == c:
                a_, b_ = PIECES[p_][1], PIECES[p_][2]
                ast.add_patch(Rectangle((a_ / MB, 0), (b_ - a_) / MB, 1, facecolor="white", edgecolor="#9A9890",
                                        hatch="....", lw=0.3))
        for a_, b_, nm in arms:
            ast.text((a_ + b_) / 2 / MB, 0.5, nm, ha="center", va="center", fontsize=6.5, color=INK)
        if BLOCKS.get(c) or arms or any(PIECES[p_][0] == c for p_ in NOHOM):
            ast.add_patch(Rectangle((0, 0), L2n[c] / MB, 1, facecolor="none", edgecolor="#9A9890", lw=0.4))
        else:
            ast.spines["left"].set_visible(False)
        for a_ in (ax, az, ast):
            for sp in ("top", "right"):
                a_.spines[sp].set_visible(False)
            a_.tick_params(labelsize=7)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(az.get_xticklabels(), visible=False)
        if rr == nrow - 1 or i == len(REFC) - 1:
            ast.set_xlabel("Mb along the chromosome copy", fontsize=8)
        else:
            plt.setp(ast.get_xticklabels(), visible=False)
    hs = [Line2D([], [], color=INK, lw=1.2, label="crossovers, cM/Mb (each call spread over its interval)"),
          Patch(color="#B4B2A9", alpha=0.55, label="95% band, nuclei resampled")]
    if DEAD is not None and len(DEAD["rows"]):
        hs.append(Patch(color=RED, alpha=0.15, label="dead spot (no crossover, longer than chance)"))
    if any(CUTS.values()):
        hs += [Patch(facecolor="none", edgecolor="#9A9890", hatch="////", label="near a cut of the reference: not tested"),
               Line2D([], [], color=INK, ls="--", lw=0.9, label="cut (a join of the assembly)")]
    if SPANS:
        hs.append(Patch(color="#534AB7", alpha=0.15, label="coupling interval (pollen)"))
    hs += [Line2D([], [], marker="v", ls="", color=CCOL[I45], label="45S rDNA >= 10 kb"),
           Line2D([], [], marker="v", ls="", color=CCOL[I5], label="5S rDNA >= 10 kb"),
           Patch(color=CCOL[ISAT], label="satellite array >= 50 kb")]
    hz = [Line2D([], [], color=col, ls=ls, lw=lw, label=v + (" (strongest repeat association)" if v in tv else ""))
          for v, col, ls, lw in zvars]
    if HAVE_ST:
        hz += [Patch(color=KIND_COL[k], label="other haplotype " + k) for k in ("collinear", "inverted", "translocated")]
        hz += [Patch(facecolor="white", edgecolor="#9A9890", label="other haplotype: not aligned")]
    if NOHOM:
        hz += [Patch(facecolor="white", edgecolor="#9A9890", hatch="....", label="homolog in the same haplotype assembly")]
    fig.legend(handles=hs, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.5 / H), ncol=5, fontsize=7, frameon=False)
    fig.legend(handles=hz, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.95 / H), ncol=4, fontsize=7, frameon=False,
               title="middle row: each variable in standard deviations from its mean over all windows (z)%s" % (
                   "; bottom strip: where the other haplotype aligns" if HAVE_ST or NOHOM else ""), title_fontsize=7)
    fig.suptitle("%s (%s): crossover density against the repeats, %g Mb windows; %d calls in %d nuclei, median "
                 "interval %.1f Mb" % (SP, S, W / MB, len(co), NN, WIDTH_MB), x=0.005, y=1 - 0.12 / H, ha="left",
                 fontsize=11)
    fig.text(0.005, 0.12 / H, "What to look for: whether the crossover line rises and falls with the coloured lines "
             "(the repeat classes or lineages most associated with it), or with genes and markers instead, or with "
             "the structure strip (crossovers need a collinear homolog). Dead spots and the coupling interval are "
             "the deepest 'downs'.", fontsize=7.5, color=MUTED)
    p = os.path.join(A.out, "recombination_vs_repeats_%s.png" % S)
    fig.savefig(p, dpi=120)
    plt.close(fig)
    return p


# ============================================================ figure 2: the tests
def fig_tests():
    fig = plt.figure(figsize=(17, max(9.5, 0.27 * MT + 3)))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[1.15, 1], left=0.25, right=0.985, top=0.9, bottom=0.07,
                  hspace=0.55, wspace=0.42)
    ax = fig.add_subplot(gs[:, 0])
    if TA is not None:
        t = TA.copy()
        order = {"position": 0, "markers": 1, "genes": 2, "structure": 3, "repeat class": 4, "LTR lineage": 5}
        t["o"] = t.group.map(order)
        t = t.sort_values(["o", "variable"]).reset_index(drop=True)
        yy = np.arange(len(t))[::-1]
        for y_, (_, r) in zip(yy, t.iterrows()):
            if np.isfinite(r.get("null_lo", np.nan)):
                ax.plot([r.null_lo, r.null_hi], [y_ + 0.12] * 2, color="#D3D1C7", lw=5, solid_capstyle="butt")
            if np.isfinite(r.get("pnull_lo", np.nan)):
                ax.plot([r.pnull_lo, r.pnull_hi], [y_ - 0.12] * 2, color="#9A9890", lw=2, solid_capstyle="butt")
            for c_, v in r["_pc"] if isinstance(r["_pc"], list) else []:
                ax.plot(v, y_ + 0.12, "o", ms=2.2, color="#9A9890", zorder=3)
            ax.plot(r.rho, y_ + 0.12, "o", mfc="white", mec=INK, ms=6, zorder=4)
            if np.isfinite(r.rho_partial):
                ax.plot(r.rho_partial, y_ - 0.12, "o", color=RED if r.beyond_chance else INK, ms=6, zorder=5)
            elif r.beyond_chance:
                ax.plot(r.rho, y_ + 0.12, "o", color=RED, ms=6, zorder=5)
        ax.set_yticks(yy)
        ax.set_yticklabels(["%s%s" % (v, "  *" if b else "") for v, b in zip(t.variable, t.beyond_chance)],
                           fontsize=7.5)
        for tl, b in zip(ax.get_yticklabels(), t.beyond_chance):
            if b:
                tl.set_color(RED)
                tl.set_fontweight("bold")
        prev = None
        for y_, g in zip(yy, t.group):
            if prev is not None and g != prev:
                ax.axhline(y_ + 0.5, color="#E3E1DA", lw=0.8)
            prev = g
        ax.axvline(0, color=MUTED, lw=0.8)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-0.7, len(t) - 0.3)
        ax.set_xlabel("Spearman rho with crossover density (windows)", fontsize=8)
        ax.set_title("a. One variable at a time: open = raw, filled = after distance to the end and markers;\n"
                     "bars = 95%% of %d rotations of the variable along its chromosome; grey dots = each chromosome;\n"
                     "red * = beyond chance after Bonferroni (p < %.1e, %d variables)" % (A.nperm, ALPHA, MT),
                     fontsize=8, loc="left")
    ax.tick_params(labelsize=7)
    ab = fig.add_subplot(gs[0, 1])
    if GLM is not None:
        R = GLM["R"]
        names = [k for k in ("position", "markers", "genes", "structure", "repeats", "position + markers",
                             "position + markers + repeats", "everything but the repeats", "everything") if k in R]
        yy = np.arange(len(names))[::-1]
        ab.barh(yy, [100 * R[k] for k in names], color=["#185FA5" if "repeats" in k and "but" not in k else "#B4B2A9"
                                                         for k in names], height=0.6)
        nl = GLM["nul"]
        ref = {"repeats": 0.0, "position + markers + repeats": R["position + markers"],
               "everything": R["everything but the repeats"]}
        for y_, k in zip(yy, names):
            ab.text(100 * R[k] + 0.5, y_, "%.1f%%" % (100 * R[k]), va="center", fontsize=7)
            if k in ref:
                j = list(ref).index(k)
                ab.plot([100 * (ref[k] + np.percentile(nl[:, j], 95))] * 2, [y_ - 0.35, y_ + 0.35], color=RED, lw=1.5)
        ab.set_yticks(yy)
        ab.set_yticklabels(names, fontsize=7.5)
        ab.set_xlabel("crossover deviance explained, %% (Poisson model of %d windows)" % GLM["n"], fontsize=8)
        ab.set_title("b. How much of the ups and downs each set explains;\nred tick = 95th percentile with the repeat "
                     "tracks rotated (%d rotations)" % GLM["nn"], fontsize=8, loc="left")
        ab.set_xlim(0, max(5, 100 * max(R.values()) * 1.25))
    ab.tick_params(labelsize=7)
    ac = fig.add_subplot(gs[1, 1])
    if TARR is not None and len(TARR) and "ratio" in TARR:
        t = TARR[np.isfinite(TARR.get("ratio", pd.Series(dtype=float)))].reset_index(drop=True)
        yy = np.arange(len(t))[::-1]
        LO, HI = 0.1, 10.0
        cl = lambda v: min(max(v, LO * 1.06), HI / 1.06)
        for y_, (_, r) in zip(yy, t.iterrows()):
            ac.plot([cl(r.null_lo), cl(r.null_hi)], [y_, y_], color="#D3D1C7", lw=5, solid_capstyle="butt")
            ac.plot(cl(r.ratio), y_, "<" if r.ratio <= LO else "o", color=RED if r.p < 0.05 else INK, ms=6)
            if r.ratio <= LO:
                ac.text(LO * 1.25, y_ + 0.18, "%.2f" % r.ratio, fontsize=6.5, color=INK)
        ac.set_yticks(yy)
        ac.set_yticklabels(["%s (%d, %.0f Mb)  p %s" % (r.set, r.n, r.zone_mb, ptxt(r.p, r.nn))
                            for _, r in t.iterrows()], fontsize=7)
        ac.axvline(1, color=MUTED, lw=0.8)
        ac.set_xscale("log")
        ac.set_xlim(LO, HI)
        tk = [0.1, 0.25, 0.5, 1, 2, 4, 10]
        ac.set_xticks(tk)
        ac.set_xticklabels(["%g" % v for v in tk])
        ac.minorticks_off()
        ac.set_xlabel("crossovers within %g Mb of the arrays: observed / expected from the chromosome" % A.array_mb,
                      fontsize=8)
    ac.set_title("c. At the big tandem arrays; bar = 95% with the arrays rotated along their chromosome", fontsize=8,
                 loc="left")
    ac.tick_params(labelsize=7)
    ad = fig.add_subplot(gs[2, 1])
    tv = top_repeat_vars(1)
    if tv:
        v = tv[0]
        m = AN & np.isfinite(WT[v].values)
        cols = plt.get_cmap("tab20")
        for i, c in enumerate(REFC):
            mm = m & (WT.chrom.values == c)
            ad.scatter(100 * WT[v].values[mm], WT.cM_per_Mb.values[mm], s=7, color=cols(i % 20), label=c, alpha=0.8,
                       lw=0)
        xv, yv = 100 * WT[v].values[m], WT.cM_per_Mb.values[m]
        qs = np.unique(np.percentile(xv, np.linspace(0, 100, 11)))
        if len(qs) > 2:
            bi = np.clip(np.searchsorted(qs, xv, side="right") - 1, 0, len(qs) - 2)
            bx = [xv[bi == k].mean() for k in range(len(qs) - 1) if (bi == k).any()]
            by = [yv[bi == k].mean() for k in range(len(qs) - 1) if (bi == k).any()]
            ad.plot(bx, by, color=INK, lw=1.6, label="mean per tenth")
        r = TA.set_index("variable").loc[v]
        ad.set_xlabel("%s, %% of the window" % v, fontsize=8)
        ad.set_ylabel("cM/Mb", fontsize=8)
        ad.set_title("d. The repeat variable with the strongest association: rho %+.2f (p %s);\n"
                     "after distance to the end and markers %+.2f (p %s)" % (
                         r.rho, ptxt(r.p_rot, A.nperm), r.rho_partial, ptxt(r.p_partial, A.nperm)), fontsize=8,
                     loc="left")
        ad.legend(fontsize=6, ncol=2 if len(REFC) <= 8 else 4, frameon=False, markerscale=1.5)
    ad.tick_params(labelsize=7)
    fig.suptitle("%s (%s): do the repeats explain the crossover landscape? %d windows of %g Mb tested" % (
        SP, S, int(AN.sum()), W / MB), x=0.01, ha="left", fontsize=11)
    fig.text(0.01, 0.012, "What to look for: a filled dot outside its bar (an association beyond what the variable's "
             "own patchiness gives by chance), the blue bars against their red ticks (repeats explaining more than "
             "rotated repeat tracks would), and dots off 1 in c (crossovers depleted or enriched at the arrays).",
             fontsize=7.5, color=MUTED)
    p = os.path.join(A.out, "recombination_tests_%s.png" % S)
    fig.savefig(p, dpi=120)
    plt.close(fig)
    return p


# ============================================================ figure 3: the breakpoints
def fig_breaks():
    BT, TT = BPR["BT"], BPR["tests"]
    nb = len(BT)
    fig = plt.figure(figsize=(17, max(7.5, 0.38 * max(nb, NK) + 4.3)))
    H = fig.get_figheight()
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 1.3], left=0.11, right=0.985, top=1 - 1.25 / H,
                  bottom=2.35 / H, wspace=0.62)
    ax = fig.add_subplot(gs[0, 0])
    t = TT[TT.test.str.startswith("breakpoints: share")]
    sets_ = [s_ for s_ in ("all joins", "haplotype breakpoints") if s_ in set(t.set)]
    mk_ = {"all joins": ("o", INK), "haplotype breakpoints": ("D", "#1D9E75")}
    yy = np.arange(NK)[::-1]
    fin = t.ratio[np.isfinite(t.ratio) & (t.ratio > 0)]
    LO, HI = 0.05, float(min(1e4, max(50.0, 2.0 * fin.max()))) if len(fin) else 50.0
    cl = lambda v: min(max(v, LO * 1.08), HI / 1.08)
    for si, s_ in enumerate(sets_):
        tt = t[t.set == s_].set_index("what")
        dy = 0.15 * (1 if si == 0 else -1) if len(sets_) > 1 else 0
        for y_, k in zip(yy, CLASSES):
            if k not in tt.index:
                continue
            r = tt.loc[k]
            if not (np.isfinite(r.null_mean) and r.null_mean > 0):
                continue
            ax.plot([cl(r.null_lo / r.null_mean), cl(r.null_hi / r.null_mean)], [y_ + dy] * 2, color="#D3D1C7", lw=4,
                    solid_capstyle="butt")
            col = RED if r.p < 0.05 else mk_[s_][1]
            ax.plot(cl(r.ratio), y_ + dy, "<" if r.ratio <= LO else mk_[s_][0], color=col, mfc=col, ms=6)
            if r.ratio >= HI / 1.08:
                ax.text(HI / 1.3, y_ + dy + 0.22, "%.0fx" % r.ratio, fontsize=6.5, color=col, ha="right")
    ax.set_yticks(yy)
    ax.set_yticklabels(CLASSES, fontsize=7.5)
    ax.axvline(1, color=MUTED, lw=0.8)
    ax.set_xscale("log")
    ax.set_xlim(LO, HI)
    ax.set_xlabel("share within %g Mb of the breakpoints / share at random places on the same copies" % A.near_mb,
                  fontsize=8)
    ax.legend(handles=[Line2D([], [], marker=mk_[s_][0], ls="", color=mk_[s_][1],
                              label="%s (%d)" % (s_, int(t[t.set == s_].n.iloc[0]))) for s_ in sets_] +
              [Patch(color="#D3D1C7", label="95%% of %d random placements" % BPR["nn"]),
               Line2D([], [], marker="o", ls="", color=RED, label="p < 0.05 (two-sided, not corrected)"),
               Line2D([], [], marker="<", ls="", color=INK, label="none at all (share 0)")],
              fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.075), ncol=3)
    ax.set_title("a. What the breakpoints are made of", fontsize=8.5, loc="left")
    ax.tick_params(labelsize=7)
    ab = fig.add_subplot(gs[0, 1])
    yy = np.arange(nb)[::-1]
    sym = {"45S": ("v", CCOL[I45], 0.2), "5S": ("v", CCOL[I5], 0.0), "satellite": ("s", CCOL[ISAT], -0.2)}
    for y_, (_, r) in zip(yy, BT.iterrows()):
        for lab, (m_, col, dy) in sym.items():
            d = r["nearest_%s_mb" % lab]
            if np.isfinite(d):
                ab.plot(max(d, 0.01), y_ + dy, m_, color=col, ms=6)
            da = r.get("nearest_%s_any_mb" % lab, np.nan)
            if lab != "satellite" and np.isfinite(da) and not (np.isfinite(d) and abs(da - d) < 1e-9):
                ab.plot(max(da, 0.01), y_ + dy, m_, mfc="none", mec=col, ms=6)
    ab.axvline(A.near_mb, color=MUTED, lw=0.8, ls=":")
    ab.set_yticks(yy)
    ab.set_yticklabels(["%s%s" % (r["name"], "  [%s, phi %+.2f]" % (r.verdict, r.phi) if r.set == "join" else "")
                        for _, r in BT.iterrows()], fontsize=7)
    for tl, (_, r) in zip(ab.get_yticklabels(), BT.iterrows()):
        if r.set == "join" and not r.joined:
            tl.set_color(RED)
    ab.set_xscale("log")
    ab.set_xlim(0.008, 500)
    tk = [0.01, 0.1, 1, 10, 100]
    ab.set_xticks(tk)
    ab.set_xticklabels(["0 (on it)", "0.1", "1", "10", "100"])
    ab.minorticks_off()
    ab.set_xlabel("Mb from the breakpoint to the nearest array on the same copy (dotted: %g Mb)" % A.near_mb,
                  fontsize=8)
    tj = TT[(TT.test == "breakpoints: with %s within %g Mb" % ("rDNA", A.near_mb)) & (TT.set == "all joins")]
    extra = ""
    if len(tj):
        r = tj.iloc[0]
        extra = "\njoins with rDNA (>= %d kb) within %g Mb: %d of %d; at random places %.1f (p %s)" % (
            RDNA_MIN // 1000, A.near_mb, r.obs, r.n, r.null_mean, ptxt(r.p, BPR["nn"]))
    ab.set_title("b. Nearest 45S (red) and 5S (orange) rDNA array >= %d kb (filled; open = any size) and satellite "
                 "array >= %d kb (purple)%s" % (RDNA_MIN // 1000, SAT_TICK // 1000, extra), fontsize=8.5, loc="left")
    ab.tick_params(labelsize=7)
    fig.suptitle("%s (%s): the repeats at the joins of the assembly (pollen verdict in brackets; red = the pollen do "
                 "not inherit the join)\nand at the inversion / translocation breakpoints between the haplotypes" % (
                     SP, S), x=0.01, y=1 - 0.1 / H, ha="left", va="top", fontsize=10.5)
    fig.text(0.01, 0.12 / H, "\n".join(textwrap.wrap(
        "What to look for: a class far right of its grey bar (enriched at the breakpoints), and rDNA or one satellite "
        "family next to several joins. Two cautions: joins that share an arm end are not independent (the same arm end "
        "is seen on both homologs), so the joins' p-values are optimistic; and a join of the assembly is where the "
        "scaffolder joined two contigs, and contigs tend to end at long arrays the assembler cannot cross (rDNA, "
        "satellites), so rDNA at a join can be the assembly's doing. The pollen verdict separates the two: a join the "
        "pollen inherit is a real junction.", 230)), fontsize=7.5, color=MUTED, va="bottom")
    p = os.path.join(A.out, "breakpoints_vs_repeats_%s.png" % S)
    fig.savefig(p, dpi=120)
    plt.close(fig)
    return p


for nm, fn in (("figure: tracks", fig_tracks), ("figure: tests", fig_tests)):
    p = section(nm, fn)
    if p:
        OUTS.append(p)
if BPR is not None:
    p = section("figure: breakpoints", fig_breaks)
    if p:
        OUTS.append(p)

# ============================================================ summary
L_ = []


def P(s=""):
    L_.append(s)


P("REPEATS x RECOMBINATION  %s (%s), CARP %s, %s" % (S, SP, carp_version or "?", time.strftime("%Y-%m-%d %H:%M")))
P("crossover reference: " + ", ".join("%s = %s" % (c, " + ".join(a_[2] for a_ in arms_of(c))) if arms_of(c) else c
                                      for c in REFC))
for j in JOINS:
    P("  join %s at %.2f Mb: pollen phi %+.2f (%s)" % (j["name"], j["pos"] / MB, j["phi"],
                                                    "inherited as joined" if j["joined"] else "NOT inherited as joined"))
P("crossovers: %d calls in %d nuclei; median call interval %.2f Mb (the resolution of everything below)" % (
    len(co), NN, WIDTH_MB))
P("windows: %g Mb; %d tested of %d (%d within %.1f Mb of a cut left out)" % (W / MB, int(AN.sum()), NW,
                                                                           int(WT.at_cut.sum()), W / 2 / MB))
tot = sum(L2n[c] for c in REFC)
cls_bp = sum(CSK[OFF[c] + NF[c]] - CSK[OFF[c]] for c in REFC)
P("repeats on the crossover reference (%.2f Gb): %s" % (tot / 1e9, ", ".join(
    "%s %.1f%%" % (CLASSES[k], 100 * cls_bp[k] / tot) for k in np.argsort(-cls_bp) if cls_bp[k] / tot >= 0.001)))
v = WT.loc[AN, "all repeats"].values
if len(v):
    P("  'all repeats' per window: median %.0f%%, middle 90%% %.0f-%.0f%%; crossovers per window: median %.2f cM/Mb, "
      "middle 90%% %.2f-%.2f" % (100 * np.median(v), 100 * np.percentile(v, 5), 100 * np.percentile(v, 95),
                                 np.median(WT.cM_per_Mb.values[AN]), np.percentile(WT.cM_per_Mb.values[AN], 5),
                                 np.percentile(WT.cM_per_Mb.values[AN], 95)))
P()
P("Q1a ONE VARIABLE AT A TIME (Spearman rho across windows; p from %d rotations along each chromosome;" % A.nperm)
P("    Bonferroni: p < %.1e over %d variables; 'after' = partial, after distance to the end and markers)" % (ALPHA, MT))
if TA is not None:
    t = TA.copy()
    t["key"] = np.where(np.isfinite(t.rho_partial), t.rho_partial.abs(), t.rho.abs())
    oc = "outside_coupling_rho" in t
    if oc:
        P("    'outside': the same with the windows inside the coupling interval left out (partial; raw for the "
          "distance to the end and the markers)")
    for _, r in t.sort_values(["beyond_chance", "key"], ascending=[False, False]).iterrows():
        if not np.isfinite(r.rho):
            P("    %-40s not tested: %d windows with a value, or no variation" % (r.variable[:40], r.n_windows))
            continue
        P("  %s %-40s rho %+.2f (p %-7s) after %s (p %-7s)%s same sign on %s chr" % (
            "*" if r.beyond_chance else " ", r.variable[:40], r.rho, ptxt(r.p_rot, A.nperm),
            "%+.2f" % r.rho_partial if np.isfinite(r.rho_partial) else "  -  ",
            ptxt(r.p_partial, A.nperm) if np.isfinite(r.p_partial) else "-",
            (" outside %s (p %-7s)" % ("%+.2f" % r.outside_coupling_rho, ptxt(r.outside_coupling_p, A.nperm))
             if np.isfinite(r.get("outside_coupling_rho", np.nan)) else " outside   -   (p -      )") if oc else "",
            r.same_sign))
P()
if GLM is not None:
    R = GLM["R"]
    P("Q1b HOW MUCH EACH SET EXPLAINS (Poisson deviance explained, %d windows):" % GLM["n"])
    P("  alone: " + ", ".join("%s %.1f%%" % (k, 100 * R[k]) for k in ("position", "markers", "genes", "structure",
                                                                    "repeats") if k in R))
    P("  position + markers %.1f%%; + repeats %.1f%%; everything but the repeats %.1f%%; everything %.1f%%" % (
        100 * R["position + markers"], 100 * R["position + markers + repeats"], 100 * R["everything but the repeats"],
        100 * R["everything"]))
    for lab, o_, nl, p in zip(GLM["labels"], GLM["obs"], GLM["nul"].T, GLM["p"]):
        P("  %-34s %+5.1f%% (rotated repeat tracks: median %.1f%%, 95th percentile %.1f%%; p %s)" % (
            lab, 100 * o_, 100 * np.median(nl), 100 * np.percentile(nl, 95), ptxt(p, GLM["nn"])))
    P("  (repeat classes in the model: %s; structure = %s)" % (", ".join(GLM["classes"]), " + ".join(
        (["the synteny blocks"] if HAVE_ST else []) + (["the coupling interval"] if COUP.any() else [])) or "none"))
P()
if TARR is not None and len(TARR):
    P("Q1c CROSSOVERS WITHIN %g Mb OF THE BIG TANDEM ARRAYS (observed / expected from the chromosome's own density;"
      % A.array_mb)
    P("    null: the arrays rotated along their chromosome)")
    for _, r in TARR.iterrows():
        if np.isfinite(r.get("ratio", np.nan)):
            P("  %-40s %3d arrays, %5.0f Mb of zone: %.2f x (rotated: %.2f, 95%% %.2f-%.2f; p %s)" % (
                r.set[:40], r.n, r.zone_mb, r.ratio, r.null_median, r.null_lo, r.null_hi, ptxt(r.p, r.nn)))
        else:
            P("  %-40s %s" % (r.set[:40], r.get("note", "")))
    P()
if DEAD is not None:
    d = DEAD["rows"]
    P("Q1d DEAD SPOTS: runs >= %d Mb without a crossover midpoint (chance: %.1f such runs, 95%% %.0f-%.0f); "
      "repeats inside vs the whole copy" % ((DEAD["thr"],) + DEAD["chance"]))
    for _, r in d.iterrows():
        P("  %-10s %5.0f-%5.0f Mb: repeats %.0f%% (copy %.0f%%), Ty3 %.0f%% (%.0f%%), Ty1 %.0f%% (%.0f%%), satellite %.0f kb, "
          "rDNA arrays %d; markers %.2fx the copy; other haplotype collinear %s, inverted %s, translocated %s, "
          "not aligned %s; coupling %.0f%%%s" % (
              r.chrom, r.start_mb, r.end_mb, 100 * r.repeats, 100 * r.repeats_chrom, 100 * r.gypsy, 100 * r.gypsy_chrom,
              100 * r.copia, 100 * r.copia_chrom, r.satellite_kb, r.rdna_arrays, r.markers_vs_chrom,
              ff(r.collinear), ff(r.inverted), ff(r.translocated), ff(r.not_aligned), 100 * r.coupling_share,
              "; a cut inside" if r.cut_inside else ""))
    if not len(d):
        P("  none")
    P()
def near_txt(r, lab):
    d = r["nearest_%s_mb" % lab]
    t = "%.2f Mb %s" % (d, r["nearest_%s_side" % lab]) if np.isfinite(d) else "none >= %d kb" % (RDNA_MIN // 1000)
    da = r.get("nearest_%s_any_mb" % lab, np.nan)
    if np.isfinite(da) and not (np.isfinite(d) and abs(da - d) < 1e-9):
        t += " (any size: %.2f Mb, %.1f kb)" % (da, r["nearest_%s_any_kb" % lab])
    return t


if BPR is not None:
    BT, TT = BPR["BT"], BPR["tests"]
    P("Q2 BREAKPOINTS: repeats within %g Mb, against %d random places on the same copy" % (A.near_mb, BPR["nn"]))
    P("    left / right = toward the start / the end of the copy: at a join 'L1|P', left is the L1 side")
    for _, r in BT.iterrows():
        P("  %-46s%s" % (r["name"][:46], "  [%s, phi %+.2f]" % (r.verdict, r.phi) if r.set == "join" else ""))
        P("      nearest 45S %s, 5S %s, satellite %s; within %g Mb: repeats %.0f%% (Ty3 %.0f%%, Ty1 %.0f%%), %s" % (
            near_txt(r, "45S"), near_txt(r, "5S"),
            "%.2f Mb %s (%s)" % (r.nearest_satellite_mb, r.nearest_satellite_side, r.get("nearest_satellite_family", ""))
            if np.isfinite(r.nearest_satellite_mb) else "none", A.near_mb,
            100 * np.nansum([r[k] for k in CLASSES]), 100 * r["Ty3/gypsy"], 100 * r["Ty1/copia"],
            "satellites: " + r.satellite_families_within if r.satellite_families_within else "no satellite array"))
    for s_ in ("all joins", "joins the pollen inherit", "joins the pollen do not inherit", "haplotype breakpoints"):
        tt = TT[TT.set == s_]
        if not len(tt):
            continue
        cnt = tt[tt.test.str.startswith("breakpoints: with")]
        P("  %s (%d): with rDNA nearby: %s" % (s_, int(tt.n.iloc[0]), "; ".join(
            "%s within %g Mb %d (random %.1f, p %s)" % (r["what"], r.R_mb, r.obs, r.null_mean, ptxt(r.p, BPR["nn"]))
            for _, r in cnt[cnt.what.isin(["45S", "5S", "satellite"])].iterrows() if r.R_mb == A.near_mb)))
        P("      within 0.5 / 1 / 2.5 Mb, any rDNA: %s" % ", ".join(
            "%d (random %.1f, p %s)" % (r.obs, r.null_mean, ptxt(r.p, BPR["nn"])) for _, r in
            cnt[cnt.what == "rDNA"].sort_values("R_mb").iterrows()))
        sh = tt[tt.test.str.startswith("breakpoints: share") & (tt.p < 0.05)]
        P("      classes off chance (p < 0.05, uncorrected): %s" % (", ".join(
            "%s %.2fx (p %s)" % (r["what"], r.ratio, ptxt(r.p, BPR["nn"])) for _, r in sh.iterrows()) or "none"))
    P()
for n_ in NOTES:
    P("note: " + n_)
for f_ in FAILED:
    P("PART FAILED: " + f_)
P("files read:\n  " + "\n  ".join(USED))
sp_ = os.path.join(A.out, "summary_%s.txt" % S)
OUTS.append(sp_)
P("wrote " + ", ".join(OUTS))
P("(%.0f s)" % (time.time() - T0))
open(sp_, "w").write("\n".join(L_) + "\n")
print("\n".join(L_))
if FAILED:
    raise SystemExit(1)
