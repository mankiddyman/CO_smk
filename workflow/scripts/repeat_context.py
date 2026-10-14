#!/usr/bin/env python3
"""repeat_context.py -- the repeats (CARP) where the pollen and the assembly disagree, and FISH probe candidates.

Reads the CARP repeat annotation of the species' two-haplotype assembly, as reproducible_phd publishes it
(refs/{species}_repeats/: Repeat_Annotation_Unified.gff3 = one classification per base; TideCluster's family table),
and for a cut reference the karyotype test of cosegregation_2n.py (reports/cosegregation_2n/karyotype_test_SAMPLE.tsv,
the coupling point's 95% interval from logs/cosegregation_2n_SAMPLE.log, hic_scan_SAMPLE.tsv) and the breakpoints
(qc/reference/SAMPLE/breakpoints.csv).

Writes to OUT/:
  repeat_landscape_SAMPLE.png  every chromosome copy at 1 Mb: share of each repeat class, 45S / 5S rDNA arrays, the
                               biggest tandem-repeat families; arms, joins and the coupling interval for a cut reference
  repeat_joins_SAMPLE.png      (cut reference) +-3 Mb around each join of the two-haplotype assembly, with the
                               pollen's verdict on the join and where the two haplotypes stop aligning
  repeat_coupling_SAMPLE.png   (when the pollen couple two chromosomes) both copies of each, the coupling point and
                               its interval, the Hi-C scan's strongest blocks, the tandem-repeat families in both intervals
  rdna_loci_SAMPLE.csv         45S and 5S rDNA arrays (pieces within 20 kb merged)
  fish_probes_SAMPLE.csv       every tandem-repeat family: size, arrays, the chromosome carrying most of it, at a join,
                               in the coupling interval; ranked for FISH (chromosome-specific and large first)
Run from the CO_smk root. About 1-3 min, most of it reading the annotation.
"""
import argparse
import csv
import os
import re
import sys
from urllib.parse import unquote

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("--out", default="reports/repeat_context")
ap.add_argument("--refs", default="/netscratch/dep_mercier/grp_marques/Aaryan/refs")
ap.add_argument("--coseg", default="reports/cosegregation_2n")
ap.add_argument("--zoom_mb", type=float, default=3.0, help="Mb either side of a join")
ap.add_argument("--bin_mb", type=float, default=1.0, help="bin of the landscape")
A = ap.parse_args()
S = A.sample
MB = 1_000_000
BIN = int(A.bin_mb * MB)
MERGE = 20_000             # pieces of one family (or of rDNA) closer than this are one array
ARRAY_MIN = 5_000          # arrays shorter than this are not listed as arrays
os.makedirs(A.out, exist_ok=True)
USED, NOTES = [], []

# repeat classes, bottom to top of the stacked share; rDNA last so it shows on top
CLASSES = ["Ty3/gypsy", "Ty1/copia", "other LTR", "LINE", "DNA transposon", "simple / low complexity",
           "other repeat", "satellite", "5S rDNA", "45S rDNA"]
CCOL = ["#185FA5", "#0F6E56", "#7FA7D4", "#8DA33B", "#B07A10", "#C9C7BE", "#E3E1DA", "#534AB7", "#E07B39", "#A32D2D"]
I45, I5, ISAT = CLASSES.index("45S rDNA"), CLASSES.index("5S rDNA"), CLASSES.index("satellite")
ARM_COL = {"L1": "#5DCAA5", "L2": "#F0997B", "P": "#AFA9EC", "Q": "#FAC775", "X": "#85B7EB", "core": "#D3D1C7"}
FAM_COL = ["#7B2D8E", "#D4A017", "#2E86AB", "#C0392B", "#16A085", "#8E6E53", "#E84393", "#2C3E50", "#7F8C2D",
           "#5D6D7E"]
INK, MUTED = "#2C2C2A", "#6B6A64"


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


def natural(c):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", c)]


def num(c):
    m = re.match(r"chr(\d+)", c)
    return int(m.group(1)) if m else 0


def hap(c):
    return c[-1]


def fai_lens(p):
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(p) if l.strip()}


# ------------------------------------------------------------------ the sample, its species, the two-haplotype assembly
rows = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}
if S not in rows:
    raise SystemExit("%s is not in config/samples.csv" % S)
row = rows[S]
SPECIES = row["species"]
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

REP = os.path.join(A.refs, "%s_repeats" % SPECIES)
GFF = os.path.join(REP, "Repeat_Annotation_Unified.gff3")
if not os.path.exists(GFF):
    raise SystemExit("no CARP annotation published for %s yet (%s); reproducible_phd publishes it with "
                     "repeats.smk (publish_repeats)" % (SPECIES, GFF))
USED.append(GFF)
seqf = os.path.join(REP, "sequences.tsv")
if os.path.exists(seqf):
    LA = {r["seqid"]: int(r["length"]) for r in csv.DictReader(open(seqf), delimiter="\t")}
    bad = [c for c in COPIES if LA.get(c) != L2n[c]]
    if bad:
        raise SystemExit("the CARP annotation was made on a different assembly (lengths differ for %s): %s vs %s"
                         % (", ".join(bad[:4]), seqf, "the assembly in MANIFEST.txt"))
    USED.append(seqf)
carp_version = ""
mf = os.path.join(REP, "MANIFEST.txt")
if os.path.exists(mf):
    m = re.search(r"^carp_version:\s*(\S+)", open(mf).read(), re.M)
    carp_version = m.group(1) if m else ""
    USED.append(mf)

# ------------------------------------------------------------------ the annotation: one row per top-level feature
fams = {}                       # TRC name -> index
cs, ss, es, ks, fs = [], [], [], [], []
n_lines = 0
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
        if "ID=UA_L2_" in at:          # nested members of a tandem LTR container: counted in the container
            continue
        n_lines += 1
        k = at.find("classification=")
        cls = at[k + 15:].split(";", 1)[0].rstrip() if k >= 0 else "Unknown"
        if "%" in cls:
            cls = unquote(cls)
        kc = coarse(cls)
        fi = -1
        if kc in (ISAT, I45, I5):
            k2 = at.find("Name=")
            nm = at[k2 + 5:].split(";", 1)[0].rstrip() if k2 >= 0 else ""
            if nm.startswith("TRC_"):
                fi = fams.setdefault(nm, len(fams))
        cs.append(ci)
        ss.append(int(f[3]) - 1)
        es.append(int(f[4]))
        ks.append(kc)
        fs.append(fi)
F = pd.DataFrame({"c": np.array(cs, dtype=np.int16), "s": np.array(ss, dtype=np.int64),
                  "e": np.array(es, dtype=np.int64), "k": np.array(ks, dtype=np.int8),
                  "f": np.array(fs, dtype=np.int32)})
del cs, ss, es, ks, fs
FAMS = {i: n for n, i in fams.items()}
print("%s (%s): %d features of the CARP annotation%s on %d chromosome copies" % (
    S, SPECIES, len(F), " (CARP %s)" % carp_version if carp_version else "", len(COPIES)))

# TideCluster's family table: monomer length and rDNA calls
MONO, RDNA_T = {}, {}
for tab, pref in (("trc_table_default.tsv", "TRC_"), ("trc_table_short.tsv", "TRC_S_")):
    p = os.path.join(REP, tab)
    if not os.path.exists(p):
        continue
    USED.append(p)
    for r in csv.DictReader(open(p), delimiter="\t"):
        tid = r.get("TRC_ID", "")
        name = pref + tid.split("TRC_", 1)[1] if tid.startswith("TRC_") else tid
        for col in ("monomer_tarean", "monomer_kite", "prevalent_founder"):
            v = (r.get(col) or "").strip()
            if v and v.replace(".", "", 1).isdigit() and float(v) > 0:
                MONO[name] = int(float(v))
                break
        if (r.get("rDNA_type") or "").strip():
            RDNA_T[name] = r["rDNA_type"].strip()


def bin_cov(s, e, nb, B):
    """bp of the intervals [s, e) in each bin of B bp."""
    cov = np.zeros(nb)
    if not len(s):
        return cov
    b0, b1 = s // B, (e - 1) // B
    one = b0 == b1
    np.add.at(cov, b0[one], (e - s)[one])
    for a, b in zip(s[~one], e[~one]):
        for k in range(a // B, min((b - 1) // B, nb - 1) + 1):
            cov[k] += min(b, (k + 1) * B) - max(a, k * B)
    return np.minimum(cov, B)


def merge(s, e, gap=MERGE):
    """Merge sorted intervals closer than gap: returns starts, ends, number of pieces."""
    if not len(s):
        return [], [], []
    o = np.argsort(s)
    s, e = s[o], e[o]
    S_, E_, N_ = [s[0]], [e[0]], [1]
    for a, b in zip(s[1:], e[1:]):
        if a - E_[-1] <= gap:
            E_[-1] = max(E_[-1], b)
            N_[-1] += 1
        else:
            S_.append(a), E_.append(b), N_.append(1)
    return S_, E_, N_


BYC = {c: F[F.c == i] for c, i in CIDX.items()}
NB = {c: -(-L2n[c] // BIN) for c in COPIES}
COV = {c: np.stack([bin_cov(BYC[c].s.values[BYC[c].k.values == k], BYC[c].e.values[BYC[c].k.values == k], NB[c], BIN)
                    for k in range(len(CLASSES))], axis=1) for c in COPIES}

# rDNA loci (all pieces of a type, TideCluster arrays and RepeatMasker subunits, merged within 20 kb)
RLOCI = []
for c in COPIES:
    d = BYC[c]
    for k, nm in ((I45, "45S"), (I5, "5S")):
        x = d[d.k == k]
        for a, b, n in zip(*merge(x.s.values, x.e.values)):
            RLOCI.append(dict(type=nm, copy=c, start=int(a), end=int(b), kb=(b - a) / 1e3, pieces=n))

# tandem-repeat families: arrays per copy
SAT = F[F.f >= 0]
FAMROWS = []
ARRAYS = {}                      # family -> list of (copy, start, end)
for fi, g in SAT.groupby("f"):
    name = FAMS[fi]
    arrs = []
    per_copy = {}
    for ci, gc in g.groupby("c"):
        c = COPIES[ci]
        per_copy[c] = int((gc.e - gc.s).sum())
        for a, b, n in zip(*merge(gc.s.values, gc.e.values)):
            if b - a >= ARRAY_MIN:
                arrs.append((c, int(a), int(b)))
    ARRAYS[name] = arrs
    tot = sum(per_copy.values())
    by_num = {}
    for c, v in per_copy.items():
        by_num[num(c)] = by_num.get(num(c), 0) + v
    topn = max(by_num, key=by_num.get)
    kinds = g.k.value_counts()
    kind = "45S rDNA" if kinds.get(I45, 0) >= kinds.max() else ("5S rDNA" if kinds.get(I5, 0) >= kinds.max()
                                                               else "satellite")
    if name in RDNA_T and kind == "satellite":
        kind = RDNA_T[name] + " rDNA"
    FAMROWS.append(dict(family=name, kind=kind, monomer_bp=MONO.get(name, ""), total_kb=round(tot / 1e3, 1),
                        arrays=len(arrs), largest_array_kb=round(max([b - a for _, a, b in arrs], default=0) / 1e3, 1),
                        top_chromosome="chr%d" % topn, share_on_top=round(by_num[topn] / tot, 3) if tot else 0,
                        copies="; ".join("%s %.0f kb" % (c, v / 1e3) for c, v in
                                         sorted(per_copy.items(), key=lambda kv: -kv[1]) if v >= 0.1 * tot)))
FT = pd.DataFrame(FAMROWS, columns=["family", "kind", "monomer_bp", "total_kb", "arrays", "largest_array_kb",
                                   "top_chromosome", "share_on_top", "copies"])
if len(FT):
    FT = FT.sort_values("total_kb", ascending=False).reset_index(drop=True)
else:
    FT.to_csv(os.path.join(A.out, "fish_probes_%s.csv" % S), index=False)
TOPF = [r.family for r in FT.itertuples() if r.kind == "satellite"][:len(FAM_COL)]
FCOL = {f: FAM_COL[i] for i, f in enumerate(TOPF)}

# ------------------------------------------------------------------ the pollen's verdicts (cut reference)
JOINS, COUPLES, SCAN = [], [], []
kt = os.path.join(A.coseg, "karyotype_test_%s.tsv" % S)
if os.path.exists(kt):
    USED.append(kt)
    for r in csv.DictReader(open(kt), delimiter="\t"):
        if r["kind"] == "join":
            m = re.match(r"^(chr\d+_hap[12])(?: \(([A-D])\))?: (\S+) \| (\S+) at ([\d.]+) Mb", r["what"])
            if m and m.group(1) in CIDX:
                JOINS.append(dict(copy=m.group(1), letter=m.group(2) or "", left=m.group(3), right=m.group(4),
                                  pos=float(m.group(5)) * MB, phi=float(r["phi"]), lo=float(r["ci_lo"] or "nan"),
                                  hi=float(r["ci_hi"] or "nan"), detail=r["detail"]))
        elif r["kind"] == "couple":
            m = re.match(r"^(\S+) ([\d.]+) Mb x (\S+) ([\d.]+) Mb", r["what"])
            if m:
                COUPLES.append(dict(p=m.group(1), pp=float(m.group(2)) * MB, q=m.group(3), qp=float(m.group(4)) * MB,
                                    phi=float(r["phi"]), ci_p=None, ci_q=None))
pieces = {}
pp_ = os.path.join(refdir, "pieces.tsv")
if os.path.exists(pp_):
    for r in csv.DictReader(open(pp_), delimiter="\t"):
        pieces[r["piece"]] = (r["source"], int(r["start"]) - 1)
lg = "logs/cosegregation_2n_%s.log" % S
if COUPLES and os.path.exists(lg):
    USED.append(lg)
    for l in open(lg):
        m = re.search(r"COUPLING (\S+) x (\S+): point at ([\d.]+) Mb \(95% ([\d.]+)-([\d.]+)\) x ([\d.]+) Mb "
                      r"\(95% ([\d.]+)-([\d.]+)\)", l)
        if m:
            for cp in COUPLES:
                if (cp["p"], cp["q"]) == (m.group(1), m.group(2)):
                    cp["ci_p"] = (float(m.group(4)) * MB, float(m.group(5)) * MB)
                    cp["ci_q"] = (float(m.group(7)) * MB, float(m.group(8)) * MB)
for cp in COUPLES:            # piece coordinates -> copy coordinates
    for k in ("p", "q"):
        src, off = pieces.get(cp[k], (cp[k], 0))
        cp[k + "_copy"] = src
        cp[k + "p_copy"] = off + cp[k + "p"]
        ci = cp["ci_" + k]
        cp["ci_%s_copy" % k] = (off + ci[0], off + ci[1]) if ci else None
    if cp["ci_p"] is None:
        NOTES.append("no 95%% interval for the coupling %s x %s (not in %s): the point only" % (cp["p"], cp["q"], lg))
sc = os.path.join(A.coseg, "hic_scan_%s.tsv" % S)
if os.path.exists(sc):
    USED.append(sc)
    SCAN = list(csv.DictReader(open(sc), delimiter="\t"))
BP = {}
bpf = "qc/reference/%s/breakpoints.csv" % S
if os.path.exists(bpf):
    USED.append(bpf)
    for r in csv.DictReader(open(bpf)):
        BP.setdefault(r["chrom"], []).append(r)
for j in JOINS:
    j["gap"] = None
    for r in BP.get(j["copy"], []):
        if abs(int(r["cut"]) - j["pos"]) < 2 * MB:
            a = int(r["left_last_base"]) if (r.get("left_last_base") or "").strip() else None
            b = int(r["right_first_base"]) if (r.get("right_first_base") or "").strip() else None
            j["gap"] = (a if a is not None else int(r["cut"]), b if b is not None else int(r["cut"]))
JOINS.sort(key=lambda j: j["phi"])

# hap1 <-> hap2 positions through the synteny blocks of the haplotype alignment
SYN = None
for r in rows.values():
    if r["species"] == SPECIES and r["haplotype"] == "hap1":
        p = "qc/translocations/%s/%s_synteny_blocks.tsv" % (r["sample_id"], r["sample_id"])
        if os.path.exists(p):
            SYN = pd.read_csv(p, sep="\t")
            USED.append(p)
        break


def map_pos(c_from, x, c_to):
    """Position on c_to aligned to x on c_from (None without synteny blocks between them)."""
    if SYN is None or hap(c_from) == hap(c_to):
        return None
    if hap(c_from) == "1":
        b = SYN[(SYN.hap1 == c_from) & (SYN.hap2 == c_to)]
        s0, s1, d0, d1 = "t_start", "t_end", "q_start", "q_end"
    else:
        b = SYN[(SYN.hap2 == c_from) & (SYN.hap1 == c_to)]
        s0, s1, d0, d1 = "q_start", "q_end", "t_start", "t_end"
    if not len(b):
        return None
    inside = b[(b[s0] <= x) & (b[s1] >= x)]
    if not len(inside):
        dist = np.minimum((b[s0] - x).abs(), (b[s1] - x).abs())
        inside = b.loc[[dist.idxmin()]]
    r = inside.loc[inside["aligned"].idxmax()] if "aligned" in inside else inside.iloc[0]
    f_ = min(max((x - r[s0]) / max(r[s1] - r[s0], 1), 0), 1)
    return r[d0] + f_ * (r[d1] - r[d0]) if r["strand"] == "+" else r[d1] - f_ * (r[d1] - r[d0])


# coupling intervals on all four copies, and the families inside them
for cp in COUPLES:
    cp["spans"] = {}
    for k in ("p", "q"):
        c = cp[k + "_copy"]
        ci = cp["ci_%s_copy" % k] or (cp[k + "p_copy"] - 10 * MB, cp[k + "p_copy"] + 10 * MB)
        cp["spans"][c] = (ci[0], ci[1], cp[k + "p_copy"], "pollen")
        other = re.sub(r"_hap\d$", "_hap%s" % ("2" if hap(c) == "1" else "1"), c)
        if other in CIDX:
            ys = [map_pos(c, v, other) for v in (ci[0], ci[1], cp[k + "p_copy"])]
            if None not in ys:
                cp["spans"][other] = (min(ys[0], ys[1]), max(ys[0], ys[1]), ys[2], "mapped")
    inside = {}
    for c, (a, b, _, _) in cp["spans"].items():
        d = BYC[c]
        d = d[(d.f >= 0) & (d.e > a) & (d.s < b)]
        for fi, v in (np.minimum(d.e, b) - np.maximum(d.s, a)).groupby(d.f).sum().items():
            inside.setdefault(FAMS[fi], {})[c] = int(v)
    cp["families"] = inside
    # a family 'shared' by the two intervals: >= 1 kb inside the interval on each chromosome; its enrichment =
    # (share of the family's bp inside the intervals) / (share of the genome inside the intervals)
    span_bp = sum(b - a for a, b, _, _ in cp["spans"].values())
    cp["genome_share"] = span_bp / sum(L2n[c] for c in COPIES)
    ftot = SAT.assign(l=SAT.e - SAT.s).groupby("f").l.sum()
    cp["shared"] = []
    for f, v in inside.items():
        bp_p = sum(x for c, x in v.items() if num(c) == num(cp["p_copy"]))
        bp_q = sum(x for c, x in v.items() if num(c) == num(cp["q_copy"]))
        if bp_p >= 1000 and bp_q >= 1000:
            tot_f = float(ftot.get(fams[f], 0))
            share = sum(v.values()) / tot_f if tot_f else 0
            cp["shared"].append(dict(family=f, bp_p=bp_p, bp_q=bp_q, per_copy=v, total=tot_f, share=share,
                                     enrich=share / cp["genome_share"] if cp["genome_share"] else 0))
    cp["shared"].sort(key=lambda r: -min(r["bp_p"], r["bp_q"]))
    cp["scan"] = [d_ for d_ in SCAN if {num(d_["a"]), num(d_["b"])} == {num(cp["p_copy"]), num(cp["q_copy"])}][:3]


# ------------------------------------------------------------------ tables
def near_join(c, a, b, reach=2 * MB):
    out = []
    for j in JOINS:
        if j["copy"] == c and a - reach <= j["pos"] <= b + reach:
            out.append("%s%s %.1f Mb" % (j["copy"], " (%s)" % j["letter"] if j["letter"] else "", j["pos"] / MB))
    return out


for r in RLOCI:
    r["at_join"] = "; ".join(near_join(r["copy"], r["start"], r["end"]))
    js = [abs((r["start"] + r["end"]) / 2 - j["pos"]) / MB for j in JOINS if j["copy"] == r["copy"]]
    r["Mb_to_nearest_join"] = round(min(js), 2) if js else ""
RL = pd.DataFrame(RLOCI, columns=["type", "copy", "start", "end", "kb", "pieces", "at_join", "Mb_to_nearest_join"])
RL["start"] = RL["start"] + 1                          # 1-based, inclusive, as in GFF3
RL.to_csv(os.path.join(A.out, "rdna_loci_%s.csv" % S), index=False, float_format="%.1f")
if len(RL):
    RL["_n"], RL["_h"] = RL["copy"].map(num), RL["copy"].str[-1]
    RL = RL.sort_values(["type", "_n", "_h", "start"]).drop(columns=["_n", "_h"])
    RL.to_csv(os.path.join(A.out, "rdna_loci_%s.csv" % S), index=False, float_format="%.1f")
if len(FT):
    FT["at_join"] = ["; ".join(sorted(set(sum((near_join(c, a, b, 1 * MB) for c, a, b in ARRAYS[f]), []))))
                     for f in FT.family]
    FT["in_coupling_interval"] = ["; ".join("%s %.0f kb" % (c, v / 1e3) for cp in COUPLES
                                            for c, v in cp["families"].get(f, {}).items() if v >= 1000)
                                  for f in FT.family]
    FT["fish_score"] = (FT.share_on_top * np.log10(1 + FT.largest_array_kb)).round(2)
    rk = {"45S rDNA": 0, "5S rDNA": 1}
    FT["_r"] = FT.kind.map(lambda k: rk.get(k, 2))
    FT = FT.sort_values(["_r", "fish_score"], ascending=[True, False]).drop(columns="_r")
    FT.to_csv(os.path.join(A.out, "fish_probes_%s.csv" % S), index=False)


# ------------------------------------------------------------------ figure helpers
def stacked(ax, c, a=0, b=None, B=BIN, cov=None):
    """Stacked share of each repeat class per bin over [a, b) of copy c."""
    b = L2n[c] if b is None else b
    if cov is None:
        d = BYC[c]
        d = d[(d.e > a) & (d.s < b)]
        nb = -(-(b - a) // B)
        cov = np.stack([bin_cov(np.clip(d.s.values[d.k.values == k] - a, 0, None),
                                np.clip(d.e.values[d.k.values == k] - a, None, b - a), nb, B)
                        for k in range(len(CLASSES))], axis=1)
        x0 = a
    else:
        x0 = 0
    nb = cov.shape[0]
    x = (x0 + (np.arange(nb) + 0.5) * B) / MB
    den = np.full(nb, float(B))
    den[-1] = max((b - x0) - (nb - 1) * B, 1)             # the last bin is shorter
    ax.stackplot(x, (cov / den[:, None]).T, colors=CCOL, linewidth=0)
    return cov


def rdna_marks(ax, c, y, a=0, b=None, size=7):
    b = L2n[c] if b is None else b
    for r in RLOCI:
        if r["copy"] == c and r["end"] > a and r["start"] < b and r["kb"] >= 10:
            ax.plot((r["start"] + r["end"]) / 2 / MB, y, marker="v", ms=size, color=CCOL[I45 if r["type"] == "45S"
                                                                                        else I5],
                    mec="white", mew=0.5, clip_on=False, zorder=5)


def family_ticks(ax, c, y0, h, a=0, b=None, fams=None, minbp=20_000):
    b = L2n[c] if b is None else b
    for f in (fams or TOPF):
        for c2, s_, e_ in ARRAYS.get(f, []):
            if c2 == c and e_ > a and s_ < b and e_ - s_ >= minbp:
                w = max((e_ - s_) / MB, (b - a) / MB / 400)
                ax.add_patch(Rectangle((s_ / MB, y0), w, h, color=FCOL.get(f, "#534AB7"), lw=0, clip_on=False))


def arm_strip(ax, c, y0, h):
    js = [j for j in JOINS if j["copy"] == c]
    if not js:
        return
    j = js[0]
    for a, b, arm in ((0, j["pos"], j["left"]), (j["pos"], L2n[c], j["right"])):
        ax.add_patch(Rectangle((a / MB, y0), (b - a) / MB, h, color=ARM_COL.get(arm, "#EFEDE6"), lw=0, clip_on=False))
        ax.text((a + b) / 2 / MB, y0 + h / 2, arm, ha="center", va="center", fontsize=6, color=INK, clip_on=False)


def famlabel(f):
    m = MONO.get(f)
    tot = FT.set_index("family").total_kb.get(f, 0) if len(FT) else 0
    return "%s%s, %.1f Mb" % (f, " (%d bp)" % m if m else "", tot / 1e3)


class_handles = [Patch(color=CCOL[k], label=CLASSES[k]) for k in range(len(CLASSES) - 1, -1, -1)]

# ------------------------------------------------------------------ figure 1: the landscape
H = 0.62 * len(COPIES) + 1.9
fig, axes = plt.subplots(len(COPIES), 1, figsize=(16, H), sharex=True)
axes = np.atleast_1d(axes)
xmax = max(L2n[c] for c in COPIES) / MB
for ax, c in zip(axes, COPIES):
    stacked(ax, c, cov=COV[c])
    ax.set_ylim(0, 1)
    ax.set_xlim(0, xmax)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["0", "1"], fontsize=6)
    lab = c + ("  (%s)" % next((j["letter"] for j in JOINS if j["copy"] == c and j["letter"]), "")
               if any(j["copy"] == c and j["letter"] for j in JOINS) else "")
    ax.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=8)
    ax.axvline(L2n[c] / MB, color=MUTED, lw=0.6)
    rdna_marks(ax, c, 1.08)
    family_ticks(ax, c, -0.30, 0.18)
    arm_strip(ax, c, -0.55, 0.2)
    for j in JOINS:
        if j["copy"] == c:
            ax.axvline(j["pos"] / MB, color=INK if j["phi"] >= 0.5 else "#A32D2D", lw=1.0, ls="--")
    for cp in COUPLES:
        sp_ = cp["spans"].get(c)
        if sp_:
            ax.axvspan(sp_[0] / MB, sp_[1] / MB, color="#534AB7", alpha=0.12, lw=0)
            ax.axvline(sp_[2] / MB, color="#534AB7", lw=1.2, ls="-" if sp_[3] == "pollen" else ":")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
axes[-1].set_xlabel("Mb along each chromosome copy")
fh = [Line2D([], [], marker="v", ls="", color=CCOL[I45], label="45S rDNA array >= 10 kb"),
      Line2D([], [], marker="v", ls="", color=CCOL[I5], label="5S rDNA array >= 10 kb")]
fh += [Patch(color=FCOL[f], label=famlabel(f)) for f in TOPF]
leg1 = fig.legend(handles=class_handles, loc="upper left", bbox_to_anchor=(0.005, 1 - 0.75 / H), ncol=5, fontsize=7,
                  frameon=False, title="share of each %g Mb that is (stacked):" % A.bin_mb, title_fontsize=7.5)
fig.legend(handles=fh, loc="upper left", bbox_to_anchor=(0.50, 1 - 0.75 / H), ncol=3, fontsize=7, frameon=False,
           title="markers above / ticks below each copy (largest tandem-repeat families, arrays >= 20 kb):",
           title_fontsize=7.5)
fig.suptitle("%s (%s): where the repeats are on every chromosome copy (CARP%s)" % (
    SPECIES.replace("_", " "), S, " " + carp_version if carp_version else ""), x=0.01, y=1 - 0.12 / H, ha="left",
    fontsize=11)
what = ("What to look for: rDNA (triangles) and big satellite arrays (ticks) at the joins (dashed; red = the pollen "
        "reject the join) and inside the coupling interval (purple band; solid line = the coupling point, dotted = "
        "the same point on the other haplotype, through the alignment). The strip under a copy names its arms."
        if JOINS or COUPLES else
        "What to look for: where the 45S and 5S rDNA (triangles) and the big satellite families (ticks) sit on each "
        "chromosome (FISH landmarks), and whether both haplotypes of a chromosome carry them at the same place.")
fig.text(0.01, 0.12 / H, what, fontsize=7.5, color=MUTED, va="bottom", wrap=True)
fig.subplots_adjust(left=0.09, right=0.99, top=1 - 1.75 / H, bottom=0.75 / H, hspace=0.75)
p1 = os.path.join(A.out, "repeat_landscape_%s.png" % S)
fig.savefig(p1, dpi=130)
plt.close(fig)
OUTS = [p1]

# ------------------------------------------------------------------ figure 2: the joins up close
def lanes(ax, c, a, b, fams, cols, minw):
    """One row per family: its arrays on copy c within [a, b)."""
    ax.set_xlim(a / MB, b / MB)
    ax.set_ylim(max(len(fams), 1), 0)
    ax.set_yticks([])
    for k_, f in enumerate(fams):
        for c2, s_, e_ in ARRAYS.get(f, []):
            if c2 == c and e_ > a and s_ < b:
                ax.add_patch(Rectangle((s_ / MB, k_ + 0.15), max((e_ - s_) / MB, minw), 0.7, color=cols[f], lw=0))
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=7)


def gap_classes(c, g0, g1):
    """Share of [g0, g1) covered by each class (union of each class's features)."""
    d = BYC[c]
    d = d[(d.e > g0) & (d.s < g1)]
    out = {}
    for k, gk in d.groupby("k"):
        ms, me, _ = merge(np.clip(gk.s.values, g0, g1), np.clip(gk.e.values, g0, g1), gap=0)
        out[k] = sum(e_ - s_ for s_, e_ in zip(ms, me)) / max(g1 - g0, 1)
    return out


Z = int(A.zoom_mb * MB)
ZB = 50_000
if JOINS:
    import textwrap
    PER = 2.05
    H2 = PER * len(JOINS) + 1.35
    fig = plt.figure(figsize=(16, H2))
    for i, j in enumerate(JOINS):
        c = j["copy"]
        a, b = max(0, int(j["pos"] - Z)), min(L2n[c], int(j["pos"] + Z))
        top = H2 - 1.0 - i * PER                    # inches from the bottom
        ax = fig.add_axes([0.06, (top - 1.05) / H2, 0.58, 0.8 / H2])
        stacked(ax, c, a, b, ZB)
        ax.set_xlim(a / MB, b / MB)
        ax.set_ylim(0, 1)
        ax.set_xticklabels([])
        ax.set_ylabel("share of\n50 kb", fontsize=7)
        bad = j["phi"] < 0.5
        for x_ in (ax,):
            if j["gap"] and j["gap"][1] > j["gap"][0]:
                x_.axvspan(j["gap"][0] / MB, j["gap"][1] / MB, color="#2C2C2A", alpha=0.13, lw=0)
            x_.axvline(j["pos"] / MB, color="#A32D2D" if bad else INK, lw=1.4, ls="--")
        rdna_marks(ax, c, 1.09, a, b, size=8)
        d = BYC[c]
        d = d[(d.f >= 0) & (d.e > a) & (d.s < b)]
        near = (np.minimum(d.e, b) - np.maximum(d.s, a)).groupby(d.f).sum().sort_values(ascending=False)
        nf = [FAMS[fi] for fi in near.index[:4]]
        cols = {f: FCOL.get(f, FAM_COL[(k_ + 5) % len(FAM_COL)]) for k_, f in enumerate(nf)}
        lx = fig.add_axes([0.06, (top - 1.05 - 0.42) / H2, 0.58, 0.38 / H2])
        lanes(lx, c, a, b, nf, cols, 0.004)
        if j["gap"] and j["gap"][1] > j["gap"][0]:
            lx.axvspan(j["gap"][0] / MB, j["gap"][1] / MB, color="#2C2C2A", alpha=0.13, lw=0)
        lx.axvline(j["pos"] / MB, color="#A32D2D" if bad else INK, lw=1.0, ls="--")
        lx.set_xlabel("Mb on %s" % c, fontsize=7, labelpad=1)
        for k_, f in enumerate(nf):
            lx.text(b / MB, k_ + 0.5, "  " + famlabel(f), fontsize=6.5, va="center", color=cols[f], clip_on=False)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(labelsize=7)
        fig.text(0.06, (top - 0.18) / H2, "%s%s: %s | %s at %.1f Mb -- pollen phi %+.2f (%+.2f to %+.2f): %s" % (
            c, " (%s)" % j["letter"] if j["letter"] else "", j["left"], j["right"], j["pos"] / MB, j["phi"], j["lo"],
            j["hi"], "NOT joined in the pollen" if bad else "joined in the pollen"), fontsize=9,
            color="#A32D2D" if bad else INK, va="bottom")
        lines = []
        for t in ("45S", "5S"):
            dl = [(abs((r["start"] + r["end"]) / 2 - j["pos"]), r) for r in RLOCI if r["copy"] == c and r["type"] == t]
            if dl:
                dd, r = min(dl, key=lambda x: x[0])
                lines.append("nearest %s rDNA: %.2f Mb from the join (%.0f kb at %.1f Mb)" % (
                    t, dd / MB, r["kb"], r["start"] / MB))
            else:
                lines.append("no %s rDNA on %s" % (t, c))
        if j["gap"] and j["gap"][1] > j["gap"][0]:
            g0, g1 = j["gap"]
            gc = gap_classes(c, g0, g1)
            lines += textwrap.wrap("the haplotypes stop aligning over %.0f kb (grey): %s" % (
                (g1 - g0) / 1e3, ", ".join("%s %.0f%%" % (CLASSES[k], 100 * v) for k, v in
                                          sorted(gc.items(), key=lambda kv: -kv[1]) if v >= 0.02) or "no repeat"),
                74, subsequent_indent="   ")
        for part in j["detail"].split(" | "):
            lines += textwrap.wrap(part, 74, subsequent_indent="   ")
        fig.text(0.735, (top - 0.25) / H2, "\n".join(lines[:10]), fontsize=6.9, va="top", color=INK)
    fig.legend(handles=class_handles, loc="lower left", bbox_to_anchor=(0.06, 0.05 / H2), ncol=10, fontsize=7,
               frameon=False)
    fig.suptitle("%s (%s): the repeats at every join of the two-haplotype assembly, +-%g Mb" % (
        SPECIES.replace("_", " "), S, A.zoom_mb), x=0.01, y=1 - 0.12 / H2, ha="left", fontsize=11)
    fig.text(0.01, 1 - 0.42 / H2, "\n".join(textwrap.wrap(
        "Dashed line = the join (red: the pollen inherit its two sides independently). Grey band = where the two "
        "haplotypes stop aligning (breakpoints.csv). Triangles = rDNA arrays; the rows under each plot = the four "
        "tandem-repeat families with most bp in the window. What to look for: rDNA or one satellite family at the "
        "joins the pollen reject (A, D) that is absent from the joins they confirm on the same arms (C, B).", 230)),
        fontsize=7.5, color=MUTED, va="top")
    p2 = os.path.join(A.out, "repeat_joins_%s.png" % S)
    fig.savefig(p2, dpi=130)
    plt.close(fig)
    OUTS.append(p2)

# ------------------------------------------------------------------ figure 3: the coupling
for n_cp, cp in enumerate(COUPLES):
    import textwrap
    cps = sorted({c for c in COPIES if num(c) in (num(cp["p_copy"]), num(cp["q_copy"]))}, key=lambda c: (num(c), c))
    sh = [r["family"] for r in cp["shared"][:8]]
    scol = {f: FCOL.get(f, FAM_COL[(i + 3) % len(FAM_COL)]) for i, f in enumerate(sh)}
    PER = 1.25 + 0.11 * max(len(sh), 1)
    H3 = PER * len(cps) + 1.45
    fig = plt.figure(figsize=(16, H3))
    xm = max(L2n[c] for c in cps) / MB
    for i, c in enumerate(cps):
        top = H3 - 1.0 - i * PER
        ax = fig.add_axes([0.08, (top - 0.85) / H3, 0.58, 0.8 / H3])
        stacked(ax, c, cov=COV[c])
        ax.set_xlim(0, xm)
        ax.set_ylim(0, 1)
        ax.set_xticklabels([])
        ax.set_ylabel(c, rotation=0, ha="right", va="center", fontsize=9)
        rdna_marks(ax, c, 1.09)
        lx = fig.add_axes([0.08, (top - 0.87 - 0.11 * max(len(sh), 1)) / H3, 0.58, 0.11 * max(len(sh), 1) / H3])
        lanes(lx, c, 0, int(xm * MB), sh, scol, xm / 400)
        lx.set_xlim(0, xm)
        for k_, f in enumerate(sh):
            lx.text(xm, k_ + 0.5, " " + f, fontsize=5.5, va="center", color=scol[f], clip_on=False)
        sp_ = cp["spans"].get(c)
        for x_ in (ax, lx):
            if sp_:
                x_.axvspan(sp_[0] / MB, sp_[1] / MB, color="#534AB7", alpha=0.12, lw=0)
                x_.axvline(sp_[2] / MB, color="#534AB7", lw=1.6, ls="-" if sp_[3] == "pollen" else ":")
        if i < len(cps) - 1:
            lx.set_xticklabels([])
        k_ = 0
        for d_ in cp["scan"]:
            for cc, pos in ((d_["a"], d_["pa"]), (d_["b"], d_["pb"])):
                if cc == c:
                    ax.plot(float(pos), 0.93, marker="X", ms=8, color=INK, zorder=6)
                    ax.text(float(pos) + xm / 150, 0.93 - 0.17 * (k_ % 3), "Hi-C %.1fx (#%s, %s x %s)" % (
                        float(d_["x"]), d_["rank"], d_["a"], d_["b"]), fontsize=6.3, va="center", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
                    k_ += 1
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(labelsize=7)
    lx.set_xlabel("Mb along each copy", fontsize=8)
    yy = H3 - 1.0
    fig.text(0.70, yy / H3, "Tandem-repeat families inside the interval on BOTH chromosomes\n(>= 1 kb on each; the "
             "rows under the copies, same colours). 'xN' = their share of\nthe family's bp over the intervals' share "
             "of the genome (%.1f%%):" % (100 * cp["genome_share"]), fontsize=7.3, va="top", color=INK)
    yy -= 0.55
    for r in cp["shared"][:8]:
        fig.text(0.70, yy / H3, "%s  [%s]   %.0f%% of its bp inside, x%.1f" % (
            famlabel(r["family"]), FT.set_index("family").kind.get(r["family"], "") if len(FT) else "",
            100 * r["share"], r["enrich"]), fontsize=7.3, va="top", color=scol[r["family"]], weight="bold")
        fig.text(0.70, (yy - 0.14) / H3, "\n".join(textwrap.wrap("inside: " + ", ".join(
            "%s %.0f kb" % (c, v / 1e3) for c, v in sorted(r["per_copy"].items(),
                                                        key=lambda kv: (num(kv[0]), kv[0]))), 80,
            subsequent_indent="   ")), fontsize=7, va="top", color=INK)
        yy -= 0.42
    if not cp["shared"]:
        fig.text(0.70, yy / H3, "none", fontsize=7.3, va="top", color=INK)
    ci_txt = ("%s %.0f Mb (95%% %.0f-%.0f) x %s %.0f Mb (95%% %.0f-%.0f)" % (
        cp["p_copy"], cp["pp_copy"] / MB, cp["ci_p_copy"][0] / MB, cp["ci_p_copy"][1] / MB, cp["q_copy"],
        cp["qp_copy"] / MB, cp["ci_q_copy"][0] / MB, cp["ci_q_copy"][1] / MB) if cp["ci_p_copy"] else
              "%s %.0f Mb x %s %.0f Mb" % (cp["p_copy"], cp["pp_copy"] / MB, cp["q_copy"], cp["qp_copy"] / MB))
    fig.suptitle("%s (%s): the repeats where the pollen couple chr%d and chr%d: coupling point %s" % (
        SPECIES.replace("_", " "), S, num(cp["p_copy"]), num(cp["q_copy"]), ci_txt), x=0.01, y=1 - 0.12 / H3,
        ha="left", fontsize=10.5)
    fig.text(0.01, 1 - 0.42 / H3, "\n".join(textwrap.wrap(
        "Purple band = the coupling point's 95% interval from the pollen (solid line = the point; on the other "
        "haplotype the same stretch through the alignment, dotted). X = the strongest Hi-C blocks between copies of "
        "the two chromosomes (hic_scan, x the median between chromosomes). What to look for: rDNA or a tandem-repeat "
        "family inside BOTH bands and enriched there (a shared repeat that could pair or exchange), or none (then "
        "the repeats do not explain the coupling).", 230)), fontsize=7.5, color=MUTED, va="top")
    fig.legend(handles=class_handles, loc="lower left", bbox_to_anchor=(0.08, 0.03 / H3), ncol=10, fontsize=7,
               frameon=False)
    p3 = os.path.join(A.out, "repeat_coupling_%s%s.png" % (S, "" if n_cp == 0 else "_%d" % (n_cp + 1)))
    fig.savefig(p3, dpi=130)
    plt.close(fig)
    OUTS.append(p3)

# ------------------------------------------------------------------ report
tot = sum(L2n[c] for c in COPIES)
cls_bp = sum(COV[c].sum(axis=0) for c in COPIES)
print("REPEATS (top-level features on the %d chromosome copies, %.2f Gb): %s" % (
    len(COPIES), tot / 1e9, ", ".join("%s %.1f%%" % (CLASSES[k], 100 * cls_bp[k] / tot)
                                      for k in np.argsort(-cls_bp) if cls_bp[k] > 0)))
for t in ("45S", "5S"):
    rl = [r for r in RLOCI if r["type"] == t and r["kb"] >= 10]
    print("%s rDNA arrays >= 10 kb: %d%s" % (t, len(rl), "" if not rl else ": " + "; ".join(
        "%s %.1f Mb (%.0f kb)" % (r["copy"], r["start"] / MB, r["kb"]) for r in rl[:16])))
if JOINS:
    print("JOINS (pollen phi; nearest rDNA on the same copy; families with most bp within %g Mb):" % A.zoom_mb)
    for j in JOINS:
        c = j["copy"]
        dl = {t: min([abs((r["start"] + r["end"]) / 2 - j["pos"]) / MB for r in RLOCI if r["copy"] == c
                      and r["type"] == t] or [float("nan")]) for t in ("45S", "5S")}
        d = BYC[c]
        d = d[(d.f >= 0) & (d.e > j["pos"] - Z) & (d.s < j["pos"] + Z)]
        near = (d.e - d.s).groupby(d.f).sum().sort_values(ascending=False)
        print("  %-32s phi %+.2f   45S %s, 5S %s   %s" % (
            "%s%s %s|%s %.1f Mb" % (c, " (%s)" % j["letter"] if j["letter"] else "", j["left"], j["right"],
                                    j["pos"] / MB), j["phi"],
            "%.2f Mb" % dl["45S"] if dl["45S"] == dl["45S"] else "none", "%.2f Mb" % dl["5S"] if dl["5S"] == dl["5S"]
            else "none", ", ".join("%s %.0f kb" % (FAMS[fi], v / 1e3) for fi, v in list(near.items())[:3])))
for cp in COUPLES:
    print("COUPLING %s x %s: families inside both intervals (the intervals hold %.1f%% of the genome): %s" % (
        cp["p_copy"], cp["q_copy"], 100 * cp["genome_share"], "; ".join(
            "%s %.0f kb + %.0f kb = %.0f%% of the family (x%.1f)" % (r["family"], r["bp_p"] / 1e3, r["bp_q"] / 1e3,
                                                                   100 * r["share"], r["enrich"])
            for r in cp["shared"][:6]) or "none"))
if len(FT):
    print("FISH CANDIDATES (tandem-repeat families: share on one chromosome x size of the largest array):")
    for r in FT.head(12).itertuples():
        print("  %-10s %-9s %6s bp  %8.0f kb in %3d arrays (largest %6.0f kb)  %.0f%% on %s%s%s" % (
            r.family, r.kind, r.monomer_bp if r.monomer_bp != "" else "?", r.total_kb, r.arrays, r.largest_array_kb,
            100 * r.share_on_top, r.top_chromosome, "  at join " + r.at_join if r.at_join else "",
            "  in coupling " + r.in_coupling_interval if r.in_coupling_interval else ""))
for n_ in NOTES:
    print("note:", n_)
print("files read:\n  " + "\n  ".join(USED))
print("wrote " + ", ".join(OUTS + [os.path.join(A.out, "%s_%s.csv" % (t, S)) for t in ("rdna_loci", "fish_probes")]))
