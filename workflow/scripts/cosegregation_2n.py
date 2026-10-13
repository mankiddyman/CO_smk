#!/usr/bin/env python3
"""cosegregation_2n.py (version 2) -- what the pollen inherit together, drawn on both haplotypes like a Hi-C map,
and the assembly's joins tested against the pollen.

THE MAP. For every pair of WIN_MB windows of the crossover reference, phi = the correlation of the two windows'
genotypes across pollen nuclei (+1 = always inherited together, 0 = independently, -1 = never together). In each
window a nucleus carries either the reference's allele (REF) or the other copy's (ALT). The REF copy is the reference
sequence itself; the ALT copy is its homolog in the other haplotype or, for an arm whose other copy sits in the same
haplotype, that copy. So every window sits twice on the two-haplotype map: on its REF copy, and on its ALT copy with
the sign of phi flipped. The two copies of a chromosome come out at -1 (a nucleus has one or the other), each copy is
red along its own diagonal, and copies of different chromosomes that travel together through meiosis share a red
block. Axes: every chromosome copy, hap1 then hap2 for each number. Coloured strips: the arms named by the pieces of
a cut reference (D. paradoxa: L1, L2, P, Q of chr1/chr2 in the supervisor report's colours; X and the core of chr4);
A-D name the chr1/chr2 copies by their arms (A = L1+P, B = L2+P, C = L1+Q, D = L2+Q). Beside it, the tissue Hi-C on
the same axes (contacts_all.npz).

THE KARYOTYPE TEST (cut references only). A join is a point where one arm meets another on a chromosome copy of the
two-haplotype assembly. If the join is real, the 20 Mb on either side sit on one chromosome and the pollen inherit
them together (phi near +1, lowered only by crossovers between them); if the two arms are on different chromosomes,
phi is near 0. Per join: phi with a 95% interval (nuclei resampled 1,000 times) and, for a join the pollen do not
confirm, the arm end the pollen do inherit with each side (its pollen partner). Pieces that no join connects but that
the pollen inherit together (median phi over their window pairs beyond COUPLE_MIN) are fitted with one coupling
point: phi(w1, w2) = c * phi(w1, a1) * phi(w2, a2), complete co-inheritance at a1 and a2 diluted by each chromosome's
own recombination between a window and its point. A junction between the two predicts this, as does anything else
that holds a1 and a2 together (a selected gene pair, a pairing or segregation determinant); a fusion is the special
case with a1 and a2 at chromosome ends. Reported: the best (a1, a2) by least squares, its 95% interval from resamples
of the nuclei, and R2 = share of the block's variance explained. Control: the same prediction for a between-
chromosome join the assembly has and the pollen confirm, made at the join itself.

HI-C CHECKS (with the maps). Contacts per Mb2 between the 10 Mb on either side of each join, and between each side
of a join the pollen reject and its pollen partner, against the median between chromosome copies of different
numbers ("x trans"; uniquely mapped reads when contacts_unique.npz exists). Cis neighbours are many times trans,
unrelated chromosomes about 1x. And a scan of every pair of copies of different numbers for the strongest 10 x 10 Mb
block (2 Mb steps): a junction the assembly lacks shows there as a hotspot.

Usage: cosegregation_2n.py SAMPLE [--out reports/cosegregation_2n] [--win_mb 10] [--bin_mb 10]
                           [--phd_root ...] [--no_hic] [--boot 100]
Writes OUT/cosegregation_2n_SAMPLE.{png,tsv}; for a cut reference also OUT/karyotype_test_SAMPLE.{png,tsv}; with
Hi-C also OUT/hic_scan_SAMPLE.tsv. Run from the CO_smk root.
"""
import argparse
import csv
import math
import os
import re
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle

warnings.filterwarnings("ignore", category=RuntimeWarning)
ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("--out", default="reports/cosegregation_2n")
ap.add_argument("--win_mb", type=float, default=10)
ap.add_argument("--bin_mb", type=float, default=10)
ap.add_argument("--phd_root", default="/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd")
ap.add_argument("--no_hic", action="store_true")
ap.add_argument("--boot", type=int, default=100, help="resamples of the nuclei for the coupling point's interval")
A = ap.parse_args()
S, W, BIN = A.sample, int(A.win_mb * 1e6), int(A.bin_mb * 1e6)
MIN_MOL, GT_MAJ, PHI_N = 3, 0.8, 20          # as in supervisor_report.py: >= 3 molecules, >= 80% one allele
END_K = 2                 # windows at a piece end that stand for one side of a join (20 Mb at 10 Mb windows)
CLEAR = 0.3               # a nucleus's call over those windows needs an ALT share below 0.2 or above 0.8
COUPLE_MIN = 0.3          # median phi over two unjoined pieces' window pairs above which they count as coupled
FLANK = 10_000_000        # bp either side of a join for the Hi-C check
SCAN_BIN, SCAN_BOX = 2_000_000, 5             # Hi-C scan: 2 Mb bins, 5 x 5 bins = 10 x 10 Mb blocks
ARM_COL = {"L1": "#5DCAA5", "L2": "#F0997B", "P": "#AFA9EC", "Q": "#FAC775",   # the report's arm colours
           "X": "#85B7EB", "core": "#D3D1C7"}
WHOLE_COL = "#EFEDE6"
LETTER = {frozenset(("L1", "P")): "A", frozenset(("L2", "P")): "B", frozenset(("L1", "Q")): "C",
          frozenset(("L2", "Q")): "D"}
RNG = np.random.default_rng(20261013)
os.makedirs(A.out, exist_ok=True)
USED, NOTES = [], []


def natural(c):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", c)]


def fai_lens(p):
    USED.append(p)
    return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(p) if l.strip()}


def num(c):
    m = re.match(r"chr(\d+)", c)
    return int(m.group(1)) if m else 10 ** 6


def hap(c):
    m = re.search(r"_hap(\d)$", c)
    return m.group(1) if m else None


# ------------------------------------------------------------------ the reference and its two-haplotype source
row = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}[S]
L = fai_lens(row["assembly_fasta"] + ".fai")
pp = os.path.join(os.path.dirname(row["assembly_fasta"]), "pieces.tsv")
CUT = os.path.exists(pp)
if CUT:
    USED.append(pp)
    pcs = list(csv.DictReader(open(pp), delimiter="\t"))
    seqs = [r["piece"] for r in pcs]
    src = {r["piece"]: (r["source"], int(r["start"]) - 1) for r in pcs}
else:
    n = int(row["chr_number_2n"]) // 2
    seqs = sorted(sorted(L, key=lambda c: -L[c])[:n], key=natural)
    src = {c: (c, 0) for c in seqs}
Ls = {s: L[s] for s in seqs}
L2n = None
for man in (os.path.join(os.path.dirname(row["assembly_fasta"]), "MANIFEST.txt"),
            "/netscratch/dep_mercier/grp_marques/Aaryan/refs/%s_hap1/MANIFEST.txt" % row["species"]):
    if os.path.exists(man):
        USED.append(man)
        m = re.search(r"^source_fasta:\s*(\S+)", open(man).read(), re.M)
        if m and os.path.exists(m.group(1) + ".fai"):
            L2n = fai_lens(m.group(1) + ".fai")
            break
if L2n is None:
    raise SystemExit("no two-haplotype assembly index found via MANIFEST.txt for %s" % S)
COPIES = sorted([c for c in L2n if re.match(r"^chr\d+_hap[12]$", c)], key=lambda c: (num(c), c))

# ------------------------------------------------------------------ where each piece's other copy sits
how = {}
plan = "config/split_plans/%s.csv" % S
if os.path.exists(plan):
    USED.append(plan)
    how = {r["piece"]: r.get("how", "") for r in csv.DictReader(l for l in open(plan) if not l.startswith("#"))}
ARMS = {}           # copy -> [(start, end, homolog or 'none')]
bpf = "qc/reference/%s/breakpoints.csv" % S
if os.path.exists(bpf):
    USED.append(bpf)
    bp = {}
    for r in csv.DictReader(open(bpf)):
        bp.setdefault(r["chrom"], []).append((int(r["cut"]), r["left_homolog"], r["right_homolog"]))
    for c, cuts in bp.items():
        cuts.sort()
        edges = [0] + [x for x, _, _ in cuts] + [L2n.get(c, cuts[-1][0])]
        homs = [cuts[0][1]] + [rh for _, _, rh in cuts]
        ARMS[c] = [(edges[i], edges[i + 1], homs[i]) for i in range(len(homs))]
SYN, sid = None, None
for r in csv.DictReader(open("config/samples.csv")):
    if r["species"] == row["species"] and r["haplotype"] == "hap1":
        sid = r["sample_id"]
        break
syf = "qc/translocations/%s/%s_synteny_blocks.tsv" % (sid, sid)
if sid and os.path.exists(syf):
    USED.append(syf)
    SYN = pd.read_csv(syf, sep="\t")


def arms_of(c):
    return ARMS.get(c, [(0, L2n[c], None)])


def strand(a, b, a0, a1):
    """Main strand of the alignment between copies a and b over a's stretch a0-a1 (None: no alignment)."""
    if SYN is None or hap(a) == hap(b):
        return None
    x = SYN[(SYN.hap1 == a) & (SYN.hap2 == b)] if hap(a) == "1" else SYN[(SYN.hap2 == a) & (SYN.hap1 == b)]
    s0, s1 = (x.t_start, x.t_end) if hap(a) == "1" else (x.q_start, x.q_end)
    x = x[(s1 > a0) & (s0 < a1)]
    if not len(x):
        return None
    by = x.groupby("strand").aligned.sum()
    return by.idxmax()


def other_copy(p):
    """(copy, start, end, strand, how it was found) of the other copy of piece p."""
    s, off = src[p]
    note = how.get(p, "")
    m = re.search(r"homolog (chr\d+_hap\d)", note)
    if m:
        H = m.group(1)
    elif "homolog none" in note:
        cand = [c for c, arms in ARMS.items() if hap(c) == hap(s) and c != s and any(h == "none" for _, _, h in arms)]
        if len(cand) != 1:
            return None
        H = cand[0]
    else:
        H = re.sub(r"_hap\d$", "_hap%s" % ("2" if hap(s) == "1" else "1"), s)
    if H not in L2n:
        return None
    want = "none" if "homolog none" in note else s
    match = [(a, b) for a, b, h in arms_of(H) if h == want or h is None]
    if not match:
        return None
    a, b = min(x for x, _ in match), max(y for _, y in match)
    st = strand(s, H, off, off + Ls[p])
    return H, a, b, st or "+", "alignment" if st else ("same haplotype" if hap(H) == hap(s) else "no blocks")


OTHER = {p: other_copy(p) for p in seqs}
miss = [p for p, o in OTHER.items() if o is None]
same = [p for p, o in OTHER.items() if o and o[4] == "same haplotype"]
noal = [p for p, o in OTHER.items() if o and o[4] == "no blocks"]
if miss:
    NOTES.append("other copy not found for %s (drawn on their own copy only)" % ", ".join(miss))
if same:
    NOTES.append("other copy in the same haplotype, no alignment to it, same orientation assumed: %s" % "; ".join(
        "%s on %s %.1f-%.1f Mb" % (p, OTHER[p][0], OTHER[p][1] / 1e6, OTHER[p][2] / 1e6) for p in same))
if noal:
    NOTES.append("no synteny blocks for %s, same orientation assumed%s" % (
        "these" if len(noal) < len(seqs) else "any chromosome", (": " + ", ".join(noal)) if len(noal) < len(seqs) else ""))

# ------------------------------------------------------------------ window genotypes of every nucleus
cells = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()]
USED.append("results/cell_qc/%s/good_cells.tsv" % S)
masked = set()
mp = "results/crossovers/%s/masks.tsv" % S
if os.path.exists(mp):
    USED.append(mp)
    mk = pd.read_csv(mp, sep="\t", dtype={"barcode": str, "chrom": str})
    masked = set(zip(mk.barcode[mk.masked == "yes"], mk.chrom[mk.masked == "yes"]))
nwin = {s: max(1, int(round(Ls[s] / W))) for s in seqs}
wins = [(s, i) for s in seqs for i in range(nwin[s])]
widx = {w: k for k, w in enumerate(wins)}
G = np.full((len(cells), len(wins)), np.nan)
for ci, bc in enumerate(cells):
    p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        continue
    m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"], dtype={0: str})
    m = m[(m.rc != m.ac) & m.chrom.isin(Ls)].copy()
    if not len(m):
        continue
    m["alt"] = (m.ac > m.rc).astype(int)
    m["w"] = [min(int(q // W), nwin[c] - 1) for c, q in zip(m.chrom, m.pos)]
    for (c, w), r in m.groupby(["chrom", "w"]).alt.agg(["mean", "size"]).iterrows():
        if (bc, c) in masked or r["size"] < MIN_MOL:
            continue
        if r["mean"] >= GT_MAJ - 1e-9 or 1 - r["mean"] >= GT_MAJ - 1e-9:
            G[ci, widx[(c, w)]] = 1 if r["mean"] >= GT_MAJ - 1e-9 else 0
USED.append("results/cell_data_mol/%s/<barcode>.tsv (%d nuclei)" % (S, len(cells)))
G = G[~np.all(np.isnan(G), axis=1)]


def phis(Ga, Gb):
    """phi between the columns of Ga and of Gb over the nuclei called in both (NaN below PHI_N nuclei); also the
    discordant and total nuclei per pair."""
    X1a, X0a, X1b, X0b = (Ga == 1) * 1.0, (Ga == 0) * 1.0, (Gb == 1) * 1.0, (Gb == 0) * 1.0
    n11, n00, n10, n01 = X1a.T @ X1b, X0a.T @ X0b, X1a.T @ X0b, X0a.T @ X1b
    den = np.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
    ph = (n11 * n00 - n10 * n01) / den
    nn = n11 + n00 + n10 + n01
    ph[(den == 0) | (nn < PHI_N)] = np.nan
    return ph, n10 + n01, nn


PHI = phis(G, G)[0]
np.fill_diagonal(PHI, 1.0)

# ------------------------------------------------------------------ each window on both copies, binned like Hi-C
nb = {c: int(math.ceil(L2n[c] / BIN)) for c in COPIES}
boff, o = {}, 0
for c in COPIES:
    boff[c], o = o, o + nb[c]
NB = o
Amat = np.zeros((len(wins), NB))              # +1 where a window's REF copy sits, -1 where its ALT copy sits
for k, (s, i) in enumerate(wins):
    x = (i + 0.5) * W if i < nwin[s] - 1 else ((nwin[s] - 1) * W + Ls[s]) / 2.0     # window middle on the piece
    c, off = src[s]
    if c in boff:
        Amat[k, boff[c] + min(int((off + x) // BIN), nb[c] - 1)] += 1
    ot = OTHER.get(s)
    if ot:
        H, a, b, st, _ = ot
        f = x / Ls[s]
        y = a + f * (b - a) if st == "+" else b - f * (b - a)
        Amat[k, boff[H] + min(int(y // BIN), nb[H] - 1)] -= 1
F = np.isfinite(PHI).astype(float)
P0 = np.nan_to_num(PHI)
M = (Amat.T @ P0 @ Amat) / (np.abs(Amat).T @ F @ np.abs(Amat))
NPAIR = np.abs(Amat).T @ F @ np.abs(Amat)
M[NPAIR == 0] = np.nan
with open(os.path.join(A.out, "cosegregation_2n_%s.tsv" % S), "w") as fo:
    fo.write("copy_a\tMb_a\tcopy_b\tMb_b\tsigned_phi\twindow_pairs\n")
    lab = [(c, j) for c in COPIES for j in range(nb[c])]
    for i in range(NB):
        for j in range(i, NB):
            if NPAIR[i, j] > 0:
                fo.write("%s\t%d\t%s\t%d\t%.3f\t%d\n" % (lab[i][0], lab[i][1] * A.bin_mb, lab[j][0], lab[j][1] * A.bin_mb,
                                                         M[i, j], NPAIR[i, j]))

# ------------------------------------------------------------------ arms of every copy
def arm_name(p):
    s = src[p][0]
    return p[len(s) + 1:] if p.startswith(s + "_") and len(p) > len(s) + 1 else None


SEGS = {c: [] for c in COPIES}      # copy -> its pieces, REF (the piece itself) or ALT (where its other copy sits)
for p in seqs:
    c, off = src[p]
    if c in SEGS:
        SEGS[c].append(dict(a=off, b=off + Ls[p], arm=arm_name(p), piece=p, role="REF", strand="+"))
    o_ = OTHER.get(p)
    if o_ and o_[0] in SEGS:
        SEGS[o_[0]].append(dict(a=o_[1], b=o_[2], arm=arm_name(p), piece=p, role="ALT", strand=o_[3]))
for c in SEGS:
    SEGS[c].sort(key=lambda s_: s_["a"])
LET = {c: LETTER[frozenset(s_["arm"] for s_ in ss if s_["arm"])] for c, ss in SEGS.items()
       if frozenset(s_["arm"] for s_ in ss if s_["arm"]) in LETTER}
HAS_ARMS = any(s_["arm"] for ss in SEGS.values() for s_ in ss)


def cname(c):
    return "%s (%s)" % (c, LET[c]) if c in LET else c


def win_idx(p):
    return np.array([widx[(p, i)] for i in range(nwin[p])])


def wmid(p, i):
    """Middle of window i of piece p, Mb on the piece."""
    return ((i + 0.5) * W if i < nwin[p] - 1 else ((nwin[p] - 1) * W + Ls[p]) / 2.0) / 1e6


# ------------------------------------------------------------------ the karyotype test
def piece_end(seg, side):
    """Which end of the piece a segment's left ('a') or right ('b') side is."""
    fwd = seg["role"] == "REF" or seg["strand"] == "+"
    return ("start" if side == "a" else "end") if fwd else ("end" if side == "a" else "start")


def endcall(p, end, k=END_K, centre=None):
    """Per nucleus: 1 = ALT, 0 = REF over the k windows at one end of piece p (or around window `centre`), NaN if
    uncalled or mixed."""
    ix = win_idx(p)
    if centre is not None:
        ix = ix[max(0, centre - 1):centre + 2]
    else:
        ix = ix[:k] if end == "start" else ix[-k:]
    v = np.nanmean(G[:, ix], axis=1)
    out = np.full(len(G), np.nan)
    ok = np.isfinite(v) & (np.abs(v - 0.5) > CLEAR)
    out[ok] = (v[ok] > 0.5).astype(float)
    return out


def sgn(seg):
    """+1 if a nucleus carries this copy's allele when it is ALT for the piece, -1 when REF."""
    return -1.0 if seg["role"] == "REF" else 1.0


def phi_pair(x, y, idx=None):
    if idx is not None:
        x, y = x[idx], y[idx]
    ok = np.isfinite(x) & np.isfinite(y)
    X, Y = x[ok] > 0.5, y[ok] > 0.5
    n11, n00, n10, n01 = np.sum(X & Y), np.sum(~X & ~Y), np.sum(X & ~Y), np.sum(~X & Y)
    den = math.sqrt(float((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00)))
    return ((n11 * n00 - n10 * n01) / den if den else float("nan")), int(ok.sum()), int(n11 + n00)


def boot_ci(x, y, s, nb_=1000):
    v = [s * phi_pair(x, y, RNG.integers(0, len(x), len(x)))[0] for _ in range(nb_)]
    v = [t for t in v if np.isfinite(t)]
    return (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))) if v else (float("nan"),) * 2


def fit_point(Gs, p, q):
    """One coupling point between pieces p (rows) and q (columns): mean squared error of c * outer(phi(., a1),
    phi(., a2)) against the observed phi for every candidate (a1, a2), c by least squares."""
    ip, iq = win_idx(p), win_idx(q)
    O, disc, nn = phis(Gs[:, ip], Gs[:, iq])
    Pp, Pq = phis(Gs[:, ip], Gs[:, ip])[0], phis(Gs[:, iq], Gs[:, iq])[0]
    np.fill_diagonal(Pp, 1.0)
    np.fill_diagonal(Pq, 1.0)
    fin = np.isfinite(O)
    var = float(np.mean((O[fin] - O[fin].mean()) ** 2)) if fin.any() else float("nan")
    mse, cc = np.full(O.shape, np.nan), np.full(O.shape, np.nan)
    for i in range(len(ip)):
        for j in range(len(iq)):
            Pr = np.outer(Pp[:, i], Pq[:, j])
            m_ = fin & np.isfinite(Pr)
            if m_.sum() < 0.5 * fin.sum() or not (Pr[m_] ** 2).sum():
                continue
            c_ = (Pr[m_] * O[m_]).sum() / (Pr[m_] ** 2).sum()
            mse[i, j], cc[i, j] = np.mean((O[m_] - c_ * Pr[m_]) ** 2), c_
    return dict(O=O, Pp=Pp, Pq=Pq, mse=mse, c=cc, var=var, disc=disc, nn=nn)


def near_end(f_, a0, b0, reach=3):
    """The candidate point pair nearest (a0, b0), within `reach` windows of each, that has enough data to fit."""
    best = None
    for i in range(reach):
        for j in range(reach):
            a_, b_ = a0 + (i if a0 == 0 else -i), b0 + (j if b0 == 0 else -j)
            if 0 <= a_ < f_["mse"].shape[0] and 0 <= b_ < f_["mse"].shape[1] and np.isfinite(f_["mse"][a_, b_]):
                key = (i + j, f_["mse"][a_, b_])
                if best is None or key < best[0]:
                    best = (key, (a_, b_))
    return best[1] if best else None


JOINS, ENDS, COUPLES, CONTROL = [], [], [], None
if CUT and HAS_ARMS:
    for c in COPIES:
        ss = SEGS[c]
        for s1, s2 in zip(ss, ss[1:]):
            if abs(s2["a"] - s1["b"]) > 2e6:
                continue
            e1, e2 = (s1["piece"], piece_end(s1, "b")), (s2["piece"], piece_end(s2, "a"))
            x, y = endcall(*e1), endcall(*e2)
            sg = sgn(s1) * sgn(s2)
            ph, n_, agree = phi_pair(x, y)
            JOINS.append(dict(copy=c, pos=s1["b"], left=s1, right=s2, e1=e1, e2=e2, s1=sgn(s1), s2=sgn(s2),
                              phi=sg * ph, n=n_, agree=agree if sg > 0 else n_ - agree, ci=boot_ci(x, y, sg),
                              between=num(src[s1["piece"]][0]) != num(src[s2["piece"]][0])))
    for j in JOINS:
        ENDS.append(dict(join=j, side="left", seg=j["left"], e=j["e1"], s=j["s1"]))
        ENDS.append(dict(join=j, side="right", seg=j["right"], e=j["e2"], s=j["s2"]))

    def partner(end):
        x, best = endcall(*end["e"]), None
        for o_ in ENDS:
            if o_["join"] is end["join"] or o_["e"] == end["e"]:
                continue
            ph = end["s"] * o_["s"] * phi_pair(x, endcall(*o_["e"]))[0]
            if np.isfinite(ph) and (best is None or ph > best[1]):
                best = (o_, ph)
        return best
    for j in JOINS:
        j["partners"] = [partner(e) for e in ENDS if e["join"] is j] if j["phi"] < 0.5 else []
    # pieces no join connects that the pollen inherit together
    par = {p: p for p in seqs}

    def find(p):
        while par[p] != p:
            p = par[p]
        return p
    for j in JOINS:
        par[find(j["left"]["piece"])] = find(j["right"]["piece"])
    for i_, p in enumerate(seqs):
        for q in seqs[i_ + 1:]:
            if find(p) == find(q) or num(src[p][0]) == num(src[q][0]):
                continue
            med = float(np.nanmedian(PHI[np.ix_(win_idx(p), win_idx(q))]))
            if np.isfinite(med) and abs(med) > COUPLE_MIN:
                COUPLES.append(dict(p=p, q=q, med=med))
    COUPLES.sort(key=lambda d_: -abs(d_["med"]))
    for cp in COUPLES:
        f_ = fit_point(G, cp["p"], cp["q"])
        k = np.nanargmin(f_["mse"])
        i_, j_ = np.unravel_index(k, f_["mse"].shape)
        last_p, last_q = nwin[cp["p"]] - 1, nwin[cp["q"]] - 1
        corners = [near_end(f_, a_, b_) for a_ in (0, last_p) for b_ in (0, last_q)]
        corners = {k_: f_["mse"][k_] for k_ in corners if k_}
        ca = min(corners, key=corners.get) if corners else None
        bi, bj = [], []
        for _ in range(A.boot):
            g_ = fit_point(G[RNG.integers(0, len(G), len(G))], cp["p"], cp["q"])
            if np.isfinite(g_["mse"]).any():
                a_, b_ = np.unravel_index(np.nanargmin(g_["mse"]), g_["mse"].shape)
                bi.append(wmid(cp["p"], a_))
                bj.append(wmid(cp["q"], b_))
        x, y = endcall(cp["p"], None, centre=i_), endcall(cp["q"], None, centre=j_)
        sg = 1.0 if cp["med"] > 0 else -1.0
        ph, n_, agree = phi_pair(x, y)
        cp.update(fit=f_, a=(i_, j_), cfit=float(f_["c"][i_, j_]), r2=1 - f_["mse"][i_, j_] / f_["var"],
                  corner=ca, r2_corner=(1 - corners[ca] / f_["var"]) if ca else float("nan"),
                  ci_p=(float(np.percentile(bi, 2.5)), float(np.percentile(bi, 97.5))) if bi else (float("nan"),) * 2,
                  ci_q=(float(np.percentile(bj, 2.5)), float(np.percentile(bj, 97.5))) if bj else (float("nan"),) * 2,
                  phi=sg * ph, n=n_, agree=agree if sg > 0 else n_ - agree, ci=boot_ci(x, y, sg))
    # control: the between-chromosome join the pollen confirm best, predicted at the join itself
    ok_ = [j for j in JOINS if j["between"] and j["phi"] >= 0.5]
    if ok_:
        j = max(ok_, key=lambda d_: d_["phi"])
        p, q = j["e1"][0], j["e2"][0]
        f_ = fit_point(G, p, q)
        at = near_end(f_, 0 if j["e1"][1] == "start" else nwin[p] - 1, 0 if j["e2"][1] == "start" else nwin[q] - 1)
        if at:
            CONTROL = dict(join=j, p=p, q=q, fit=f_, a=at, r2=1 - f_["mse"][at] / f_["var"], cfit=float(f_["c"][at]))

# ------------------------------------------------------------------ Hi-C
HIC_ALL, CHK = None, None


class HiC:
    """A contact map of the two-haplotype assembly (npz: M, chroms, bin), read in bp coordinates."""

    def __init__(self, path):
        z = np.load(path)
        self.M, self.ch, self.bin, self.path = z["M"], [str(x) for x in z["chroms"]], int(z["bin"]), path
        unknown = [c for c in self.ch if c not in L2n]
        if unknown:
            raise KeyError("%s has sequences the assembly index lacks: %s" % (path, ", ".join(unknown[:5])))
        self.off, o = {}, 0
        for c in self.ch:
            self.off[c], o = o, o + -(-L2n[c] // self.bin)

    def idx(self, c, a, b):
        """Bins wholly inside a-b (a bin straddling a join belongs to neither side), at least one."""
        a, b = max(0, int(a)), min(L2n[c], int(b))
        i0, i1 = -(-a // self.bin), b // self.bin
        if i1 <= i0:
            i0, i1 = a // self.bin, max(a // self.bin + 1, -(-b // self.bin))
        return self.off[c] + i0, self.off[c] + i1, (i1 - i0) * self.bin / 1e6

    def block(self, c1, a1, b1, c2, a2, b2):
        """Contacts per Mb2 between two stretches."""
        i0, i1, l1 = self.idx(c1, a1, b1)
        j0, j1, l2 = self.idx(c2, a2, b2)
        return float(self.M[i0:i1, j0:j1].sum()) / (l1 * l2)

    def agg(self, order, step):
        """Contacts summed into `step`-bp bins over the copies in `order`: matrix, bins per copy."""
        f = max(1, int(round(step / self.bin)))
        sel, grp, per, k = [], [], [], 0
        for c in order:
            nbc = -(-L2n[c] // self.bin)
            na = -(-nbc // f)
            per.append(na)
            for i in range(na):
                b0, b1 = i * f, min((i + 1) * f, nbc)
                sel += list(range(self.off[c] + b0, self.off[c] + b1))
                grp += [k] * (b1 - b0)
                k += 1
        sel, grp = np.array(sel), np.array(grp)
        starts = np.r_[0, np.where(np.diff(grp) != 0)[0] + 1]
        Ms = self.M[np.ix_(sel, sel)].astype(np.float32)
        return np.add.reduceat(np.add.reduceat(Ms, starts, axis=0), starts, axis=1), per, f * self.bin


if not A.no_hic:
    hdir = os.path.join(A.phd_root, "results", row["species"], "qc/hic_remap")
    pa, pu = os.path.join(hdir, "contacts_all.npz"), os.path.join(hdir, "contacts_unique.npz")
    try:
        h_ = HiC(pa) if os.path.exists(pa) else None
    except KeyError as e:
        h_ = None
        NOTES.append(str(e))
    if h_ is not None:
        if all(c in h_.off for c in COPIES):
            USED.append(pa)
            Hm, per_, step_ = h_.agg(COPIES, BIN)
            if per_ == [nb[c] for c in COPIES]:
                HIC_ALL = Hm
            else:
                NOTES.append("Hi-C bins do not line up with the %g Mb map: map panel left out" % A.bin_mb)
            if not os.path.exists(pu):
                CHK = h_
        else:
            NOTES.append("Hi-C map: chromosome copies missing from %s" % pa)
        if CHK is None:
            del h_
    elif not os.path.exists(pa):
        NOTES.append("no Hi-C map at %s" % pa)
    try:
        h_ = HiC(pu) if os.path.exists(pu) else None
    except KeyError as e:
        h_ = None
        NOTES.append(str(e))
    if h_ is not None:
        if all(c in h_.off for c in COPIES):
            USED.append(pu)
            CHK = h_
HIC = {}
if CHK is not None:
    tr = [CHK.block(a, 0, L2n[a], b, 0, L2n[b]) for i_, a in enumerate(COPIES) for b in COPIES[i_ + 1:] if num(a) != num(b)]
    BASE = float(np.median(tr))
    HIC["src"], HIC["base"] = os.path.basename(CHK.path).replace(".npz", ""), BASE

    def xt(c1, a1, b1, c2, a2, b2):
        return CHK.block(c1, a1, b1, c2, a2, b2) / BASE

    def side_span(end_seg_side, copy, pos):
        """The FLANK next to a join on its side: left of pos for a left arm, right of it for a right arm."""
        return (copy, pos - FLANK, pos) if end_seg_side == "left" else (copy, pos, pos + FLANK)
    for j in JOINS:
        j["hic"] = xt(j["copy"], j["pos"] - FLANK, j["pos"], j["copy"], j["pos"], j["pos"] + FLANK)
        for e, pr in zip([e for e in ENDS if e["join"] is j], j.get("partners", [])):
            if pr is None:
                continue
            o_ = pr[0]
            sp1 = side_span(e["side"], j["copy"], j["pos"])
            sp2 = side_span(o_["side"], o_["join"]["copy"], o_["join"]["pos"])
            e["hic_partner"] = xt(*sp1, *sp2)
    # scan every pair of copies of different numbers for the strongest 10 x 10 Mb block
    H2, per2, st2 = CHK.agg(COPIES, SCAN_BIN)
    e2 = np.cumsum([0] + per2)
    mb2 = st2 / 1e6
    scan = []
    for i_, a in enumerate(COPIES):
        for j_, b in enumerate(COPIES):
            if j_ <= i_ or num(a) == num(b):
                continue
            sub = H2[e2[i_]:e2[i_ + 1], e2[j_]:e2[j_ + 1]].astype(float)
            k = min(SCAN_BOX, sub.shape[0], sub.shape[1])
            cs = np.pad(sub.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
            box = cs[k:, k:] - cs[:-k, k:] - cs[k:, :-k] + cs[:-k, :-k]
            v = box / ((k * mb2) ** 2) / BASE
            r_, c_ = np.unravel_index(np.argmax(v), v.shape)
            scan.append(dict(a=a, b=b, x=float(v[r_, c_]), pa=(r_ + k / 2) * mb2, pb=(c_ + k / 2) * mb2,
                             med=float(np.median(v))))
    scan.sort(key=lambda d_: -d_["x"])
    for k, d_ in enumerate(scan):
        d_["rank"] = k + 1
    HIC["scan"] = scan
    for cp in COUPLES:
        n1, n2 = num(src[cp["p"]][0]), num(src[cp["q"]][0])
        cp["hic_scan"] = [d_ for d_ in scan if {num(d_["a"]), num(d_["b"])} == {n1, n2}]
    pd.DataFrame(scan).to_csv(os.path.join(A.out, "hic_scan_%s.tsv" % S), sep="\t", index=False,
                              float_format="%.3f")
    del H2

# ------------------------------------------------------------------ figure 1: the map
def strip(ax, order, tot, where):
    """The arms of every copy as a coloured strip along the bottom ('x') or the left ('y') of a map."""
    sa = ax.inset_axes([0, -0.032, 1, 0.024] if where == "x" else [-0.032, 0, 0.024, 1])
    o_ = 0.0
    for c in order:
        for s_ in SEGS[c] or [dict(a=0, b=L2n[c], arm=None)]:
            a0, b0 = o_ + s_["a"] / 1e6, o_ + min(s_["b"], nb[c] * BIN) / 1e6
            col = ARM_COL.get(s_["arm"], "#B4B2A9") if s_["arm"] else WHOLE_COL
            sa.add_patch(Rectangle((a0, 0), b0 - a0, 1, color=col, lw=0) if where == "x" else
                         Rectangle((0, a0), 1, b0 - a0, color=col, lw=0))
            if s_["arm"] and b0 - a0 > 0.015 * tot:
                if where == "x":
                    sa.text((a0 + b0) / 2, 0.5, s_["arm"], ha="center", va="center", fontsize=5.5)
                else:
                    sa.text(0.5, (a0 + b0) / 2, s_["arm"], ha="center", va="center", fontsize=5.5, rotation=90)
        o_ += nb[c] * A.bin_mb
        if where == "x":
            sa.axvline(o_, color="white", lw=1.2)
        else:
            sa.axhline(o_, color="white", lw=1.2)
    if where == "x":
        sa.set_xlim(0, tot)
        sa.set_ylim(0, 1)
    else:
        sa.set_ylim(tot, 0)
        sa.set_xlim(0, 1)
    sa.axis("off")


def frame(ax, order, tot):
    edges = np.cumsum([0] + [nb[c] for c in order]) * A.bin_mb
    for e in edges[1:-1]:
        ax.axhline(e, color="#888780", lw=0.4)
        ax.axvline(e, color="#888780", lw=0.4)
    ax.set_xticks((edges[:-1] + edges[1:]) / 2)
    ax.set_xticklabels([cname(c) for c in order], rotation=90, fontsize=6)
    ax.set_yticks((edges[:-1] + edges[1:]) / 2)
    ax.set_yticklabels([cname(c) for c in order], fontsize=6)
    for i, c in enumerate(order):        # arm boundaries of every copy (breakpoints.csv)
        for a, b, _ in ARMS.get(c, [])[1:]:
            p_ = edges[i] + a / 1e6
            ax.plot([p_, p_], [0, tot], color="#2C2C2A", lw=0.6, ls=":")
            ax.plot([0, tot], [p_, p_], color="#2C2C2A", lw=0.6, ls=":")
    ax.set_xlim(0, tot)
    ax.set_ylim(tot, 0)
    if HAS_ARMS:
        strip(ax, order, tot, "x")
        strip(ax, order, tot, "y")
        ax.tick_params(axis="x", pad=13)
        ax.tick_params(axis="y", pad=13)


tot = NB * A.bin_mb
ncol = 2 if HIC_ALL is not None else 1
fig, axs = plt.subplots(1, ncol, figsize=(7.8 * ncol, 8.4))
axs = np.atleast_1d(axs)
fig.subplots_adjust(left=0.1 if ncol == 2 else 0.16, right=0.95 if ncol == 2 else 0.88, top=0.86, bottom=0.17,
                    wspace=0.32)
im = axs[0].imshow(np.ma.masked_invalid(M), cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest",
                   extent=(0, tot, tot, 0))
frame(axs[0], COPIES, tot)
axs[0].set_title("pollen genotypes, %s nuclei, %g Mb bins: red = inherited together, blue = never together, "
                 "white = independent" % ("{:,}".format(len(G)), A.bin_mb), loc="left", fontsize=8)
cb = fig.colorbar(im, ax=axs[0], fraction=0.035, pad=0.02)
cb.set_label("signed phi between the copies' windows", fontsize=7)
if HIC_ALL is not None:
    lab_ = np.repeat(np.arange(len(COPIES)), [nb[c] for c in COPIES])
    inter = lab_[:, None] != lab_[None, :]
    Xh = np.log10(HIC_ALL.astype(float) + 1)
    lo, hi = np.percentile(Xh[inter], 75), np.percentile(Xh[~inter], 99.5)
    im2 = axs[1].imshow(Xh, cmap="Reds", vmin=lo, vmax=max(hi, lo + 1e-6), interpolation="nearest",
                        extent=(0, tot, tot, 0))
    frame(axs[1], COPIES, tot)
    axs[1].set_title("tissue Hi-C, all reads (contacts_all), same bins: log10 contacts, white = typical between "
                     "chromosomes", loc="left", fontsize=8)
    cb2 = fig.colorbar(im2, ax=axs[1], fraction=0.035, pad=0.02)
    cb2.set_label("log10 contacts", fontsize=7)
title = "%s (%s): which chromosome copies the pollen inherit together%s" % (
    row["species"].replace("_", " "), S, ", beside the tissue Hi-C" if HIC_ALL is not None else "")
fig.suptitle(title, x=0.01, y=0.985, ha="left", fontsize=10)
if HAS_ARMS:
    present = [a for a in ARM_COL if any(s_["arm"] == a for ss in SEGS.values() for s_ in ss)]
    ARM_TXT = {"L1": "arm L1", "L2": "arm L2", "P": "arm P", "Q": "arm Q", "X": "X: on chr4 in hap1, chr3 in hap2",
               "core": "chr4 without X"}
    fig.legend([Patch(color=ARM_COL[a]) for a in present] + [Patch(color=WHOLE_COL)],
               [ARM_TXT.get(a, "arm " + a) for a in present] + ["whole chromosome"], loc="upper left",
               ncol=len(present) + 1, fontsize=7, frameon=False, bbox_to_anchor=(0.01, 0.955))
txt = ("Axes: every chromosome copy of the two-haplotype assembly, hap1 then hap2 for each number%s. Each reference "
       "window is drawn twice: on its own copy, and on its other copy with the sign flipped (placed by relative "
       "position along the matching arm)%s. %s" % (
           "; strips = arms, dotted = arm boundaries (breakpoints.csv); A = chr1_hap1 (L1+P), B = chr2_hap1 (L2+P), "
           "C = chr1_hap2 (L1+Q), D = chr2_hap2 (L2+Q)" if LET else "",
           "" if not HAS_ARMS else "", "; ".join(NOTES)))
fig.text(0.01, 0.01, txt, fontsize=6.5, ha="left", va="bottom", wrap=True)
png = os.path.join(A.out, "cosegregation_2n_%s.png" % S)
fig.savefig(png, dpi=150)
plt.close(fig)


# ------------------------------------------------------------------ figure 2: the karyotype test
def segname(seg):
    """An arm by its name, a whole chromosome by its number."""
    return seg["arm"] or src[seg["piece"]][0].split("_")[0]


def joinlabel(j):
    c = j["copy"]
    return "%s%s: %s | %s at %.1f Mb" % (c, " (%s)" % LET[c] if c in LET else "", segname(j["left"]),
                                          segname(j["right"]), j["pos"] / 1e6)


def endlabel(e):
    j = e["join"]
    return "%s's %s (%s %s %.1f Mb)" % (LET.get(j["copy"], j["copy"]), segname(e["seg"]), j["copy"],
                                       "up to" if e["side"] == "left" else "from", j["pos"] / 1e6)


KT = None
if JOINS:
    rows = []
    for j in sorted(JOINS, key=lambda d_: -d_["phi"]):
        note = "%d of %d nuclei agree" % (j["agree"], j["n"])
        if "hic" in j:
            note += "; Hi-C across the join %.1fx trans" % j["hic"]
        for e, pr in zip([e for e in ENDS if e["join"] is j], j.get("partners", [])):
            if pr:
                note += "\n%s goes with %s (phi %+.2f%s)" % (
                    endlabel(e), endlabel(pr[0]), pr[1], "; Hi-C %.1fx" % e["hic_partner"] if "hic_partner" in e else "")
        rows.append(dict(label=joinlabel(j), phi=j["phi"], ci=j["ci"], note=note, kind="join"))
    for cp in COUPLES:
        a1, a2 = wmid(cp["p"], cp["a"][0]), wmid(cp["q"], cp["a"][1])
        note = ("%d of %d nuclei agree around the point; one coupling point explains R2 %.2f of the block, "
                "ends joined %.2f" % (cp["agree"], cp["n"], cp["r2"], cp["r2_corner"]))
        if cp.get("hic_scan"):
            best = cp["hic_scan"][0]
            note += "\nHi-C: strongest %s x %s block %.1fx trans at %.0f / %.0f Mb (rank %d of %d copy pairs)" % (
                best["a"], best["b"], best["x"], best["pa"], best["pb"], best["rank"], len(HIC["scan"]))
        rows.append(dict(label="%s %.0f Mb x %s %.0f Mb: no join in the assembly" % (cp["p"], a1, cp["q"], a2),
                         phi=cp["phi"], ci=cp["ci"], note=note, kind="couple"))
    KT = rows
    nmap = (1 if COUPLES else 0) + (1 if CONTROL is not None else 0)
    head, top_in, below, gap, rect_in = 0.95, 0.55 * len(rows) + 0.7, 0.45, 1.5, 3.3   # inches
    H = head + top_in + below + nmap * (gap + rect_in) + 0.45
    fig = plt.figure(figsize=(16, H))

    def fy(inch):
        """Figure fraction of a height given in inches from the top."""
        return 1 - inch / H
    ax = fig.add_axes([0.25, fy(head + top_in), 0.3, top_in / H])
    y = np.arange(len(rows))[::-1]
    for yy, r in zip(y, rows):
        col = "#534AB7" if r["kind"] == "couple" else ("#2C2C2A" if r["phi"] >= 0.5 else "#A32D2D")
        ax.plot(r["ci"], [yy, yy], color=col, lw=1.6)
        ax.plot(r["phi"], yy, "o", color=col, ms=6)
        ax.text(1.03, yy, r["note"], fontsize=6.8, va="center", ha="left", color=col if col != "#2C2C2A" else
                "#444441", transform=ax.get_yaxis_transform())
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=7.5)
    for t, r in zip(ax.get_yticklabels(), rows):
        t.set_color("#534AB7" if r["kind"] == "couple" else ("#2C2C2A" if r["phi"] >= 0.5 else "#A32D2D"))
    ax.axvline(0, color="#B4B2A9", lw=0.8)
    ax.axvline(1, color="#B4B2A9", lw=0.8, ls="--")
    ax.set_xlim(-0.6, 1.1)
    ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel("phi between the %g Mb on either side, across nuclei (95%% interval):\n+1 = one chromosome, "
                  "0 = inherited independently" % (END_K * A.win_mb), fontsize=7.5)
    ax.tick_params(axis="x", labelsize=7)
    fig.text(0.01, fy(head - 0.3), "a  Every join of the two-haplotype assembly, tested by the pollen. Black: the "
             "pollen confirm it; red: they do not (and what they inherit with each side instead); purple: a coupling "
             "for which the assembly has no join.", fontsize=9, va="bottom")

    def mapax(rect, Mx, p, q, title, anchor=None, ci=None, zd=None):
        ax_ = fig.add_axes(rect)
        ext = (0, nwin[q] * A.win_mb, nwin[p] * A.win_mb, 0)
        im_ = ax_.imshow(np.ma.masked_invalid(Mx), cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest",
                         extent=ext, aspect="auto")
        if zd is not None and zd.any():
            ii, jj = np.nonzero(zd)
            ax_.plot((jj + 0.5) * A.win_mb, (ii + 0.5) * A.win_mb, ".", color="#2C2C2A", ms=2.2)
        if anchor is not None:
            ax_.plot((anchor[1] + 0.5) * A.win_mb, (anchor[0] + 0.5) * A.win_mb, "x", color="#2C2C2A", ms=9, mew=2)
        if ci is not None:
            yy_, xx_ = (anchor[0] + 0.5) * A.win_mb, (anchor[1] + 0.5) * A.win_mb
            ax_.plot([ci[1][0], ci[1][1]], [yy_, yy_], color="#2C2C2A", lw=1)
            ax_.plot([xx_, xx_], [ci[0][0], ci[0][1]], color="#2C2C2A", lw=1)
        ax_.set_xlabel("%s, Mb" % q, fontsize=7.5)
        ax_.set_ylabel("%s, Mb" % p, fontsize=7.5)
        ax_.set_title(title, loc="left", fontsize=8.3)
        ax_.tick_params(labelsize=7)
        return ax_, im_
    cols = [0.06, 0.385, 0.71]
    wd = 0.25
    cur = head + top_in + below
    im_ = None
    if COUPLES:
        cp = COUPLES[0]
        f_ = cp["fit"]
        i_, j_ = cp["a"]
        top_ = cur + gap
        yb, rect_h = fy(top_ + rect_in), rect_in / H
        zd = (f_["disc"] == 0) & (f_["nn"] >= PHI_N) if cp["med"] > 0 else None
        axb, im_ = mapax([cols[0], yb, wd, rect_h], f_["O"], cp["p"], cp["q"],
                         "b  observed in the pollen (dots: none of the >= %d nuclei\ncalled in both windows disagree)"
                         % PHI_N, (i_, j_), (cp["ci_p"], cp["ci_q"]), zd)
        cbx = axb.inset_axes([1.02, 0.0, 0.03, 1.0])
        fig.colorbar(im_, cax=cbx).ax.tick_params(labelsize=6)
        pred = cp["cfit"] * np.outer(f_["Pp"][:, i_], f_["Pq"][:, j_])
        mapax([cols[1], yb, wd, rect_h], pred, cp["p"], cp["q"],
              "c  predicted: complete co-inheritance at x only (lines: 95%%\ninterval), c %.2f; explains R2 %.2f of b" % (
                  cp["cfit"], cp["r2"]), (i_, j_), (cp["ci_p"], cp["ci_q"]))
        if cp["corner"]:
            a_, b_ = cp["corner"]
            predc = f_["c"][a_, b_] * np.outer(f_["Pp"][:, a_], f_["Pq"][:, b_])
            mapax([cols[2], yb, wd, rect_h], predc, cp["p"], cp["q"],
                  "d  predicted: ends joined (a fusion; best of the 4 end pairs:\n%s %s, %s %s), R2 %.2f" % (
                      cp["p"], "start" if a_ < nwin[cp["p"]] / 2 else "end", cp["q"],
                      "start" if b_ < nwin[cp["q"]] / 2 else "end", cp["r2_corner"]), (a_, b_))
        fig.text(0.01, fy(top_ - 0.75), "%s x %s: where the pollen couple them, and two models of why" % (
            cp["p"], cp["q"]), fontsize=9, va="bottom")
        cur = top_ + rect_in
    if CONTROL is not None:
        f_ = CONTROL["fit"]
        a_, b_ = CONTROL["a"]
        top_ = cur + gap
        yb, rect_h = fy(top_ + rect_in), rect_in / H
        j = CONTROL["join"]
        _, im2_ = mapax([cols[0], yb, wd, rect_h], f_["O"], CONTROL["p"], CONTROL["q"],
                        "e  observed in the pollen", (a_, b_))
        im_ = im_ or im2_
        pred = CONTROL["cfit"] * np.outer(f_["Pp"][:, a_], f_["Pq"][:, b_])
        mapax([cols[1], yb, wd, rect_h], pred, CONTROL["p"], CONTROL["q"],
              "f  predicted from the join alone (x: the windows next to it),\nc %.2f; explains R2 %.2f of e" % (
                  CONTROL["cfit"], CONTROL["r2"]), (a_, b_))
        fig.text(0.01, fy(top_ - 0.75), "Control, the same model for a join the assembly has and the pollen "
                 "confirm: %s" % joinlabel(j), fontsize=9, va="bottom")
    if COUPLES:
        cp = COUPLES[0]
        i_, j_ = cp["a"]
        f_ = cp["fit"]
        if CONTROL is None:
            top_ = cur + gap
            yb, rect_h = fy(top_ + rect_in), rect_in / H
        axp = fig.add_axes([cols[2], yb, wd, rect_h])
        for vec, p, a_, col in ((f_["Pp"][:, i_], cp["p"], i_, "#0F6E56"), (f_["Pq"][:, j_], cp["q"], j_, "#534AB7")):
            axp.plot([wmid(p, k) for k in range(nwin[p])], vec, color=col, lw=1.6, marker="o", ms=2.5,
                     label="%s: phi of each window with %.0f Mb" % (p, wmid(p, a_)))
        axp.axhline(0, color="#B4B2A9", lw=0.8)
        axp.set_ylim(-0.3, 1.05)
        axp.set_xlabel("Mb along the chromosome", fontsize=7.5)
        axp.set_ylabel("phi within the chromosome", fontsize=7.5)
        axp.legend(fontsize=6.5, frameon=False, loc="lower center")
        axp.set_title("g  the two profiles c multiplies: how often each window\nrecombines away from the point (1 = "
                      "never)", loc="left", fontsize=8.3)
        axp.tick_params(labelsize=7)
    fig.suptitle("%s (%s): the two-haplotype assembly's karyotype tested against %d pollen nuclei" % (
        row["species"].replace("_", " "), S, len(G)), x=0.01, y=fy(0.2), ha="left", va="top", fontsize=10.5)
    png2 = os.path.join(A.out, "karyotype_test_%s.png" % S)
    fig.savefig(png2, dpi=150)
    plt.close(fig)
    with open(os.path.join(A.out, "karyotype_test_%s.tsv" % S), "w") as fo:
        fo.write("kind\twhat\tphi\tci_lo\tci_hi\tdetail\n")
        for r in rows:
            fo.write("%s\t%s\t%.3f\t%.3f\t%.3f\t%s\n" % (r["kind"], r["label"], r["phi"], r["ci"][0], r["ci"][1],
                                                       r["note"].replace("\n", " | ")))
        if CONTROL is not None:
            fo.write("control\t%s\t%.3f\t\t\tpredicted from the join alone: c %.2f, R2 %.2f\n" % (
                joinlabel(CONTROL["join"]), CONTROL["join"]["phi"], CONTROL["cfit"], CONTROL["r2"]))

# ------------------------------------------------------------------ report
print("%s: %d nuclei, %d windows, %d copies; other copies: %s" % (
    S, len(G), len(wins), len(COPIES), "; ".join("%s -> %s %.1f-%.1f Mb (%s)" % (p, o[0], o[1] / 1e6, o[2] / 1e6, o[3])
                                                 if o else "%s -> none" % p for p, o in OTHER.items())))
for n_ in NOTES:
    print("note:", n_)
if KT:
    print("KARYOTYPE TEST (phi of the %d Mb either side across nuclei; 95%% interval; +1 = one chromosome, 0 = "
          "independent)" % (END_K * A.win_mb))
    for r in KT:
        print("  %-48s phi %+.2f (%+.2f to %+.2f)\n      %s" % (r["label"], r["phi"], r["ci"][0], r["ci"][1],
                                                               r["note"].replace("\n", "\n      ")))
    for cp in COUPLES:
        print("  COUPLING %s x %s: point at %.0f Mb (95%% %.0f-%.0f) x %.0f Mb (95%% %.0f-%.0f), c %.2f, R2 %.2f; "
              "ends joined R2 %.2f" % (cp["p"], cp["q"], wmid(cp["p"], cp["a"][0]), cp["ci_p"][0], cp["ci_p"][1],
                                      wmid(cp["q"], cp["a"][1]), cp["ci_q"][0], cp["ci_q"][1], cp["cfit"], cp["r2"],
                                      cp["r2_corner"]))
    if CONTROL is not None:
        print("  CONTROL %s: predicted from the join alone c %.2f, R2 %.2f" % (
            joinlabel(CONTROL["join"]), CONTROL["cfit"], CONTROL["r2"]))
if HIC.get("scan"):
    print("HI-C SCAN (%s; strongest 10 x 10 Mb block per pair of copies of different numbers, x the median between "
          "them = %.3g contacts per Mb2):" % (HIC["src"], HIC["base"]))
    for d_ in HIC["scan"][:10]:
        print("  %2d. %-10s x %-10s %6.1fx at %.0f / %.0f Mb" % (d_["rank"], d_["a"], d_["b"], d_["x"], d_["pa"],
                                                                 d_["pb"]))
    for cp in COUPLES:
        for d_ in cp.get("hic_scan", []):
            print("  %s x %s: strongest block %.1fx at %.0f / %.0f Mb, rank %d of %d" % (
                d_["a"], d_["b"], d_["x"], d_["pa"], d_["pb"], d_["rank"], len(HIC["scan"])))
print("files read:\n  " + "\n  ".join(USED))
print("wrote %s and %s" % (png, os.path.join(A.out, "cosegregation_2n_%s.tsv" % S)))
if KT:
    print("wrote %s and %s" % (os.path.join(A.out, "karyotype_test_%s.png" % S),
                               os.path.join(A.out, "karyotype_test_%s.tsv" % S)))
if HIC.get("scan"):
    print("wrote %s" % os.path.join(A.out, "hic_scan_%s.tsv" % S))
