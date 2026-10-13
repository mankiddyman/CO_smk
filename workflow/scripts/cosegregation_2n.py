#!/usr/bin/env python3
"""cosegregation_2n.py -- what the pollen inherit together, drawn on both haplotypes like a Hi-C map.

For every pair of WIN_MB windows of the crossover reference, phi = the correlation of the two windows' genotypes
across pollen nuclei (+1 = always inherited together, 0 = independently, -1 = never together). In each window a
nucleus carries either the reference's allele (REF) or the other copy's (ALT). The REF copy is the reference
sequence itself; the ALT copy is its homolog in the other haplotype, or, for an arm whose other copy sits in the
same haplotype (an arm with no homolog in the other one), that copy. So every window sits twice on the
two-haplotype map: on its REF copy, and on its ALT copy with the sign of phi flipped. Homologous copies then come
out at -1 (a nucleus has one or the other: blue lines where Hi-C has its homolog lines), each copy is red along
its own diagonal, and two copies of different chromosomes that travel together through meiosis share a red block.
Next to it, the Hi-C map of the tissue on the same axes (contacts_all.npz), as in the supervisor report.

Where the other copy sits, all read from files: the split plan (config/split_plans/SAMPLE.csv, 'homolog
chrN_hapM'), the breakpoints of every chromosome of both haplotypes (qc/reference/SAMPLE/breakpoints.csv, from
breakpoints_from_paf.py: which homolog each arm has, 'none' = its other copy is not in the other haplotype),
and the hap2-on-hap1 synteny blocks (qc/translocations/<hap1 sample>/<hap1 sample>_synteny_blocks.tsv) for the
orientation of each arm against its other copy. Without breakpoints.csv (a reference not cut, e.g. D. binata) a
chromosome's other copy is the other haplotype's chromosome of the same number, end to end. A window is placed
on its other copy by its relative position along the matching arm. An arm whose other copy is in the same
haplotype has no alignment to it and is placed in the same orientation; the caption says which arms those are.

Usage: cosegregation_2n.py SAMPLE [--out reports/cosegregation_2n] [--win_mb 10] [--bin_mb 10]
                           [--phd_root ...] [--no_hic]
Writes OUT/cosegregation_2n_SAMPLE.png and OUT/cosegregation_2n_SAMPLE.tsv (bin pairs: copy, Mb, copy, Mb,
mean signed phi, window pairs). Run from the CO_smk root.
"""
import argparse
import csv
import math
import os
import re

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("--out", default="reports/cosegregation_2n")
ap.add_argument("--win_mb", type=float, default=10)
ap.add_argument("--bin_mb", type=float, default=10)
ap.add_argument("--phd_root", default="/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd")
ap.add_argument("--no_hic", action="store_true")
A = ap.parse_args()
S, W, BIN = A.sample, int(A.win_mb * 1e6), int(A.bin_mb * 1e6)
MIN_MOL, GT_MAJ, PHI_N = 3, 0.8, 20          # as in supervisor_report.py: >= 3 molecules, >= 80% one allele
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
if os.path.exists(pp):
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
X1, X0 = (G == 1).astype(float), (G == 0).astype(float)
n11, n00, n10, n01 = X1.T @ X1, X0.T @ X0, X1.T @ X0, X0.T @ X1
den = np.sqrt((n11 + n10) * (n01 + n00) * (n11 + n01) * (n10 + n00))
with np.errstate(divide="ignore", invalid="ignore"):
    PHI = (n11 * n00 - n10 * n01) / den
PHI[(den == 0) | (n11 + n00 + n10 + n01 < PHI_N)] = np.nan
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
with np.errstate(divide="ignore", invalid="ignore"):
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


# ------------------------------------------------------------------ figure
def frame(ax, order, tot):
    edges = np.cumsum([0] + [nb[c] for c in order]) * A.bin_mb
    for e in edges[1:-1]:
        ax.axhline(e, color="#888780", lw=0.4)
        ax.axvline(e, color="#888780", lw=0.4)
    ax.set_xticks((edges[:-1] + edges[1:]) / 2)
    ax.set_xticklabels([c.replace("_hap", " h") for c in order], rotation=90, fontsize=6)
    ax.set_yticks((edges[:-1] + edges[1:]) / 2)
    ax.set_yticklabels([c.replace("_hap", " h") for c in order], fontsize=6)
    for i, c in enumerate(order):        # arm boundaries of every copy (breakpoints.csv)
        for a, b, _ in ARMS.get(c, [])[1:]:
            p_ = edges[i] + a / 1e6
            ax.plot([p_, p_], [0, tot], color="#2C2C2A", lw=0.6, ls=":")
            ax.plot([0, tot], [p_, p_], color="#2C2C2A", lw=0.6, ls=":")
    ax.set_xlim(0, tot)
    ax.set_ylim(tot, 0)


hic = None
if not A.no_hic:
    hp = os.path.join(A.phd_root, "results", row["species"], "qc/hic_remap", "contacts_all.npz")
    if os.path.exists(hp):
        USED.append(hp)
        z = np.load(hp)
        HM, hch, hbin = z["M"], [str(x) for x in z["chroms"]], int(z["bin"])
        if all(c in hch for c in COPIES):
            f = max(1, int(round(BIN / hbin)))
            hnb = {c: -(-L2n[c] // hbin) for c in hch}
            hoff, o = {}, 0
            for c in hch:
                hoff[c], o = o, o + hnb[c]
            sel, grp = [], []
            for c in COPIES:
                for j in range(nb[c]):
                    b0, b1 = j * f, min((j + 1) * f, hnb[c])
                    if b1 <= b0:            # the Hi-C bins are a touch coarser at the end of a copy
                        b0, b1 = hnb[c] - 1, hnb[c]
                    sel += list(range(hoff[c] + b0, hoff[c] + b1))
                    grp += [boff[c] + j] * (b1 - b0)
            sel, grp = np.array(sel), np.array(grp)
            starts = np.r_[0, np.where(np.diff(grp) != 0)[0] + 1]
            Ms = HM[np.ix_(sel, sel)].astype(np.float32)
            hic = np.add.reduceat(np.add.reduceat(Ms, starts, axis=0), starts, axis=1)
            del Ms
        else:
            NOTES.append("Hi-C map: chromosome copies missing from %s" % hp)
        del HM
    else:
        NOTES.append("no Hi-C map at %s" % hp)
tot = NB * A.bin_mb
ncol = 2 if hic is not None else 1
fig, axs = plt.subplots(1, ncol, figsize=(7.4 * ncol, 7.6))
axs = np.atleast_1d(axs)
im = axs[0].imshow(np.ma.masked_invalid(M), cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest",
                   extent=(0, tot, tot, 0))
frame(axs[0], COPIES, tot)
axs[0].set_title("pollen: inherited together (red) or never together (blue), %s nuclei, %g Mb bins"
                 % ("{:,}".format(len(G)), A.bin_mb), loc="left", fontsize=8.5)
cb = fig.colorbar(im, ax=axs[0], fraction=0.035, pad=0.02)
cb.set_label("signed phi between the copies' windows", fontsize=7)
if hic is not None:
    lab_ = np.repeat(np.arange(len(COPIES)), [nb[c] for c in COPIES])
    inter = lab_[:, None] != lab_[None, :]
    Xh = np.log10(hic.astype(float) + 1)
    lo, hi = np.percentile(Xh[inter], 75), np.percentile(Xh[~inter], 99.5)
    im2 = axs[1].imshow(Xh, cmap="Reds", vmin=lo, vmax=max(hi, lo + 1e-6), interpolation="nearest",
                        extent=(0, tot, tot, 0))
    frame(axs[1], COPIES, tot)
    axs[1].set_title("tissue Hi-C, all reads (contacts_all), same bins: log10 contacts, white = typical between "
                     "chromosomes", loc="left", fontsize=8.5)
    cb2 = fig.colorbar(im2, ax=axs[1], fraction=0.035, pad=0.02)
    cb2.set_label("log10 contacts", fontsize=7)
title = "%s (%s), both haplotypes: what the pollen inherit together" % (row["species"].replace("_", " "), S)
if hic is not None:
    title += ", next to the tissue's Hi-C"
if ARMS:
    title += "\ndotted = arm boundaries of every copy (breakpoints.csv)"
fig.suptitle(title, x=0.01, ha="left", fontsize=9.5)
txt = ("Each reference window is drawn twice: on its own copy, and on its other copy with the sign flipped "
       "(placed by relative position along the matching arm). " + ("; ".join(NOTES) if NOTES else ""))
fig.text(0.01, 0.005, txt, fontsize=6.5, ha="left", va="bottom", wrap=True)
fig.tight_layout(rect=(0, 0.04, 1, 0.95))
png = os.path.join(A.out, "cosegregation_2n_%s.png" % S)
fig.savefig(png, dpi=150)
plt.close(fig)
print("%s: %d nuclei, %d windows, %d copies; other copies: %s" % (
    S, len(G), len(wins), len(COPIES), "; ".join("%s -> %s %.1f-%.1f Mb (%s)" % (p, o[0], o[1] / 1e6, o[2] / 1e6, o[3])
                                                 if o else "%s -> none" % p for p, o in OTHER.items())))
for n_ in NOTES:
    print("note:", n_)
print("files read:\n  " + "\n  ".join(USED))
print("wrote %s and %s" % (png, os.path.join(A.out, "cosegregation_2n_%s.tsv" % S)))
