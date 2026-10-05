#!/usr/bin/env python3
"""chr5_chr6_check.py -- why are chr5 and chr6 inherited together in the paradoxa pollen?

Two dry-lab looks, both from files that already exist:
  A. POLLEN: where along chr5 and chr6 the linkage sits (linkage_scan r matrix on the
     crossover reference). One character per pair of 10 Mb windows:
       '#' r <= 0.15 (inherited together)   '+' r <= 0.30   '.' r > 0.30 (independent)   ' ' too few cells
     A big reciprocal exchange links long stretches; a small exchanged piece, kept by
     pollen that die without it, links a patch that fades with distance.
  B. ASSEMBLY, fine scale: pieces of hap2 chr5/chr6 (>= 100 kb, MAPQ >= 20) whose best
     match is on the OTHER hap1 chromosome (chr6/chr5), from the hap2-on-hap1 alignment.
     'moved' = hap2 has no copy of that stretch at its usual place (a small translocation,
     too small for the 2 Mb Hi-C map); 'extra copy' = it is also at its usual place.

Run from the CO_smk root.
Usage: chr5_chr6_check.py R_MATRIX.npy WINDOWS.tsv HAP2_ON_HAP1.paf OUT.txt
"""
import sys

import numpy as np
import pandas as pd

RM, WT, PAF, OUT = sys.argv[1:5]
out = []

# ---- A. pollen linkage between chr5 and chr6 windows -------------------------------------
R = np.load(RM)
w = pd.read_csv(WT, sep="\t")
c5 = w.index[w.chrom == "chr5_hap1"].to_numpy()
c6 = w.index[w.chrom == "chr6_hap1"].to_numpy()
out += ["A. POLLEN: r between chr5_hap1 windows (rows) and chr6_hap1 windows (columns), 10 Mb each",
        "   '#' r <= 0.15   '+' <= 0.30   '.' > 0.30   ' ' too few cells", ""]
head = "".join(str((k * 10) // 100 % 10) if (k * 10) % 100 == 0 else " " for k in range(len(c6)))
out.append("   chr6 Mb (x100) ->   " + head)
for i in c5:
    row = ""
    for j in c6:
        r = R[i, j]
        row += " " if np.isnan(r) else ("#" if r <= 0.15 else ("+" if r <= 0.30 else "."))
    vals = R[i, c6]
    best = np.nanargmin(vals) if np.isfinite(vals).any() else None
    tag = "  min r %.2f at chr6 %d-%d Mb" % (vals[best], best * 10, best * 10 + 10) if best is not None else ""
    out.append("   chr5 %3d-%3d Mb   %s%s" % (w.window[i] * 10, w.window[i] * 10 + 10, row, tag))
blk = R[np.ix_(c5, c6)]
ctl = R[np.ix_(c5, w.index[w.chrom == "chr1_hap1"].to_numpy())]
out += ["", "   chr5 x chr6: median r %.2f over %d window pairs; control chr5 x chr1_hap1: median r %.2f" % (
    np.nanmedian(blk), np.isfinite(blk).sum(), np.nanmedian(ctl)), ""]

# ---- B. small moved pieces between chr5 and chr6 in the assembly -------------------------
p = pd.read_csv(PAF, sep="\t", header=None, usecols=range(12),
                names=["q", "ql", "qs", "qe", "st", "t", "tl", "ts", "te", "nm", "al", "mq"])
if not p.q.isin(["chr5_hap2", "chr6_hap2"]).any():
    sys.exit("no chr5_hap2/chr6_hap2 queries in %s; query names look like: %s" % (PAF, ", ".join(p.q.unique()[:6])))
p = p[(p.mq >= 20) & (p.qe - p.qs >= 20000)]
BIN = 10000
cov = {}
for c in ("chr5", "chr6"):
    hits = p[(p.q == c + "_hap2") & (p.t == c + "_hap1")]
    tl = int(p[p.t == c + "_hap1"].tl.max()) if (p.t == c + "_hap1").any() else 0
    v = np.zeros(tl // BIN + 1, dtype=bool)
    for a, b in zip(hits.ts, hits.te):
        v[a // BIN:b // BIN + 1] = True
    cov[c] = v
out += ["B. ASSEMBLY: pieces of hap2 chr5/chr6 (>= 100 kb merged, MAPQ >= 20) that match the OTHER chromosome on hap1"]
rows = []
for qc, tc in (("chr5_hap2", "chr6_hap1"), ("chr6_hap2", "chr5_hap1")):
    h = p[(p.q == qc) & (p.t == tc)].sort_values("qs")
    seg = None
    for r in h.itertuples():
        if seg and r.qs - seg["qe"] <= 200000 and abs(r.ts - seg["te"]) <= 500000 + abs(r.te - r.ts):
            seg["qe"] = max(seg["qe"], r.qe)
            seg["ts"], seg["te"] = min(seg["ts"], r.ts), max(seg["te"], r.te)
            seg["nm"] += r.nm
            seg["al"] += r.al
            seg["len"] += r.qe - r.qs
        else:
            if seg:
                rows.append(seg)
            seg = dict(q=qc, qs=r.qs, qe=r.qe, t=tc, ts=r.ts, te=r.te, nm=r.nm, al=r.al, len=r.qe - r.qs)
    if seg:
        rows.append(seg)
rows = [s for s in rows if s["len"] >= 100000]
if not rows:
    out.append("   none: no stretch >= 100 kb of hap2 chr5/chr6 matches the other chromosome on hap1")
for s in sorted(rows, key=lambda s: -s["len"])[:25]:
    usual = s["t"].split("_")[0]          # the hap1 chromosome it matched; is it ALSO on hap2's own copy of that chromosome?
    v = cov.get(usual)
    frac = float(v[s["ts"] // BIN:s["te"] // BIN + 1].mean()) if v is not None and len(v) else float("nan")
    kind = "moved (no copy at the usual place on hap2)" if frac < 0.2 else (
        "extra copy (also at the usual place)" if frac > 0.6 else "partly")
    out.append("   %s %7.2f-%7.2f Mb  ->  %s %7.2f-%7.2f Mb   %5.0f kb aligned, identity %.3f   %s" % (
        s["q"], s["qs"] / 1e6, s["qe"] / 1e6, s["t"], s["ts"] / 1e6, s["te"] / 1e6, s["len"] / 1e3,
        s["nm"] / max(s["al"], 1), kind))
open(OUT, "w").write("\n".join(out) + "\n")
print("\n".join(out))
