#!/usr/bin/env python3
"""score_verdicts.py -- every hapCO setting scored against your crossover verdicts from the review page.

The evidence (committed under docs/co_review/):
  SAMPLE_verdicts.csv   exported from review.html with 'e': one row per event with your verdict (real / false /
                        unsure; source user or auto), missed crossovers, masked sequences, bad cells
  SAMPLE_events.tsv.gz  the page's events: which of the sweep's settings call each one (a 0/1 string)
  SAMPLE_settings.tsv   the sweep's settings, in the order of that string
Scoring: real = your 'real' plus auto-real (unless --no_auto); 'unsure' is left out; so are events on masked
sequences and in bad cells; a missed crossover counts as a miss for every setting.
  TP  real events the setting calls   FP  false events it calls   FN  real events (and missed ones) it does not
  precision TP/(TP+FP), recall TP/(TP+FN), F1 = 2PR/(P+R)
Uncertainty: the cells are resampled with replacement (--boot): 95% interval of F1, and the share of resamples
in which a setting is (tied) best.
--extra b0.3_mn6_bs2_t5,...  settings outside the sweep, re-called here with hapco_port.py from the pipeline's
             per-cell tables and matched to the page's events (same direction, overlapping possible position).
             Before that, the re-call is checked on the sweep's own columns (the current setting, and --expect):
             they must come out identical, else exit 1. Calls of an extra setting that match no reviewed event
             are reported as 'unreviewed' (not scored).
--expect S   exit 1 unless S ties for the best F1 (the config stamp is written only if it does).
Usage: score_verdicts.py SAMPLE [--dir docs/co_review] [--extra ...] [--expect ...] [--out FILE]
Run from the CO_smk directory.
"""
import argparse
import csv
import gzip
import os
import re
import sys
from multiprocessing import get_context

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hapco_port as H  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("--dir", default="docs/co_review")
ap.add_argument("--extra", default="")
ap.add_argument("--expect", default="")
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--no_auto", action="store_true")
ap.add_argument("--out")
A = ap.parse_args()
S = A.sample
OUT = A.out or "qc/co_review/%s/scores.tsv" % S
GN = {"REF": 1, "ALT": 0, "mixed": -1}

SET = list(csv.DictReader(open("%s/%s_settings.tsv" % (A.dir, S)), delimiter="\t"))
IDS = [r["id"] for r in SET]
NS = len(IDS)
CUR = next((r["id"] for r in SET if r.get("label") == "current"), None)
with gzip.open("%s/%s_events.tsv.gz" % (A.dir, S), "rt") as fh:
    EV = list(csv.DictReader(fh, delimiter="\t"))
BYID = {int(e["event"]): e for e in EV}
for e in EV:
    if len(e["settings"]) != NS:
        sys.exit("events and settings do not belong together: %d settings, an event has %d" % (NS, len(e["settings"])))
    e["pre"], e["post"] = (GN[x] for x in e["direction"].split(">"))
    for k in ("start", "end", "cstart", "cend"):
        e[k] = int(e[k])
    e["bits"] = np.frombuffer(e["settings"].encode(), dtype=np.uint8) == 49
    e["v"] = None

rows = [l for l in open("%s/%s_verdicts.csv" % (A.dir, S)) if not l.startswith("#")]
MASK, BAD, MISSED, bad_rows = set(), set(), [], []
for r in csv.DictReader(rows):
    if r["type"] == "event":
        e = BYID.get(int(r["event"]))
        if (e is None or e["barcode"] != r["barcode"] or e["chrom"] != r["chrom"] or e["start"] != int(r["start"])
                or e["end"] != int(r["end"]) or e["direction"] != r["direction"]):
            bad_rows.append(r["event"])
            continue
        if r["verdict"] in ("real", "false", "unsure") and (r["source"] == "user" or (r["source"] == "auto" and not A.no_auto)):
            e["v"] = r["verdict"]
    elif r["type"] == "mask" and r["verdict"] == "masked":
        MASK.add((r["barcode"], r["chrom"]))
    elif r["type"] == "cell" and r["verdict"] == "bad":
        BAD.add(r["barcode"])
    elif r["type"] == "missed":
        MISSED.append((r["barcode"], r["chrom"], int(r["start"])))
if bad_rows:
    sys.exit("%d verdict rows do not match the events (first: event %s): verdicts and events are from different pages"
             % (len(bad_rows), bad_rows[0]))

cf = "qc/co_sweep/%s/cells.tsv" % S
if os.path.exists(cf):
    CELLS = [r["barcode"] for r in csv.DictReader(open(cf), delimiter="\t") if r["called"] == "1"]
else:
    CELLS = sorted(set(e["barcode"] for e in EV))
CI = {b: i for i, b in enumerate(CELLS)}
for e in EV:
    if e["barcode"] not in CI:
        sys.exit("event %s is in a cell the sweep did not call (%s)" % (e["event"], e["barcode"]))
scored = [e for e in EV if (e["barcode"], e["chrom"]) not in MASK and e["barcode"] not in BAD]
miss = [m for m in MISSED if (m[0], m[1]) not in MASK and m[0] not in BAD]


def per_cell(called):
    """called: list (per event, in EV order) of booleans for one setting -> per-cell TP, FP, FN columns."""
    tp, fp, fn = (np.zeros(len(CELLS)) for _ in range(3))
    for e, c in zip(EV, called):
        if (e["barcode"], e["chrom"]) in MASK or e["barcode"] in BAD:
            continue
        i = CI[e["barcode"]]
        if e["v"] == "real":
            if c:
                tp[i] += 1
            else:
                fn[i] += 1
        elif e["v"] == "false" and c:
            fp[i] += 1
    for b, _, _ in miss:
        fn[CI[b]] += 1
    return tp, fp, fn


BITS = np.array([e["bits"] for e in EV])                     # events x settings
cols = [per_cell(BITS[:, k]) for k in range(NS)]
names = list(IDS)
calls_per_cell = [BITS[:, k].sum() / len(CELLS) for k in range(NS)]
unreviewed = [0] * NS


# ---------------------------------------------------------------- settings outside the sweep: re-call
def parse(sid):
    m = re.fullmatch(r"b([\d.]+)_mn(\d+)_bs([\d.]+)_t(\d+)", sid)
    if not m:
        sys.exit("cannot read setting %s (form b0.3_mn6_bs2_t4)" % sid)
    return float(m.group(1)), int(m.group(2)), float(m.group(3)) * 1e6, int(m.group(4))


EXTRA = [x for x in A.extra.split(",") if x]
if EXTRA:
    d = {}
    for l in open("results/crossovers/%s/co_calling_params.txt" % S):
        m = re.match(r"\s*([A-Za-z_ ]+?)(?:\s*\(-\w\))?:\s*(.*)", l)
        if m:
            d[m.group(1).strip()] = m.group(2).strip()
    INDIR = d.get("input tables") or "results/cell_data/%s" % S
    FIX = dict(window_af=float(d["windowAF"]), genotype=float(d["genotype"]), cell_markers=int(float(d["cell_markers"])))
    CMAP = H.read_chrom_map("results/cell_data/%s/chrom_map.tsv" % S)
    EVC = {}
    for j, e in enumerate(EV):
        EVC.setdefault((e["barcode"], e["chrom"]), []).append(j)
    CHECK = [x for x in dict.fromkeys([CUR, A.expect]) if x and x in IDS]
    RECALL = CHECK + [x for x in EXTRA if x not in CHECK]

    def recall_cell(bc):
        rows_ = H.read_cell("%s/%s.tsv" % (INDIR, bc))
        prep = H.prepare(rows_, CMAP, FIX["cell_markers"])
        out = {}
        for sid in RECALL:
            b, mn, bs, t = parse(sid)
            hit, un = [], 0
            if prep is not None:
                sm = H.smooth_cell(prep, b, FIX["window_af"])
                for c, size in CMAP:
                    if c not in sm:
                        continue
                    pos, smt = sm[c]
                    filtered = H.filter_blocks(H.genotype_blocks(pos, smt, FIX["genotype"]), bs, mn, t)
                    if not filtered:
                        continue
                    bps = H.get_breakpoints(filtered)
                    bpu = H.fine_breakpoints(bps, pos, smt) if bps else []
                    for co, fi in zip(bps, bpu):
                        st, en, pre, post, cs, ce = int(fi[0]), int(fi[1]), int(fi[2]), int(fi[3]), int(co[0]), int(co[1])
                        best, bd = None, None
                        for j in EVC.get((bc, c), []):
                            e = EV[j]
                            if e["pre"] != pre or e["post"] != post:
                                continue
                            if (cs <= e["end"] and e["start"] <= ce) or (e["cstart"] <= en and st <= e["cend"]):
                                dd = abs((st + en) - (e["start"] + e["end"]))
                                if bd is None or dd < bd:
                                    best, bd = j, dd
                        if best is None:
                            un += 1
                        else:
                            hit.append(best)
            out[sid] = (hit, un)
        return out

    with get_context("fork").Pool(A.threads) as pool:
        res = pool.map(recall_cell, CELLS)
    for sid in RECALL:
        called = np.zeros(len(EV), dtype=bool)
        un = 0
        for r in res:
            called[r[sid][0]] = True
            un += r[sid][1]
        if sid in IDS:
            k = IDS.index(sid)
            same = bool((called == BITS[:, k]).all()) and un == 0
            print("RE-CALL CHECK  %s from the pipeline's tables (%s): %s" % (
                sid, INDIR, "identical to the sweep's column" if same else
                "DIFFERS from the sweep (%d events, %d unmatched calls)" % (int((called != BITS[:, k]).sum()), un)))
            if not same:
                sys.exit(1)
        else:
            cols.append(per_cell(called))
            names.append(sid)
            calls_per_cell.append(called.sum() / len(CELLS))
            unreviewed.append(un)

# ---------------------------------------------------------------- scores and bootstrap over cells
TPm = np.array([c[0] for c in cols]).T                       # cells x settings
FPm = np.array([c[1] for c in cols]).T
FNm = np.array([c[2] for c in cols]).T


def f1(tp, fp, fn):
    with np.errstate(divide="ignore", invalid="ignore"):
        p = tp / (tp + fp)
        r = tp / (tp + fn)
        return p, r, 2 * p * r / (p + r)


tp, fp, fn = TPm.sum(0), FPm.sum(0), FNm.sum(0)
P, R, F = f1(tp, fp, fn)
rng = np.random.default_rng(1)
W = rng.multinomial(len(CELLS), [1.0 / len(CELLS)] * len(CELLS), size=A.boot).astype(float)
_, _, FB = f1(W @ TPm, W @ FPm, W @ FNm)
FB = np.nan_to_num(FB, nan=-1)
lo, hi = np.percentile(FB, 2.5, axis=0), np.percentile(FB, 97.5, axis=0)
best_share = (FB >= FB.max(axis=1, keepdims=True) - 1e-12).mean(axis=0)

nv = {v: sum(1 for e in scored if e["v"] == v) for v in ("real", "false", "unsure")}
n_auto = sum(1 for e in scored if e["v"] == "real" and e["triage"] == "auto")
print("VERDICTS  %s: %d events in %d cells; scored: %d real (%d of them auto-real), %d false; %d unsure and %d "
      "without a verdict left out; %d on masked sequences or in bad cells left out; %d missed crossovers" % (
          S, len(EV), len(CELLS), nv["real"], n_auto, nv["false"], nv["unsure"],
          sum(1 for e in scored if e["v"] is None), len(EV) - len(scored), len(miss)))
order = sorted(range(len(names)), key=lambda k: (-np.nan_to_num(F[k], nan=-1), names[k]))
show = order[:12] + [k for k in range(len(names)) if names[k] in (CUR, A.expect) or names[k] in EXTRA]
show = list(dict.fromkeys(show))
print("%-20s %9s %5s %4s %4s %9s %7s %6s %15s %6s" % ("setting", "COs/cell", "TP", "FP", "FN", "precision", "recall",
                                                     "F1", "F1 95% (cells)", "best"))
for k in show:
    tag = " <- current when swept" if names[k] == CUR else (" <- expected" if names[k] == A.expect else "")
    tag += " (re-called; %d calls outside the reviewed events)" % unreviewed[k] if names[k] in EXTRA else ""
    print("%-20s %9.2f %5d %4d %4d %9.3f %7.3f %6.3f %7.3f-%.3f %5.0f%%%s" % (
        names[k], calls_per_cell[k], tp[k], fp[k], fn[k], P[k], R[k], F[k], lo[k], hi[k], 100 * best_share[k], tag))
print("best F1 at each value (all settings scored):")
for key, i in (("base_af", 0), ("marker_num", 1), ("block_size Mb", 2), ("terminal", 3)):
    vals = {}
    for k, sid in enumerate(names):
        v = parse(sid)[i]
        v = v / 1e6 if i == 2 else v
        if not np.isnan(F[k]) and F[k] > vals.get(v, (-1, ""))[0]:
            vals[v] = (F[k], sid)
    print("  %-14s %s" % (key, "   ".join("%g: %.3f" % (v, vals[v][0]) for v in sorted(vals))))
os.makedirs(os.path.dirname(OUT) or ".", exist_ok=True)
with open(OUT, "w") as fo:
    fo.write("setting\tcos_per_cell\tTP\tFP\tFN\tprecision\trecall\tF1\tF1_lo\tF1_hi\tbest_share\tunreviewed_calls\tsource\n")
    for k in order:
        fo.write("%s\t%.3f\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.3f\t%d\t%s\n" % (
            names[k], calls_per_cell[k], tp[k], fp[k], fn[k], P[k], R[k], F[k], lo[k], hi[k], best_share[k],
            unreviewed[k], "re-called" if names[k] in EXTRA else "sweep"))
print("table: %s" % OUT)
if A.expect:
    if A.expect not in names:
        sys.exit("--expect %s was not scored" % A.expect)
    k = names.index(A.expect)
    top = np.nanmax(F)
    if F[k] < top - 1e-12:
        print("EXPECT FAILED  %s F1 %.3f, best %.3f (%s)" % (A.expect, F[k], top, names[order[0]]))
        sys.exit(1)
    print("EXPECT OK  %s ties for the best F1 (%.3f)" % (A.expect, F[k]))
