#!/usr/bin/env python3
"""review_page.py -- every crossover call of a settings sweep, cell by cell, on one HTML page you review with keys.

What it does (run from the CO_smk directory, after co_sweep.py):
  1. Pools the calls of all settings into EVENTS, per cell and sequence: calls of the same direction (REF->ALT
     or ALT->REF) whose possible positions overlap are one event (at most one call per setting; a call whose
     range touches two events joins the nearer). Each event knows which settings call it.
  2. Measures each event: how many molecules on each side agree with the genotype it switches between, up to
     the nearest event called by at least as many settings ("left 23 molecules, 96% REF").
  3. Triage. AUTO (pre-marked real, can be overridden): called by >= --auto_frac of the settings with >=
     --auto_support molecules at >= --auto_purity on each side. REVIEW: everything else. MASKED: the event's
     sequence is doubled, missing or mixed in that cell (--mask_classes; classes as calling_overview.py:
     too few = expected < 8 molecules; missing < 0.25 x expected; doubled = switches between neighbouring
     molecules >= 0.25 at >= 1.2 x expected; mixed = switches >= 0.25 at normal coverage). Masked events are
     not scored; 'g' on the page unmasks (or masks) a sequence in a cell.
  4. Pairs: the two edges of one short block (opposite directions, called by the same settings) get the same
     verdict together.
  5. Writes OUTDIR/review.html (self-contained, data inside) and OUTDIR/events.tsv.gz (the events, with the
     settings that call each as a 0/1 string in settings.tsv order).

On the page: r real, f false, u unsure (then it jumps to the next event that needs you), x clear, n/p next or
previous event that needs you, arrows: events (left/right) and cells (up/down), b bad cell, g mask this
sequence, m missed-crossover mode (click where a crossover is that no setting calls), i Mb or molecule axis,
+/- zoom, h help, e export. A scoreboard ranks the settings against your verdicts (precision, recall, F1).
Verdicts stay in the browser (if it allows) and go out as verdicts.csv with 'e'; 'import' reads one back.

Usage: review_page.py SAMPLE [--sweep qc/co_sweep/S] [--outdir qc/co_review/S] [--cells FILE | --sample_n N]
                       [--display current,b0.2_mn3_bs1_t4,b0.3_mn8_bs5_t12] [--auto_frac 0.9]
                       [--auto_support 8] [--auto_purity 0.8] [--mask_classes doubled,missing,mixed]
"""
import argparse
import bisect
import csv
import datetime
import gzip
import hashlib
import json
import os
import re
import statistics
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hapco_port as H  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("sample")
ap.add_argument("--sweep")
ap.add_argument("--outdir")
ap.add_argument("--cells")
ap.add_argument("--sample_n", type=int, default=0)
ap.add_argument("--display", default="current,b0.2_mn3_bs1_t4,b0.3_mn8_bs5_t12")
ap.add_argument("--auto_frac", type=float, default=0.9)
ap.add_argument("--auto_support", type=int, default=8)
ap.add_argument("--auto_purity", type=float, default=0.8)
ap.add_argument("--mask_classes", default="doubled,missing,mixed")
A = ap.parse_args()
S = A.sample
SW = A.sweep or "qc/co_sweep/%s" % S
OUT = A.outdir or "qc/co_review/%s" % S
T0 = time.time()
for f in ("settings.tsv", "cells.tsv", "calls.tsv.gz", "blocks.tsv.gz"):
    if not os.path.exists(os.path.join(SW, f)):
        sys.exit("review_page: %s/%s missing -- run co_sweep.py %s first" % (SW, f, S))
os.makedirs(OUT, exist_ok=True)


def params_as_run():
    d = {}
    for l in open("results/crossovers/%s/co_calling_params.txt" % S):
        m = re.match(r"\s*([A-Za-z_ ]+?)(?:\s*\(-\w\))?:\s*(.*)", l)
        if m:
            d[m.group(1).strip()] = m.group(2).strip()
    return d


RUNP = params_as_run()
INDIR = RUNP.get("input tables") or "results/cell_data/%s" % S
WINDOW_AF = float(RUNP["windowAF"])
CMAP = H.read_chrom_map("results/cell_data/%s/chrom_map.tsv" % S)
CHROMS = [c for c, _ in CMAP]
CLEN = dict(CMAP)

# ---------------------------------------------------------------- settings and cells
SET = list(csv.DictReader(open(os.path.join(SW, "settings.tsv")), delimiter="\t"))
NS = len(SET)
SIDX = {r["id"]: k for k, r in enumerate(SET)}
CURRENT = [r["id"] for r in SET if r.get("label") == "current"]
CURRENT = CURRENT[0] if CURRENT else SET[0]["id"]
DISPLAY, DNAMES = [], []
for tok in A.display.split(","):
    sid = CURRENT if tok == "current" else tok
    if sid in SIDX and SIDX[sid] not in DISPLAY:
        DISPLAY.append(SIDX[sid])
        DNAMES.append("current" if sid == CURRENT else sid)
    elif sid not in SIDX:
        print("WARNING  display setting %s is not in the sweep" % sid)
BAFS = sorted(set(float(r["base_af"]) for r in SET))

swc = list(csv.DictReader(open(os.path.join(SW, "cells.tsv")), delimiter="\t"))
GOOD = [r["barcode"] for r in swc]
CALLED = [r["barcode"] for r in swc if r["called"] == "1"]
NMOL = {r["barcode"]: int(r["molecules"]) for r in swc}

# crossovers per cell at the current setting (to spread a sample over them); the calls themselves are read
# below for the shown cells only
ncur = {bc: 0 for bc in CALLED}
with gzip.open(os.path.join(SW, "calls.tsv.gz"), "rt") as fh:
    for l in fh:
        if l.startswith(CURRENT + "\t"):
            bc = l.split("\t", 2)[1]
            if bc in ncur:
                ncur[bc] += 1

if A.cells:
    want = [l.split()[0] for l in open(A.cells) if l.strip() and not l.startswith("#")]
    SHOW = [bc for bc in want if bc in set(CALLED)]
    print("cells: %d listed in %s, %d of them called by hapCO" % (len(want), A.cells, len(SHOW)))
elif A.sample_n and A.sample_n < len(CALLED):
    ranked = sorted(CALLED, key=lambda b: (ncur[b], b))
    step = (len(ranked) - 1) / (A.sample_n - 1) if A.sample_n > 1 else 0
    idx = sorted(set(int(round(i * step)) for i in range(A.sample_n)))
    pick = set(ranked[i] for i in idx)
    SHOW = [bc for bc in CALLED if bc in pick]
    print("cells: %d of %d, spread evenly over the crossover counts at %s (%d to %d)" % (
        len(SHOW), len(CALLED), CURRENT, ncur[ranked[0]], ncur[ranked[-1]]))
else:
    SHOW = list(CALLED)
SHOWSET = set(SHOW)
if not SHOW:
    sys.exit("review_page: no cells to show")
calls_by = {}
with gzip.open(os.path.join(SW, "calls.tsv.gz"), "rt") as fh:
    rd = csv.reader(fh, delimiter="\t")
    next(rd)
    for s, bc, c, st, en, pre, post, cs, ce in rd:
        if bc in SHOWSET:
            calls_by.setdefault((bc, c), []).append((SIDX[s], int(st), int(en), int(pre), int(post), int(cs), int(ce)))

# ---------------------------------------------------------------- molecule tables, dosage (as calling_overview.py)
tables, dos_rows = {}, []
for bc in GOOD:                                      # every selected cell: the dosage classes need the medians
    p = "%s/%s.tsv" % (INDIR, bc)
    if not os.path.exists(p):
        continue
    rows = H.read_cell(p)
    if bc in SHOWSET:
        tables[bc] = rows
    by = {c: [] for c in CHROMS}
    tot = 0
    for c, p, rc, ac in rows:
        if rc != ac and c in by:
            by[c].append((p, 1 if ac > rc else 0))
            tot += 1
    for c in CHROMS:
        a = [g for _, g in sorted(by[c])]
        n = len(a)
        sw = sum(1 for x, y in zip(a, a[1:]) if x != y) / (n - 1) if n > 1 else float("nan")
        dos_rows.append([bc, c, n, tot, sw])
share = {c: statistics.median([r[2] / max(r[3], 1) for r in dos_rows if r[1] == c]) for c in CHROMS} if dos_rows else {}
DOS = {}
for bc, c, n, tot, sw in dos_rows:
    exp = tot * share[c]
    cov = n / exp if exp else float("nan")
    if exp < 8:
        cl = "too few"
    elif cov < 0.25:
        cl = "missing"
    elif sw >= 0.25:
        cl = "doubled" if cov >= 1.2 else "mixed"
    else:
        cl = "haploid"
    DOS[(bc, c)] = (cl, cov, sw, exp)
MASKC = set(x for x in A.mask_classes.split(",") if x)
ovp = "qc/calling_overview/%s_dosage.tsv" % S
if os.path.exists(ovp):
    same = tot = 0
    for r in csv.DictReader(open(ovp), delimiter="\t"):
        k = (r["barcode"], r["chrom"])
        if k in DOS:
            tot += 1
            same += DOS[k][0] == r["class"]
    print("dosage classes: %d of %d cell x sequence the same as in %s" % (same, tot, ovp))

hap = {}
ht = "qc/haplotypes/%s/haplotype_tracks.tsv.gz" % S
if os.path.exists(ht):
    with gzip.open(ht, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                hap[r["barcode"]] = float(r["haploidness"])
            except (ValueError, KeyError, TypeError):
                pass
SRC = {}
pf = "results/reference/%s/pieces.tsv" % S
if os.path.exists(pf):
    for r in csv.DictReader(open(pf), delimiter="\t"):
        SRC[r["piece"]] = "%s %.1f-%.1f Mb" % (r["source"], int(r["start"]) / 1e6, int(r["end"]) / 1e6)

blocks = {}
with gzip.open(os.path.join(SW, "blocks.tsv.gz"), "rt") as fh:
    rd = csv.reader(fh, delimiter="\t")
    next(rd)
    for s, bc, c, st, en, g in rd:
        if bc in SHOWSET and s in SIDX and SIDX[s] in DISPLAY:
            blocks.setdefault((bc, c, SIDX[s]), []).append([int(st), int(en), int(g)])
for k, sid in zip(DISPLAY, DNAMES):
    if not any(key[2] == k for key in blocks):
        print("WARNING  no blocks for display setting %s in blocks.tsv.gz (co_sweep.py --show)" % SET[k]["id"])


# ---------------------------------------------------------------- events
def pool(calls):
    """Calls of one cell and sequence -> events (see the docstring)."""
    sites = {}
    for k, st, en, pre, post, cs, ce in calls:
        x = sites.setdefault((st, en, pre, post), [set(), cs, ce])
        x[0].add(k)
        x[1], x[2] = min(x[1], cs), max(x[2], ce)
    ev = []
    for (st, en, pre, post), (ks, cs, ce) in sorted(sites.items(), key=lambda kv: (-len(kv[1][0]), kv[0][0])):
        best, bd = None, None
        for e in ev:
            if e["pre"] != pre or e["post"] != post or (e["set"] & ks):
                continue
            if (cs <= e["en"] and e["st"] <= ce) or (e["cs"] <= en and st <= e["ce"]):
                d = abs((st + en) - (e["st"] + e["en"]))
                if bd is None or d < bd:
                    best, bd = e, d
        if best is None:
            ev.append(dict(st=st, en=en, pre=pre, post=post, cs=cs, ce=ce, set=set(ks)))
        else:
            best["set"] |= ks
            best["cs"], best["ce"] = min(best["cs"], cs), max(best["ce"], ce)
    ev.sort(key=lambda e: (e["st"], e["en"]))
    return ev


def agrees(g, s):
    return s > 0.5 if g == 1 else (s < 0.5 if g == 0 else 0.25 <= s <= 0.75)


def measure(ev, pos, shr):
    n = len(pos)
    for e in ev:
        e["i"] = bisect.bisect_left(pos, e["en"])
    for e in ev:
        lvl = len(e["set"])
        lft = [x["i"] for x in ev if x is not e and len(x["set"]) >= lvl and x["i"] < e["i"]]
        rgt = [x["i"] for x in ev if x is not e and len(x["set"]) >= lvl and x["i"] > e["i"]]
        lb, rb = (max(lft) if lft else 0), (min(rgt) if rgt else n)
        for side, a, b, g in (("l", lb, e["i"], e["pre"]), ("r", e["i"], rb, e["post"])):
            seg = [s for s in shr[a:b] if s != 0.5]
            m = sum(1 for s in seg if agrees(g, s))
            e["s" + side], e["p" + side] = m, (m / len(seg) if seg else 0.0)
        e["auto"] = (len(e["set"]) / NS >= A.auto_frac and min(e["sl"], e["sr"]) >= A.auto_support
                     and min(e["pl"], e["pr"]) >= A.auto_purity)
    ev.sort(key=lambda e: (e["i"], e["st"]))
    for e in ev:
        e["pair"] = None
    for j, e in enumerate(ev):                       # pairs: the two edges of one short block
        if e["pair"] is not None or e["auto"]:
            continue
        for f in ev[j + 1:]:
            if len(f["set"]) < 0.5 * len(e["set"]):
                continue                             # a much weaker event inside the block: look past it
            if (f["pre"] == e["post"] and f["post"] == e["pre"] and f["pair"] is None and not f["auto"]
                    and len(e["set"] & f["set"]) / len(e["set"] | f["set"]) >= 0.9):
                e["pair"], f["pair"] = f, e
            break
    return ev


def sm_chars(v):
    return "1" if v == 1.0 else ("0" if v == 0.0 else "m")


EVENTS, CELLS_OUT = [], []
counts = dict(auto=0, review=0, masked=0, fringe=0, contested=0, weak=0)
for ci, bc in enumerate(SHOW):
    rows = tables.get(bc, [])
    per = {c: [] for c in CHROMS}
    for c, p, rc, ac in rows:
        if (rc > 0 or ac > 0) and c in per:
            per[c].append((p, rc / (rc + ac)))
    cell = dict(bc=bc, mol=NMOL.get(bc, 0), hap=round(hap[bc], 3) if bc in hap else None, ncur=ncur.get(bc, 0),
                chroms=[])
    num = 0
    for chi, c in enumerate(CHROMS):
        x = per[c]
        pos = [p for p, _ in x]
        shr = [s for _, s in x]
        cl, cov, sw, exp = DOS.get((bc, c), ("too few", float("nan"), float("nan"), 0))
        msk = cl in MASKC
        ev = measure(pool(calls_by.get((bc, c), [])), pos, shr) if pos else []
        sm = []
        for b in BAFS:
            sm.append("".join(sm_chars(v) for v in H.smooth(shr, b, WINDOW_AF)) if len(shr) >= 2 else
                      "".join(sm_chars(v) for v in shr))
        ids = []
        for e in ev:
            num += 1
            e["id"], e["ci"], e["ch"], e["num"] = len(EVENTS), ci, chi, num
            EVENTS.append(e)
            ids.append(e["id"])
            fr = len(e["set"]) / NS
            if msk:
                counts["masked"] += 1
            elif e["auto"]:
                counts["auto"] += 1
            else:
                counts["review"] += 1
                counts["fringe" if fr < 0.1 else ("weak" if fr >= A.auto_frac else "contested")] += 1
        deltas, prev = [], 0
        for p in pos:
            deltas.append(p - prev)
            prev = p
        cell["chroms"].append(dict(
            c=c, len=CLEN[c], src=SRC.get(c, ""), cls=cl, cov=None if cov != cov else round(cov, 2),
            sw=None if sw != sw else round(sw, 2), mask=msk, p=deltas, s=[int(round(100 * s)) for s in shr], sm=sm,
            blk=[blocks.get((bc, c, k), []) for k in DISPLAY], ev=ids))
    CELLS_OUT.append(cell)

cpc = [0] * NS
for (bc, c), cl in calls_by.items():
    if bc in SHOWSET:
        for r in cl:
            cpc[r[0]] += 1
ev_out = []
for e in EVENTS:
    bits = "".join("1" if k in e["set"] else "0" for k in range(NS))
    ev_out.append(dict(id=e["id"], ci=e["ci"], ch=e["ch"], num=e["num"], st=e["st"], en=e["en"], cs=e["cs"],
                       ce=e["ce"], pre=e["pre"], post=e["post"], n=len(e["set"]), bits=bits, sl=e["sl"], sr=e["sr"],
                       pl=round(e["pl"], 3), pr=round(e["pr"], 3), auto=bool(e["auto"]),
                       pair=e["pair"]["id"] if e["pair"] is not None else -1, i=e["i"]))
GN = {1: "REF", 0: "ALT", -1: "mixed"}
tsv = ["event\tbarcode\tchrom\tstart\tend\tcstart\tcend\tdirection\tn_settings\tfraction\ttriage\tleft_molecules\t"
       "right_molecules\tleft_purity\tright_purity\tpair\tsettings"]
for e in ev_out:
    c = CELLS_OUT[e["ci"]]
    ch = c["chroms"][e["ch"]]
    tri = "masked" if ch["mask"] else ("auto" if e["auto"] else "review")
    tsv.append("\t".join(map(str, [e["id"], c["bc"], ch["c"], e["st"], e["en"], e["cs"], e["ce"],
                                   GN.get(e["pre"], "?") + ">" + GN.get(e["post"], "?"), e["n"],
                                   "%.3f" % (e["n"] / NS), tri, e["sl"], e["sr"], e["pl"], e["pr"],
                                   e["pair"] if e["pair"] >= 0 else "", e["bits"]])))
body = "\n".join(tsv) + "\n"
with gzip.open(os.path.join(OUT, "events.tsv.gz"), "wt") as fo:
    fo.write(body)
page_id = hashlib.md5(body.encode()).hexdigest()[:10]
try:
    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE,
                            stderr=subprocess.DEVNULL, universal_newlines=True).stdout.strip()
    mine = ["workflow/scripts/%s" % f for f in ("hapco_port.py", "co_sweep.py", "review_page.py", "review_page.html")]
    if commit and subprocess.run(["git", "status", "--porcelain", "--"] + mine, stdout=subprocess.PIPE,
                                 stderr=subprocess.DEVNULL, universal_newlines=True).stdout.strip():
        commit += " + uncommitted sweep/review scripts"
except OSError:
    commit = ""
data = dict(sample=S, page_id=page_id, built=datetime.datetime.now().isoformat(timespec="minutes"), commit=commit,
            settings=[dict(id=r["id"], b=float(r["base_af"]), mn=int(r["marker_num"]), bs=int(r["block_size"]),
                           t=int(r["terminal_marker_num"]), cur=r["id"] == CURRENT,
                           cpc=round(cpc[k] / max(len(SHOW), 1), 2)) for k, r in enumerate(SET)],
            display=DISPLAY, dnames=DNAMES, baf=BAFS,
            rule=dict(frac=A.auto_frac, support=A.auto_support, purity=A.auto_purity),
            mask_classes=sorted(MASKC), cells=CELLS_OUT, events=ev_out,
            n_good=len(GOOD), n_called=len(CALLED), sweep=SW)
blob = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
HTML = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "review_page.html")).read() \
    if os.path.exists(os.path.join(os.path.dirname(os.path.abspath(__file__)), "review_page.html")) else None
if HTML is None:
    sys.exit("review_page: review_page.html (the page template) must sit next to this script")
page = HTML.replace("/*DATA*/", blob)
with open(os.path.join(OUT, "review.html"), "w") as fo:
    fo.write(page)
n_ev = len(ev_out)
print("REVIEW PAGE  %s: %d cells (of %d hapCO calls; %d selected), %d settings, %d events -> %s/review.html (%.1f MB) in %.0f s"
      % (S, len(SHOW), len(CALLED), len(GOOD), NS, n_ev, OUT, len(page) / 1e6, time.time() - T0))
print("  auto-accepted (>= %.0f%% of settings, >= %d molecules at >= %.0f%% each side): %d" % (
    100 * A.auto_frac, A.auto_support, 100 * A.auto_purity, counts["auto"]))
print("  need your verdict: %d  (called by < 10%% of settings: %d; 10-%.0f%%: %d; >= %.0f%% but thin support: %d)" % (
    counts["review"], counts["fringe"], 100 * A.auto_frac, counts["contested"], 100 * A.auto_frac, counts["weak"]))
nmask = sum(1 for c in CELLS_OUT for ch in c["chroms"] if ch["mask"])
print("  masked: %d events on %d cell x sequence (%s)" % (counts["masked"], nmask, ", ".join(sorted(MASKC))))
print("  pairs (both edges of one short block, one verdict): %d" % (sum(1 for e in ev_out if e["pair"] >= 0) // 2))
