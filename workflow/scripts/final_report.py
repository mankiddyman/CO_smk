#!/usr/bin/env python3
"""final_report.py -- the final evidence report on D. paradoxa's karyotype as its own pollen inherit it, against its
two-haplotype assembly (the chr1/chr2 arms A-D, the X arm of chr3/chr4, the chr5 x chr6 coupling, the inversions
between the haplotypes); the crossover landscapes of D. paradoxa and D. binata; and what the CARP repeats do and do
not explain.

Made for a presentation: every figure is also saved as its own numbered PNG in OUT/fig/ (new diagrams and the
pipeline's own figures, in the order they appear), each section opens with a one-line slide message, and
OUT/slides_outline.md lists the section messages and the figures in order, with what to look for in each.

Reads only what the pipeline wrote (run from the CO_smk root; every input is optional, a missing one drops its part):
  reports/supervisor_v5/{numbers.txt, fig/*.png}                         supervisor_report.py v5
  reports/cosegregation_2n/{karyotype_test,cosegregation_2n,hic_scan}_P.tsv and PNGs, logs/cosegregation_2n_P.log
  reports/x_arm/x_arm_P.{tsv,png}                                        x_arm_dosage.py
  reports/repeat_context/{rdna_loci,fish_probes}_S.csv and PNGs          repeat_context.py
  reports/repeat_recombination{,_near5}/{windows,assoc,tests,breakpoints}_S.csv, summary_S.txt, PNGs
  qc/translocations/HAP1/HAP1_synteny_blocks.tsv                          translocation_map.py
  qc/reference/P/breakpoints.csv, config/split_plans/P.csv, results/reference/P/{pieces.tsv, MANIFEST.txt},
  config/samples.csv
Writes OUT/report.html (self-contained), OUT/numbers.txt, OUT/fig/NN_*.png, OUT/slides_outline.md.
Usage: python workflow/scripts/final_report.py [--out reports/final_report]. About a minute.
"""
import argparse
import csv
import html
import math
import os
import re
import shutil
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon, Rectangle

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from report_kit import Report, binom_p, fi, ff, ok   # noqa: E402

ap = argparse.ArgumentParser(description="Final report: karyotype vs assembly, crossover landscapes, repeats")
ap.add_argument("--paradoxa", default="Dparadoxa_split", help="the D. paradoxa crossover sample (cut reference)")
ap.add_argument("--binata", default="Dbinata_hap1", help="the D. binata crossover sample")
ap.add_argument("--hap1", default="Dparadoxa_hap1", help="synteny blocks: qc/translocations/HAP1/HAP1_synteny_blocks.tsv")
ap.add_argument("--v5", default="reports/supervisor_v5", help="the supervisor report v5 folder (numbers.txt, fig/)")
ap.add_argument("--out", default="reports/final_report")
A = ap.parse_args()
P, B = A.paradoxa, A.binata
MB = 1e6
SPS = [("binata", B), ("paradoxa", P)]

# ------------------------------------------------------------------ one palette for the whole report
SP_COL = {"binata": "#0F6E56", "paradoxa": "#534AB7"}
SP_NAME = {"binata": "D. binata", "paradoxa": "D. paradoxa"}
HAP_COL = {"hap1": "#378ADD", "hap2": "#D85A30"}
ARM_COL = {"L1": "#5DCAA5", "L2": "#F0997B", "P": "#AFA9EC", "Q": "#FAC775", "X": "#85B7EB", "core": "#D3D1C7"}
CHR_COL = {"chr3": "#E6D3A3", "chr5": "#C2D69B", "chr6": "#E7B7C8"}
ST_COL = {"collinear": "#B4B2A9", "inverted": "#C2185B", "translocated": "#1D9E75"}
CPL_COL = "#7A2E8E"            # the chr5 x chr6 coupling
INK, MUTED, BAD, LINE = "#2C2C2A", "#6B6A64", "#A32D2D", "#B4B2A9"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})

R = Report("<i>D. paradoxa</i>: the karyotype its pollen inherit against its assembly, the crossover landscapes, "
           "and what the repeats explain (final report)", A.out, "workflow/scripts/final_report.py")
R.parts.append("<style>.slide{background:#F1EFE8;border-left:4px solid #534AB7;padding:6px 10px;font-weight:bold;"
               "margin:10px 0}.small{font-size:11px;color:#5F5E5A}ul{margin:4px 0 10px}li{margin:2px 0}</style>")
TOC, SLIDES, FIGS, SECTION_MSGS = [], [], {}, []
NFIG = [0]


# ------------------------------------------------------------------ helpers
def fnum(x):
    try:
        v = float(str(x).strip())
        return v if math.isfinite(v) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def H(s):
    return html.escape(str(s), quote=False)


def pf(p):
    """p-value text (plain: escape with H() inside HTML)."""
    if not ok(p):
        return "p –"
    if p < 0.001:
        return "p < 0.001"
    return "p %.3f" % p if p < 0.01 else "p %.2f" % p


def pe(p):
    """An exact p-value (binomial tests) with two significant digits."""
    return "p –" if not ok(p) else ("p %.2g" % p)


def pc(x, d=0):
    return "–" if not ok(x) else "%.*f%%" % (d, 100.0 * x)


def sg(x, d=2):
    return "–" if not ok(x) else "%+.*f" % (d, x)


def num_of(c):
    m = re.match(r"chr(\d+)", str(c))
    return int(m.group(1)) if m else 0


def hap_of(c):
    m = re.search(r"_(hap[12])", str(c))
    return m.group(1) if m else ""


def wmedian(v, w):
    v, w = np.asarray(v, float), np.asarray(w, float)
    if not len(v):
        return 0.0
    o = np.argsort(v)
    c = np.cumsum(w[o])
    return float(v[o][np.searchsorted(c, c[-1] / 2.0)])


def h2(id_, t):
    TOC.append((id_, t))
    R.parts.append("<h2 id='%s'>%s</h2>" % (id_, t))


def slide(msg):
    R.parts.append("<p class='slide'>Slide message: %s</p>" % msg)
    SECTION_MSGS.append((TOC[-1][1] if TOC else "", strip_tags(msg)))
    return msg


def look(t):
    R.look(t)
    if SLIDES and SLIDES[-1][3] is None:
        SLIDES[-1][3] = strip_tags(t)


def bullets(items):
    items = [i for i in items if i]
    if items:
        R.parts.append("<ul>%s</ul>" % "".join("<li>%s</li>" % i for i in items))


def missing_note(what):
    R.parts.append("<p class='missing'>%s</p>" % H(what))


def need(*paths):
    miss = [p for p in paths if not R.have(p)]
    if miss:
        raise FileNotFoundError("missing input: " + ", ".join(miss))


def strip_tags(s):
    return re.sub(r"<[^>]+>", "", html.unescape(s)) if s else s


def emb(key, name, title, message, width=100):
    """Embed a figure: a matplotlib figure built earlier (FIGS[key]) or an existing PNG (key = a path). Saved as
    fig/NN_name.png in the order of appearance, and listed for the slide outline."""
    is_path = isinstance(key, str) and key.endswith(".png")
    src = key if is_path else FIGS.get(key)
    if src is None or (is_path and not R.have(src)):
        R.img(src, width)
        return None
    NFIG[0] += 1
    dst = os.path.join(R.fig_dir, "%02d_%s.png" % (NFIG[0], name))
    if is_path:
        shutil.copyfile(src, dst)
    else:
        src.savefig(dst, dpi=170, bbox_inches="tight")
        plt.close(src)
    R.img(dst, width)
    SLIDES.append([os.path.basename(dst), title, strip_tags(message), None])
    return dst


# ------------------------------------------------------------------ inputs
SAMPLES = {}
if R.have("config/samples.csv"):
    SAMPLES = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}


def two_hap_lengths(sample):
    """Lengths of the two-haplotype assembly behind a crossover reference (its MANIFEST's source_fasta index)."""
    row = SAMPLES.get(sample)
    cands = [os.path.join("results/reference", sample, "MANIFEST.txt")]
    if row:
        cands.insert(0, os.path.join(os.path.dirname(row["assembly_fasta"]), "MANIFEST.txt"))
    for man in cands:
        if os.path.exists(man):
            m = re.search(r"^source_fasta:\s*(\S+)", open(man).read(), re.M)
            if m and os.path.exists(m.group(1) + ".fai"):
                R.have(man)
                R.have(m.group(1) + ".fai")
                return {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(m.group(1) + ".fai") if l.strip()}
    R.have(cands[0])
    return {}


L2N = R.guarded("two-haplotype lengths", lambda: two_hap_lengths(P)) or {}
COPIES = sorted([c for c in L2N if re.match(r"^chr\d+_hap[12]$", c)], key=lambda c: (num_of(c), c))


def read_pieces():
    p = os.path.join("results/reference", P, "pieces.tsv")
    need(p)
    return [dict(piece=r["piece"], source=r["source"], start=int(r["start"]) - 1, end=int(r["end"]))
            for r in csv.DictReader(open(p), delimiter="\t")]


PIECES = R.guarded("pieces", read_pieces) or []
REFC = []
for pc_ in PIECES:
    if pc_["source"] not in REFC:
        REFC.append(pc_["source"])
CUTS = {c: [pc_["start"] / MB for pc_ in PIECES if pc_["source"] == c and pc_["start"] > 0] for c in REFC}


def read_plan():
    p = "config/split_plans/%s.csv" % P
    if not R.have(p):
        return {}
    return {r["piece"]: r for r in csv.DictReader(l for l in open(p) if not l.startswith("#"))}


PLAN = R.guarded("split plan", read_plan) or {}

V5 = {"binata": {}, "paradoxa": {}}
if R.have(os.path.join(A.v5, "numbers.txt")):
    for l in open(os.path.join(A.v5, "numbers.txt")):
        if "\t" in l:
            k, v = l.rstrip("\n").split("\t", 1)
            sp, _, key = k.partition(".")
            if sp in V5:
                V5[sp][key] = v


def v5n(sp, key):
    return fnum(V5[sp].get(key))


def v5s(sp, key):
    return V5[sp].get(key, "–")


def v5c(sp, key):
    """A v5 'chance' string without nested brackets: '0.6 (95% 0–2)' -> '0.6, 95% 0–2'."""
    return re.sub(r" \((95%[^)]*)\)", r", \1", v5s(sp, key))


def v5fig(name):
    return os.path.join(A.v5, "fig", name)


def shape_word(sp):
    s = v5s(sp, "shape")
    for k, w in (("U-shaped", "packed at the chromosome ends"), ("flat", "spread along the chromosomes, no end bias"),
                 ("ends depleted", "scarce at the chromosome ends"),
                 ("mild end bias", "slightly more at the chromosome ends")):
        if s.startswith(k):
            return w
    return "of unknown shape"


# ---- the karyotype test (cosegregation_2n.py)
END_RE = r"(\S+)'s (\S+) \((\S+) (up to|from) ([\d.]+) Mb\)"
PART_RE = re.compile(END_RE + r" goes with " + END_RE + r" \(phi ([+-]?[\d.]+)(?:; Hi-C ([\d.]+)x)?\)")
KT_TSV = "reports/cosegregation_2n/karyotype_test_%s.tsv" % P
COS_TSV = "reports/cosegregation_2n/cosegregation_2n_%s.tsv" % P
SCAN_TSV = "reports/cosegregation_2n/hic_scan_%s.tsv" % P
COS_LOG = "logs/cosegregation_2n_%s.log" % P


def parse_kt():
    need(KT_TSV)
    J, C, CT = [], None, None
    for r in csv.DictReader(open(KT_TSV), delimiter="\t"):
        d = dict(kind=r["kind"], what=r["what"], phi=fnum(r["phi"]), lo=fnum(r["ci_lo"]), hi=fnum(r["ci_hi"]),
                 detail=r["detail"])
        m = re.search(r"(\d+) of (\d+) nuclei agree", r["detail"])
        d["agree"], d["n"] = (int(m.group(1)), int(m.group(2))) if m else (None, None)
        if r["kind"] in ("join", "control"):
            m = re.match(r"^(\S+?)(?: \(([A-Z])\))?: (\S+) \| (\S+) at ([\d.]+) Mb", r["what"])
            if not m:
                continue
            d.update(copy=m.group(1), letter=m.group(2) or "", left=m.group(3), right=m.group(4),
                     pos=float(m.group(5)))
            m = re.search(r"Hi-C across the join ([\d.]+)x", r["detail"])
            d["hic"] = float(m.group(1)) if m else float("nan")
            d["partners"] = [dict(own=g[0], arm=g[1], copy=g[2], side=g[3], pos=float(g[4]), o_own=g[5], o_arm=g[6],
                                  o_copy=g[7], o_side=g[8], o_pos=float(g[9]), phi=float(g[10]), hic=fnum(g[11]))
                             for g in PART_RE.findall(r["detail"])]
            if r["kind"] == "join":
                J.append(d)
            else:
                m = re.search(r"c ([\d.]+), R2 (-?[\d.]+)", r["detail"])
                d["c"], d["r2"] = (float(m.group(1)), float(m.group(2))) if m else (float("nan"), float("nan"))
                CT = d
        elif r["kind"] == "couple":
            m = re.match(r"^(\S+) ([\d.]+) Mb x (\S+) ([\d.]+) Mb", r["what"])
            if m:
                d.update(p=m.group(1), pt_p=float(m.group(2)), q=m.group(3), pt_q=float(m.group(4)))
            m = re.search(r"explains R2 (-?[\d.]+) of the block, ends joined (-?[\d.]+)", r["detail"])
            d["r2"], d["r2_ends"] = (float(m.group(1)), float(m.group(2))) if m else (float("nan"), float("nan"))
            m = re.search(r"Hi-C: strongest (\S+) x (\S+) block ([\d.]+)x trans at ([\d.]+) / ([\d.]+) Mb "
                          r"\(rank (\d+) of (\d+) copy pairs\)", r["detail"])
            d["hic"] = dict(a=m.group(1), b=m.group(2), x=float(m.group(3)), pa=float(m.group(4)),
                            pb=float(m.group(5)), rank=int(m.group(6)), of=int(m.group(7))) if m else None
            C = d
    return J, C, CT


JOINS, COUPLE, CONTROL = R.guarded("karyotype test", parse_kt) or ([], None, None)
JBY = {j["copy"]: j for j in JOINS}
LETTERS = {j["letter"]: j for j in JOINS if j["letter"]}
ABCD = all(k in LETTERS for k in "ABCD")


def parse_coupling_log():
    need(COS_LOG)
    out = dict(nuclei=None)
    for l in open(COS_LOG):
        m = re.search(r"COUPLING (\S+) x (\S+): point at ([\d.]+) Mb \(95% ([\d.]+)-([\d.]+)\) x ([\d.]+) Mb "
                      r"\(95% ([\d.]+)-([\d.]+)\), c ([\d.]+), R2 (-?[\d.]+); ends joined R2 (-?[\d.]+)", l)
        if m:
            g = m.groups()
            out.update(p=g[0], q=g[1], pt_p=float(g[2]), lo_p=float(g[3]), hi_p=float(g[4]), pt_q=float(g[5]),
                       lo_q=float(g[6]), hi_q=float(g[7]), c=float(g[8]), r2=float(g[9]), r2_ends=float(g[10]))
        m = re.match(r"^%s: (\d+) nuclei" % re.escape(P), l)
        if m:
            out["nuclei"] = int(m.group(1))
    return out


CPL = R.guarded("coupling log", parse_coupling_log) or {}


def read_bps():
    p = "qc/reference/%s/breakpoints.csv" % P
    need(p)
    return {r["chrom"]: r for r in csv.DictReader(open(p))}


BPS = R.guarded("breakpoints.csv", read_bps) or {}


def own_name(c):
    """A copy without a join is named after the arm its joined homolog keeps (chr4_hap2 = 'core'), else its number."""
    for j in JOINS:
        if num_of(j["copy"]) == num_of(c) and j["copy"] != c:
            b = BPS.get(j["copy"], {})
            for side, arm in (("left", j["left"]), ("right", j["right"])):
                if num_of(b.get("%s_homolog" % side, "")) == num_of(c):
                    return arm
    return "chr%d" % num_of(c)


def segs_of(c):
    """The arms of a copy of the two-haplotype assembly, in Mb: [(start, end, name)]."""
    L = L2N.get(c, 0) / MB
    j = JBY.get(c)
    if j:
        return [(0.0, j["pos"], j["left"]), (j["pos"], L, j["right"])]
    return [(0.0, L, own_name(c))]


def col_of(name):
    return ARM_COL.get(name) or CHR_COL.get(name) or "#EFEDE6"


def arm_copies(arm):
    return [c for c in COPIES if any(nm == arm for _, _, nm in segs_of(c))]


# ---- synteny blocks between the haplotypes (translocation_map.py)
SYN_TSV = "qc/translocations/%s/%s_synteny_blocks.tsv" % (A.hap1, A.hap1)


def read_syn():
    need(SYN_TSV)
    s = pd.read_csv(SYN_TSV, sep="\t", keep_default_na=False)
    s["kind"] = s["kind"].astype(str).str.strip().replace({"": "collinear", "nan": "collinear"}).str.lower()
    for c in ("q_start", "q_end", "t_start", "t_end", "aligned"):
        s[c] = pd.to_numeric(s[c], errors="coerce")
    return s


SYN = R.guarded("synteny blocks", read_syn)


def main_strand(t, q):
    if SYN is None:
        return "+"
    x = SYN[(SYN.hap1 == t) & (SYN.hap2 == q)]
    if not len(x):
        return "+"
    return x.groupby("strand").aligned.sum().idxmax()


def merged_inversions(copy, min_mb=0.5):
    """Inverted blocks on one copy, merged when closer than 1 Mb: [(start, end)] in Mb of that copy."""
    if SYN is None:
        return []
    h = hap_of(copy)
    x = SYN[(SYN.kind == "inverted") & (SYN.hap1 == copy if h == "hap1" else SYN.hap2 == copy)]
    iv = sorted((r.t_start / MB, r.t_end / MB) if h == "hap1" else (r.q_start / MB, r.q_end / MB)
                for r in x.itertuples())
    out = []
    for a, b in iv:
        if out and a - out[-1][1] <= 1.0:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return [(a, b) for a, b in out if b - a >= min_mb]


# ---- X arm copies (x_arm_dosage.py)
XARM_TSV, XARM_PNG = "reports/x_arm/x_arm_%s.tsv" % P, "reports/x_arm/x_arm_%s.png" % P


# ---- repeats
def rr_path(kind, sp_sample, near5=False):
    d = "reports/repeat_recombination_near5" if near5 else "reports/repeat_recombination"
    ext = "txt" if kind == "summary" else ("png" if kind in ("recombination_vs_repeats", "recombination_tests",
                                                           "breakpoints_vs_repeats") else "csv")
    return os.path.join(d, "%s_%s.%s" % (kind, sp_sample, ext))


def rc_path(kind, sp_sample):
    ext = "csv" if kind in ("rdna_loci", "fish_probes") else "png"
    return os.path.join("reports/repeat_context", "%s_%s.%s" % (kind, sp_sample, ext))


def load_rr(sample, near5=False):
    out = {}
    for k in ("windows", "assoc", "tests", "breakpoints"):
        p = rr_path(k, sample, near5)
        if k == "breakpoints" and sample != P:        # written only where there are joins or synteny blocks
            if not os.path.exists(p):
                continue
        if R.have(p):
            out[k] = pd.read_csv(p, keep_default_na=True)
    p = rr_path("summary", sample, near5)
    if R.have(p):
        out["summary"] = s = open(p).read()
        m = re.search(r"crossovers: (\d+) calls in (\d+) nuclei; median call interval ([\d.]+) Mb", s)
        if m:
            out["calls"], out["nuclei"], out["interval_mb"] = int(m.group(1)), int(m.group(2)), float(m.group(3))
        m = re.search(r"CARP ([\d.]+)", s)
        out["carp"] = m.group(1) if m else "?"
        m = re.search(r"windows: ([\d.]+) Mb", s)
        out["win_mb"] = float(m.group(1)) if m else float("nan")
        m = re.search(r"DEAD SPOTS: runs >= ([\d.]+) Mb without a crossover midpoint \(chance: ([\d.]+) such runs, "
                      r"95% ([\d.]+)-([\d.]+)\)", s)
        if m:
            out["dead_mb"], out["dead_chance"] = float(m.group(1)), "%s, 95%% %s–%s" % m.group(2, 3, 4)
    if "windows" in out:
        w = out["windows"]
        for c in ("coupling", "dead", "at_cut"):
            if c in w:
                w[c] = w[c].astype(str).str.lower().isin(["true", "1"])
    return out


RRD = {sp: R.guarded("repeat_recombination %s" % sp, lambda s=s: load_rr(s)) or {} for sp, s in SPS}
RR5 = R.guarded("repeat_recombination near5", lambda: load_rr(P, True)) or {}
DEAD_MB = RRD.get("paradoxa", {}).get("dead_mb", v5n("paradoxa", "dead_spots_min_mb"))


def load_rc(sample):
    out = {}
    for k in ("rdna_loci", "fish_probes"):
        p = rc_path(k, sample)
        if R.have(p):
            out[k] = pd.read_csv(p, keep_default_na=False)
    return out


RCD = {sp: R.guarded("repeat_context %s" % sp, lambda s=s: load_rc(s)) or {} for sp, s in SPS}


def assoc_row(sp, var):
    a = RRD.get(sp, {}).get("assoc")
    if a is None:
        return None
    x = a[a.variable == var]
    return x.iloc[0] if len(x) else None


def test_row(sp, test_prefix, what, src=None):
    t = (src or RRD.get(sp, {})).get("tests")
    if t is None:
        return None
    x = t[t.test.astype(str).str.startswith(test_prefix) & (t.what == what)]
    return x.iloc[0] if len(x) else None


def array_mb(sp):
    t = RRD.get(sp, {}).get("tests")
    if t is None:
        return float("nan")
    m = t.test.astype(str).str.extract(r"crossover density within ([\d.]+) Mb")[0].dropna()
    return float(m.iloc[0]) if len(m) else float("nan")


def rep_verdict(sp):
    """Do the repeats explain the crossover landscape, against repeat tracks shifted along their chromosome?"""
    alone = test_row(sp, "deviance explained: rotation null", "repeats alone")
    after = test_row(sp, "deviance explained: rotation null", "repeats after position + markers")
    sat = test_row(sp, "crossover density within", "satellite arrays >= 50 kb")
    a = assoc_row(sp, "all repeats")
    v = dict(alone=alone, after=after, sat=sat, rep=a,
             alone_sig=alone is not None and ok(alone.p) and alone.p < 0.05,
             after_sig=after is not None and ok(after.p) and after.p < 0.05,
             sat_avoid=sat is not None and ok(sat.obs) and ok(sat.p) and sat.obs < 1 and sat.p < 0.05)
    if v["alone_sig"] and v["after_sig"]:
        v["word"] = "yes"
    elif v["alone_sig"] or v["sat_avoid"]:
        v["word"] = "yes, but tangled with position"
    else:
        v["word"] = "no"
    return v


RV = {sp: R.guarded("repeat verdict %s" % sp, lambda sp=sp: rep_verdict(sp)) or {} for sp, _ in SPS}


# =================================================================== derived results
# ---- joins
def joins_summary():
    good = [j for j in JOINS if ok(j["phi"]) and j["phi"] >= 0.5]
    bad = [j for j in JOINS if ok(j["phi"]) and j["phi"] < 0.5]
    hic_good = [j["hic"] for j in good if ok(j["hic"]) and j["hic"] > 0.05]
    out = dict(good=good, bad=bad, hic_lo=min(hic_good) if hic_good else float("nan"),
               hic_hi=max(hic_good) if hic_good else float("nan"))
    pairs = {}
    for j in bad:
        for pr in j["partners"]:
            key = frozenset([(pr["own"], pr["arm"]), (pr["o_own"], pr["o_arm"])])
            pairs.setdefault(key, pr)
    out["pairs"] = list(pairs.values())
    for j in JOINS:
        R.put("join.%s.phi" % (j["letter"] or j["copy"]), j["phi"])
        R.put("join.%s.hic_x" % (j["letter"] or j["copy"]), j["hic"])
    return out


JS = R.guarded("joins", joins_summary) or dict(good=[], bad=[], pairs=[], hic_lo=float("nan"), hic_hi=float("nan"))


def hic_supports(x):
    """True: as high as the confirmed joins; False: about unrelated chromosomes; None: unmappable or unclear."""
    if not ok(x) or x < 0.05:
        return None
    if ok(JS["hic_lo"]) and x >= 0.5 * JS["hic_lo"]:
        return True
    return False if x <= 3 else None


def hic_word(x):
    s = hic_supports(x)
    if not ok(x):
        return "not measured"
    if x < 0.05:
        return "%.1f× (no usable reads: an unmappable stretch reads 0, so this says nothing)" % x
    if s:
        return "%.1f× (as high as the confirmed joins, %.0f–%.0f×: the tissue reads it as one chromosome)" % (
            x, JS["hic_lo"], JS["hic_hi"])
    if s is False:
        return "%.1f× (the level of unrelated chromosomes, ~1×)" % x
    return "%.1f×" % x


def pair_of(a1, a2):
    return next((p for p in JS["pairs"] if {p["arm"], p["o_arm"]} == {a1, a2}), None)


# ---- clopper-pearson
def binom_cdf(k, n, p):
    if k < 0:
        return 0.0
    if k >= n:
        return 1.0
    return sum(math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) +
                        (i * math.log(p) if i else 0.0) + ((n - i) * math.log1p(-p) if n - i else 0.0))
               for i in range(k + 1))


def cp_ci(k, n, level=0.95):
    if n == 0:
        return float("nan"), float("nan")
    a = (1 - level) / 2

    def bis(f):
        lo, hi = 1e-12, 1 - 1e-12
        for _ in range(80):
            mid = (lo + hi) / 2
            if f(mid):
                hi = mid
            else:
                lo = mid
        return (lo + hi) / 2
    lo = 0.0 if k == 0 else bis(lambda p: 1 - binom_cdf(k - 1, n, p) >= a)
    hi = 1.0 if k == n else bis(lambda p: binom_cdf(k, n, p) <= a)
    return lo, hi


# ---- X arm
def x_stats():
    need(XARM_TSV)
    x = pd.read_csv(XARM_TSV, sep="\t", keep_default_na=False)
    for c in ("X_coverage", "X_alt_share", "X_switch"):
        x[c] = pd.to_numeric(x[c], errors="coerce")
    s = dict(n_all=len(x))
    both = x[x.chr4_next_to_X.isin(["REF", "ALT"]) & x.chr3_end.isin(["REF", "ALT"])]
    s["N"] = len(both)
    s["RR"] = int(((both.chr4_next_to_X == "REF") & (both.chr3_end == "REF")).sum())
    s["AA"] = int(((both.chr4_next_to_X == "ALT") & (both.chr3_end == "ALT")).sum())
    s["RA"] = int(((both.chr4_next_to_X == "REF") & (both.chr3_end == "ALT")).sum())   # X twice
    s["AR"] = int(((both.chr4_next_to_X == "ALT") & (both.chr3_end == "REF")).sum())   # X never
    s["table"] = pd.crosstab(x.chr4_next_to_X, x.chr3_end)
    t, nv, N = s["RA"], s["AR"], s["N"]
    r_ref, c_ref = (s["RR"] + s["RA"]) / max(N, 1), (s["RR"] + s["AR"]) / max(N, 1)
    s["exp"] = dict(RR=N * r_ref * c_ref, AA=N * (1 - r_ref) * (1 - c_ref), RA=N * r_ref * (1 - c_ref),
                    AR=N * (1 - r_ref) * c_ref)
    called = both[both.agree.isin(["yes", "no"])]
    s["agree"], s["agree_n"] = int((called.agree == "yes").sum()), len(called)
    tw = x[x.observed_X == "2: both"]
    one = x[x.observed_X.isin(["1: hap1's X", "1: hap2's X", "1: crossover in X"])]
    s["twice_obs_all"] = len(tw)
    s["twice_from_RA"] = int(((tw.chr4_next_to_X == "REF") & (tw.chr3_end == "ALT")).sum())
    s["twice_cov"] = (tw.X_coverage.min(), tw.X_coverage.median(), tw.X_coverage.max()) if len(tw) else None
    s["twice_switch"] = (tw.X_switch.min(), tw.X_switch.max()) if len(tw) else None
    s["one_switch_max"] = one.X_switch.max() if len(one) else float("nan")
    s["one_cov"] = (one.X_coverage.min(), one.X_coverage.max()) if len(one) else None
    s["n_one"] = len(one)
    s["n_xco"] = int((x.observed_X == "1: crossover in X").sum())
    unres = x[x.predicted_X == "unresolved"]
    s["unres"] = len(unres)
    s["unres_one"] = int(unres.observed_X.isin(["1: hap1's X", "1: hap2's X", "1: crossover in X"]).sum())
    s["doubled_masks"] = int((x.co_masks_class == "doubled").sum())
    s["p_twice_never"] = binom_p(min(t, nv), t + nv, 0.5) if t + nv else float("nan")
    s["lost"] = nv < t and ok(s["p_twice_never"]) and s["p_twice_never"] < 0.05
    s["share"] = t / N if N else float("nan")
    s["share_ci"] = cp_ci(t, N)
    s["p_indep"] = binom_p(t, N, 1 / 3.0) if N else float("nan")
    s["not_indep"] = ok(s["p_indep"]) and s["p_indep"] < 0.05 and s["share"] < 1 / 3.0
    s["xco_share"] = s["n_xco"] / s["n_one"] if s["n_one"] else float("nan")
    s["xco_meioses"] = min(1.0, 2 * s["xco_share"]) if ok(s["xco_share"]) else float("nan")
    a2 = (1 - s["xco_meioses"]) / 2 if ok(s["xco_meioses"]) else float("nan")
    s["chain_only_share"] = a2 / (2 - a2) if ok(a2) else float("nan")
    s["chain_only_adjacent"] = a2
    s["p_chain_only"] = binom_p(t, N, s["chain_only_share"]) if N and ok(a2) else float("nan")

    def adj(sh):
        return 2 * sh / (1 + sh) if ok(sh) else float("nan")      # share = a / (2 - a)  ->  a = 2 share / (1 + share)
    s["adjacent"] = adj(s["share"])
    s["adjacent_ci"] = (adj(s["share_ci"][0]), adj(s["share_ci"][1]))
    s["abort"] = s["adjacent"] / 2
    for k in ("N", "RR", "AA", "RA", "AR", "agree", "agree_n", "n_xco", "n_one", "unres", "unres_one",
              "doubled_masks", "p_twice_never", "share", "p_indep", "xco_meioses", "chain_only_share",
              "p_chain_only", "adjacent", "abort"):
        R.put("x.%s" % k, s[k])
    R.put("x.adjacent_ci", "%.3f–%.3f" % s["adjacent_ci"])
    return s


XS = R.guarded("X arm", x_stats)


# ---- chr5 x chr6
def cpl_stats():
    s = dict(CPL)
    if COUPLE:
        s.update(agree=COUPLE["agree"], n=COUPLE["n"], phi=COUPLE["phi"], hic=COUPLE["hic"])
        for k in ("p", "q", "pt_p", "pt_q", "r2", "r2_ends"):
            if s.get(k) is None and COUPLE.get(k) is not None:
                s[k] = COUPLE[k]
    if R.have(SCAN_TSV):
        sc = pd.read_csv(SCAN_TSV, sep="\t")
        s["scan_n"] = len(sc)
        s["scan_median"] = float(sc.x.median())
        s["scan_p90"] = float(np.percentile(sc.x, 90))
    hb = s.get("hic")
    s["junction"] = None if not hb or not ok(JS["hic_lo"]) else hb["x"] >= 0.5 * JS["hic_lo"]
    s["joined"] = any(num_of(j["copy"]) in (num_of(s.get("p")), num_of(s.get("q"))) for j in JOINS) \
        if s.get("p") else False
    for k in ("agree", "n", "pt_p", "lo_p", "hi_p", "pt_q", "lo_q", "hi_q", "c", "r2", "r2_ends"):
        if k in s:
            R.put("chr5x6.%s" % k, s[k])
    return s


CS = R.guarded("chr5 x chr6", cpl_stats) or {}


# ---- deserts and structure (repeat_recombination windows, crossover reference coordinates)
def runs_of(w, col):
    """Merge consecutive windows with w[col] True per chromosome: [(chrom, start_mb, end_mb, rows)]."""
    out = []
    for c in list(dict.fromkeys(w.chrom)):
        x = w[w.chrom == c].sort_values("start")
        cur = None
        for _, r in x.iterrows():
            if r[col]:
                if cur and abs(r.start - cur[2]) < 1:
                    cur[2] = r.end
                    cur[3].append(r)
                else:
                    if cur:
                        out.append(cur)
                    cur = [c, r.start, r.end, [r]]
            elif cur:
                out.append(cur)
                cur = None
        if cur:
            out.append(cur)
    return [(c, a / MB, b / MB, pd.DataFrame(rows)) for c, a, b, rows in out]


DESERT_RULES = dict(inverted=0.25, translocated=0.25, coupling=0.5, none=0.5, markers=0.6)


def desert_stats():
    w = RRD.get("paradoxa", {}).get("windows")
    if w is None:
        raise FileNotFoundError("missing input: %s" % rr_path("windows", P))
    copy_mk, copy_rep = {}, {}
    for c in w.chrom.unique():
        x = w[w.chrom == c]
        m = x.markers_per_mb.notna().values
        copy_mk[c] = float(np.average(x.markers_per_mb.values[m], weights=x.len_mb.values[m])) if m.any() else np.nan
        copy_rep[c] = float(np.average(x["all repeats"].values, weights=x.len_mb.values))
    rows = []
    for c, a, b, x in runs_of(w, "dead"):
        wt = x.len_mb.values.astype(float)

        def avg(col):
            v = x[col].values.astype(float)
            m = np.isfinite(v)
            return float(np.average(v[m], weights=wt[m])) if m.any() else float("nan")
        r = dict(chrom=c, start=a, end=b, mb=b - a, inverted=avg("struct_inverted"),
                 translocated=avg("struct_translocated"), collinear=avg("struct_collinear"),
                 coupling=float(np.average(x.coupling.astype(float), weights=wt)),
                 markers=avg("markers_per_mb") / copy_mk.get(c, np.nan),
                 cut=any(a < u < b for u in CUTS.get(c, [])), repeats=avg("all repeats"),
                 repeats_copy=copy_rep.get(c, np.nan))
        known = np.isfinite(x.struct_collinear.values.astype(float))
        r["none"] = float(np.average((~known).astype(float) + np.where(known, x.struct_none.values.astype(float), 0),
                                     weights=wt))
        lab = []
        if ok(r["inverted"]) and r["inverted"] >= DESERT_RULES["inverted"]:
            lab.append("an inversion between the haplotypes (%s of it)" % pc(r["inverted"]))
        if ok(r["translocated"]) and r["translocated"] >= DESERT_RULES["translocated"]:
            lab.append("the X arm (%s translocated)" % pc(r["translocated"]))
        if r["coupling"] >= DESERT_RULES["coupling"]:
            lab.append("the coupled chr5 × chr6 middle")
        r["kind"] = "structure" if lab else None
        if r["none"] >= DESERT_RULES["none"]:
            lab.append("structure unknown: not aligned to the other copy (%s)" % pc(r["none"]))
            r["kind"] = r["kind"] or "unknown"
        if ok(r["markers"]) and r["markers"] < DESERT_RULES["markers"]:
            lab.append("few markers (%.2f× the copy)" % r["markers"])
            r["kind"] = r["kind"] or "weak"
        if r["cut"]:
            lab.append("a cut of the crossover reference inside")
            r["kind"] = r["kind"] or "weak"
        r["kind"] = r["kind"] or "unexplained"
        r["labels"] = lab or ["nothing in these data"]
        rows.append(r)
    D = pd.DataFrame(rows)
    out = dict(rows=D, n=len(D), mb=float(D.mb.sum()) if len(D) else 0.0)
    for k in ("structure", "weak", "unknown", "unexplained"):
        out["n_" + k] = int((D.kind == k).sum()) if len(D) else 0
    st = D[(D.inverted >= DESERT_RULES["inverted"]) | (D.coupling >= DESERT_RULES["coupling"])] if len(D) else D
    out["struct_markers"] = (float(st.markers.min()), float(st.markers.max())) if len(st) else None
    # crossover density by structure class (window majority)
    an = w[(~w.at_cut) & (w.n_cells >= 10)].copy()

    def cls(r):
        if r.coupling:
            return "coupled chr5 × chr6 middle"
        if ok(r.struct_translocated) and r.struct_translocated >= 0.5:
            return "X arm (translocated)"
        if ok(r.struct_inverted) and r.struct_inverted >= 0.5:
            return "inverted between the haplotypes"
        if not ok(r.struct_collinear) or (ok(r.struct_none) and r.struct_none >= 0.5):
            return "not aligned (P, Q: structure unknown)"
        return "collinear"
    an["cls"] = an.apply(cls, axis=1)
    tot_co, tot_mb = an.co_mid.sum(), an.len_mb.sum()
    rr = []
    for k in ("collinear", "inverted between the haplotypes", "coupled chr5 × chr6 middle", "X arm (translocated)",
              "not aligned (P, Q: structure unknown)"):
        x = an[an.cls == k]
        if not len(x):
            continue
        exp = tot_co * x.len_mb.sum() / tot_mb
        rr.append(dict(cls=k, windows=len(x), mb=x.len_mb.sum(), co=int(x.co_mid.sum()), exp=exp,
                       ratio=x.co_mid.sum() / exp if exp else np.nan,
                       dead=float(np.average(x.dead.astype(float), weights=x.len_mb))))
    out["by_class"] = pd.DataFrame(rr)
    for r in rr:
        R.put("structure.%s.co_ratio" % r["cls"], r["ratio"])
    for k in ("n", "n_structure", "n_weak", "n_unknown", "n_unexplained"):
        R.put("deserts.%s" % k, out[k])
    return out


DS = R.guarded("deserts", desert_stats)


def deserts_sentence():
    if not DS:
        return ""
    bits = ["%d sit on structure (an inversion between the haplotypes, the coupled chr5 × chr6 middle or X)"
            % DS["n_structure"]]
    if DS["n_weak"]:
        bits.append("%d on few markers or a cut" % DS["n_weak"])
    if DS["n_unknown"]:
        bits.append("%d where the structure is unknown (no alignment)" % DS["n_unknown"])
    bits.append("%d on nothing these data show" % DS["n_unexplained"])
    return "%d stretches of ≥ %s Mb without a crossover (chance %s): %s." % (
        DS["n"], ff(DEAD_MB, 0), RRD.get("paradoxa", {}).get("dead_chance", "–"), ", ".join(bits))


# =================================================================== figures
def bar(ax, x0, x1, y, h, col, ec="none", lw=0.0, z=2, hatch=None, alpha=1.0):
    ax.add_patch(Rectangle((x0, y - h / 2), x1 - x0, h, fc=col, ec=ec, lw=lw, zorder=z, hatch=hatch, alpha=alpha))


def chrom(ax, c, y, h=0.42, x0=0.0, rev=False, label=True, fs=8, minw=18, segs=None):
    """A copy of the two-haplotype assembly to scale, coloured by arm; rev draws it right to left."""
    segs = segs or segs_of(c)
    L = max(b for _, b, _ in segs)
    for a, b, nm in segs:
        xa, xb = (x0 + a, x0 + b) if not rev else (x0 + L - b, x0 + L - a)
        bar(ax, xa, xb, y, h, col_of(nm))
        if label and xb - xa >= minw:
            ax.text((xa + xb) / 2, y, nm, ha="center", va="center", fontsize=fs, fontweight="bold", color=INK, zorder=4)
    bar(ax, x0, x0 + L, y, h, "none", ec=INK, lw=0.7, z=3)
    return L


def join_mark(ax, x, y, h, j, fs=7.5, text=True):
    good = j["phi"] >= 0.5
    col = INK if good else BAD
    ax.plot([x, x], [y - h / 2 - 0.05, y + h / 2 + 0.05], color=col, lw=1.6, zorder=5)
    yy = y + h / 2 + 0.14
    ax.plot(x, yy, marker="o" if good else "X", ms=7, color=col, mec="white", mew=0.8, zorder=6)
    if text:
        ax.text(x + 4, yy, "pollen φ %s" % sg(j["phi"]), fontsize=fs, color=col, va="center", ha="left", zorder=6)


def fig_overview():
    """D0: the copies of the two-haplotype assembly, coloured by arm, with the pollen's verdicts."""
    if not COPIES:
        raise FileNotFoundError("no two-haplotype assembly index")
    nums = sorted({num_of(c) for c in COPIES})
    Lmax = max(L2N[c] for c in COPIES) / MB
    fig, ax = plt.subplots(figsize=(13, 1.19 * len(nums) + 1.6))
    ys, y = {}, 0.0
    for n in nums:
        for h in ("hap1", "hap2"):
            c = "chr%d_%s" % (n, h)
            if c in L2N:
                ys[c] = y
                y -= 0.62
        y -= 0.5
    for c, yy in ys.items():
        chrom(ax, c, yy, h=0.36, fs=7.5)
        let = JBY.get(c, {}).get("letter", "")
        ax.text(-8, yy, ("%s  " % let if let else "") + c, ha="right", va="center", fontsize=8.5,
                color=HAP_COL.get(hap_of(c), INK), fontweight="bold" if let else "normal")
        for a, b in merged_inversions(c):
            bar(ax, a, b, yy - 0.25, 0.07, ST_COL["inverted"], z=3)
    for c, j in JBY.items():
        if c in ys:
            join_mark(ax, j["pos"], ys[c], 0.36, j, fs=7)
    if CS.get("p") in ys and CS.get("q") in ys and ok(CS.get("pt_p")):
        x1, y1, x2, y2 = CS["pt_p"], ys[CS["p"]], CS["pt_q"], ys[CS["q"]]
        ax.annotate("", xy=(x2, y2 + 0.2), xytext=(x1, y1 - 0.2),
                    arrowprops=dict(arrowstyle="<->", color=CPL_COL, lw=1.6, ls="--",
                                    connectionstyle="arc3,rad=-0.3"), zorder=7)
        xr = max(L2N.get(CS["p"], 0), L2N.get(CS["q"], 0)) / MB + 12
        ax.text(xr, (y1 + y2) / 2, "%s and %s: inherited together\naround %.0f / %.0f Mb in %s of %s nuclei;\n%s" % (
            CS["p"], CS["q"], x1, x2, CS.get("agree", "–"), CS.get("n", "–"),
            "no join in the assembly" if not CS.get("joined") else "joined in the assembly"), color=CPL_COL,
                fontsize=7.8, va="center")
    ax.set_xlim(-95, Lmax + 120)
    ax.set_ylim(y + 0.3, 0.75)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("Mb along each copy")
    LEG = {"L1": "arm L1", "L2": "arm L2", "P": "arm P", "Q": "arm Q", "X": "X: on chr4 in hap1, on chr3 in hap2",
           "core": "chr4 without X (core)"}
    present = [a for a in list(ARM_COL) + list(CHR_COL) if any(nm == a for c in COPIES for _, _, nm in segs_of(c))]
    hd = [Patch(fc=col_of(a), ec=INK, lw=0.5) for a in present] + \
         [Patch(fc=ST_COL["inverted"]), Line2D([], [], marker="o", color=INK, ls="", ms=7),
          Line2D([], [], marker="X", color=BAD, ls="", ms=7), Line2D([], [], color=CPL_COL, ls="--", lw=1.6)]
    lb = [LEG.get(a, a) for a in present] + ["inverted between the haplotypes", "join the pollen confirm",
                                              "join the pollen reject", "coupled in the pollen"]
    ax.legend(hd, lb, loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=5, fontsize=7.5, frameon=False)
    ax.set_title("The two-haplotype assembly of D. paradoxa (hap1 blue, hap2 orange), with what the pollen say about "
                 "each join", loc="left", fontsize=10)
    return fig


def fig_arms():
    """D1: chr1/chr2 arms. a: the assembly's four copies; b: what the pollen inherit; c, d: meiosis I pictures."""
    need(KT_TSV)
    if not ABCD:
        raise ValueError("the karyotype test does not name the four copies A-D")
    fig = plt.figure(figsize=(13.5, 8.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 1], hspace=0.32, wspace=0.05)
    axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    axc, axd = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    order = ["A", "C", "B", "D"]
    Lmax = max(L2N[LETTERS[k]["copy"]] for k in order) / MB
    for i, k in enumerate(order):
        j = LETTERS[k]
        c, y = j["copy"], -i * 1.0
        chrom(axa, c, y, h=0.36)
        join_mark(axa, j["pos"], y, 0.36, j, text=False)
        axa.text(-10, y, "%s  %s" % (k, c), ha="right", va="center", fontsize=9, fontweight="bold",
                 color=HAP_COL[hap_of(c)])
        axa.text(0, y - 0.33, "%s|%s join at %.1f Mb: pollen φ %s (95%% %s to %s) → %s" % (
            j["left"], j["right"], j["pos"], sg(j["phi"]), sg(j["lo"]), sg(j["hi"]),
            "one piece" if j["phi"] >= 0.5 else "NOT one piece"), fontsize=7.4, va="center",
                 color=INK if j["phi"] >= 0.5 else BAD)
        axa.text(0, y - 0.5, "%s of %s nuclei agree; tissue Hi-C across the join %s" % (
            j["agree"], j["n"], ("%.1f×" % j["hic"]) if ok(j["hic"]) else "–"), fontsize=7.1, va="center",
                 color=MUTED)
    axa.set_xlim(-95, Lmax + 20)
    axa.set_ylim(-3.65, 0.45)
    axa.axis("off")
    axa.set_title("a  The assembly (tissue): four copies, A–D, each two arms joined", loc="left", fontsize=10)
    rows = []
    for k in order:
        j = LETTERS[k]
        if j["phi"] >= 0.5:
            rows.append(dict(name="%s (as assembled)" % k, segs=segs_of(j["copy"]), phi=j["phi"], src=[k, k]))
    for pr in JS["pairs"]:
        ends = [(pr["arm"], pr["pos"] if pr["side"] == "up to" else L2N[pr["copy"]] / MB - pr["pos"], pr["own"],
                 pr["side"]),
                (pr["o_arm"], pr["o_pos"] if pr["o_side"] == "up to" else L2N[pr["o_copy"]] / MB - pr["o_pos"],
                 pr["o_own"], pr["o_side"])]
        ends.sort(key=lambda e: 0 if e[3] == "up to" else 1)      # the arm that starts its copy goes first
        first, second = ends
        rows.append(dict(name="%s's %s + %s's %s" % (first[2], first[0], second[2], second[0]),
                         segs=[(0.0, first[1], first[0]), (first[1], first[1] + second[1], second[0])], phi=pr["phi"],
                         src=[first[2], second[2]]))
    rows.sort(key=lambda r: (r["segs"][0][2], r["name"]))
    for i, r in enumerate(rows):
        y = -i * 1.0
        chrom(axb, None, y, h=0.36, segs=r["segs"])
        for (a, b, nm), s in zip(r["segs"], r["src"]):
            if b - a > 40:
                axb.text((a + b) / 2, y - 0.33, "from %s" % s, ha="center", va="center", fontsize=7, color=MUTED)
        axb.text(-10, y, r["name"], ha="right", va="center", fontsize=8.5, fontweight="bold")
        axb.text(r["segs"][-1][1] + 10, y, "φ %s" % sg(r["phi"]), fontsize=8, va="center", color=INK)
    axb.set_xlim(-140, Lmax + 60)
    axb.set_ylim(-3.65, 0.45)
    axb.axis("off")
    combos = sorted({"·".join(s[2] for s in r["segs"]) for r in rows})
    axb.set_title("b  What the pollen inherit together: %s" % ", ".join(combos), loc="left", fontsize=10)
    arm_len = {}
    for k in order:
        for a, b, nm in segs_of(LETTERS[k]["copy"]):
            arm_len.setdefault(nm, []).append(b - a)
    arm_len = {k: float(np.mean(v)) for k, v in arm_len.items()}
    dirs = {"L1": (-1, 0), "P": (0, 1), "L2": (1, 0), "Q": (0, -1)}
    sc, d = 1 / 300.0, 0.045
    comp = {k: (LETTERS[k]["left"], LETTERS[k]["right"]) for k in order}
    ok_cross = all(set(v) <= set(dirs) and dirs[v[0]][0] * dirs[v[1]][0] + dirs[v[0]][1] * dirs[v[1]][1] == 0
                   for v in comp.values())
    if ok_cross:
        done = set()
        for k, (a1, a2) in comp.items():
            cx = (dirs[a1][0] + dirs[a2][0]) * d
            cy = (dirs[a1][1] + dirs[a2][1]) * d
            for arm in (a1, a2):
                dx, dy = dirs[arm]
                ln = arm_len[arm] * sc
                axc.plot([cx, cx + dx * ln], [cy, cy + dy * ln], color=ARM_COL[arm], lw=6.5, solid_capstyle="butt",
                         zorder=2)
                if arm not in done:
                    done.add(arm)
                    axc.text(dx * (ln + 0.1), dy * (ln + 0.1), arm, ha="center", va="center", fontsize=8.5,
                             fontweight="bold", color=INK)
            axc.text(cx * 4.2, cy * 4.2, k, ha="center", va="center", fontsize=10, fontweight="bold",
                     color=HAP_COL[hap_of(LETTERS[k]["copy"])])
        axc.set_xlim(-1.15, 1.15)
        axc.set_ylim(-0.75, 0.75)
        axc.set_aspect("equal")
    axc.axis("off")
    ring, arm = [order[0]], comp[order[0]][1]          # walk the ring: each copy shares one arm with the next
    while len(ring) < len(order):
        nxt = next((k for k in order if k not in ring and arm in comp[k]), None)
        if nxt is None:
            break
        ring.append(nxt)
        arm = [a for a in comp[nxt] if a != arm][0]
    axc.set_title("c  Meiosis I if the assembly is right: one cross of four\n(%s, each pairing one arm with the "
                  "next)" % "–".join(ring + [ring[0]]), loc="left", fontsize=9.5)
    for i, r in enumerate(rows):
        y = 0.9 - i * 0.55 - (0.25 if i >= 2 else 0)
        x = 0.0
        for a, b, nm in r["segs"]:
            axd.add_patch(Rectangle((x, y - 0.11), (b - a) * sc, 0.22, fc=col_of(nm), ec="none"))
            axd.text(x + (b - a) * sc / 2, y, nm, ha="center", va="center", fontsize=8, fontweight="bold")
            x += (b - a) * sc
        axd.add_patch(Rectangle((0, y - 0.11), x, 0.22, fc="none", ec=INK, lw=0.6))
        axd.text(x + 0.04, y, r["name"], va="center", fontsize=7.6, color=MUTED)
    axd.set_xlim(-0.05, 2.4)
    axd.set_ylim(-1.35, 1.2)
    axd.axis("off")
    axd.set_title("d  Meiosis I if the pollen are right: ordinary pairs (bivalents)", loc="left", fontsize=9.5)
    nbad = [j["letter"] or j["copy"] for j in JS["bad"] if j["copy"] in [LETTERS[k]["copy"] for k in order]]
    fig.suptitle("chr1 / chr2: %s" % ("the pollen do not inherit the assembly's joins %s" % " and ".join(nbad)
                                      if nbad else "the pollen inherit all four joins of the assembly"),
                 x=0.01, y=0.995, ha="left", fontsize=11)
    return fig


def fig_x():
    """D2: chr3 / chr4 / X."""
    need(XARM_TSV)
    s = XS
    fig = plt.figure(figsize=(13.5, 9.8))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1.1], hspace=0.3, wspace=0.12)
    axa, axb = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    axc, axd = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    cs = ["chr3_hap1", "chr3_hap2", "chr4_hap1", "chr4_hap2"]
    Lmax = max(L2N[c] for c in cs) / MB
    tl = SYN[SYN.kind == "translocated"] if SYN is not None else None
    for i, c in enumerate(cs):
        y = -i * 1.0
        chrom(axa, c, y, h=0.42)
        axa.text(-10, y, c, ha="right", va="center", fontsize=9, color=HAP_COL[hap_of(c)], fontweight="bold")
        j = JBY.get(c)
        if j:
            join_mark(axa, j["pos"], y, 0.42, j)
        for a, b, nm in segs_of(c):     # X's orientation: from its free (0 Mb) end
            if nm == "X":
                rev = False
                if tl is not None and len(tl):
                    x_ = tl[(tl.hap2 == c) | (tl.hap1 == c)]
                    rev = bool(len(x_)) and x_.iloc[0].strand == "-" and hap_of(c) == "hap2"
                xa, xb = (b, a) if rev else (a, b)
                axa.annotate("", xy=(xb - (xb - xa) * 0.12, y - 0.33), xytext=(xa + (xb - xa) * 0.12, y - 0.33),
                             arrowprops=dict(arrowstyle="->", color=INK, lw=0.9))
                axa.text(xa, y - 0.33, " 0 " if rev else "0 ", ha="right" if not rev else "left", va="center",
                         fontsize=6.5, color=MUTED)
    axa.set_xlim(-80, Lmax + 20)
    axa.set_ylim(-3.6, 0.75)
    axa.axis("off")
    axa.set_title("a  The assembly: X is the start of chr4 in hap1 and the end of chr3 in hap2\n"
                  "(arrow: X from its free 0 Mb end to the end that is joined)", loc="left", fontsize=9.5)
    L = {c: L2N[c] / MB for c in cs}
    j3, j4 = JBY.get("chr3_hap2"), JBY.get("chr4_hap1")
    x3 = L["chr3_hap2"] - j3["pos"] if j3 else 100.0      # hap2's X length
    x4 = j4["pos"] if j4 else 85.0                          # hap1's X length
    rows = [("chr3_hap1", x3 + j3["pos"] - L["chr3_hap1"] if j3 else 0, True),
            ("chr3_hap2", 0.0, True), ("chr4_hap1", 0.0, False), ("chr4_hap2", x4, False)]
    for i, (c, x0, rev) in enumerate(rows):
        y = -i * 0.8
        chrom(axb, c, y, h=0.34, x0=x0, rev=rev, fs=7.5)
        axb.text(x0 - 8, y, c, ha="right", va="center", fontsize=8.5, color=HAP_COL[hap_of(c)], fontweight="bold")
    spans = [(0, 1, x3, min(x3 + j3["pos"], x3 + L["chr3_hap1"]) if j3 else 100, "chr3 pairs with chr3"),
             (1, 2, 0, min(x3, x4), "X pairs with X"),
             (2, 3, x4, min(L["chr4_hap1"], x4 + L["chr4_hap2"]), "chr4 core pairs with core")]
    for r1, r2, a, b, t in spans:
        y1, y2 = -r1 * 0.8 - 0.17, -r2 * 0.8 + 0.17
        for xx in np.linspace(a + 3, b - 3, max(4, int((b - a) / 9))):
            axb.plot([xx, xx], [y1, y2], color=LINE, lw=0.7, zorder=1)
        axb.text((a + b) / 2 if r1 != 1 else b + 6, (y1 + y2) / 2, t, fontsize=7.5, ha="center" if r1 != 1 else "left",
                 va="center", color=INK, bbox=dict(fc="white", ec="none", pad=0.6), zorder=3)
    axb.set_xlim(-95, max(L["chr4_hap1"], x4 + L["chr4_hap2"]) + 10)
    axb.set_ylim(-2.85, 0.6)
    axb.axis("off")
    axb.set_title("b  At meiosis I the four copies can pair as one chain\n(chr3_hap1 – chr3_hap2 – chr4_hap1 – "
                  "chr4_hap2; chr3 drawn right to left)", loc="left", fontsize=9.5)
    cells = [("RR", "chr3_hap1", "chr4_hap1", "X once (hap1's)"), ("RA", "chr3_hap2", "chr4_hap1", "X twice"),
             ("AR", "chr3_hap1", "chr4_hap2", "X never"), ("AA", "chr3_hap2", "chr4_hap2", "X once (hap2's)")]
    scl = 1 / 720.0
    for k, (key, c3, c4, lab) in enumerate(cells):
        cx, cy = (k % 2) * 1.05, -(k // 2) * 1.0
        unbalanced = key in ("RA", "AR")
        axc.add_patch(Rectangle((cx, cy - 0.92), 0.98, 0.92, fc="#FBEFEF" if unbalanced else "#FBFAF7",
                                ec=BAD if unbalanced else LINE, lw=1.0))
        axc.text(cx + 0.04, cy - 0.12, lab, fontsize=9.5, fontweight="bold", color=BAD if unbalanced else INK,
                 va="center")
        for m_, cc in enumerate((c3, c4)):
            yy = cy - 0.3 - m_ * 0.15
            axc.text(cx + 0.04, yy, cc, fontsize=7, va="center", color=HAP_COL[hap_of(cc)])
            x0 = cx + 0.27
            for a, b, nm in segs_of(cc):
                axc.add_patch(Rectangle((x0 + a * scl, yy - 0.045), (b - a) * scl, 0.09, fc=col_of(nm), ec="none"))
            axc.add_patch(Rectangle((x0, yy - 0.045), L2N[cc] / MB * scl, 0.09, fc="none", ec=INK, lw=0.5))
        axc.text(cx + 0.04, cy - 0.62, "%d nuclei" % s[key], fontsize=12, fontweight="bold", va="center",
                 color=BAD if unbalanced else INK)
        axc.text(cx + 0.04, cy - 0.78, "expected if chr3, chr4 independent: %.1f" % s["exp"][key],
                 fontsize=7.2, va="center", color=MUTED)
    axc.set_xlim(-0.02, 2.06)
    axc.set_ylim(-1.96, 0.03)
    axc.axis("off")
    axc.set_title("c  The four spores a nucleus can come from (%d nuclei with chr3 and chr4 both read)" % s["N"],
                  loc="left", fontsize=9.5)
    sh, (lo, hi) = s["share"], s["share_ci"]
    ms = [("observed: %d of %d (95%% %s–%s)" % (s["RA"], s["N"], pc(lo), pc(hi)), sh, (lo, hi), None),
          ("chr3 and chr4 independent, X-less spores die", 1 / 3.0, None, s["p_indep"]),
          ("held together only when X has a crossover (%s of meioses), X-less spores die" % pc(s["xco_meioses"]),
           s["chain_only_share"], None, s["p_chain_only"])]
    for i, (nm, v, ci, p) in enumerate(ms):
        y = -i * 1.25
        axd.barh(y, v, height=0.5, color=SP_COL["paradoxa"] if i == 0 else "#D3D1C7", ec=INK, lw=0.6)
        if ci:
            axd.errorbar(v, y, xerr=[[v - ci[0]], [ci[1] - v]], fmt="none", ecolor=INK, capsize=4, lw=1.2)
        axd.text(0.003, y + 0.43, nm, fontsize=8.2, va="center")
        axd.text(max(v, ci[1] if ci else v) + 0.012, y, pc(v) + ("   against observed: %s" % pe(p) if p is not None
                                                                else ""), va="center", fontsize=8.2)
    axd.set_yticks([])
    axd.spines["left"].set_visible(False)
    axd.set_xlim(0, 0.62)
    axd.set_ylim(-2.85, 0.85)
    axd.set_xlabel("share of the surviving nuclei that carry X twice")
    axd.set_title("d  How often X twice? Observed against two models\n(X never: %d of %d; twice against never, if both "
                  "survived alike: %s)" % (s["AR"], s["N"], pe(s["p_twice_never"])), loc="left", fontsize=9.5)
    fig.suptitle("chr3 / chr4 / X: X sits on chr4 in hap1 and on chr3 in hap2%s" % (
        "; pollen without X are missing" if s["lost"] else ""), x=0.01, y=0.995, ha="left", fontsize=11)
    return fig


def cos_matrix(a, b):
    c = pd.read_csv(COS_TSV, sep="\t")
    x = c[(c.copy_a == a) & (c.copy_b == b)]
    if not len(x):
        x = c[(c.copy_a == b) & (c.copy_b == a)].rename(columns={"Mb_a": "Mb_b", "Mb_b": "Mb_a"})
    if not len(x):
        return None
    allm = np.sort(c[c.copy_a == c.copy_b].Mb_a.unique())
    step = float(np.median(np.diff(allm))) if len(allm) > 1 else 10.0
    na, nb = int(math.ceil(L2N[a] / MB / step)), int(math.ceil(L2N[b] / MB / step))
    M = np.full((na, nb), np.nan)
    for r in x.itertuples():
        i, j = int(round(r.Mb_a / step)), int(round(r.Mb_b / step))
        if 0 <= i < na and 0 <= j < nb:
            M[i, j] = r.signed_phi
    return M, step


def fig_coupling():
    """D3: chr5 x chr6."""
    need(COS_TSV)
    p, q = CS.get("p", "chr5_hap1"), CS.get("q", "chr6_hap1")
    fig = plt.figure(figsize=(13.5, 9.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.3, 1], height_ratios=[1, 1.1], hspace=0.3, wspace=0.16)
    axa, axb, axc = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[:, 1]), fig.add_subplot(gs[1, 0])
    cs = [c for c in ("%s_hap1" % p.split("_")[0], "%s_hap2" % p.split("_")[0], "%s_hap1" % q.split("_")[0],
                      "%s_hap2" % q.split("_")[0]) if c in L2N]
    w = RRD.get("paradoxa", {}).get("windows")
    Lmax = max(L2N[c] for c in cs) / MB
    for i, c in enumerate(cs):
        y = -i * 1.0 - (0.3 if i >= 2 else 0)
        chrom(axa, c, y, h=0.38)
        axa.text(-8, y, c, ha="right", va="center", fontsize=9, color=HAP_COL[hap_of(c)], fontweight="bold")
        for a, b in merged_inversions(c):
            bar(axa, a, b, y - 0.26, 0.08, ST_COL["inverted"], z=3)
        if w is not None and c in set(w.chrom):
            for _, a, b, _ in runs_of(w[w.chrom == c], "dead"):
                bar(axa, a, b, y + 0.3, 0.07, INK, z=3)
        for k in ("p", "q"):
            if CS.get(k) == c and ok(CS.get("lo_%s" % k)):
                lo_, hi_, pt = CS["lo_%s" % k], CS["hi_%s" % k], CS["pt_%s" % k]
                axa.plot([lo_, hi_], [y + 0.45, y + 0.45], color=CPL_COL, lw=2.2, solid_capstyle="butt")
                axa.plot(pt, y + 0.45, marker="D", color=CPL_COL, ms=6)
                axa.text(lo_, y + 0.62, "coupling point %.0f Mb (95%% %.0f–%.0f)" % (pt, lo_, hi_), fontsize=7.5,
                         color=CPL_COL, va="center", ha="left")
    axa.set_xlim(-75, Lmax + 15)
    axa.set_ylim(-3.75, 0.8)
    axa.axis("off")
    axa.set_title("a  The assembly: %s and %s are separate chromosomes in both haplotypes\n(purple: coupling point, "
                  "95%% interval; black: desert ≥ %s Mb; pink: inverted)" % (
                      p.split("_")[0], q.split("_")[0], ff(DEAD_MB, 0)), loc="left", fontsize=9.3)
    mm = cos_matrix(p, q)
    if mm is not None:
        M, step = mm
        im = axb.imshow(np.ma.masked_invalid(M), cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest",
                        extent=(0, M.shape[1] * step, M.shape[0] * step, 0), aspect="auto")
        if ok(CS.get("pt_p")):
            axb.plot(CS["pt_q"], CS["pt_p"], "x", color=INK, ms=11, mew=2.2)
            if ok(CS.get("lo_p")):
                axb.plot([CS["lo_q"], CS["hi_q"]], [CS["pt_p"], CS["pt_p"]], color=INK, lw=1)
                axb.plot([CS["pt_q"], CS["pt_q"]], [CS["lo_p"], CS["hi_p"]], color=INK, lw=1)
        axb.set_xlabel("%s, Mb" % q)
        axb.set_ylabel("%s, Mb" % p)
        cb = fig.colorbar(im, ax=axb, fraction=0.05, pad=0.03, shrink=0.55)
        cb.set_label("φ across nuclei: red = inherited together, white = independent", fontsize=8)
    axb.set_title("b  The pollen: every window of %s against\nevery window of %s (x: coupling point; lines: 95%%)" % (
        p, q), loc="left", fontsize=9.3)
    axc.axis("off")
    axc.set_xlim(0, 2.1)
    axc.set_ylim(-1.55, 0.25)
    sc_ = 1 / 520.0
    L5, L6 = L2N.get(p, 290e6) / MB, L2N.get(q, 360e6) / MB
    b5, b6 = CS.get("pt_p", L5 / 2), CS.get("pt_q", L6 / 2)
    n5, n6 = p.split("_")[0], q.split("_")[0]

    def mini(x0, y, parts, lab):
        x = x0
        for ln, colr in parts:
            axc.add_patch(Rectangle((x, y - 0.045), ln * sc_, 0.09, fc=colr, ec="none"))
            x += ln * sc_
        axc.add_patch(Rectangle((x0, y - 0.045), x - x0, 0.09, fc="none", ec=INK, lw=0.5))
        axc.text(x + 0.02, y, lab, fontsize=7, va="center", color=MUTED)
    c5, c6 = col_of(n5), col_of(n6)
    axc.text(0.0, 0.17, "1  A hidden reciprocal swap between %s and %s (in one haplotype)" % (n5, n6), fontsize=8.8,
             fontweight="bold")
    mini(0.02, 0.02, [(L5, c5)], "%s (normal)" % n5)
    mini(0.02, -0.1, [(L6, c6)], "%s (normal)" % n6)
    mini(1.05, 0.02, [(b5, c5), (L6 - b6, c6)], "%s·%s (swapped)" % (n5[3:], n6[3:]))
    mini(1.05, -0.1, [(b6, c6), (L5 - b5, c5)], "%s·%s (swapped)" % (n6[3:], n5[3:]))
    axc.text(0.02, -0.25, "a spore with one swapped and one normal chromosome lacks a segment and dies, so the middles "
             "go together;\npredicts a %s–%s junction in HiFi reads and Hi-C, and a ring or chain of four at meiosis I"
             % (n5, n6), fontsize=7.6, va="top", color=INK)
    axc.text(0.0, -0.72, "2  No physical link: mixed spores are lost (selection) or never made (directed segregation)",
             fontsize=8.8, fontweight="bold")
    for k, (lab, keep) in enumerate((("%s hap1 + %s hap1" % (n5, n6), True), ("%s hap2 + %s hap2" % (n5, n6), True),
                                     ("%s hap1 + %s hap2" % (n5, n6), False), ("%s hap2 + %s hap1" % (n5, n6), False))):
        x0, y = 0.02 + (k % 2) * 1.03, -0.87 - (k // 2) * 0.14
        axc.text(x0, y, ("kept:  " if keep else "lost:  ") + lab, fontsize=7.8, color=INK if keep else BAD,
                 va="center")
    axc.text(0.02, -1.2, "predicts no junction anywhere, two ordinary pairs at meiosis I, and either ~half the pollen "
             "aborted\n(selection) or normal pollen (directed segregation)", fontsize=7.6, va="top", color=INK)
    axc.set_title("c  Two ways to get it (the pollen cannot tell them apart)", loc="left", fontsize=9.5)
    fig.suptitle("%s × %s: the middles of the two chromosomes are inherited as one unit%s" % (
        n5, n6, ", with no join in the assembly" if not CS.get("joined") else ""), x=0.01, y=0.995, ha="left",
                 fontsize=11)
    return fig


def fig_structure():
    """D4: per chromosome, the crossover reference copy (crossovers above it), the other copy, and the synteny
    blocks between them: collinear, inverted, translocated."""
    w = RRD.get("paradoxa", {}).get("windows")
    if SYN is None or w is None or not REFC:
        raise FileNotFoundError("missing synteny blocks, windows or pieces")
    rows = []
    for c in REFC:
        other = "chr%d_%s" % (num_of(c), "hap2" if hap_of(c) == "hap1" else "hap1")
        if other in L2N and c in L2N:
            rows.append((c, other))
    cm_max = float(np.nanpercentile(w.cM_per_Mb, 99)) or 1.0
    geo = {}
    for ref, oth in rows:
        t, q = (ref, oth) if hap_of(ref) == "hap1" else (oth, ref)
        flip = {ref: False, oth: main_strand(t, q) == "-"}
        x0 = {ref: 0.0, oth: 0.0}

        def drawn(c, v):
            return (L2N[c] - v) / MB if flip[c] else v / MB
        blk = SYN[(SYN.hap1 == t) & (SYN.hap2 == q) & (SYN.kind == "collinear")]
        if len(blk):
            tm = [drawn(t, (r.t_start + r.t_end) / 2) for r in blk.itertuples()]
            qm = [drawn(q, (r.q_start + r.q_end) / 2) for r in blk.itertuples()]
            dlt = [(a - b) if t == ref else (b - a) for a, b in zip(tm, qm)]      # ref minus other
            x0[oth] = wmedian(dlt, blk.aligned.values)
        geo[ref] = (t, q, flip, x0)
    lo_x = min([0.0] + [geo[r][3][o] for r, o in rows])
    hi_x = max(max(L2N[r] / MB, L2N[o] / MB + geo[r][3][o]) for r, o in rows)
    fig, axs = plt.subplots(len(rows), 1, figsize=(13.5, 2.1 * len(rows) + 1.0), sharex=True)
    axs = np.atleast_1d(axs)
    for ax, (ref, oth) in zip(axs, rows):
        t, q, flip, x0 = geo[ref]
        yr, yo, hb, base = 1.45, 0.3, 0.3, 1.92
        x = w[w.chrom == ref].sort_values("start")
        if len(x):
            xs = np.r_[x.start.values, x.end.values[-1]] / MB
            ys = base + 0.95 * np.clip(x.cM_per_Mb.values / cm_max, 0, 1.2)
            ax.fill_between(xs, base, np.r_[ys, ys[-1]], step="post", color=SP_COL["paradoxa"], alpha=0.55, lw=0)
        ax.plot([0, L2N[ref] / MB], [base, base], color=LINE, lw=0.6)
        for _, a, b, _ in runs_of(x, "dead"):
            bar(ax, a, b, 1.75, 0.07, INK, z=3)
        partner = CS.get("q") if ref == CS.get("p") else (CS.get("p") if ref == CS.get("q") else None)
        for _, a, b, _ in runs_of(x, "coupling"):
            ax.plot([a, b], [3.12, 3.12], color=CPL_COL, lw=2.4, solid_capstyle="butt")
            ax.text(b + 3, 3.12, "coupled to %s in the pollen" % partner if partner else "coupled", color=CPL_COL,
                    fontsize=7.4, va="center")
        for u in CUTS.get(ref, []):
            ax.plot([u, u], [yr - hb / 2, 3.0], color=INK, lw=0.8, ls=":")
        lab_x = lo_x - 6
        for c, y in ((ref, yr), (oth, yo)):
            chrom(ax, c, y, h=hb, x0=x0[c], rev=flip[c], fs=7, minw=22)
            ax.text(lab_x, y, c + (" (drawn reversed)" if flip[c] else ""), ha="right", va="center", fontsize=8.3,
                    color=HAP_COL[hap_of(c)], fontweight="bold" if c == ref else "normal")

        def px(c, v):
            return x0[c] + ((L2N[c] - v) / MB if flip[c] else v / MB)
        blk = SYN[(SYN.hap1 == t) & (SYN.hap2 == q) & SYN.kind.isin(["collinear", "inverted"])]
        y_t, y_q = (yr - hb / 2, yo + hb / 2) if t == ref else (yo + hb / 2, yr - hb / 2)
        for r in blk.sort_values("aligned").itertuples():
            q_at_t0, q_at_t1 = (r.q_start, r.q_end) if r.strand == "+" else (r.q_end, r.q_start)
            poly = [(px(t, r.t_start), y_t), (px(t, r.t_end), y_t), (px(q, q_at_t1), y_q), (px(q, q_at_t0), y_q)]
            ax.add_patch(Polygon(poly, closed=True, fc=ST_COL[r.kind], ec="none",
                                 alpha=0.85 if r.kind == "inverted" else 0.5, zorder=1))
        nlab = {ref: 0, oth: 0}
        for r in SYN[SYN.kind == "translocated"].sort_values("aligned", ascending=False).itertuples():
            for c, a, b, other_c in ((r.hap1, r.t_start, r.t_end, r.hap2), (r.hap2, r.q_start, r.q_end, r.hap1)):
                if c in (ref, oth) and (b - a) >= 2e6:
                    y = yr if c == ref else yo
                    xa, xb = sorted((px(c, a), px(c, b)))
                    bar(ax, xa, xb, y, hb + 0.08, "none", ec=ST_COL["translocated"], lw=2.2, z=5)
                    if (b - a) >= 10e6 and nlab[c] < 2:
                        ax.text(lab_x, y - 0.25 - 0.17 * nlab[c], "%.0f–%.0f Mb: homolog on %s" % (
                            a / MB, b / MB, other_c), ha="right", va="center", fontsize=6.8,
                                color=ST_COL["translocated"])
                        nlab[c] += 1
        ax.set_ylim(-0.05, 3.25)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlim(lo_x - 95, hi_x + 50)
        ax.set_xticks(np.arange(0, hi_x + 1, 100))
    axs[-1].set_xlabel("Mb along each copy (the second copy shifted to line up with the first)")
    hd = [Patch(fc=SP_COL["paradoxa"], alpha=0.55), Patch(fc=INK), Patch(fc=ST_COL["collinear"], alpha=0.5),
          Patch(fc=ST_COL["inverted"], alpha=0.85), Patch(fc="none", ec=ST_COL["translocated"], lw=2),
          Line2D([], [], color=CPL_COL, lw=2.4), Line2D([], [], color=INK, ls=":", lw=0.8)]
    lb = ["crossovers, cM/Mb (top of scale: %.2f)" % cm_max, "no crossover over ≥ %s Mb (a desert)" % ff(DEAD_MB, 0),
          "same order in both haplotypes", "inverted between the haplotypes", "translocated (X)",
          "coupled chr5 × chr6 middle", "cut of the crossover reference"]
    fig.legend(hd, lb, loc="lower center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle("Crossovers (purple) on the copy they were called on, above the two copies of each chromosome and how "
                 "they align", x=0.01, y=0.998, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0.04, 1, 0.985))
    return fig


def fig_repeats_summary():
    """R1: repeats against crossovers in both species, and what the repeats add beyond chance."""
    have = [sp for sp, _ in SPS if "windows" in RRD.get(sp, {})]
    if not have:
        raise FileNotFoundError("missing repeat_recombination windows")
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for ax, sp in zip(axs[:2], ("binata", "paradoxa")):
        d = RRD.get(sp, {})
        if "windows" not in d:
            ax.axis("off")
            continue
        w = d["windows"]
        w = w[(~w.at_cut) & (w.n_cells >= 10)]
        ax.scatter(100 * w["all repeats"], w.cM_per_Mb, s=6, color=SP_COL[sp], alpha=0.5, lw=0)
        ar = assoc_row(sp, "all repeats")
        ttl = "%s, %g Mb windows" % (SP_NAME[sp], d.get("win_mb", float("nan")))
        if ar is not None:
            ttl += "\nall repeats vs crossovers: ρ %s (%s); after position + markers %s" % (
                sg(ar.rho), pf(ar.p_rot), sg(ar.rho_partial))
        ax.set_title(ttl, loc="left", fontsize=8.8, color=SP_COL[sp])
        ax.set_xlabel("all repeats, % of the window (CARP)")
        ax.set_ylabel("crossovers, cM/Mb")
    ax = axs[2]
    labs = ["repeats alone", "repeats after position + markers", "repeats after everything else"]
    for k, sp in enumerate(("binata", "paradoxa")):
        for i, lab in enumerate(labs):
            r = test_row(sp, "deviance explained: rotation null", lab)
            if r is None:
                continue
            y = -(i * 2.6 + k)
            ax.barh(y, 100 * r.obs, height=0.8, color=SP_COL[sp], alpha=0.85)
            ax.plot([100 * r.null_hi] * 2, [y - 0.45, y + 0.45], color=INK, lw=1.6)
            ax.text(max(100 * r.obs, 100 * r.null_hi) + 0.8, y, "%.1f%% (%s)" % (100 * r.obs, pf(r.p)), va="center",
                    fontsize=7.5)
    ax.set_yticks([-(i * 2.6 + 0.5) for i in range(len(labs))])
    ax.set_yticklabels(labs, fontsize=8)
    ax.set_xlabel("crossover ups and downs explained (Poisson deviance, %)")
    ax.set_title("What the repeats explain; black tick: the 95th percentile\nof the same repeat tracks shifted at "
                 "random along their chromosome", loc="left", fontsize=8.8)
    ax.legend([Patch(color=SP_COL["binata"]), Patch(color=SP_COL["paradoxa"])],
              [SP_NAME["binata"], SP_NAME["paradoxa"]], fontsize=7.5, frameon=False, loc="lower right")
    fig.tight_layout()
    return fig


for key, fn in (("D0", fig_overview), ("D1", fig_arms), ("D2", fig_x), ("D3", fig_coupling), ("D4", fig_structure),
                ("R1", fig_repeats_summary)):
    f_ = R.guarded("figure %s" % key, fn)
    if f_ is not None:
        FIGS[key] = f_


# =================================================================== the report
def nom_arms():
    return ("A = chr1_hap1 (L1 + P), B = chr2_hap1 (L2 + P), C = chr1_hap2 (L1 + Q), D = chr2_hap2 (L2 + Q); "
            "L1, L2, P, Q are the four arms of chr1/chr2; hap1 and hap2 are the two haplotypes of the assembly")


def nom_x():
    xl = JBY["chr4_hap1"]["pos"] if "chr4_hap1" in JBY else float("nan")
    return ("X = the %.1f Mb arm that starts chr4 in hap1 (chr4_hap1 = X + core) and ends chr3 in hap2 (chr3_hap2 = "
            "chr3 + X); core = chr4 without X" % xl)


def tissue_vs_pollen():
    """What the tissue Hi-C says about A and D: (verdict, sentence with the numbers, short sentence)."""
    if not ABCD:
        return None, "", ""
    sup = [hic_supports(LETTERS[k]["hic"]) for k in "AD"]
    hx = "A %s and D %s" % tuple(("%.1f×" % LETTERS[k]["hic"]) if ok(LETTERS[k]["hic"]) else "–" for k in "AD")
    ref = "the confirmed joins read %.0f–%.0f×" % (JS["hic_lo"], JS["hic_hi"]) if ok(JS["hic_lo"]) else ""
    if all(s is True for s in sup):
        return True, "The tissue Hi-C reads A and D as real joins (%s; %s): tissue and pollen disagree." % (hx, ref), \
            "Tissue and pollen disagree."
    if all(s is False for s in sup):
        return False, "The tissue Hi-C agrees with the pollen (%s; %s): no contacts across A or D beyond unrelated " \
                      "chromosomes." % (hx, ref), "Tissue and pollen agree."
    return None, "The tissue Hi-C is mixed or unclear at A and D (%s; %s)." % (hx, ref), \
        "The tissue Hi-C is unclear here."


# ---- 0 summary
def section_summary():
    h2("s0", "0. Summary")
    sP, sB = "paradoxa", "binata"
    nP, nB = v5n(sP, "cells_selected"), v5n(sB, "cells_selected")
    coP, coB = v5n(sP, "co_total"), v5n(sB, "co_total")
    R.p("Pollen nuclei were genotyped along the assemblies: %s <i>D. paradoxa</i> nuclei (%s crossovers) and %s "
        "<i>D. binata</i> nuclei (%s crossovers). Each nucleus is one meiotic product, so whatever it inherits together "
        "tells us which pieces sit on one chromosome in the plant that made the pollen." % (
            fi(nP), fi(coP), fi(nB), fi(coB)))
    tv, tv_txt, _ = tissue_vs_pollen()
    if ABCD:
        p1, p2 = pair_of("L1", "Q"), pair_of("L2", "P")
        rej = [k for k in "ABCD" if LETTERS[k]["phi"] < 0.5]
        acc = [k for k in "ABCD" if LETTERS[k]["phi"] >= 0.5]
        R.p("<b>chr1/chr2 arms</b> (%s): the pollen reject the assembly's join%s %s and confirm %s%s.%s %s" % (
            nom_arms(), "s" if len(rej) != 1 else "", " and ".join("%s (φ %s)" % (k, sg(LETTERS[k]["phi"])) for k in rej)
            or "–", " and ".join("%s (%s)" % (k, sg(LETTERS[k]["phi"])) for k in acc) or "none",
            "; instead A's L1 goes with D's Q (φ %s) and D's L2 with A's P (φ %s)" % (sg(p1["phi"]), sg(p2["phi"]))
            if p1 and p2 else "",
            " So the pollen carry L1·Q and L2·P twice each: two ordinary chromosome pairs, not the ring of four the "
            "assembly implies." if p1 and p2 else "", tv_txt))
    if XS:
        s = XS
        R.p("<b>X arm</b> (%s): of %d nuclei whose chr3 and chr4 could both be read, %d carry X twice and %d never%s. "
            "%s" % (nom_x(), s["N"], s["RA"], s["AR"], ": pollen without X are lost" if s["lost"] else "",
                    "X twice is rarer than if chr3 and chr4 assorted independently (%s of survivors against 1/3; %s), so "
                    "the chain mostly segregates the balanced way; about %s of meioses go the unbalanced (adjacent) way "
                    "(95%% %s–%s)." % (pc(s["share"]), H(pe(s["p_indep"])), pc(s["adjacent"]),
                                       pc(s["adjacent_ci"][0]), pc(s["adjacent_ci"][1])) if s["not_indep"] else
                    "X twice is as common as independent assortment of chr3 and chr4 predicts (%s; %s)." % (
                        pc(s["share"]), H(pe(s["p_indep"])))))
    if CS.get("agree") is not None:
        hb = CS.get("hic")
        R.p("<b>chr5 × chr6</b>: around one coupling point (%s ~%.0f Mb × %s ~%.0f Mb) %s of %s nuclei carry matching "
            "copies (both hap1 or both hap2), although the assembly has %s between them%s." % (
                CS.get("p"), CS.get("pt_p", float("nan")), CS.get("q"), CS.get("pt_q", float("nan")), CS["agree"],
                CS["n"], "no join" if not CS.get("joined") else "a join",
                "" if not hb else (" and the tissue Hi-C shows no junction at the level of a real join (strongest "
                                   "block %.1f×, rank %d of %d pairs)" if CS.get("junction") is False else
                                   " (tissue Hi-C, strongest block %.1f×, rank %d of %d pairs)") % (
                    hb["x"], hb["rank"], hb["of"])))
    if DS:
        R.p("<b>Crossover deserts</b> (<i>D. paradoxa</i>): " + deserts_sentence())
    vp, vb = RV.get("paradoxa", {}), RV.get("binata", {})
    if vp.get("after") is not None:
        R.p("<b>Repeats</b>: in <i>D. paradoxa</i> the repeats %s the crossover ups and downs (they add %.1f%% after "
            "position and markers; repeat tracks shifted at random add up to %.1f%%; %s). In <i>D. binata</i> %s." % (
                "do not explain" if vp["word"] == "no" else "partly explain", 100 * vp["after"].obs,
                100 * vp["after"].null_hi, H(pf(vp["after"].p)),
                ("they %s: repeats alone explain %.0f%% (%s)%s" % (
                    "matter" if vb["word"] != "no" else "do not explain them either", 100 * vb["alone"].obs,
                    H(pf(vb["alone"].p)),
                    ", and crossovers within %s Mb of satellite arrays run at %.2f× the expected rate (%s)" % (
                        ff(array_mb("binata"), 1), vb["sat"].obs, H(pf(vb["sat"].p))) if vb.get("sat") is not None and
                    ok(vb["sat"].obs) else "")) if vb.get("alone") is not None else "–"))
    rows = [["pollen nuclei analysed", fi(nB), fi(nP)],
            ["crossovers called (after the dosage masks)", fi(coB), fi(coP)],
            ["crossovers per nucleus, mean", ff(v5n(sB, "co_mean"), 2), ff(v5n(sP, "co_mean"), 2)],
            ["genetic map, cM", fi(v5n(sB, "map_cm")), fi(v5n(sP, "map_cm"))],
            ["landscape: chromosome ends against middles", v5s(sB, "shape"), v5s(sP, "shape")],
            ["crossover deserts (runs without a crossover, longer than chance allows)"] +
            ["%s ≥ %s Mb (chance %s)" % (fi(v5n(sp, "dead_spots")), fi(v5n(sp, "dead_spots_min_mb")),
                                         v5c(sp, "dead_spots_chance")) for sp in (sB, sP)]]
    rr_ = []
    for sp in ("binata", "paradoxa"):
        w = RRD.get(sp, {}).get("windows")
        rr_.append(pc(np.average(w["all repeats"], weights=w.len_mb)) if w is not None else "–")
    rows.append(["repeats, share of the crossover reference (CARP)"] + rr_)
    for lab, test, what in (("repeats add, after position + markers (p against shifted tracks)",
                             "deviance explained: rotation null", "repeats after position + markers"),
                            ("crossovers near satellite arrays ≥ 50 kb, × expected",
                             "crossover density within", "satellite arrays >= 50 kb")):
        cells = []
        for sp in ("binata", "paradoxa"):
            r = test_row(sp, test, what)
            if r is None or not ok(r.obs):
                cells.append("–")
            elif test.startswith("deviance"):
                cells.append("%+.1f%% (%s)" % (100 * r.obs, pf(r.p)))
            else:
                cells.append("%.2f× within %s Mb (%s; %d arrays)" % (r.obs, ff(array_mb(sp), 1), pf(r.p), r.n))
        rows.append([lab] + cells)
    rows.append(["do the repeats explain where crossovers fall?", vb.get("word", "–"), vp.get("word", "–")])
    R.table(["", "D. binata", "D. paradoxa"], rows)
    st = []
    if ABCD:
        st.append(["chr1/chr2: joins %s rejected by the pollen" % ", ".join(j["letter"] for j in JS["bad"]
                                                                          if j["letter"]),
                   "solid in the pollen (φ near 0 against near 1 for the joins they keep)",
                   {True: "the tissue Hi-C reads them as real: cytology decides (a cross of four or two pairs at "
                          "meiosis I)",
                    False: "the tissue Hi-C agrees: most likely scaffolding errors",
                    None: "the tissue Hi-C is unclear: cytology decides"}[tv]])
    if XS:
        st.append(["X on chr4 (hap1) and chr3 (hap2)%s" % ("; X-less pollen lost" if XS["lost"] else ""),
                   "solid (%d of %d twice, %d never)" % (XS["RA"], XS["N"], XS["AR"]),
                   "the adjacent share is provisional: X-twice nuclei are probably undercounted (Limitations)"])
    if CS.get("agree") is not None:
        st.append(["chr5 × chr6 middles inherited as one unit", "solid as a pollen pattern (%s of %s nuclei)" % (
            CS["agree"], CS["n"]), "cause open: HiFi split reads, chromosome painting, pollen viability"])
    if DS:
        st.append(["deserts on structure", "%d of %d deserts" % (DS["n_structure"], DS["n"]),
                   "the inversions come from the hap2-on-hap1 alignment and are not yet confirmed by reads"])
    if vp:
        st.append(["repeats and crossovers: D. paradoxa %s, D. binata %s" % (vp.get("word", "–"), vb.get("word", "–")),
                   "tested against shifted repeat tracks",
                   "D. paradoxa has few crossovers (power); in D. binata repeats and position go together"])
    R.p("<b>Status</b>")
    R.table(["item", "how sure", "what remains"], st)
    bad = [j["letter"] or j["copy"] for j in JS["bad"]]
    slide("The pollen test the assembly: %s, %s; %s%s." % (
        "%d join%s rejected (%s)" % (len(bad), "s" if len(bad) != 1 else "", ", ".join(bad)) if bad else
        "every join confirmed",
        "the X move confirmed" if "chr4_hap1" in JBY and JBY["chr4_hap1"]["phi"] >= 0.5 else "X unclear",
        "a hidden chr5–chr6 coupling found" if CS.get("agree") is not None and not CS.get("joined") else
        "no hidden coupling",
        "; %d of %d crossover deserts sit on structure; repeats explain %s crossovers%s" % (
            DS["n_structure"], DS["n"],
            "<i>D. binata</i>'s" if vb.get("word", "no") != "no" else "neither species'",
            ", not <i>D. paradoxa</i>'s" if vp.get("word") == "no" and vb.get("word", "no") != "no" else "")
        if DS else ""))
    emb("D0", "overview_assembly_vs_pollen", "The assembly and what the pollen say",
        "hap1 and hap2 copies to scale; the pollen confirm %d joins%s%s." % (
            len(JS["good"]), ", reject %s" % " and ".join(bad) if bad else "",
            ", and couple %s with %s" % (CS.get("p"), CS.get("q")) if CS.get("agree") is not None else ""))
    look("black circles are joins the pollen inherit as one piece, red crosses joins they do not; the purple arrow "
         "links two chromosomes the pollen inherit together; pink strips are stretches inverted between the "
         "haplotypes.")


# ---- 1 material and pipeline
def section_material():
    h2("s1", "1. Material and pipeline")
    slide("Pollen from one plant per species: nuclei sequenced one by one (10x single-cell RNA), genotyped at the "
          "plant's own heterozygous sites, crossovers called per nucleus.")
    rows = []
    for lab, key, d in (("barcodes tested", "barcodes_tested", 0), ("nuclei called", "cells_called", 0),
                        ("haploid nuclei selected (analysed)", "cells_selected", 0),
                        ("informative markers (good)", "good_n", 0), ("markers per Mb", "good_per_mb", 1),
                        ("two-haplotype assembly, 2C Mb", "assembly_2C_chromosomes_mb", 0),
                        ("crossovers called", "co_total", 0),
                        ("median crossover interval, kb (resolution)", "width_median_kb", 0)):
        rows.append([lab] + [(fi(v5n(sp, key)) if d == 0 else ff(v5n(sp, key), d)) for sp in ("binata", "paradoxa")])
    rows.insert(3, ["near-identical pairs (two nuclei of one grain): observed (chance)"] +
                ["%s (chance %s)" % (fi(v5n(sp, "dup_pairs")), v5c(sp, "dup_pairs_chance"))
                 for sp in ("binata", "paradoxa")])
    rows.insert(6, ["crossover reference, 1C Mb (sequences)"] +
                ["%s (%s)" % (fi(v5n(sp, "reference_1C_mb")), fi(v5n(sp, "reference_sequences")))
                 for sp in ("binata", "paradoxa")])
    rows.append(["crossover-calling settings (from the eye review)"] + [v5s(sp, "score_chosen")
                                                                        for sp in ("binata", "paradoxa")])
    R.table(["", "D. binata (%s)" % B, "D. paradoxa (%s)" % P], rows)
    if PIECES:
        R.p("The <i>D. paradoxa</i> crossover reference (%s) takes one copy of each chromosome from the two-haplotype "
            "assembly and cuts it where two arms meet, so that a join the pollen might not inherit cannot fake a "
            "crossover: %s. Every <i>D. paradoxa</i> position in this report is in the coordinates of these copies."
            % (P, "; ".join("%s = %s" % (c, " + ".join(p_["piece"].replace(c + "_", "") for p_ in PIECES
                                                       if p_["source"] == c)
                                         if any(p_["source"] == c and p_["piece"] != c for p_ in PIECES) else "whole")
                            for c in REFC)))
    R.p("<b>Names used throughout.</b> hap1 / hap2 = the two haplotypes of the assembly (blue / orange in every "
        "figure). %s. %s. φ (phi) = how often two windows are inherited together across nuclei, corrected for allele "
        "frequencies: +1 = always together (one chromosome), 0 = independent, −1 = never together." % (
            nom_arms(), nom_x()))
    R.table(["step", "script", "what it gives"], [
        ["select haploid pollen nuclei, call markers", "Snakefile (cell_qc, markers)", "%s / %s nuclei" % (
            fi(v5n("binata", "cells_selected")), fi(v5n("paradoxa", "cells_selected")))],
        ["call crossovers per nucleus, settings chosen from an eye review", "co_sweep.py, review_page.py",
         "%s / %s crossovers" % (fi(v5n("binata", "co_total")), fi(v5n("paradoxa", "co_total")))],
        ["landscapes, deserts, interference", "supervisor_report.py (v5)", "section 2"],
        ["which pieces the pollen inherit together; every join tested; Hi-C at each join",
         "cosegregation_2n.py", "sections 3a, 3c"],
        ["X copies per nucleus", "x_arm_dosage.py", "section 3b"],
        ["how the two haplotypes align (synteny blocks)", "translocation_map.py", "section 3d"],
        ["CARP repeat annotation; repeats at joins, rDNA, FISH candidates", "repeats.smk (CARP), repeat_context.py",
         "section 4"],
        ["repeats against crossovers and breakpoints, with chance baselines", "repeat_recombination.py", "section 4"],
        ["this report", "final_report.py", "every number below is read from the files listed in the appendix"]])


# ---- 2 landscapes
def section_landscapes():
    h2("s2", "2. Crossover landscapes")
    sB, sP = "binata", "paradoxa"
    slide("<i>D. paradoxa</i>: %s crossovers per nucleus, %s, with %s deserts (chance %s); <i>D. binata</i>: %s, %s." % (
        ff(v5n(sP, "co_mean"), 1), shape_word(sP), fi(v5n(sP, "dead_spots")), v5c(sP, "dead_spots_chance"),
        ff(v5n(sB, "co_mean"), 1), shape_word(sB)))
    emb(v5fig("fig_crossovers_per_grain.png"), "crossovers_per_nucleus", "Crossovers per nucleus",
        "Mean %s (D. binata) and %s (D. paradoxa) crossovers per nucleus." % (ff(v5n(sB, "co_mean"), 1),
                                                                              ff(v5n(sP, "co_mean"), 1)))
    look("a nucleus carries half of its meiosis's crossovers on average; one crossover per chromosome pair gives n / 2 "
         "per nucleus.")
    emb(v5fig("fig_landscape_relative.png"), "landscape_relative", "Where crossovers fall along chromosomes",
        "D. binata: %s; D. paradoxa: %s." % (v5s(sB, "shape"), v5s(sP, "shape")))
    look("crossover density along the chromosome (all chromosomes pooled) next to the markers and genes that make "
         "crossovers detectable: a curve that only follows the markers is not a finding.")
    R.table(["", "D. binata", "D. paradoxa"], [
        ["crossovers per nucleus, mean (SD)", "%s (%s)" % (ff(v5n(sB, "co_mean"), 2), ff(v5n(sB, "co_sd"), 2)),
         "%s (%s)" % (ff(v5n(sP, "co_mean"), 2), ff(v5n(sP, "co_sd"), 2))],
        ["crossovers per chromosome pair per meiosis (2 × per nucleus / chromosomes)",
         ff(2 * v5n(sB, "co_mean") / v5n(sB, "chromosomes"), 2), ff(2 * v5n(sP, "co_mean") / v5n(sP, "chromosomes"), 2)],
        ["genetic map, cM", fi(v5n(sB, "map_cm")), fi(v5n(sP, "map_cm"))],
        ["ends against middles", v5s(sB, "shape"), v5s(sP, "shape")],
        ["close crossover pairs, observed / shuffled (interference < 1)", v5s(sB, "interference_short_ratio"),
         v5s(sP, "interference_short_ratio")],
        ["deserts (runs without a crossover, longer than chance allows)"] +
        ["%s ≥ %s Mb (chance %s)" % (fi(v5n(sp, "dead_spots")), fi(v5n(sp, "dead_spots_min_mb")),
                                     v5c(sp, "dead_spots_chance")) for sp in (sB, sP)],
        ["median crossover interval, kb", fi(v5n(sB, "width_median_kb")), fi(v5n(sP, "width_median_kb"))]])
    bullets(["<i>D. binata</i>: crossovers %s. Section 4 compares this with the repeats." % shape_word(sB),
             "<i>D. paradoxa</i>: %s crossovers per chromosome pair per meiosis, %s; its %s deserts (chance %s) are "
             "what sections 3d and 4 take apart." % (
                 ff(2 * v5n(sP, "co_mean") / v5n(sP, "chromosomes"), 1), shape_word(sP), fi(v5n(sP, "dead_spots")),
                 v5c(sP, "dead_spots_chance")),
             "The <i>D. paradoxa</i> calls are placed less precisely (median interval %s kb against %s kb), so its "
             "maps are at ~%s Mb resolution." % (fi(v5n(sP, "width_median_kb")), fi(v5n(sB, "width_median_kb")),
                                                 ff(v5n(sP, "width_median_kb") / 1000, 0))])


# ---- 3 karyotype vs assembly
def diff_rows():
    """The differences from the assembly, region by region: [region, assembly, pollen, tissue Hi-C, meaning]."""
    rows = []
    if ABCD:
        tv, _, _ = tissue_vs_pollen()
        p1, p2 = pair_of("L1", "Q"), pair_of("L2", "P")
        ok_ = [k for k in "ABCD" if LETTERS[k]["phi"] >= 0.5]
        no_ = [k for k in "ABCD" if LETTERS[k]["phi"] < 0.5]
        rows.append(["chr1/chr2 arms",
                     ", ".join("%s = %s·%s" % (k, LETTERS[k]["left"], LETTERS[k]["right"]) for k in "ABCD"),
                     "%s%s%s" % ("%s confirmed (φ %s)" % (", ".join(ok_), ", ".join(sg(LETTERS[k]["phi"]) for k in ok_))
                                 if ok_ else "",
                                 "; %s not (φ %s)" % (", ".join(no_), ", ".join(sg(LETTERS[k]["phi"]) for k in no_))
                                 if no_ else "",
                                 "; L1 goes with Q (φ %s), L2 with P (φ %s)" % (sg(p1["phi"]), sg(p2["phi"]))
                                 if p1 and p2 else ""),
                     "A %s, D %s%s" % ("%.1f×" % LETTERS["A"]["hic"] if ok(LETTERS["A"]["hic"]) else "–",
                                       "%.1f×" % LETTERS["D"]["hic"] if ok(LETTERS["D"]["hic"]) else "–",
                                       "; confirmed joins %.0f–%.0f×" % (JS["hic_lo"], JS["hic_hi"])
                                       if ok(JS["hic_lo"]) else ""),
                     {True: "tissue and pollen disagree: cytology decides",
                      False: "A and D are assembly errors" if no_ else "the assembly is right",
                      None: "the pollen reject %s; the Hi-C is unclear" % " and ".join(no_) if no_ else
                      "the assembly is right"}[tv]])
    if XS and "chr4_hap1" in JBY:
        s = XS
        jx = [JBY[c] for c in ("chr4_hap1", "chr3_hap2") if c in JBY]
        rows.append(["chr3/chr4/X", "X (%.1f Mb) starts chr4 in hap1 and ends chr3 in hap2" % JBY["chr4_hap1"]["pos"],
                     "both X joins %s (φ %s); X twice in %d of %d nuclei, never in %d" % (
                         "hold" if all(j["phi"] >= 0.5 for j in jx) else "tested",
                         ", ".join(sg(j["phi"]) for j in jx), s["RA"], s["N"], s["AR"]),
                     ", ".join("%s|%s %s" % (j["left"], j["right"], ("%.1f×" % j["hic"]) if ok(j["hic"]) else "–")
                               for j in jx) + (" (0.0× = unmappable)" if any(ok(j["hic"]) and j["hic"] < 0.05
                                                                            for j in jx) else ""),
                     "the assembly is right%s%s" % ("; X-less spores die" if s["lost"] else "",
                                                    "; ~%s of meioses adjacent" % pc(s["adjacent"])
                                                    if s["not_indep"] else "")])
    if CS.get("agree") is not None:
        hb = CS.get("hic")
        rows.append(["%s × %s" % (str(CS.get("p")).split("_")[0], str(CS.get("q")).split("_")[0]),
                     "separate chromosomes in both haplotypes" if not CS.get("joined") else "joined",
                     "middles inherited together: %s of %s nuclei around %.0f × %.0f Mb" % (
                         CS["agree"], CS["n"], CS.get("pt_p", float("nan")), CS.get("pt_q", float("nan"))),
                     "strongest block %.1f× (rank %d of %d)%s" % (
                         hb["x"], hb["rank"], hb["of"], ": no junction" if CS.get("junction") is False else "")
                     if hb else "–",
                     "a link the assembly misses: a hidden swap or selection" if not CS.get("joined") else "–"])
    if DS and SYN is not None:
        inv = [(c, a, b) for c in REFC for a, b in merged_inversions(c)]
        bc = DS["by_class"]
        ir = bc[bc.cls == "inverted between the haplotypes"] if len(bc) else bc
        nd = int((DS["rows"].inverted >= DESERT_RULES["inverted"]).sum()) if len(DS["rows"]) else 0
        rows.append(["inversions between the haplotypes",
                     "%d stretches ≥ 0.5 Mb on the crossover copies (%.0f Mb)" % (len(inv), sum(b - a for _, a, b in inv)),
                     "crossovers at %s× the genome mean inside them; %d of %d deserts on them" % (
                         ff(ir.ratio.iloc[0], 2) if len(ir) else "–", nd, DS["n"]),
                     "not checked here", "fits heterozygous inversions; confirm with reads across the breakpoints"])
    return rows


DIFF_HEAD = ["region", "the assembly says", "the pollen say", "tissue Hi-C (× unrelated chromosomes)", "so"]


def section_karyotype():
    h2("s3", "3. The karyotype the pollen inherit, against the assembly")
    R.p("The assembly is built from somatic tissue (HiFi reads, scaffolded with Hi-C). The pollen are the plant's own "
        "meiotic products: two pieces that sit on one chromosome are inherited together (φ ≈ +1) unless a crossover "
        "separates them; pieces on different chromosomes are inherited independently (φ ≈ 0). Every join of the "
        "assembly where one arm meets another was tested this way (cosegregation_2n.py, the 20 Mb on either side of "
        "the join, %s nuclei)." % (CPL.get("nuclei") or "–"))
    dr = diff_rows()
    if dr:
        R.p("<b>Differences from the assembly, region by region</b> (details in 3a–3d):")
        R.table(DIFF_HEAD, dr)
    if JOINS:
        R.table(["join in the assembly", "pollen φ (95%)", "nuclei agreeing", "pollen verdict",
                 "tissue Hi-C across the join (× unrelated chromosomes)"],
                [["%s%s: %s | %s at %.1f Mb" % (j["letter"] + " = " if j["letter"] else "", j["copy"], j["left"],
                                                j["right"], j["pos"]),
                  "%s (%s to %s)" % (sg(j["phi"]), sg(j["lo"]), sg(j["hi"])), "%s of %s" % (j["agree"], j["n"]),
                  "confirmed" if j["phi"] >= 0.5 else "REJECTED", "%.1f×" % j["hic"] if ok(j["hic"]) else "–"]
                 for j in JOINS])
        R.p("<span class='small'>Hi-C: contacts in the 10 Mb either side of the join, over the median between copies "
            "of different chromosomes (cosegregation_2n.py). 0.0× means no usable reads there (an unmappable "
            "stretch), not proof of no contact.</span>")

    # 3a
    h2("s3a", "3a. chr1 / chr2: the arms L1, L2, P, Q")
    if ABCD:
        La, Lb, Lc, Ld = (LETTERS[k] for k in "ABCD")
        p1, p2 = pair_of("L1", "Q"), pair_of("L2", "P")
        tv, _, tv_short = tissue_vs_pollen()
        if p1 and p2:
            slide("The pollen inherit L1 with Q and L2 with P (φ %s, %s): two ordinary chromosome pairs. The assembly's "
                  "A (L1·P, φ %s) and D (L2·Q, φ %s) are not inherited as one piece." % (
                      sg(p1["phi"]), sg(p2["phi"]), sg(La["phi"]), sg(Ld["phi"])))
        else:
            slide("The pollen test the four chr1/chr2 copies A–D: φ A %s, B %s, C %s, D %s." % tuple(
                sg(LETTERS[k]["phi"]) for k in "ABCD"))
        emb("D1", "chr1_chr2_arms", "chr1/chr2: assembly against pollen",
            "A φ %s, D φ %s; the pollen pair L1 with Q and L2 with P; a cross of four against two pairs." % (
                sg(La["phi"]), sg(Ld["phi"])))
        look("in a, red joins are the ones the pollen do not inherit; in b, each bar is a chromosome as the pollen "
             "inherit it, with the copy each arm comes from. If the assembly were right, the pollen would show A and D "
             "as one piece each (φ near +1), as they do for B and C.")
        pc_ = [c for c in arm_copies("P")]
        qc_ = [c for c in arm_copies("Q")]
        same_hap = len({hap_of(c) for c in pc_}) == 1 and len({hap_of(c) for c in qc_}) == 1 and pc_ and qc_
        plan_none = [p_ for p_, r in PLAN.items() if "homolog none" in r.get("how", "")]
        bullets([
            "Names: %s." % nom_arms(),
            "Assembly: L1 sits with P on A but with Q on C; L2 sits with P on B but with Q on D. At meiosis that "
            "needs a cross of four chromosomes pairing arm by arm (panel c).",
            "Pollen: B (φ %s) and C (φ %s) hold; A (φ %s) and D (φ %s) are coin flips (%s of %s and %s of %s nuclei "
            "agree).%s" % (sg(Lb["phi"]), sg(Lc["phi"]), sg(La["phi"]), sg(Ld["phi"]), La["agree"], La["n"],
                           Ld["agree"], Ld["n"],
                           " Instead A's L1 goes with D's Q (φ %s) and D's L2 with A's P (φ %s): the pollen behave as "
                           "if the plant had L1·Q twice and L2·P twice (panel d)." % (sg(p1["phi"]), sg(p2["phi"]))
                           if p1 and p2 else ""),
            "Tissue Hi-C across A: %s; across D: %s. The confirmed joins B and C read %s and %s." % (
                hic_word(La["hic"]), hic_word(Ld["hic"]), "%.1f×" % Lb["hic"] if ok(Lb["hic"]) else "–",
                "%.1f×" % Lc["hic"] if ok(Lc["hic"]) else "–"),
            "Not a crossover hotspot: the windows next to A's join are the L1 end that C also carries and the P start "
            "that B also carries; a crossover there would break C or B too (they hold at φ %s and %s), and a "
            "crossover cannot give an arm a partner on another chromosome." % (sg(Lc["phi"]), sg(Lb["phi"]))
            if Lb["phi"] >= 0.5 and Lc["phi"] >= 0.5 else "",
            "How the assembly could get there: both P copies sit in %s (%s) and both Q copies in %s (%s)%s, so "
            "scaffolding each haplotype on its own could only join L1 to P in one haplotype and L2 to Q in the "
            "other." % (hap_of(pc_[0]), ", ".join(pc_), hap_of(qc_[0]), ", ".join(qc_),
                        "; the split plan lists %s as 'homolog none'" % " and ".join(plan_none) if plan_none else "")
            if same_hap else "",
            tv_short + (" Readings: the Hi-C is misled near rDNA (nucleolar clustering), the germline karyotype "
                        "differs from the tissue's, or meiosis rewires the arms. Meiosis I cytology decides: a cross "
                        "of four (assembly) or two pairs (pollen)." if tv else
                        " The joins A and D are then most likely scaffolding errors; meiosis I cytology would confirm "
                        "two ordinary pairs." if tv is False else " Meiosis I cytology decides.")])
    else:
        missing_note("The karyotype test does not name A–D; section 3a skipped.")

    # 3b
    h2("s3b", "3b. chr3 / chr4 / X: an arm that moved")
    if XS:
        s = XS
        slide("X sits on chr4 in hap1 and on chr3 in hap2. %d of %d nuclei carry X twice, %d never%s%s." % (
            s["RA"], s["N"], s["AR"], ": pollen without X are lost" if s["lost"] else "",
            ", and ~%s of meioses go the unbalanced way" % pc(s["adjacent"]) if s["not_indep"] else ""))
        emb("D2", "chr3_chr4_X", "chr3/chr4/X: the chain and its spores",
            "X twice in %d of %d, never in %d; independent assortment %s." % (
                s["RA"], s["N"], s["AR"], "rejected (%s)" % pe(s["p_indep"]) if s["not_indep"] else "not rejected"))
        look("panel c is the test: a nucleus with chr4_hap1 and chr3_hap2 should carry X twice, one with chr4_hap2 and "
             "chr3_hap1 none. If X-less spores survived, the 'X never' box would hold about as many nuclei as "
             "'X twice'.")
        ct = s["table"]
        lab4 = {"REF": "chr4_hap1 (X + core)", "ALT": "chr4_hap2 (core only)", "mixed": "chr4 unclear",
                "few": "chr4 too few molecules"}
        lab3 = {"REF": "chr3_hap1 (no X)", "ALT": "chr3_hap2 (+ X)", "mixed": "chr3 unclear", "few": "too few"}
        cols = [c for c in ("REF", "ALT", "mixed", "few") if c in ct.columns]
        R.table(["which chr4 (next to X)"] + [lab3[c] for c in cols],
                [[lab4.get(r, r)] + [int(ct.loc[r, c]) for c in cols] for r in ("REF", "ALT", "mixed", "few")
                 if r in ct.index])
        jx = [JBY.get("chr4_hap1"), JBY.get("chr3_hap2")]
        bullets([
            "Names: %s. Which chr4 and chr3 a nucleus got is read next to the joins (chr4's core next to X, chr3's "
            "end); how many X copies it has is read from X itself (molecules over the expected number, and whether "
            "both alleles are interleaved)." % nom_x(),
            "The X joins in the pollen: %s." % "; ".join("%s %s|%s φ %s" % (j["copy"], j["left"], j["right"],
                                                                           sg(j["phi"])) for j in jx if j),
            "X copies match the prediction from chr3 and chr4 in %d of %d nuclei with both read; the %d nuclei with "
            "an unclear side carry X once in %d cases." % (s["agree"], s["agree_n"], s["unres"], s["unres_one"]),
            "X twice: %d nuclei, %s chr4_hap1 + chr3_hap2, coverage %s× (median %s; one-copy nuclei %s–%s×), allele "
            "switching %s–%s against at most %s for one copy." % (
                s["twice_obs_all"], "all" if s["twice_from_RA"] == s["twice_obs_all"] else "%d of them" %
                s["twice_from_RA"],
                "%.2f–%.2f" % (s["twice_cov"][0], s["twice_cov"][2]) if s["twice_cov"] else "–",
                ff(s["twice_cov"][1], 2) if s["twice_cov"] else "–",
                ff(s["one_cov"][0], 2) if s["one_cov"] else "–", ff(s["one_cov"][1], 2) if s["one_cov"] else "–",
                ff(s["twice_switch"][0], 2) if s["twice_switch"] else "–",
                ff(s["twice_switch"][1], 2) if s["twice_switch"] else "–", ff(s["one_switch_max"], 2)),
            "X never: %d. If X-twice and X-less spores were made and survived alike they would be equally common "
            "(%d against %d: %s)%s." % (
                s["AR"], s["RA"], s["AR"], H(pe(s["p_twice_never"])),
                ": the X-less ones are lost, most likely as dead microspores (a haploid missing %.1f Mb)" % (
                    JBY["chr4_hap1"]["pos"] if "chr4_hap1" in JBY else float("nan")) if s["lost"] else ""),
            "X is a one-way move (a terminal piece sits on chr4 in one haplotype and on chr3 in the other), not a "
            "swap: so one unbalanced spore type gains X and the other loses it."])
        R.table(["model (X-less spores die in each)", "X twice among survivors", "observed %d of %d" % (s["RA"], s["N"])],
                [["chr3 and chr4 assort independently", "1 in 3", pe(s["p_indep"])],
                 ["held together only when X has a crossover (%s of meioses), else independent" % pc(s["xco_meioses"]),
                  pc(s["chain_only_share"]), pe(s["p_chain_only"])],
                 ["a chain of four, balanced (alternate) in most meioses", "fits: %s of meioses adjacent (95%% %s–%s)"
                  % (pc(s["adjacent"]), pc(s["adjacent_ci"][0]), pc(s["adjacent_ci"][1])), "–"]])
        bullets([
            "Adjacent share: if only X-less spores die, X-twice nuclei are a/(2 − a) of survivors, where a is the "
            "share of meioses that segregate the chain the unbalanced way; %d of %d gives a = %s (95%% %s–%s). Then "
            "~%s of all microspores abort from this exchange alone." % (
                s["RA"], s["N"], pc(s["adjacent"]), pc(s["adjacent_ci"][0]), pc(s["adjacent_ci"][1]), pc(s["abort"])),
            "X recombines: %d of %d one-copy nuclei switch haplotype inside X (%s), so X has a crossover in ~%s of "
            "meioses (a crossover shows in half of a meiosis's products)." % (
                s["n_xco"], s["n_one"], pc(s["xco_share"]), pc(s["xco_meioses"])),
            "Open: a chain needs a chiasma in X, seen in only ~%s of meioses; if that were the only link, %s of "
            "survivors would carry X twice, against %s observed (%s). %s" % (
                pc(s["xco_meioses"]), pc(s["chain_only_share"]), pc(s["share"]), H(pe(s["p_chain_only"])),
                "Either X crossovers are undercounted, or something besides an X chiasma holds chr3 and chr4 "
                "together." if ok(s["p_chain_only"]) and s["p_chain_only"] < 0.05 else
                "The data fit chains that form only with an X crossover.")])
        emb(XARM_PNG, "x_arm_per_nucleus", "X copies per nucleus (x_arm_dosage.py)",
            "Per nucleus: which chr3/chr4 it got and how many X copies it carries.")
        look("each row is one nucleus: the chr4 and chr3 calls next to the joins, then X's coverage and alleles; X twice "
             "shows as both alleles interleaved at higher coverage.")
    else:
        missing_note("x_arm_dosage.py output missing; section 3b skipped.")

    # 3c
    h2("s3c", "3c. chr5 × chr6: two chromosomes inherited as one unit")
    if CS.get("agree") is not None or CS.get("pt_p") is not None:
        hb = CS.get("hic")
        slide("The middles of %s and %s travel together: %s of %s nuclei carry matching copies around one coupling "
              "point (~%.0f × ~%.0f Mb)%s." % (
                  str(CS.get("p")).split("_")[0], str(CS.get("q")).split("_")[0], CS.get("agree", "–"), CS.get("n", "–"),
                  CS.get("pt_p", float("nan")),
                  CS.get("pt_q", float("nan")),
                  ", with no join in the assembly" + (" and no junction in the Hi-C" if CS.get("junction") is False
                                                      else "") if not CS.get("joined") else ""))
        emb("D3", "chr5_chr6_coupling", "chr5 × chr6 coupling",
            "One coupling point (R² %s); joined ends R² %s; two explanations." % (ff(CS.get("r2"), 2),
                                                                                  ff(CS.get("r2_ends"), 2)))
        look("in b, red spreads over the middles of both chromosomes and fades towards all four ends if the ends "
             "recombine away and the middles never do; a fusion of two ends would instead be reddest in one corner.")
        R.table(["", "value"], [
            ["coupling point (95% interval over resampled nuclei)", "%s %.0f Mb (%.0f–%.0f) × %s %.0f Mb (%.0f–%.0f)" % (
                CS.get("p", "–"), CS.get("pt_p", float("nan")), CS.get("lo_p", float("nan")),
                CS.get("hi_p", float("nan")), CS.get("q", "–"), CS.get("pt_q", float("nan")),
                CS.get("lo_q", float("nan")), CS.get("hi_q", float("nan")))],
            ["nuclei with matching copies around the point", "%s of %s" % (CS.get("agree", "–"), CS.get("n", "–"))],
            ["how much of the pattern one coupling point explains (R²); its strength c (1 = complete)",
             "R² %s, c %s" % (ff(CS.get("r2"), 2), ff(CS.get("c"), 2))],
            ["the same for joined ends (a fusion), best of the 4 end pairs", "R² %s" % ff(CS.get("r2_ends"), 2)],
            ["control: the %s join, predicted from the join alone" % (
                "%s|%s" % (CONTROL["left"], CONTROL["right"]) if CONTROL else "–"),
             "R² %s, c %s" % (ff(CONTROL["r2"], 2), ff(CONTROL["c"], 2)) if CONTROL else "–"],
            ["tissue Hi-C, strongest chr5 × chr6 block (× unrelated)", "%.1f× at %.0f / %.0f Mb (%s × %s), rank %d of %d"
             % (hb["x"], hb["pa"], hb["pb"], hb["a"], hb["b"], hb["rank"], hb["of"]) if hb else "–"],
            ["tissue Hi-C, strongest block of every pair of different chromosomes: median (90th percentile)",
             "%.1f× (%.1f×)" % (CS["scan_median"], CS["scan_p90"]) if "scan_median" in CS else "–"],
            ["tissue Hi-C across the confirmed joins", "%.1f–%.1f×" % (JS["hic_lo"], JS["hic_hi"])
             if ok(JS["hic_lo"]) else "–"]])
        n5, n6 = str(CS.get("p", "chr5")).split("_")[0], str(CS.get("q", "chr6")).split("_")[0]
        inv56 = ["%s %.0f–%.0f" % (c, a, b) for c in COPIES if num_of(c) in (num_of(n5), num_of(n6))
                 for a, b in merged_inversions(c)]
        bullets([
            "Assembly: %s and %s are separate chromosomes in both haplotypes%s, aligned end to end%s." % (
                n5, n6, "" if not CS.get("joined") else " (but a join involves them)",
                " apart from inversions (%s Mb)" % "; ".join(inv56) if inv56 else ""),
            "Pollen: one coupling point explains R² %s of the block (c %s); joined ends explain %s. The point cannot "
            "be narrowed further, because the middles do not recombine." % (
                ff(CS.get("r2"), 2), ff(CS.get("c"), 2), ff(CS.get("r2_ends"), 2)),
            "The desert does not by itself make the linkage: no crossovers keep each chromosome's middle in one piece, "
            "but cannot tie %s to %s. Something must hold the two middles together." % (n5, n6),
            "Explanation 1, a hidden reciprocal swap (panel c1): a spore with one swapped and one normal chromosome "
            "lacks a segment and dies, so only matching combinations reach the pollen; crossovers near the swap points "
            "are suppressed, which fits the deserts. Predicts a junction (HiFi reads split %s → %s near the point) and "
            "a ring or chain of four at meiosis I." % (n5, n6),
            "Explanation 2, no physical link (panel c2): mixed spores are lost or never made. Predicts no junction and "
            "two ordinary pairs at meiosis I.",
            "Crossover calls are unaffected (balanced spores keep each chromosome's coordinates); the genetic map is "
            "not: %s and %s behave as one linkage group." % (n5, n6)])
    else:
        missing_note("No coupling in the karyotype test; section 3c skipped.")

    # 3d
    h2("s3d", "3d. Inversions between the haplotypes, and where the crossovers are missing")
    if DS:
        inv_mb = sum(b - a for c in REFC for a, b in merged_inversions(c))
        sm = DS["struct_markers"]
        slide("%d of %d crossover deserts sit on structure: an inversion between the haplotypes, the coupled chr5 × chr6 "
              "middle or X.%s" % (DS["n_structure"], DS["n"], (
                  " In the inversion and coupling deserts markers run at %.2f–%.2f× their copy's density, so there the "
                  "crossovers are missing, not hidden." % sm) if sm and sm[0] >= DESERT_RULES["markers"] else ""))
        emb("D4", "structure_and_crossovers", "Deserts sit on structure",
            "Crossovers above each chromosome's two copies and their alignment.")
        look("black bars (deserts) over pink twisted ribbons (inverted stretches), the purple coupling bracket or the "
             "X arm; grey straight ribbons (same order) carry the crossovers. A desert over grey ribbons with normal "
             "markers would be unexplained.")
        D = DS["rows"]
        if len(D):
            R.table(["desert (crossover reference copy), Mb", "Mb", "inverted", "translocated", "not aligned",
                     "coupled", "markers × copy", "repeats (copy)", "what it overlaps"],
                    [["%s %.0f–%.0f" % (r.chrom, r.start, r.end), ff(r.mb, 0), pc(r.inverted), pc(r.translocated),
                      pc(r.none), pc(r.coupling), ff(r.markers, 2), "%s (%s)" % (pc(r.repeats), pc(r.repeats_copy)),
                      "; ".join(r.labels)] for r in D.itertuples()])
            R.p("<span class='small'>Deserts: runs of %s Mb windows without a crossover midpoint, ≥ %s Mb, so that "
                "fewer than one is expected by chance (repeat_recombination.py, the same rule as the supervisor "
                "report). Labels: inverted or translocated ≥ %d%% of the desert; coupled ≥ %d%%; not aligned ≥ %d%%; "
                "few markers below %.1f× the copy's density.</span>" % (
                    ff(RRD["paradoxa"].get("win_mb"), 0), ff(DEAD_MB, 0), 100 * DESERT_RULES["inverted"],
                    100 * DESERT_RULES["coupling"], 100 * DESERT_RULES["none"], DESERT_RULES["markers"]))
        bc = DS["by_class"]
        if len(bc):
            R.table(["windows by structure (majority)", "windows", "Mb", "crossovers (midpoints)",
                     "expected at the genome mean", "observed / expected", "share in deserts"],
                    [[r.cls, r.windows, ff(r.mb, 0), r.co, ff(r.exp, 1), ff(r.ratio, 2), pc(r.dead)]
                     for r in bc.itertuples()])
        sr = [assoc_row("paradoxa", v) for v in ("other haplotype: collinear", "other haplotype: inverted",
                                                   "in the coupling interval (pollen)")]

        def rtxt(r, partial=True):
            if r is None:
                return "–"
            v, p = (r.rho_partial, r.p_partial) if partial and ok(r.rho_partial) else (r.rho, r.p_rot)
            return "ρ %s, %s" % (sg(v), H(pf(p)))
        bullets([
            "Inverted stretches between the haplotypes total %.0f Mb on the crossover reference copies (synteny blocks "
            "of hap2 on hap1, merged within 1 Mb)." % inv_mb,
            "Window by window, after distance to the end and markers, crossovers go with 'the other haplotype is "
            "collinear here' (%s) and against inverted stretches (%s) and the coupled middles (%s); p from shifting "
            "each track along its own chromosome." % (rtxt(sr[0]), rtxt(sr[1]), rtxt(sr[2])),
            "Why structure removes crossovers: in an inversion heterozygote a crossover inside the inverted stretch "
            "makes unbalanced chromatids, so recombinant products are lost or crossovers are suppressed there; the "
            "same holds near the points of a swap.",
            "Repeats inside the deserts: %s–%s, against %s–%s for their whole copies (the 'repeats (copy)' column)." % (
                pc(D.repeats.min()), pc(D.repeats.max()), pc(D.repeats_copy.min()), pc(D.repeats_copy.max()))
            if len(D) else "",
            "P and Q have no alignment: their homologs sit in the same haplotype (P on %s, Q on %s), so the "
            "hap2-on-hap1 alignment never compared them. Aligning P with P and Q with Q is the next check." % (
                " and ".join(arm_copies("P")), " and ".join(arm_copies("Q"))) if arm_copies("P") else ""])
    else:
        missing_note("repeat_recombination windows missing; section 3d skipped.")


# ---- 4 repeats
def section_repeats():
    h2("s4", "4. Repeats: what they do and do not explain")
    vp, vb = RV.get("paradoxa", {}), RV.get("binata", {})
    bp = RRD.get("paradoxa", {}).get("breakpoints")
    t = RRD.get("paradoxa", {}).get("tests")
    ja = None
    if t is not None and "R_mb" in t:
        x = t[t.test.astype(str).str.startswith("breakpoints: with rDNA within") & (t.set == "all joins") &
              (pd.to_numeric(t.R_mb, errors="coerce") == 1)]
        ja = x.iloc[0] if len(x) else None
    slide("Repeats %s where <i>D. paradoxa</i> crossovers fall (they add %s after position and markers, %s)%s; in "
          "<i>D. binata</i> %s." % (
              "do not explain" if vp.get("word") == "no" else "partly explain",
              pc(vp["after"].obs, 1) if vp.get("after") is not None else "–",
              pf(vp["after"].p) if vp.get("after") is not None else "–",
              ", nor where the joins sit (rDNA within 1 Mb of %d of %d, random %.1f)" % (ja.obs, ja.n, ja.null_mean)
              if ja is not None and not (ok(ja.p) and ja.p < 0.05) else "",
              "they %s (repeats alone explain %s, %s)" % (
                  "do" if vb.get("word") != "no" else "do not either", pc(vb["alone"].obs, 0), pf(vb["alone"].p))
              if vb.get("alone") is not None else "–"))
    R.p("Repeats: the CARP annotation (version %s) of each haplotype assembly, read onto the crossover reference. "
        "Every comparison below is made against chance: repeat tracks shifted at random along their own chromosome "
        "(keeping their patchiness), or breakpoints moved to random places on the same copy." % (
            RRD.get("paradoxa", {}).get("carp", "?")))
    # 4a
    h2("s4a", "4a. Repeat landscapes")
    for sp, s in SPS:
        emb(rc_path("repeat_landscape", s), "repeat_landscape_%s" % sp, "Repeat landscape, %s" % SP_NAME[sp],
            "Repeat classes along every copy (repeat_context.py).")
        look("whether repeats form blocks (for example in the chromosome middles) or cover every copy at a similar "
             "level; compare with where the crossovers fall (section 4b).")
    cls_rows = []
    for c in ["Ty3/gypsy", "Ty1/copia", "other LTR", "LINE", "DNA transposon", "simple / low complexity",
              "satellite", "5S rDNA", "45S rDNA", "all repeats"]:
        row = [c]
        for sp in ("binata", "paradoxa"):
            w = RRD.get(sp, {}).get("windows")
            row.append(pc(np.average(w[c], weights=w.len_mb), 1) if w is not None and c in w else "–")
        cls_rows.append(row)
    rng_rows = []
    for sp in ("binata", "paradoxa"):
        w = RRD.get(sp, {}).get("windows")
        if w is not None:
            rng_rows.append("%s %s–%s per %g Mb window" % (
                SP_NAME[sp], pc(np.percentile(w["all repeats"], 5)), pc(np.percentile(w["all repeats"], 95)),
                RRD[sp].get("win_mb", float("nan"))))
    R.table(["share of the crossover reference", "D. binata", "D. paradoxa"], cls_rows)
    if rng_rows:
        R.p("All repeats, middle 90%% of the windows: %s." % "; ".join(rng_rows))
    # 4b
    h2("s4b", "4b. Repeats and crossovers")
    emb("R1", "repeats_vs_crossovers", "Repeats against crossovers, both species",
        "Repeats explain crossovers: D. binata %s, D. paradoxa %s (against shifted tracks)." % (
            vb.get("word", "–"), vp.get("word", "–")))
    look("bars that do not clear their black tick explain no more than repeat tracks shifted at random.")
    rows = []
    for lab in ("position", "markers", "genes", "structure", "repeats", "everything"):
        row = [lab]
        for sp in ("binata", "paradoxa"):
            r = test_row(sp, "deviance explained (Poisson", lab)
            row.append(pc(r.obs, 1) if r is not None else "–")
        rows.append(row)
    for lab in ("repeats alone", "repeats after position + markers", "repeats after everything else"):
        row = [lab + " (shifted tracks: 95th percentile; p)"]
        for sp in ("binata", "paradoxa"):
            r = test_row(sp, "deviance explained: rotation null", lab)
            row.append("%s (%s; %s)" % (pc(r.obs, 1), pc(r.null_hi, 1), pf(r.p)) if r is not None else "–")
        rows.append(row)
    R.table(["share of the crossover ups and downs explained (Poisson deviance)", "D. binata", "D. paradoxa"], rows)
    R.p("<span class='small'>position = distance to the nearest chromosome end; markers = informative markers per "
        "Mb (a crossover is only seen between markers); structure = how the other haplotype aligns (D. paradoxa "
        "only) and the coupling interval.</span>")
    arr = []
    for sp in ("binata", "paradoxa"):
        tt = RRD.get(sp, {}).get("tests")
        if tt is None:
            continue
        for r in tt[tt.test.astype(str).str.startswith("crossover density within")].itertuples():
            arr.append([SP_NAME[sp], r.what, fi(r.n), "%.2f×" % r.obs if ok(r.obs) else "fewer than 2 arrays",
                        "%.2f (%.2f–%.2f)" % (r.null_mean, r.null_lo, r.null_hi) if ok(r.null_mean) else "–",
                        pf(r.p) if ok(r.p) else "–"])
    if arr:
        R.table(["species", "arrays", "n", "crossovers within %s Mb, × expected" % ff(array_mb("paradoxa"), 1),
                 "shifted arrays: mean (95%)", "p"], arr)
    for sp, s in SPS:
        emb(rr_path("recombination_vs_repeats", s), "recombination_vs_repeats_%s" % sp,
            "Crossovers and repeats along the chromosomes, %s" % SP_NAME[sp], "Tracks: crossovers, repeats, markers.")
        look("whether the crossover track rises and falls with any repeat track, chromosome by chromosome.")
    top = []
    for sp in ("binata", "paradoxa"):
        a = RRD.get(sp, {}).get("assoc")
        if a is None:
            continue
        a = a.copy()
        a["k"] = a.rho_partial.abs().fillna(a.rho.abs())
        for r in a.sort_values("k", ascending=False).head(6).itertuples():
            top.append([SP_NAME[sp], r.variable, sg(r.rho), pf(r.p_rot), sg(r.rho_partial) if ok(r.rho_partial) else "–",
                        pf(r.p_partial) if ok(r.p_partial) else "–", r.same_sign,
                        "yes" if str(r.beyond_chance) == "True" else "no"])
    if top:
        R.table(["species", "variable (per window)", "ρ with crossovers", "p (shifted)", "ρ after end distance + markers",
                 "p", "same sign on chromosomes", "beyond chance (Bonferroni)"], top)

    def species_bullet(sp):
        v = RV.get(sp, {})
        a = RRD.get(sp, {}).get("assoc")
        if a is None:
            return ""
        e, rp = assoc_row(sp, "distance to the nearest end (Mb)"), assoc_row(sp, "all repeats")
        nrep = int((a.group.isin(["repeat class", "LTR lineage"]) & (a.beyond_chance.astype(str) == "True")).sum())
        nall = int(a.group.isin(["repeat class", "LTR lineage"]).sum())
        w = RRD[sp].get("windows")
        bits = []
        if e is not None and ok(e.p_rot) and e.p_rot < 0.05:
            bits.append("crossovers are %s the chromosome ends (distance to the end ρ %s)" % (
                "densest near" if e.rho < 0 else "scarcest near", sg(e.rho)))
        if rp is not None:
            bits.append("%s where repeats are dense (all repeats ρ %s, %s)" % (
                "fewer crossovers" if rp.rho < 0 else "more crossovers", sg(rp.rho), pf(rp.p_rot))
                        if ok(rp.p_rot) and rp.p_rot < 0.05 else "no relation to the repeat share (all repeats ρ %s, %s)"
                        % (sg(rp.rho), pf(rp.p_rot)))
        if v.get("sat") is not None and ok(v["sat"].obs):
            bits.append("near satellite arrays crossovers run at %.2f× the expected rate (%s)" % (v["sat"].obs,
                                                                                               pf(v["sat"].p)))
        bits.append("%d of %d repeat variables pass the Bonferroni line after distance to the end and markers" % (
            nrep, nall))
        if w is not None:
            bits.append("repeats cover %s–%s of a window (middle 90%%)" % (
                pc(np.percentile(w["all repeats"], 5)), pc(np.percentile(w["all repeats"], 95))))
        tail = ""
        if v.get("word") == "yes, but tangled with position":
            tail = " Repeats and distance to the end go together here, so their effects cannot be pulled apart."
        return "<i>%s</i>: %s.%s" % (SP_NAME[sp], H("; ".join(bits)), tail)
    bullets([species_bullet("binata"), species_bullet("paradoxa"),
             "Power: %s <i>D. paradoxa</i> crossovers at ~%s Mb resolution; a weak repeat effect could be missed, but "
             "the structure effect is found with the same data." % (
                 fi(RRD.get("paradoxa", {}).get("calls")), ff(RRD.get("paradoxa", {}).get("interval_mb"), 1))])
    # 4c
    h2("s4c", "4c. Repeats at the joins and the inversion breakpoints")
    hits = []
    if bp is not None:
        jj = bp[bp.set == "join"]
        rows = []
        for r in jj.itertuples():
            def near(kind):
                m, side = getattr(r, "nearest_%s_mb" % kind, np.nan), getattr(r, "nearest_%s_side" % kind, "")
                if not ok(m):
                    return "none ≥ 10 kb on the copy"
                return "%.1f Mb, %s" % (m, "spans the join" if side == "spans it" else "%s side" % side)
            sat = "%.1f Mb %s (%s)" % (r.nearest_satellite_mb, r.nearest_satellite_side, r.nearest_satellite_family) \
                if ok(r.nearest_satellite_mb) else "none"
            anyk = "%.1f Mb (%.1f kb)" % (r.nearest_45S_any_mb, r.nearest_45S_any_kb) if ok(r.nearest_45S_any_mb) \
                else "none"
            rows.append([r.name, "%s (φ %s)" % (r.verdict, sg(r.phi)), near("45S"), anyk, near("5S"), sat])
            if ok(r.nearest_45S_mb) and r.nearest_45S_mb <= 5:
                jn = re.search(r"(\S+)\|(\S+)$", str(r.name))
                arm = (jn.group(1) if r.nearest_45S_side == "left" else jn.group(2) if r.nearest_45S_side == "right"
                       else "%s|%s" % jn.groups()) if jn else "?"
                hits.append(dict(name=r.name, mb=r.nearest_45S_mb, side=r.nearest_45S_side, arm=arm))
        R.table(["join", "pollen", "nearest 45S array ≥ 10 kb", "nearest 45S piece of any size", "nearest 5S ≥ 10 kb",
                 "nearest satellite array ≥ 50 kb"], rows)
        R.p("<span class='small'>left / right = toward the start / the end of the copy: at a join 'L1|P', left is the "
            "L1 side.</span>")
    rows, seen_R = [], set()
    sets = ["all joins", "joins the pollen inherit", "joins the pollen do not inherit", "haplotype breakpoints"]
    for src in (RRD.get("paradoxa", {}), RR5):
        tt = src.get("tests")
        if tt is None or "R_mb" not in tt:
            continue
        x = tt[tt.test.astype(str).str.startswith("breakpoints: with ")].copy()
        x["R"] = pd.to_numeric(x.R_mb, errors="coerce")
        x["cls"] = x.test.astype(str).str.extract(r"with (\S+) within")[0]
        for R_ in sorted(x.R.dropna().unique()):
            if R_ in seen_R:
                continue
            seen_R.add(R_)
            for st_ in sets:
                y = x[(x.R == R_) & (x.set == st_)]
                if len(y):
                    cell = {r.cls: "%d (random %.1f, %s)" % (r.obs, r.null_mean, pf(r.p)) for r in y.itertuples()}
                    rows.append([st_, R_, int(y.n.iloc[0]), cell.get("45S", "–"), cell.get("5S", "–"),
                                 cell.get("satellite", "–")])
    if rows:
        R.table(["breakpoints", "within, Mb", "n", "with 45S", "with 5S", "with a satellite array"],
                [[a, ff(b, 1), c, d, e, f] for a, b, c, d, e, f in sorted(rows, key=lambda r: (sets.index(r[0]), r[1]))])
        R.p("<span class='small'>Counts of breakpoints with an array within the distance; random = the mean over the "
            "same number of random places on the same copies. The largest distance (the near5 run) was chosen after "
            "looking at the figures, so its p-values are optimistic.</span>")
    emb(rr_path("breakpoints_vs_repeats", P), "breakpoints_vs_repeats", "Repeats at joins and breakpoints",
        "Repeats near joins and inversion breakpoints against random places.")
    look("bars near 1 are what random places on the same copies give; points left of the dotted line have an array "
         "closer than the test distance.")
    emb(rc_path("repeat_joins", P), "repeat_joins_zoom", "The joins up close (repeat_context.py)",
        "The repeat composition around each join of the assembly.")
    look("a join made at an rDNA array or a satellite would show that array right at the join line.")

    def bpc(cls, R_, set_):
        tt = RRD.get("paradoxa", {}).get("tests")
        if tt is None or "R_mb" not in tt:
            return "–"
        x = tt[tt.test.astype(str).str.startswith("breakpoints: with %s within" % cls) & (tt.set == set_) &
               (pd.to_numeric(tt.R_mb, errors="coerce") == R_)]
        return H("%d of %d (random %.1f, %s)" % (x.obs.iloc[0], x.n.iloc[0], x.null_mean.iloc[0], pf(x.p.iloc[0]))) \
            if len(x) else "–"
    arms_hit = sorted({h["arm"] for h in hits})
    bullets([
        "Within 1 Mb of the %s joins: 45S at %s, 5S at %s, a satellite array at %s." % (
            len(JOINS) or "–", bpc("45S", 1, "all joins"), bpc("5S", 1, "all joins"), bpc("satellite", 1, "all joins")),
        "45S arrays within 5 Mb of a join: %s.%s" % (
            "; ".join("%s: %.1f Mb (%s)" % (h["name"], h["mb"], "spans the join" if h["side"] == "spans it" else
                                             "%s side, the %s arm" % (h["side"], h["arm"])) for h in hits) or "none",
            (" All on the same arm end (%s) of different copies, so this is one observation; the 5 Mb distance was "
             "chosen after looking, so it is suggestive at best." % arms_hit[0]) if len(hits) > 1 and
            len(arms_hit) == 1 else (" The 5 Mb distance was chosen after looking, so this is suggestive at best."
                                     if hits else "")),
        "Inversion breakpoints between the haplotypes: rDNA within 1 Mb at %s, satellites at %s." % (
            bpc("rDNA", 1, "haplotype breakpoints"), bpc("satellite", 1, "haplotype breakpoints"))])
    # 4d
    h2("s4d", "4d. The chr5 × chr6 coupling interval")
    emb(rc_path("repeat_coupling", P), "repeat_coupling_interval", "Repeats in the coupling interval",
        "Repeats inside the chr5 and chr6 coupling intervals.")
    look("a family enriched in both intervals (and rare elsewhere) would be a candidate for what pairs them; the "
         "intervals are wide (95% intervals of the coupling point), so most families fall in them by chance.")
    fp = RCD.get("paradoxa", {}).get("fish_probes")
    if fp is not None and len(fp):
        x = fp[fp.in_coupling_interval.astype(str).str.len() > 0].copy()
        if len(x):
            x["score"] = pd.to_numeric(x.fish_score, errors="coerce")
            x = x.sort_values("score", ascending=False)
            R.table(["tandem family", "kind", "in the coupling interval (kb per copy)", "total kb"],
                    [[r.family, r.kind, r.in_coupling_interval, ff(fnum(r.total_kb), 0)] for r in x.head(10).itertuples()])
            R.p("<span class='small'>%d tandem families have arrays inside the intervals; the 10 with the highest FISH "
                "score are listed.</span>" % len(x))
    # 4e
    h2("s4e", "4e. rDNA: what the assembly holds")
    for sp, s in SPS:
        rl = RCD.get(sp, {}).get("rdna_loci")
        if rl is None:
            continue
        if not len(rl):
            R.p("%s: no rDNA array in the assembly long enough to list." % SP_NAME[sp])
            continue
        rl = rl.copy()
        rl["kbf"] = pd.to_numeric(rl.kb, errors="coerce")
        R.table(["%s: rDNA" % SP_NAME[sp], "copy", "position, Mb", "kb", "at a join", "Mb to the nearest join"],
                [[r.type, r.copy, "%.1f" % (fnum(r.start) / MB), ff(r.kbf, 1), r.at_join or "–",
                  ff(fnum(r.Mb_to_nearest_join), 1)] for r in rl.sort_values(["type", "copy"]).itertuples()])
        b45 = rl[rl.type == "45S"]
        if len(b45):
            unit, mx = 9.8, b45.kbf.max()
            R.p("%s: the largest 45S array in the assembly is %.1f kb, about %.0f units of ~%.1f kb. %s" % (
                SP_NAME[sp], mx, mx / unit, unit,
                "A nucleolus organiser holds hundreds of units, so the NORs themselves are not in the assembly: these "
                "are stray copies or the edges of unassembled arrays, and 45S FISH is needed to place the NORs."
                if mx < 200 else "Large enough to be (part of) a nucleolus organiser."))
    # 4f
    h2("s4f", "4f. FISH probe candidates")
    for sp, s in SPS:
        fp = RCD.get(sp, {}).get("fish_probes")
        if fp is None or not len(fp):
            continue
        fp = fp.copy()
        fp["score"] = pd.to_numeric(fp.fish_score, errors="coerce")
        R.table(["%s: family" % SP_NAME[sp], "kind", "monomer, bp", "total kb", "arrays", "largest, kb",
                 "share on its top chromosome", "copies"],
                [[r.family, r.kind, r.monomer_bp, ff(fnum(r.total_kb), 0), r.arrays, ff(fnum(r.largest_array_kb), 0),
                  "%s on %s" % (pc(fnum(r.share_on_top)), r.top_chromosome), r.copies]
                 for r in fp.sort_values("score", ascending=False).head(8).itertuples()])
    R.p("<span class='small'>Ranked by repeat_context.py's FISH score: the share of a family on one chromosome × the "
        "size of its largest array; a high score marks one chromosome brightly.</span>")
    # 4g
    h2("s4g", "4g. What the repeats do and do not explain")
    dd = []
    if vp.get("after") is not None:
        dd.append(["where D. paradoxa crossovers fall", vp["word"],
                   "repeats add %s after position + markers (shifted tracks: up to %s; %s)" % (
                       pc(vp["after"].obs, 1), pc(vp["after"].null_hi, 1), pf(vp["after"].p))])
    if DS and len(DS["rows"]):
        D = DS["rows"]
        dd.append(["the D. paradoxa crossover deserts",
                   "no (each within 10 points of its copy)" if (D.repeats - D.repeats_copy).abs().max() < 0.1 else
                   "not for all: see the desert table in 3d",
                   "repeats inside the deserts %s–%s against %s–%s for their copies; %d of %d sit on structure" % (
                       pc(D.repeats.min()), pc(D.repeats.max()), pc(D.repeats_copy.min()), pc(D.repeats_copy.max()),
                       DS["n_structure"], DS["n"])])
    if vb.get("alone") is not None:
        dd.append(["where D. binata crossovers fall", vb["word"],
                   "repeats alone explain %s (%s); after position + markers %s (%s)%s" % (
                       pc(vb["alone"].obs, 1), pf(vb["alone"].p),
                       pc(vb["after"].obs, 1) if vb.get("after") is not None else "–",
                       pf(vb["after"].p) if vb.get("after") is not None else "–",
                       "; near satellite arrays %.2f× expected (%s)" % (vb["sat"].obs, pf(vb["sat"].p))
                       if vb.get("sat") is not None and ok(vb["sat"].obs) else "")])
    if t is not None and "R_mb" in t:
        x = t[t.test.astype(str).str.startswith("breakpoints: with rDNA within") &
              (pd.to_numeric(t.R_mb, errors="coerce") == 1)]
        a_, h_ = x[x.set == "all joins"], x[x.set == "haplotype breakpoints"]
        if len(a_):
            sig = ok(a_.p.iloc[0]) and a_.p.iloc[0] < 0.05
            dd.append(["where the joins sit (A–D, X)", "maybe (check)" if sig else ("no (one hint)" if hits else "no"),
                       "rDNA within 1 Mb at %d of %d joins (random %.1f, %s)%s" % (
                           a_.obs.iloc[0], a_.n.iloc[0], a_.null_mean.iloc[0], pf(a_.p.iloc[0]),
                           "; 45S within 5 Mb at %s" % ", ".join(h["name"] for h in hits) if hits else "")])
        if len(h_):
            dd.append(["where the inversions between the haplotypes start and end",
                       "no" if not (ok(h_.p.iloc[0]) and h_.p.iloc[0] < 0.05) else "maybe (check)",
                       "rDNA within 1 Mb at %d of %d breakpoints (random %.1f, %s)" % (
                           h_.obs.iloc[0], h_.n.iloc[0], h_.null_mean.iloc[0], pf(h_.p.iloc[0]))])
    dd.append(["why chr5 and chr6 are coupled", "not tested as cause",
               "a desert keeps a middle intact but cannot link two chromosomes; the coupling needs a physical link or "
               "selection"])
    rl = RCD.get("paradoxa", {}).get("rdna_loci")
    if rl is not None and len(rl):
        mx = pd.to_numeric(rl[rl.type == "45S"].kb, errors="coerce").max()
        if ok(mx):
            dd.append(["where the NORs are", "cannot say" if mx < 200 else "partly",
                       "the largest 45S array assembled is %.1f kb%s" % (
                           mx, "; NORs are not assembled: FISH" if mx < 200 else "")])
    R.table(["question", "do the repeats explain it?", "evidence"], dd)


# ---- 5 open questions
def section_open():
    h2("s5", "5. Open questions and asks")
    slide("Three cheap checks settle most of it: meiosis I cytology (a cross of four? chains?), HiFi reads split at the "
          "chr5–chr6 point, and aligning P with P and Q with Q.")
    rows = [["chr1/chr2: tissue (A, D) or pollen (L1·Q, L2·P)?", "meiosis I spreads, chr1/chr2 arm paints",
             "a cross of four = assembly right; two ordinary pairs = pollen right"],
            ["chr5 × chr6: hidden swap or selection?", "HiFi reads split between %s ~%.0f Mb and %s ~%.0f Mb (known joins "
             "as positive controls); chr5/chr6 paints on somatic metaphases; pollen viability" % (
                 CS.get("p", "chr5"), CS.get("pt_p", float("nan")), CS.get("q", "chr6"), CS.get("pt_q", float("nan"))),
             "a junction or two mixed chromosomes = swap; neither = no physical link (aborted pollen = selection)"],
            ["inversions between the haplotypes: real?", "HiFi reads across each breakpoint; Hi-C",
             "reads spanning each breakpoint in one orientation per haplotype"],
            ["the end-of-P desert", "align the two P copies (%s) and the two Q copies (%s)" % (
                " vs ".join(arm_copies("P")) or "–", " vs ".join(arm_copies("Q")) or "–"),
             "an inversion between the P copies would explain it"],
            ["where are the NORs?", "45S FISH (and 5S)", "how many sites, on which chromosomes"]]
    if XS:
        s = XS
        rows.insert(2, ["chr3/chr4/X: how often adjacent?", "an X probe (oligo paint of the X arm, or an X-specific "
                        "repeat) on tetrads; pollen viability",
                        "tetrads 2-2-0-0 in ~%s (this data), ~%s if chains form only with an X crossover, 50%% if "
                        "independent; mature pollen: 1 signal in most, 2 in ~%s, 0 in %s" % (
                            pc(s["adjacent"]), pc(s["chain_only_adjacent"]), pc(s["share"]),
                            "none" if s["AR"] == 0 else pc(s["AR"] / max(s["N"], 1)))])
    R.table(["question", "test", "what each answer means"], rows)
    if XS:
        n = len(REFC)
        R.p("<b>Predicted meiosis I configurations</b> (n = %d chromosome pairs):" % n)
        R.table(["if", "meiosis I"], [
            ["the pollen are right about chr1/chr2 and chr5/chr6 is not physical",
             "1 chain of four (chr3/chr4, when X has a chiasma) + %d pairs" % (n - 2)],
            ["chr5/chr6 is also a swap", "2 multivalents of four + %d pairs" % (n - 4)],
            ["the assembly is right about chr1/chr2", "one more multivalent: a cross (ring) of four"]])
        R.p("The chr3/chr4 chain needs a chiasma in X, which the pollen see in ~%s of meioses; in the rest chr3 and chr4 "
            "may show as two pairs." % pc(XS["xco_meioses"]))


# ---- 6 limitations
def section_limits():
    h2("s6", "6. Limitations")
    items = ["Resolution: <i>D. paradoxa</i> crossovers are placed to a median %s kb and tested in %s Mb windows; "
             "nothing here is sharper than that." % (fi(v5n("paradoxa", "width_median_kb")),
                                                    ff(RRD.get("paradoxa", {}).get("win_mb"), 0)),
             "Sample sizes: %s <i>D. paradoxa</i> nuclei and %s crossovers; per-nucleus counts treat two nuclei of one "
             "pollen grain as two (near-identical pairs: %s; chance %s)." % (
                 fi(v5n("paradoxa", "cells_selected")), fi(v5n("paradoxa", "co_total")),
                 fi(v5n("paradoxa", "dup_pairs")), v5c("paradoxa", "dup_pairs_chance"))]
    if XS:
        items.append("X twice is probably undercounted: its mixed X windows lower a nucleus's haploidness, and the "
                     "selection removes nuclei below the haploidness cut. The adjacent share (%s) is therefore a lower "
                     "estimate." % pc(XS["adjacent"]))
    items += ["Hi-C at the joins reads 0.0× where the sequence is unmappable, so a missing Hi-C signal is not evidence "
              "against a join.",
              "The inversions and translocations between the haplotypes come from aligning hap2 onto hap1 "
              "(translocation_map.py); they are not yet confirmed by reads across the breakpoints. P and Q were never "
              "aligned to their homologs.",
              "Repeats: one annotation (CARP) on haplotype assemblies whose rDNA arrays are not fully assembled; in "
              "<i>D. binata</i> repeat effects and position go together."]
    nc, nv = RRD.get("paradoxa", {}).get("calls"), v5n("paradoxa", "co_total")
    if nc is not None and ok(nv):
        items.append("Consistency: the repeats analysis used %d <i>D. paradoxa</i> crossovers in %s nuclei, the "
                     "supervisor report %s in %s%s." % (nc, RRD["paradoxa"].get("nuclei"), fi(nv),
                                                        fi(v5n("paradoxa", "cells_selected")),
                                                        "" if abs(nc - nv) < 0.5 else
                                                        ": they differ, so the two were run on different call sets"))
    bullets(items)


# ---- appendix
def section_appendix():
    h2("app", "Appendix")
    R.h3("Per-chromosome landscapes (supervisor report v5)")
    emb(v5fig("fig_per_chromosome_landscapes.png"), "per_chromosome_landscapes", "Per-chromosome landscapes",
        "Crossovers along every chromosome, both species.")
    look("whether every chromosome shows its species' pattern (ends against middles), or some differ.")
    R.h3("The karyotype test and the co-inheritance map (cosegregation_2n.py)")
    emb("reports/cosegregation_2n/karyotype_test_%s.png" % P, "karyotype_test", "Every join tested",
        "Pollen φ at every join, the chr5 × chr6 coupling and its control.")
    look("black points near +1 are joins the pollen inherit; red points near 0 are joins they do not.")
    emb("reports/cosegregation_2n/cosegregation_2n_%s.png" % P, "co_inheritance_map", "Co-inheritance map",
        "Which copies the pollen inherit together, beside the tissue Hi-C.")
    look("red squares off the diagonal are pieces on different copies inherited together.")
    R.h3("Tissue Hi-C (supervisor report v5)")
    emb(v5fig("fig_hic_both_haplotypes_paradoxa.png"), "hic_both_haplotypes", "Tissue Hi-C, both haplotypes",
        "The two-haplotype Hi-C map.")
    look("contacts between copies of different chromosomes are the translocation signal in the tissue.")
    R.h3("All variables against crossovers (repeat_recombination.py)")
    for sp, s in SPS:
        a = RRD.get(sp, {}).get("assoc")
        if a is None:
            continue
        R.table(["%s: variable" % SP_NAME[sp], "group", "windows", "ρ", "p (shifted)", "ρ after end + markers", "p",
                 "same sign", "beyond chance"],
                [[r.variable, r.group, r.n_windows, sg(r.rho), pf(r.p_rot), sg(r.rho_partial) if ok(r.rho_partial)
                  else "–", pf(r.p_partial) if ok(r.p_partial) else "–", r.same_sign,
                  "yes" if str(r.beyond_chance) == "True" else "no"] for r in a.itertuples()])
    if JOINS or COUPLE:
        R.h3("The karyotype test, as written by cosegregation_2n.py")
        R.details("karyotype_test_%s.tsv" % P, path=KT_TSV)
    R.h3("Evidence files")
    for sp, s in SPS:
        if os.path.exists(rr_path("summary", s)):
            R.details("repeat_recombination summary, %s" % SP_NAME[sp], path=rr_path("summary", s))
    if RR5:
        R.details("repeat_recombination summary, near5 run", path=rr_path("summary", P, True))
    for lg in (COS_LOG, "logs/x_arm_dosage_%s.log" % P, "logs/repeat_context_%s.log" % P):
        if os.path.exists(lg):
            R.details(lg, path=lg)
    R.h3("Figures for slides, in order (fig/)")
    R.table(["file", "title", "message"], [[a, b, c] for a, b, c, _ in SLIDES])


for name, fn in (("0 summary", section_summary), ("1 material", section_material), ("2 landscapes", section_landscapes),
                 ("3 karyotype", section_karyotype), ("4 repeats", section_repeats), ("5 open", section_open),
                 ("6 limitations", section_limits), ("appendix", section_appendix)):
    R.guarded(name, fn)

R.parts.insert(1, "<p class='small'>%s</p><p class='small'>For slides: every figure is also in fig/ as its own PNG, "
                  "numbered in report order; slides_outline.md lists the section messages, the differences from the "
                  "assembly and the figures with what to look for in each.</p>" % " · ".join(
                      "<a href='#%s'>%s</a>" % (i, H(t)) for i, t in TOC))
with open(os.path.join(A.out, "slides_outline.md"), "w") as fo:
    fo.write("# Slide outline (final_report.py)\n\nBuilt with the report; the figures are in fig/, numbered in report "
             "order.\n\n## Section messages\n\n")
    for t, m in SECTION_MSGS:
        fo.write("- **%s**: %s\n" % (strip_tags(t), m))
    dr = R.guarded("outline table", diff_rows) or []
    if dr:
        fo.write("\n## Differences from the assembly\n\n| %s |\n|%s\n" % (" | ".join(DIFF_HEAD),
                                                                       "---|" * len(DIFF_HEAD)))
        for r in dr:
            fo.write("| %s |\n" % " | ".join(str(x).replace("|", "\\|") for x in r))
    fo.write("\n## Figures\n\n")
    for i, (f, t, m, lk) in enumerate(SLIDES, 1):
        fo.write("### %d. %s\n\n- figure: fig/%s\n- message: %s\n%s\n" % (
            i, t, f, m, "- what to look for: %s\n" % lk if lk else ""))
R.put("figures", len(SLIDES))
for f_ in FIGS.values():
    plt.close(f_)
R.write()
