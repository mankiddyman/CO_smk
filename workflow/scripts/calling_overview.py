#!/usr/bin/env python3
"""calling_overview.py -- from barcodes to crossover calls, with the numbers, read from the pipeline's own files.

For each sample: the markers cells are genotyped at, every step that removes cells on the way to crossover
calling (with its threshold and how many cells survive it), the hapCO settings, and what came out. Values
'as run' are read from the outputs (selection_summary.txt, co_calling_params.txt); values 'configured' from
config/config.yaml. A mismatch means the outputs are older than the config, and is flagged.

Also, per chromosome (or piece) and cell, from the molecule tables hapCO reads: is it there once (haploid),
twice (both alleles interleaved at about double the expected molecules), not at all, or mixed? A nucleus
that inherited a translocated piece twice reads like that, whatever the reference.
  expected   the cell's molecules x the median share of this sequence across cells
  switches   allele changes between neighbouring molecules / (molecules - 1): ~0 for one haplotype
             (a crossover adds one), ~0.5 when two copies are interleaved
  classes    too few (expected < 8) | missing (< 0.25 x expected) | doubled (switches >= 0.25 and
             >= 1.2 x expected) | mixed (switches >= 0.25 at normal coverage) | haploid

Usage: calling_overview.py [SAMPLE ...]   (default: Dbinata_hap1 Dparadoxa_split Dparadoxa_CB Dparadoxa_std
                                           Dparadoxa_hap1); run from the CO_smk directory
"""
import collections
import csv
import datetime
import glob
import gzip
import os
import re
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml

SAMPLES = sys.argv[1:] or ["Dbinata_hap1", "Dparadoxa_split", "Dparadoxa_CB", "Dparadoxa_std", "Dparadoxa_hap1"]
DETAIL = {"Dbinata_hap1", "Dparadoxa_split"}           # the others get the short comparison only
CFG = yaml.safe_load(open("config/config.yaml"))
ROWS = {r["sample_id"]: r for r in csv.DictReader(open("config/samples.csv"))}
MISSING = []


def say(s=""):
    print(s, flush=True)


def have(p):
    if os.path.exists(p) and os.path.getsize(p) > 0:
        return True
    MISSING.append(p)
    return False


def fi(x):
    return format(int(x), ",") if x == x and x is not None else "-"


def cfg(section, species, key):
    return ((CFG.get(section) or {}).get(species) or {}).get(key)


def kv(p):
    d = {}
    for l in open(p):
        if ":" in l:
            k, v = l.split(":", 1)
            d[k.strip()] = v.strip()
    return d


def git_commit():
    r = subprocess.run(["git", "rev-parse", "--short", "HEAD"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                       universal_newlines=True)
    return r.stdout.strip() or "?"


def num(s):
    m = re.search(r"-?[\d.]+(?:e[+-]?\d+)?", str(s))
    return float(m.group()) if m else float("nan")


def same(a, b):
    try:
        return abs(float(a) - float(b)) < 1e-9
    except (TypeError, ValueError):
        return str(a).strip().lower() == str(b).strip().lower()


def check(label, configured, as_run):
    if configured is None or as_run is None:
        return ""
    return "" if same(configured, as_run) else "   !! configured %s, outputs ran with %s" % (configured, as_run)


def chrom_lengths(S):
    p = "results/cell_data/%s/chrom_map.tsv" % S
    if not have(p):
        return collections.OrderedDict()
    t = pd.read_csv(p, sep="\t")
    return collections.OrderedDict(zip(t["name"].astype(str), t["size"].astype(int)))


def dosage(S, cells, L):
    """Per cell and sequence: molecules, switch rate, class (see the module docstring)."""
    rows = []
    for bc in cells:
        p = "results/cell_data_mol/%s/%s.tsv" % (S, bc)
        if not os.path.exists(p):
            continue
        m = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 3, 5], names=["chrom", "pos", "rc", "ac"],
                        dtype={0: str})
        m = m[(m.rc != m.ac) & m.chrom.isin(L)].sort_values(["chrom", "pos"])
        m["alt"] = (m.ac > m.rc).astype(int)
        tot = len(m)
        g = {c: x.alt.values for c, x in m.groupby("chrom")}
        for c in L:
            a = g.get(c, np.array([], dtype=int))
            n = len(a)
            sw = float((np.diff(a) != 0).sum()) / (n - 1) if n > 1 else float("nan")
            rows.append((bc, c, n, tot, sw))
    t = pd.DataFrame(rows, columns=["barcode", "chrom", "n", "total", "switch"])
    if t.empty:
        return t
    share = (t.n / t.total.clip(lower=1)).groupby(t.chrom).median()
    t["expected"] = t.total * t.chrom.map(share)
    t["coverage"] = t.n / t.expected.replace(0, np.nan)

    def cls(r):
        if r.expected < 8:
            return "too few"
        if r.coverage < 0.25:
            return "missing"
        if r.switch >= 0.25:
            return "doubled" if r.coverage >= 1.2 else "mixed"
        return "haploid"
    t["class"] = t.apply(cls, axis=1)
    return t


def step(label, value, note=""):
    say("  %-82s %10s%s" % (label, value, ("   " + note) if note else ""))


def overview(S):
    row = ROWS.get(S)
    if row is None:
        say("== %s: not in config/samples.csv" % S)
        return None
    sp = row["species"]
    L = chrom_lengths(S)
    out = {"S": S}
    detail = S in DETAIL
    say("=" * 110)
    say("%s   (%s; reference %s, %d sequences, %.0f Mb)" % (S, sp, row["assembly_fasta"], len(L), sum(L.values()) / 1e6))
    say("=" * 110)

    # ---------------------------------------------------------------- markers
    if detail:
        say("MARKERS  (heterozygous SNPs found in the plant's HiFi reads; every cell is genotyped at these)")
        p = "results/markers/%s/filter_summary.txt" % S
        if have(p):
            t = open(p).read()
            raw = re.search(r"Raw variants:\s*(\d+)", t)
            flt = re.search(r"Markers after filter:\s*(\d+)", t)
            fp = re.search(r"Filter params:\s*(.*)", t)
            say("  raw HiFi variant calls                         %12s" % fi(int(raw.group(1)) if raw else None))
            say("  after the marker filter                        %12s   (%s)" % (
                fi(int(flt.group(1)) if flt else None), fp.group(1).strip() if fp else "?"))
        p = "results/blacklist/%s/rrna_exclude.bed" % S
        if have(p):
            b = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2], comment="#")
            say("  rDNA regions left out                          %12s intervals, %.2f Mb" % (fi(len(b)),
                                                                                             (b[2] - b[1]).sum() / 1e6))
        p = "qc/markers/%s/marker_classes.tsv.gz" % S
        if have(p):
            mc = pd.read_csv(p, sep="\t", usecols=["class"])
            vc = mc["class"].value_counts()
            say("  markers seen in enough cells to be judged      %12s" % fi(len(mc)))
            say("    " + "  ".join("%s %s (%.0f%%)" % (k, fi(v), 100.0 * v / len(mc)) for k, v in vc.items()))
            say("    good = reads ALT in ~half the cells and one allele per cell; sticky = one allele in nearly every")
            say("    cell (carries no haplotype); paralog-like = both alleles inside single cells. Only GOOD markers")
            say("    feed haploidness and crossover calling.")
            out["good"] = int(vc.get("good", 0))
        say()

    # ---------------------------------------------------------------- cells
    say("CELLS  (each step keeps part of the one before; thresholds as run)")
    n_bs = n_called = None
    p = "results/cells/%s/barcode_stats.tsv" % S
    if have(p):
        with open(p) as fh:
            n_bs = sum(1 for _ in fh) - 1
    p = "results/cells/%s/barcodes_called.tsv" % S
    if have(p):
        n_called = sum(1 for l in open(p) if l.strip() and not l.lower().startswith("barcode"))
    lower, fdr, ex = (cfg("cell_calling", sp, k) for k in ("default_lower", "fdr", "exclude_rrna"))
    step("1 barcodes in STARsolo's raw matrix", fi(n_bs))
    step("2 called as cells (emptyDrops: ambient < %s UMIs, FDR %s, rDNA genes %s)" % (
        lower, fdr, "left out" if str(ex).lower() == "true" else "kept"), fi(n_called))
    sel = kv("results/cell_qc/%s/selection_summary.txt" % S) if have("results/cell_qc/%s/selection_summary.txt" % S) else {}
    h, w = num(sel.get("min_haploidness")), num(sel.get("min_windows"))
    mm, mcm = num(sel.get("min_molecules", 0)), num(sel.get("min_chrom_molecules", 0))
    ht_p = "qc/haplotypes/%s/haplotype_tracks.tsv.gz" % S
    ht = pd.read_csv(ht_p, sep="\t") if have(ht_p) else None
    good = [l.strip() for l in open("results/cell_qc/%s/good_cells.tsv" % S) if l.strip()] \
        if have("results/cell_qc/%s/good_cells.tsv" % S) else []
    out["cells"] = good
    if ht is not None and sel:
        s = ht.dropna(subset=["haploidness"])
        s1 = s[s.haploidness >= h]
        s2 = s1[s1.windows >= w]
        s3 = s2[s2.molecules >= mm]
        s4 = s3[s3.min_chrom_molecules >= mcm] if "min_chrom_molecules" in s3 else s3
        step("3 genotyped by cellsnp at the markers", fi(len(ht)))
        step("4 scored: >= 5 full windows of 15 molecules", fi(len(s)))
        step("5 haploidness >= %s" % sel.get("min_haploidness"), fi(len(s1)))
        say("      haploidness = mean over 15-molecule windows of |2f - 1| (f = ALT share): 1 = one clean")
        say("      haploid nucleus; ~0.2 = both haplotypes; a nucleus with a share a of foreign molecules ~ 1 - a")
        step("6   and >= %s full windows" % sel.get("min_windows"), fi(len(s2)))
        step("7   and >= %s molecules in the cell" % sel.get("min_molecules", "0"), fi(len(s3)))
        step("8   and >= %s molecules on its weakest chromosome%s" % (
            sel.get("min_chrom_molecules", "0"), " (pieces counted together)" if "_split" in S else ""), fi(len(s4)),
            "= good_cells.tsv" if len(s4) == len(good) else "!! good_cells.tsv has %s" % fi(len(good)))
        one = {"haploidness": (s.haploidness < h), "windows": (s.windows < w), "molecules": (s.molecules < mm),
               "weakest chromosome": (s.min_chrom_molecules < mcm) if "min_chrom_molecules" in s else (s.molecules < 0)}
        fails = pd.DataFrame(one)
        only = {k: int((fails[k] & (fails.sum(axis=1) == 1)).sum()) for k in one}
        say("    scored cells failing ONLY this one: " + ", ".join("%s %s" % (k, fi(v)) for k, v in only.items()))
        for k, cf in (("min_haploidness", "min_haploidness"), ("min_windows", "min_windows"),
                      ("min_molecules", "min_molecules"), ("min_chrom_molecules", "min_chrom_molecules")):
            msg = check(k, cfg("cell_qc", sp, cf), sel.get(k))
            if msg:
                say("    %s:%s" % (k, msg))
        g = ht[ht.barcode.isin(good)]
        if len(g):
            q = g.molecules.quantile([0.25, 0.5, 0.75])
            say("    selected cells: molecules median %s (IQR %s-%s), weakest chromosome median %s, haploidness median %.2f"
                % (fi(q[0.5]), fi(q[0.25]), fi(q[0.75]), fi(g.min_chrom_molecules.median()) if "min_chrom_molecules" in g
                   else "-", g.haploidness.median()))
            out["mol_median"] = float(q[0.5])
    log = "logs/co_calling/%s.log" % S
    skipped = sum(1 for l in open(log) if "Info: skipping" in l) if have(log) else None
    cm = cfg("co_calling", sp, "cell_markers")
    step("9 hapCO's own floor: >= %s informative rows, else no calls -- cells skipped" % cm, fi(skipped))
    say()

    # ---------------------------------------------------------------- settings
    if detail:
        say("CALLING SETTINGS  (hapCO_identification.R, one cell at a time; as run, from co_calling_params.txt)")
        p = "results/crossovers/%s/co_calling_params.txt" % S
        run = {}
        if have(p):
            for l in open(p):
                m = re.match(r"\s*([A-Za-z_ ]+?)(?:\s*\(-\w\))?:\s*(.*)", l)
                if m:
                    run[m.group(1)] = m.group(2).strip()
        expl = [("input tables", "input_rows", "molecules = good markers only, markers < 150 bp apart merged into one row"),
                ("baseAF", "base_af", "smoothing: a molecule's call is overruled when its 4 neighbours mostly disagree"),
                ("marker_num", "marker_num", "an interior block needs at least this many rows (molecules)"),
                ("block_size", "block_size", "an interior block must span at least this many bp"),
                ("terminal_marker_num", "terminal_marker_num", "a chromosome's first and last block need this many rows, "
                                                                "and no minimum length (0 = off)"),
                ("windowAF", "window_af", "inert here (changes nothing)"),
                ("genotype", "genotype", "inert here (changes nothing)"),
                ("cell_markers", "cell_markers", "a floor only (step 9)")]
        for rk, ck, text in expl:
            v = run.get(rk, "-")
            c = cfg("co_calling", sp, ck)
            if rk == "input tables":
                v_show = "molecules" if "cell_data_mol" in v else ("markers" if v != "-" else "-")
                msg = check(rk, c, v_show) if v != "-" else ""
            else:
                v_show = v
                msg = check(rk, c, v) if v != "-" else ""
            say("  %-20s %-12s %s%s" % (rk, v_show, text, msg))
        say()

    # ---------------------------------------------------------------- result
    say("RESULT  (results/crossovers/%s)" % S)
    p = "results/crossovers/%s/co_per_cell.tsv" % S
    per = pd.read_csv(p, sep="\t") if have(p) else None
    p = "results/crossovers/%s/co_intervals.bed" % S
    co = pd.read_csv(p, sep="\t", header=None, usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "barcode"],
                     dtype={0: str}) if have(p) else None
    if per is not None:
        n = per.n_cos
        out.update(n_cells=len(per), n_cos=int(n.sum()), mean=float(n.mean()))
        say("  cells %s, crossovers %s, per cell mean %.2f, median %.0f; cells with none %s" % (
            fi(len(per)), fi(n.sum()), n.mean(), n.median(), fi((n == 0).sum())))
        vc = n.value_counts().sort_index()
        cap = 15
        parts = ["%d:%d" % (k, v) for k, v in vc.items() if k < cap]
        if (n >= cap).any():
            parts.append("%d+:%d" % (cap, (n >= cap).sum()))
        say("  crossovers per cell -> cells   " + "  ".join(parts))
    if co is not None and per is not None and len(L):
        cpc = co.groupby("chrom").size().reindex(list(L), fill_value=0) / max(len(per), 1)
        say("  per %s, crossovers per cell (one per bivalent = 0.5 in a gamete):" % ("piece" if "_split" in S else "chromosome"))
        items = ["%s %.2f" % (c, v) for c, v in cpc.items()]
        for i in range(0, len(items), 6):
            say("    " + "   ".join(items[i:i + 6]))
        low = [c for c, v in cpc.items() if v < 0.5]
        say("    below 0.5: %s" % (", ".join(low) if low else "none"))
        bp = (co.end - co.start)
        say("  crossover interval (between the last molecule of one block and the first of the next): median %.2f Mb, "
            "90%% below %.2f Mb" % (bp.median() / 1e6, bp.quantile(0.9) / 1e6))
    pdfs = glob.glob("results/crossovers/%s/per_cell/*_co.pdf" % S)
    say("  per-cell plots: %s PDFs in results/crossovers/%s/per_cell/<barcode>_co.pdf" % (fi(len(pdfs)), S))
    say()

    # ---------------------------------------------------------------- dosage
    if detail and good and len(L):
        t = dosage(S, good, L)
        if not t.empty:
            unit = "piece" if "_split" in S else "chromosome"
            say("DOSAGE PER %s AND CELL  (selected cells; molecule tables; classes in the script's docstring)" % unit.upper())
            tab = t.pivot_table(index="chrom", columns="class", values="barcode", aggfunc="count", fill_value=0)
            cols = [c for c in ("haploid", "doubled", "missing", "mixed", "too few") if c in tab.columns]
            tab = tab.reindex(list(L)).fillna(0).astype(int)
            flagged = [c for c, r in tab.iterrows() if sum(r.get(k, 0) for k in cols if k != "haploid") > 0]
            if flagged:
                say("  %-16s %s" % (unit, "".join("%10s" % c for c in cols)))
                for c, r in tab.iterrows():
                    say("  %-16s %s" % (c, "".join("%10s" % fi(r.get(k, 0)) for k in cols)))
            else:
                say("  every selected cell is haploid on every %s" % unit)
            bad = t[t["class"].isin(["doubled", "missing"])]
            say("  cells with at least one %s doubled or missing: %s of %s" % (unit, fi(bad.barcode.nunique()),
                                                                            fi(t.barcode.nunique())))
            if not bad.empty:
                say("  what the flagged cells carry (most common combinations):")
                combo = collections.Counter(" + ".join("%s %s" % (c, k) for c, k in zip(x.chrom, x["class"]))
                                            for _, x in bad.groupby("barcode"))
                for k, v in combo.most_common(8):
                    say("    %4d cells   %s" % (v, k))
            os.makedirs("qc/calling_overview", exist_ok=True)
            t.to_csv("qc/calling_overview/%s_dosage.tsv" % S, sep="\t", index=False)
            say("  per cell and %s: qc/calling_overview/%s_dosage.tsv" % (unit, S))
            say()

    # ---------------------------------------------------------------- earlier tuning
    if detail:
        say("EARLIER TUNING  (what is on disk from previous rounds)")
        n_found = 0
        for lab, pat in (("review panels", "qc/review/%s/panel_seed*.txt" % S),
                         ("panel re-calls", "qc/review/%s/runs/*" % S),
                         ("block-support checks", "qc/review/%s/support*" % S),
                         ("whole-sample sweeps", "results/co_sweep/%s/*" % S),
                         ("block-size benchmark", "qc/benchmark/%s/*" % S)):
            hits = sorted(glob.glob(pat))
            if hits:
                n_found += 1
                say("  %-22s %s: %s" % (lab, os.path.dirname(pat), ", ".join(os.path.basename(h) for h in hits[:12]) +
                                        (" ... (%d)" % len(hits) if len(hits) > 12 else "")))
        if not n_found:
            say("  none found")
        say()
    return out


say("CROSSOVER CALLING OVERVIEW   %s   repo %s" % (datetime.datetime.now().isoformat(timespec="minutes"), git_commit()))
say("Every number is read from the pipeline's files; '!!' marks a setting whose outputs are older than the config.")
say()
res = [overview(S) for S in SAMPLES]
res = [r for r in res if r]
par = [r for r in res if r["S"].startswith("Dparadoxa")]
if len(par) > 1:
    say("=" * 110)
    say("PARADOXA REFERENCES COMPARED  (same reads, same settings; only the mapping reference differs)")
    say("  %-17s %8s %8s %10s %18s" % ("sample", "cells", "COs", "COs/cell", "median molecules"))
    for r in par:
        say("  %-17s %8s %8s %10s %18s" % (r["S"], fi(len(r.get("cells", []))), fi(r.get("n_cos")),
                                          "%.2f" % r["mean"] if "mean" in r else "-", fi(r.get("mol_median"))))
    sets = [set(r.get("cells", [])) for r in par if r.get("cells")]
    if len(sets) > 1:
        say("  cells selected in all of them: %s" % fi(len(set.intersection(*sets))))
if MISSING:
    say()
    say("files not found (%d): %s" % (len(MISSING), ", ".join(sorted(set(MISSING))[:20])))
