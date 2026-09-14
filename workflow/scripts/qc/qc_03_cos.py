#!/usr/bin/env python3
"""QC stage 03 -- crossover calls and the recombination landscape.

The central question: did we detect most of the real crossovers, or are we
marker-limited? On Spondias detection had NOT saturated at the achieved density
(median 132 markers/chromosome vs ~148 for 90%), which makes the reported map
length a LOWER BOUND. That is the single most important caveat on the result.

Two constraints do the heavy lifting here:

1. The obligate-CO floor. One crossover per bivalent makes 2 of 4 chromatids
   recombinant, so >=0.5 CO/gamete => >=50 cM per chromosome. This is a hard
   physical minimum, and during the Spondias audit it selected more parameters
   than any statistical criterion -- it bounded block_size from above, excluded
   base_af 0.2, and chose max_switch_rate 0.12.

2. Zero-CO chromosomes are EXPECTED, not pathological: with one CO per bivalent
   ~50% of gametes carry none for that chromosome. An earlier draft of this
   analysis called that "meiotically impossible", which was wrong.

The saturating fit (detected = true * (1 - exp(-m/k))) is an ASSUMED functional
form. On Spondias the top marker bin observed 1.03 COs/chromosome while the fit
asymptote was 0.85 -- data above the model's own ceiling means the model is
wrong and the efficiency estimate is optimistic. Both are printed so the
discrepancy stays visible.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_common import (Summary, Flags, read_fai, chrom_sort_key, save,
                       C_RAW, C_KEEP, C_REF, C_CUT)

FLOOR_CM = 50.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--co_intervals", required=True)
    ap.add_argument("--co_summary", required=True)
    ap.add_argument("--landscape_summary", required=True)
    ap.add_argument("--switches", required=True)
    ap.add_argument("--good_cells", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--block_size", type=int, required=True)
    ap.add_argument("--marker_num", type=int, required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    S = Summary(a.sample, "03_cos")
    F = Flags(a.sample, "03_cos")
    S.add("param.block_size", a.block_size)
    S.add("param.marker_num", a.marker_num)

    lens = read_fai(a.fai)

    # ---- CO intervals ------------------------------------------------------
    co_by_pair = defaultdict(int)          # (cell, chrom) -> n COs
    widths, mids = [], []                  # interval widths, (chrom, midpoint)
    for line in open(a.co_intervals):
        f = line.split()
        if len(f) < 4:
            continue
        c, s0, e0, cell = f[0], int(f[1]), int(f[2]), f[3]
        co_by_pair[(cell, c)] += 1
        widths.append(e0 - s0)
        mids.append((c, (s0 + e0) / 2))
    widths = np.array(widths)
    n_cos = len(widths)

    cells = [l.split()[0] for l in open(a.good_cells)]
    n_cells = len(cells)
    S.add("cells", n_cells)
    S.add("total_cos", n_cos)
    S.add("mean_cos_per_cell", round(n_cos / n_cells, 3) if n_cells else 0)

    # ---- per-chromosome landscape -----------------------------------------
    per_chrom = []
    total_cm = 0.0
    with open(a.landscape_summary) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for line in f:
            g = line.rstrip("\n").split("\t")
            c = g[ix["chrom"]]
            cm = float(g[ix["total_genetic_cM"]])
            per_chrom.append((c, float(g[ix["size_Mb"]]), int(g[ix["n_cos"]]),
                              cm, float(g[ix["cM_per_Mb"]])))
            total_cm += cm
    n_chrom = len(per_chrom)
    below = [p for p in per_chrom if p[3] < FLOOR_CM]
    S.add("chromosomes", n_chrom)
    S.add("total_map_cM", round(total_cm, 1))
    S.add("mean_cM_per_Mb", round(total_cm / sum(p[1] for p in per_chrom), 3))
    S.add("chromosomes_below_floor", len(below),
          f"<{FLOOR_CM:.0f} cM; obligate CO makes 2 of 4 chromatids "
          f"recombinant so this is a hard physical minimum")
    cos_per_bivalent = 2 * n_cos / n_cells / n_chrom if n_cells else 0
    S.add("cos_per_bivalent", round(cos_per_bivalent, 3),
          "2 x mean COs per gamete per chromosome; plants typically 1-3")

    blind = 2 * a.block_size / (1e6 * np.mean([p[1] for p in per_chrom]))
    S.add("co_blind_fraction_pct", round(100 * blind, 1),
          "terminal block_size at each arm cannot be called: a CO needs a "
          "supported block on BOTH sides")
    S.add("co_interval_median_kb", round(float(np.median(widths)) / 1000, 1))
    S.add("co_interval_q3_kb", round(float(np.percentile(widths, 75)) / 1000, 1))
    S.add("co_interval_p90_kb", round(float(np.percentile(widths, 90)) / 1000, 1))

    # ---- detection saturation ---------------------------------------------
    hdr = open(a.switches).readline().rstrip("\n").split("\t")
    mcol = {h.split("_")[1]: i for i, h in enumerate(hdr) if h.startswith("markers_")}
    good = set(cells)
    mk = {}
    for line in open(a.switches):
        f = line.rstrip("\n").split("\t")
        if f[0] not in good:
            continue
        for c, i in mcol.items():
            mk[(f[0], c)] = int(f[i])

    M = np.array([v for v in mk.values()])
    C = np.array([co_by_pair.get(k, 0) for k in mk])
    S.add("cell_chromosome_pairs", len(M))
    S.add("markers_per_chrom_median", int(np.median(M)))
    S.add("markers_per_chrom_q1", int(np.percentile(M, 25)))
    S.add("expected_markers_per_block",
          round(float(np.median(M)) * a.block_size /
                (1e6 * np.mean([p[1] for p in per_chrom])), 1),
          f"a typical cell's markers in one block, vs marker_num={a.marker_num}")

    edges = [30, 50, 75, 100, 150, 200, 300, 500, 10 ** 9]
    binx, biny, binn = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (M >= lo) & (M < hi)
        if s.sum() < 10:
            continue
        binx.append(float(np.median(M[s])))
        biny.append(float(C[s].mean()))
        binn.append(int(s.sum()))
    binx, biny, binn = np.array(binx), np.array(biny), np.array(binn)

    asym = k_half = eff = float("nan")
    if len(binx) >= 4:
        try:
            from scipy.optimize import curve_fit
            def sat(m, t, k):
                return t * (1 - np.exp(-m / k))
            p, _ = curve_fit(sat, binx, biny, p0=[max(biny) * 1.2, 150],
                             sigma=1 / np.sqrt(binn), maxfev=20000)
            asym, k_half = float(p[0]), float(p[1])
            eff = 100 * float(C.mean()) / asym
        except Exception as e:
            print(f"  saturation fit failed: {e}", file=sys.stderr)

    top_obs = float(biny[-1]) if len(biny) else float("nan")
    S.add("saturation_fit_asymptote", round(asym, 3))
    S.add("saturation_half_marker_count", round(k_half, 0))
    S.add("detection_efficiency_pct", round(eff, 1),
          "observed mean / fitted asymptote; MODEL-DEPENDENT")
    S.add("top_bin_observed_cos", round(top_obs, 3),
          "if this EXCEEDS the fitted asymptote, the model is wrong and the "
          "efficiency above is optimistic")
    fit_suspect = top_obs == top_obs and asym == asym and top_obs > asym
    S.add("fit_exceeded_by_data", str(fit_suspect))
    if asym == asym:
        need90 = -k_half * np.log(0.1)
        S.add("markers_needed_for_90pct", round(need90, 0))
        not_sat = float(np.median(M)) < need90
    else:
        not_sat = False

    # within-cell test: controls for cell quality and doublets entirely
    lo_r, hi_r = [], []
    for cell in good:
        d = [(mk[(cell, c)] / lens[c] * 1e6, co_by_pair.get((cell, c), 0) / lens[c] * 1e6)
             for c in mcol if (cell, c) in mk and c in lens]
        if len(d) < 8:
            continue
        d.sort()
        h = len(d) // 2
        lo_r.append(np.mean([x[1] for x in d[:h]]))
        hi_r.append(np.mean([x[1] for x in d[h:]]))
    if lo_r:
        diff = np.array(hi_r) - np.array(lo_r)
        rng = np.random.default_rng(0)
        bs = [rng.choice(diff, len(diff), replace=True).mean() for _ in range(2000)]
        ci = np.percentile(bs, [2.5, 97.5])
        S.add("within_cell_lo_density_co_per_mb", round(float(np.mean(lo_r)), 5))
        S.add("within_cell_hi_density_co_per_mb", round(float(np.mean(hi_r)), 5))
        S.add("within_cell_diff_ci", f"[{ci[0]:+.5f}, {ci[1]:+.5f}]",
              "excludes 0 => detection genuinely improves with markers, "
              "controlling for cell quality")

    # ---- flags -------------------------------------------------------------
    F.check(len(below) > 0, "FLOOR_VIOLATION",
            f"{len(below)} of {n_chrom} chromosomes below {FLOOR_CM:.0f} cM "
            f"({', '.join(f'chr{p[0]}={p[3]:.0f}' for p in below[:5])}). Either "
            f"parameters are under-calling or the chromosome count is wrong.")
    F.check(not_sat, "NOT_SATURATED",
            f"detection has not saturated: median {np.median(M):.0f} "
            f"markers/chromosome vs ~{-k_half*np.log(0.1):.0f} needed for 90%. "
            f"The reported {total_cm:.0f} cM is a LOWER BOUND.")
    F.check(not (0.5 <= cos_per_bivalent <= 4), "IMPLAUSIBLE_CO_RATE",
            f"{cos_per_bivalent:.2f} COs per bivalent is outside the normal "
            f"plant range (1-3).")
    F.check(100 * blind > 15, "BLIND_ZONE_LARGE",
            f"{100*blind:.0f}% of the genome is CO-blind at block_size="
            f"{a.block_size:,}. Consider a smaller block if marker density allows.")
    F.check(float(np.median(widths)) > a.block_size, "WIDE_INTERVALS",
            f"median CO interval {np.median(widths)/1000:.0f} kb exceeds "
            f"block_size; localisation is poor.")
    F.check(fit_suspect, "SATURATION_FIT_SUSPECT",
            f"top marker bin observes {top_obs:.2f} COs/chromosome, above the "
            f"fitted asymptote {asym:.2f}. The exponential form is wrong here; "
            f"treat detection_efficiency_pct as optimistic.")

    # ---- plots -------------------------------------------------------------
    print("\nplotting...", file=sys.stderr)

    fig, ax = plt.subplots(2, 2, figsize=(14, 9))

    A = ax[0][0]
    if len(binx):
        se = [C[(M >= lo) & (M < hi)].std() / np.sqrt(s)
              for (lo, hi), s in zip(zip(edges[:-1], edges[1:]), binn)
              if ((M >= lo) & (M < hi)).sum() >= 10]
        A.errorbar(binx, biny, yerr=se[:len(binx)], marker="o", color=C_RAW,
                   label="observed")
    if asym == asym:
        xx = np.logspace(np.log10(max(binx.min(), 10)), np.log10(binx.max() * 3), 200)
        A.plot(xx, asym * (1 - np.exp(-xx / k_half)), "k--", lw=1,
               label=f"fit: asymptote {asym:.2f}, k={k_half:.0f}")
        A.axhline(asym, color="gray", ls=":")
    A.axvline(np.median(M), color=C_KEEP, ls=":",
              label=f"median {np.median(M):.0f}")
    A.set_xscale("log")
    A.set(xlabel="markers on chromosome (per cell)", ylabel="mean COs detected",
          title="A. Detection saturation\nstill rising = marker-limited")
    A.legend(fontsize=8); A.grid(alpha=.3)

    B = ax[0][1]
    if asym == asym:
        xx = np.logspace(1, np.log10(max(M.max(), 1000)), 200)
        B.plot(xx, 100 * (1 - np.exp(-xx / k_half)), color=C_KEEP)
        for tgt in (90, 95):
            need = -k_half * np.log(1 - tgt / 100)
            B.axhline(tgt, color="gray", ls="--", lw=.7)
            B.annotate(f"{tgt}% needs {need:.0f}", (need, tgt), fontsize=8,
                       xytext=(4, -12), textcoords="offset points")
        B.axvline(np.median(M), color=C_KEEP, ls=":",
                  label=f"median {np.median(M):.0f} -> "
                        f"{100*(1-np.exp(-np.median(M)/k_half)):.0f}%")
        B.legend(fontsize=8)
    B.set_xscale("log"); B.set_ylim(0, 105)
    B.set(xlabel="markers on chromosome", ylabel="detection efficiency (%)",
          title="B. What fraction of COs do we see?")
    B.grid(alpha=.3)

    Cx = ax[1][0]
    dens = np.array([mk[k] / lens[k[1]] * 1e6 for k in mk if k[1] in lens])
    rate = np.array([co_by_pair.get(k, 0) / lens[k[1]] * 1e6 for k in mk if k[1] in lens])
    chs = np.array([k[1] for k in mk if k[1] in lens])
    band = (dens >= np.percentile(dens, 25)) & (dens <= np.percentile(dens, 75))
    rng = np.random.default_rng(0)
    rows = []
    for c in sorted(set(chs), key=chrom_sort_key):
        s = band & (chs == c)
        if s.sum() < 20:
            continue
        r = rate[s]
        bs = [rng.choice(r, len(r), replace=True).mean() for _ in range(1000)]
        rows.append((c, r.mean(), *np.percentile(bs, [2.5, 97.5])))
    rows.sort(key=lambda r: -r[1])
    if rows:
        Cx.barh(range(len(rows)), [r[1] for r in rows],
                xerr=[[r[1] - r[2] for r in rows], [r[3] - r[1] for r in rows]],
                color=C_KEEP, error_kw=dict(lw=.8))
        Cx.set_yticks(range(len(rows)))
        Cx.set_yticklabels([f"chr{r[0]}" for r in rows], fontsize=8)
        Cx.invert_yaxis()
        Cx.axvline(np.mean([r[1] for r in rows]), color=C_CUT, ls="--", lw=.8)
    Cx.set(xlabel="COs per Mb (density-matched, IQR band)",
           title="C. Per-chromosome rate\nmatched on marker density")

    D = ax[1][1]
    for c, L, n, cm, cmb in per_chrom:
        col = C_RAW if cm < FLOOR_CM else C_KEEP
        D.scatter(L, cm, color=col, s=30)
        D.annotate(c, (L, cm), fontsize=7, xytext=(4, 4), textcoords="offset points")
    D.axhline(FLOOR_CM, color=C_RAW, ls="--",
              label=f"obligate-CO floor {FLOOR_CM:.0f} cM")
    D.set(xlabel="chromosome length (Mb)", ylabel="genetic length (cM)",
          title=f"D. Map length vs physical length\ntotal {total_cm:.0f} cM")
    D.legend(fontsize=8); D.grid(alpha=.3)

    save(fig, f"{a.out_dir}/saturation.png")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    per_cell = defaultdict(int)
    for (cell, c), n in co_by_pair.items():
        per_cell[cell] += n
    counts = np.array([per_cell.get(c, 0) for c in cells])
    ax.hist(counts, bins=np.arange(-0.5, counts.max() + 1.5, 1), color=C_KEEP)
    ax.axvline(counts.mean(), color=C_CUT, ls="--",
               label=f"mean {counts.mean():.2f}")
    ax.set(xlabel="COs per cell", ylabel="cells",
           title=f"{a.sample}: COs per cell (n={len(counts)})\n"
                 f"{cos_per_bivalent:.2f} per bivalent")
    ax.legend(fontsize=8)
    save(fig, f"{a.out_dir}/co_per_cell.png")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    d2e = []
    for c, m in mids:
        if c in lens:
            d2e.append(min(m, lens[c] - m) / 1e6)
    d2e = np.array(d2e)
    ax[0].hist(d2e, bins=40, color=C_KEEP)
    ax[0].axvline(a.block_size / 1e6, color=C_RAW, ls="--",
                  label=f"block_size {a.block_size/1e6:.2g} Mb")
    ax[0].set(xlabel="distance to nearest chromosome end (Mb)", ylabel="COs",
              title="CO position relative to chromosome ends\n"
                    "a hard edge at block_size = the blind zone")
    ax[0].legend(fontsize=8)
    ax[1].hist(widths / 1000, bins=60, color=C_KEEP)
    ax[1].axvline(a.block_size / 1000, color=C_RAW, ls="--", label="block_size")
    ax[1].set_yscale("log")
    ax[1].set(xlabel="CO interval width (kb)", ylabel="COs",
              title=f"interval width (median {np.median(widths)/1000:.0f} kb)")
    ax[1].legend(fontsize=8)
    save(fig, f"{a.out_dir}/co_positions.png")

    with open(f"{a.out_dir}/per_chrom.tsv", "w") as f:
        f.write("chrom\tsize_Mb\tn_cos\tcM\tcM_per_Mb\tclears_floor\n")
        for c, L, n, cm, cmb in per_chrom:
            f.write(f"{c}\t{L:.1f}\t{n}\t{cm:.1f}\t{cmb:.2f}\t"
                    f"{'yes' if cm >= FLOOR_CM else 'NO'}\n")

    S.echo(); F.echo()
    S.write(f"{a.out_dir}/summary.tsv")
    F.write(f"{a.out_dir}/flags.txt")


if __name__ == "__main__":
    main()
