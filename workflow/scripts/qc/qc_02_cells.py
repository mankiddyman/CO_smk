#!/usr/bin/env python3
"""QC stage 02 -- alignment, cells, and per-cell genotypes.

Covers everything between alignment and the per-cell haplotype traces. This is
where the assay's fundamental limits live.

Two methods here are deliberate and were arrived at by getting them wrong first:

1. COVERAGE_HOTSPOT uses BINNED depth, not per-chromosome means. The Spondias
   rDNA array was 49 kb inside a 14.8 Mb chromosome and held 85% of the BAM;
   the chromosome mean showed 10.7x and hid it completely.

2. Contamination is measured PER MARKER, not per observation. A fixed ratio
   window ("clean = <=0.10 or >=0.90") makes deeper observations look dirtier,
   because at DP=10 a single contaminating read gives exactly 0.10. That
   produced a bogus 13%->32% "contamination rises with depth" result. A cell
   cannot be 50% contaminated at one position and clean at the next, but a bad
   marker can -- so classify markers by their aggregate behaviour across cells
   and report the rate at the well-behaved ones.

Also note: AD and DP sparse matrices have DIFFERENT non-zero entry sets (a
marker with 2 REF reads has a DP entry and no AD entry). Join on (row, col)
with AD defaulting to 0; never zip their value columns.
"""
import argparse
import os
import sys
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_common import (Summary, Flags, run, read_fai, chrom_sort_key, save,
                       C_RAW, C_KEEP, C_REF, C_CUT)


def read_mtx(path):
    """MatrixMarket coordinate file -> {(row, col): value}. Skips the header."""
    d = {}
    with open(path) as f:
        for line in f:
            if line.startswith("%"):
                continue
            break                      # this line is the dims header; discard
        for line in f:
            r, c, v = line.split()
            d[(int(r), int(c))] = int(v)
    return d


def star_log(path):
    """Parse Log.final.out into {label: value}."""
    out = {}
    for line in open(path):
        if "|" not in line:
            continue
        k, v = line.split("|", 1)
        out[k.strip()] = v.strip()
    return out


def pct(s):
    return float(s.rstrip("%")) if s.endswith("%") else float(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--star_log", required=True)
    ap.add_argument("--matrix_dir", required=True, help="Solo.out/GeneFull/raw")
    ap.add_argument("--barcodes_called", required=True)
    ap.add_argument("--cell_qc_summary", required=True, help="cells/*/qc_summary.tsv")
    ap.add_argument("--snps_dir", required=True, help="results/snps/{sample}")
    ap.add_argument("--switches", required=True)
    ap.add_argument("--bam", required=True,
                    help="RAW StarSolo BAM. Must NOT be MAPQ-filtered: "
                         "collapsed-repeat reads are MAPQ 0, so a "
                         "filtered BAM hides exactly what this looks for.")
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--min_total_markers", type=int, required=True)
    ap.add_argument("--min_per_chrom_markers", type=int, required=True)
    ap.add_argument("--max_switch_rate", type=float, required=True)
    ap.add_argument("--default_lower", type=int, required=True)
    ap.add_argument("--min_mapq", type=int, required=True)
    ap.add_argument("--bin_kb", type=int, default=50)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    S = Summary(a.sample, "02_cells")
    F = Flags(a.sample, "02_cells")
    for k in ("min_total_markers", "min_per_chrom_markers", "max_switch_rate",
              "default_lower", "min_mapq"):
        S.add(f"param.{k}", getattr(a, k))

    # ---- alignment ---------------------------------------------------------
    print("parsing STAR log...", file=sys.stderr)
    L = star_log(a.star_log)
    n_in = int(L.get("Number of input reads", 0))
    uniq = pct(L.get("Uniquely mapped reads %", "0"))
    multi = pct(L.get("% of reads mapped to multiple loci", "0"))
    many = pct(L.get("% of reads mapped to too many loci", "0"))
    short = pct(L.get("% of reads unmapped: too short", "0"))
    other = pct(L.get("% of reads unmapped: other", "0"))
    S.add("input_reads", n_in)
    S.add("uniquely_mapped_pct", uniq)
    S.add("multi_mapped_pct", multi)
    S.add("too_many_loci_pct", many)
    S.add("unmapped_too_short_pct", short)
    S.add("unmapped_other_pct", other)
    S.add("total_aligned_pct", round(uniq + multi + many, 2))

    # ---- cells -------------------------------------------------------------
    print("loading count matrix...", file=sys.stderr)
    from scipy.io import mmread
    M = mmread(f"{a.matrix_dir}/matrix.mtx").tocsc()
    bcs = [l.strip() for l in open(f"{a.matrix_dir}/barcodes.tsv")]
    umi = np.asarray(M.sum(axis=0)).ravel()
    called = set(l.split()[0] for l in open(a.barcodes_called))
    is_cell = np.array([b in called for b in bcs])

    S.add("barcodes_total", len(bcs))
    S.add("barcodes_called", int(is_cell.sum()))
    S.add("median_umi_per_cell", int(np.median(umi[is_cell])))
    genes_per_cell = np.asarray((M[:, is_cell] > 0).sum(axis=0)).ravel()
    S.add("median_genes_per_cell", int(np.median(genes_per_cell)))

    amb = umi[(~is_cell) & (umi > 0)]
    amb_frac = float(amb.mean() / np.median(umi[is_cell])) if len(amb) else 0.0
    S.add("ambient_droplets", len(amb))
    S.add("ambient_umi_share_pct", round(100 * amb.sum() / umi.sum(), 2),
          "share of ALL UMIs outside called cells; not per-cell contamination")
    S.add("ambient_per_cell_estimate_pct", round(100 * amb_frac, 2),
          "mean ambient UMI per droplet / median cell UMI")

    # ---- genotyping depth --------------------------------------------------
    print("loading cellSNP matrices...", file=sys.stderr)
    DP = read_mtx(f"{a.snps_dir}/cellSNP.tag.DP.mtx")
    AD = read_mtx(f"{a.snps_dir}/cellSNP.tag.AD.mtx")
    dp = np.fromiter(DP.values(), dtype=np.int32, count=len(DP))
    S.add("cell_marker_observations", len(dp))
    S.add("mean_dp_per_observation", round(float(dp.mean()), 3))
    S.add("pct_single_read", round(100 * float((dp == 1).mean()), 1))
    S.add("pct_dp_ge3", round(100 * float((dp >= 3).mean()), 1))
    S.add("total_umis_at_markers", int(dp.sum()))
    alt_total = sum(AD.get(k, 0) for k in DP)
    S.add("global_alt_fraction", round(alt_total / dp.sum(), 4),
          "0.5 = unbiased; below indicates reference bias via ALT dropout")

    # ---- contamination, PER MARKER ----------------------------------------
    print("estimating contamination (per-marker method)...", file=sys.stderr)
    mk_tot, mk_minor, mk_n = defaultdict(int), defaultdict(int), defaultdict(int)
    for (m, c), v in DP.items():
        if v < 5:
            continue
        alt = AD.get((m, c), 0)
        mk_tot[m] += v
        mk_minor[m] += min(alt, v - alt)
        mk_n[m] += 1
    well = [m for m in mk_n if mk_n[m] >= 5 and mk_minor[m] / mk_tot[m] <= 0.05]
    bad = [m for m in mk_n if mk_n[m] >= 5 and mk_minor[m] / mk_tot[m] > 0.15]
    n_eval = sum(1 for m in mk_n if mk_n[m] >= 5)
    S.add("markers_evaluated_for_contamination", n_eval,
          ">=5 observations at DP>=5")
    S.add("markers_well_behaved", len(well), "aggregate minor-allele rate <=5%")
    S.add("markers_persistently_intermediate", len(bad), ">15%; artefact class")
    if well:
        wset = set(well)
        num = sum(min(AD.get(k, 0), v - AD.get(k, 0))
                  for k, v in DP.items() if v >= 5 and k[0] in wset)
        den = sum(v for k, v in DP.items() if v >= 5 and k[0] in wset)
        contam = 100 * num / den
    else:
        contam = float("nan")
    S.add("contamination_pct", round(contam, 3),
          "minor-allele rate at well-behaved markers; ~sequencing error floor")

    # ---- per-cell markers and switch rate ---------------------------------
    print("loading switch diagnostics...", file=sys.stderr)
    hdr = open(a.switches).readline().rstrip("\n").split("\t")
    mcols = [i for i, h in enumerate(hdr) if h.startswith("markers_")]
    tot_mk, sw, per_chrom_mk = [], [], []
    for line in open(a.switches):
        f = line.rstrip("\n").split("\t")
        if f[0] == "barcode":
            continue
        tot_mk.append(int(f[1])); sw.append(float(f[3]))
        per_chrom_mk.extend(int(f[i]) for i in mcols)
    tot_mk = np.array(tot_mk); sw = np.array(sw)
    per_chrom_mk = np.array(per_chrom_mk)
    n_chrom = len(mcols)

    S.add("cells_in_switch_table", len(tot_mk))
    n_called = int(is_cell.sum())
    S.add("switch_table_minus_called", len(tot_mk) - n_called,
          "switch_diagnostics globs results/cell_data/, which accumulates "
          "per-cell TSVs across runs. A positive value means stale files from "
          "an earlier barcode list or marker set are present.")
    for q, lab in [(0, "min"), (25, "q1"), (50, "median"), (75, "q3"), (90, "p90")]:
        S.add(f"markers_per_cell_{lab}", int(np.percentile(tot_mk, q)))
    S.add("markers_per_cell_max", int(tot_mk.max()))
    for q, lab in [(25, "q1"), (50, "median"), (75, "q3"), (90, "p90")]:
        S.add(f"switch_rate_{lab}", round(float(np.percentile(sw, q)), 4))

    # selection funnel
    f_mk = tot_mk >= a.min_total_markers
    f_sw = sw <= a.max_switch_rate
    S.add("cells_passing_total_markers", int(f_mk.sum()))
    S.add("cells_passing_switch_rate", int(f_sw.sum()))
    S.add("cells_passing_total_and_switch", int((f_mk & f_sw).sum()))

    # full funnel: select_cells ALSO requires min_per_chrom_markers on EVERY
    # chromosome simultaneously. That last condition is much harsher than it
    # looks -- on Spondias it took 196 -> 153.
    n_final = 0
    for line in open(a.switches):
        f = line.rstrip("\n").split("\t")
        if f[0] == "barcode":
            continue
        if (int(f[1]) >= a.min_total_markers
                and float(f[3]) <= a.max_switch_rate
                and all(int(f[i]) >= a.min_per_chrom_markers for i in mcols)):
            n_final += 1
    S.add("cells_passing_all_filters", n_final,
          "adds min_per_chrom_markers on ALL chromosomes; this is the number "
          "select_cells writes to good_cells.tsv")

    # per-chromosome markers among cells that pass the total-marker filter
    sel_mk = []
    for line in open(a.switches):
        f = line.rstrip("\n").split("\t")
        if f[0] == "barcode":
            continue
        if int(f[1]) >= a.min_total_markers and float(f[3]) <= a.max_switch_rate:
            sel_mk.extend(int(f[i]) for i in mcols)
    sel_mk = np.array(sel_mk) if sel_mk else np.array([0])
    S.add("selected_markers_per_chrom_median", int(np.median(sel_mk)))
    S.add("selected_markers_per_chrom_q1", int(np.percentile(sel_mk, 25)))

    # ---- binned coverage: the rDNA catch ----------------------------------
    print(f"binned coverage ({a.bin_kb} kb windows)...", file=sys.stderr)
    lens = read_fai(a.fai)
    chroms = sorted(lens, key=chrom_sort_key)
    bed = f"{a.out_dir}/_bins.bed"
    w = a.bin_kb * 1000
    with open(bed, "w") as f:
        for c in chroms:
            for s0 in range(0, lens[c], w):
                f.write(f"{c}\t{s0}\t{min(s0 + w, lens[c])}\n")
    n_bins = sum(1 for _ in open(bed))
    out = run(f"samtools bedcov {bed} {a.bam}")
    cov = defaultdict(list)
    for line in out.splitlines():
        f = line.split("\t")
        span = int(f[2]) - int(f[1])
        if span > 0:
            cov[f[0]].append((int(f[1]), int(f[3]) / span))
    os.remove(bed)

    hotspots = []
    for c in chroms:
        if c not in cov:
            continue
        d = np.array([x[1] for x in cov[c]])
        med = np.median(d[d > 0]) if (d > 0).any() else 0
        if med <= 0:
            continue
        for (s0, dd) in cov[c]:
            if dd > 20 * med:
                hotspots.append((c, s0, dd, dd / med))
    hotspots.sort(key=lambda x: -x[3])
    S.add("coverage_bam", os.path.basename(a.bam),
          "must be the raw BAM; a MAPQ-filtered one hides collapsed repeats")
    S.add("coverage_bins", n_bins)
    S.add("coverage_hotspot_bins", len(hotspots),
          f">20x the chromosome median depth, {a.bin_kb} kb windows")
    if hotspots:
        c, s0, dd, r = hotspots[0]
        S.add("worst_hotspot", f"{c}:{s0}-{s0+w}")
        S.add("worst_hotspot_fold", round(r, 1))

    with open(f"{a.out_dir}/coverage_hotspots.tsv", "w") as f:
        f.write("chrom\tstart\tend\tmean_depth\tfold_over_chrom_median\n")
        for c, s0, dd, r in hotspots[:200]:
            f.write(f"{c}\t{s0}\t{s0+w}\t{dd:.1f}\t{r:.1f}\n")

    # ---- switch-rate bimodality -------------------------------------------
    h, edges = np.histogram(sw, bins=60)
    peak = int(np.argmax(h))
    left = h[:peak]
    bimodal = False
    if len(left) > 5:
        valley = left.min()
        second = left.max()
        bimodal = second > 3 * max(valley, 1) and second > 0.02 * h[peak]
    S.add("switch_rate_bimodal", str(bimodal),
          "a clean-haploid shoulder separated from the noise mode is what "
          "gives max_switch_rate an empirical basis")

    # ---- flags -------------------------------------------------------------
    F.check(len(hotspots) > 0, "COVERAGE_HOTSPOT",
            (f"{len(hotspots)} bins exceed 20x their chromosome median; worst "
             f"{hotspots[0][0]}:{hotspots[0][1]}-{hotspots[0][1]+w} at "
             f"{hotspots[0][3]:.0f}x. Likely a collapsed repeat (rDNA/NOR). "
             f"Check whether it dominates the library and whether markers "
             f"there should be blacklisted.") if hotspots else "")
    F.check(uniq < 50, "LOW_UNIQUE_MAPPING",
            f"uniquely mapped {uniq:.1f}%. If multi-mapping is high, check "
            f"whether it concentrates in one locus (see COVERAGE_HOTSPOT) "
            f"rather than being a genome-wide property.")
    F.check(dp.mean() < 2, "SHALLOW_GENOTYPING",
            f"mean DP per cell x marker is {dp.mean():.2f} "
            f"({100*(dp==1).mean():.1f}% single-read). Genotype calls rest on "
            f"one read; expect a background switch rate and reduced CO "
            f"detection.")
    F.check(not bimodal, "SWITCH_RATE_UNIMODAL",
            "no clear clean-haploid shoulder below the noise mode, so "
            "max_switch_rate has no empirical basis in this dataset.")
    F.check(contam == contam and contam > 5, "HIGH_CONTAMINATION",
            f"minor-allele rate at well-behaved markers is {contam:.2f}%, "
            f"well above a sequencing-error floor.")
    F.check(len(tot_mk) > n_called, "STALE_CELL_DATA",
            f"switches.tsv has {len(tot_mk):,} cells but only {n_called:,} "
            f"barcodes were called. results/cell_data/ is not cleared between "
            f"runs, so per-cell TSVs from an earlier barcode list or marker set "
            f"may be feeding into select_cells and co_calling. Verify no "
            f"SELECTED cell has a TSV older than cellSNP.tag.DP.mtx.")
    F.check(100 * amb_frac > 10, "AMBIENT_HIGH",
            f"ambient estimate {100*amb_frac:.1f}% of median cell UMI.")

    # ---- plots -------------------------------------------------------------
    print("\nplotting...", file=sys.stderr)

    fig, ax = plt.subplots(figsize=(7.5, 5))
    u = np.sort(umi[umi > 0])[::-1]
    ax.plot(np.arange(1, len(u) + 1), u, color=C_KEEP, lw=1.5)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.axhline(a.default_lower, color=C_RAW, ls="--", lw=1,
               label=f"default_lower={a.default_lower}")
    ax.axvline(is_cell.sum(), color=C_CUT, ls=":", lw=1,
               label=f"{is_cell.sum():,} called")
    ax.set(xlabel="barcode rank", ylabel="UMI count",
           title=f"{a.sample}: barcode rank\n"
                 f"median cell {np.median(umi[is_cell]):.0f} UMI")
    ax.legend(fontsize=8); ax.grid(alpha=.3, which="both")
    save(fig, f"{a.out_dir}/knee.png")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(dp[dp <= 15], bins=np.arange(0.5, 16.5, 1), color=C_RAW)
    ax.set(xlabel="reads per cell x marker (DP)", ylabel="observations",
           title=f"{a.sample}: genotyping depth\n"
                 f"mean {dp.mean():.2f}, {100*(dp==1).mean():.1f}% single-read")
    save(fig, f"{a.out_dir}/dp.png")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    r3 = np.array([AD.get(k, 0) / v for k, v in DP.items() if v >= 3])
    if len(r3):
        ax.hist(r3, bins=np.linspace(0, 1, 21), color=C_KEEP)
    ax.set(xlabel="ALT fraction (AD/DP)", ylabel="observations",
           title=f"{a.sample}: allele ratio at DP>=3 (n={len(r3):,})\n"
                 f"spikes are the DP lattice (at DP=3 only 0, 1/3, 2/3, 1), "
                 f"NOT populations")
    save(fig, f"{a.out_dir}/allele_ratio_cells.png")

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.hist(sw, bins=60, color=C_REF)
    ax.axvline(a.max_switch_rate, color=C_CUT, ls="--",
               label=f"max_switch_rate={a.max_switch_rate}")
    ax.set(xlabel="per-cell switch rate", ylabel="cells",
           title=f"{a.sample}: switch rate (n={len(sw):,})\n"
                 f"bimodal={bimodal} -- the clean shoulder justifies the cutoff")
    ax.legend(fontsize=8)
    save(fig, f"{a.out_dir}/switch_rate.png")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].hist(tot_mk, bins=60, color=C_KEEP)
    ax[0].axvline(a.min_total_markers, color=C_CUT, ls="--",
                  label=f"min_total_markers={a.min_total_markers}")
    ax[0].set_yscale("log")
    ax[0].set(xlabel="markers per cell", ylabel="cells",
              title=f"markers per cell (median {np.median(tot_mk):.0f})")
    ax[0].legend(fontsize=8)
    ax[1].hist(sel_mk, bins=60, color=C_KEEP)
    ax[1].axvline(a.min_per_chrom_markers, color=C_CUT, ls="--",
                  label=f"min_per_chrom={a.min_per_chrom_markers}")
    ax[1].set(xlabel="markers per chromosome (cells passing total+switch)",
              ylabel="cell x chromosome",
              title=f"per chromosome (median {np.median(sel_mk):.0f} "
                    f"over {n_chrom} chromosomes)")
    ax[1].legend(fontsize=8)
    save(fig, f"{a.out_dir}/markers_per_cell.png")

    ncol = 4
    nrow = int(np.ceil(len(chroms) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 2.4 * nrow),
                             squeeze=False)
    for i, c in enumerate(chroms):
        ax = axes[i // ncol][i % ncol]
        if c in cov:
            xs = np.array([x[0] for x in cov[c]]) / 1e6
            ds = np.array([x[1] for x in cov[c]])
            ax.plot(xs, ds, lw=.7, color=C_KEEP)
            m = np.median(ds[ds > 0]) if (ds > 0).any() else 0
            ax.axhline(20 * m, color=C_RAW, ls="--", lw=.8)
            ax.set_yscale("log")
        ax.set_title(f"chr{c}", fontsize=9)
        ax.tick_params(labelsize=7)
    for j in range(len(chroms), nrow * ncol):
        axes[j // ncol][j % ncol].axis("off")
    fig.suptitle(f"{a.sample}: depth in {a.bin_kb} kb bins "
                 f"(dashed = 20x chromosome median; {len(hotspots)} hotspot bins)",
                 fontsize=11)
    save(fig, f"{a.out_dir}/coverage.png", dpi=110)

    S.echo(); F.echo()
    S.write(f"{a.out_dir}/summary.tsv")
    F.write(f"{a.out_dir}/flags.txt")


if __name__ == "__main__":
    main()
