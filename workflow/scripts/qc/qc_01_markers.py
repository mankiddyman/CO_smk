#!/usr/bin/env python3
"""QC stage 01 -- marker discovery.

Characterises the het-SNP set that every downstream step depends on: how many
markers exist, how the filters attrit them, how evenly they are spaced, and
whether the depth window is placed sensibly relative to the actual coverage.

The MAX_DP_BELOW_PEAK flag exists because on Spondias max_dp=80 sat BELOW the
median het-SNP depth of 91 and silently rejected 72% of real heterozygous sites,
giving 0.13% heterozygosity against a published 0.62%. That was invisible until
CO calling failed days later.

Standalone:
  python3 qc_01_markers.py --sample S --raw_vcf ... --markers_vcf ... --fai ...
      --out_dir qc/S/01_markers --min_dp 10 --max_dp 150 --min_qual 30
      --min_alt_ratio 0.3 --max_alt_ratio 0.7 --block_size 1000000
"""
import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qc_common import (Summary, Flags, stream, read_fai, chrom_sort_key,
                       save, C_RAW, C_KEEP, C_REF, C_CUT)


def load_raw_het(vcf):
    """Biallelic het SNPs from the raw VCF -> (qual, dp, alt_fraction) arrays."""
    cmd = (f"bcftools view -m2 -M2 -v snps -g het {vcf} 2>/dev/null | "
           f"bcftools query -f '%QUAL\\t[%DP]\\t[%AD]\\n'")
    qual, dp, af = [], [], []
    for line in stream(cmd):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            d = int(f[1])
        except ValueError:
            continue
        if d == 0:
            continue
        ad = f[2].split(",")
        if len(ad) < 2:
            continue
        try:
            alt = int(ad[1])
            q = float(f[0])
        except ValueError:
            continue
        qual.append(q)
        dp.append(d)
        af.append(alt / d)
    return np.array(qual), np.array(dp, dtype=np.int64), np.array(af)


def load_marker_positions(vcf):
    """chrom -> sorted position array."""
    pos = {}
    for line in stream(f"bcftools query -f '%CHROM\\t%POS\\n' {vcf}"):
        c, p = line.split()
        pos.setdefault(c, []).append(int(p))
    return {c: np.sort(np.array(v)) for c, v in pos.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--raw_vcf", required=True)
    ap.add_argument("--markers_vcf", required=True)
    ap.add_argument("--fai", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--min_dp", type=int, required=True)
    ap.add_argument("--max_dp", type=int, required=True)
    ap.add_argument("--min_qual", type=float, required=True)
    ap.add_argument("--min_alt_ratio", type=float, required=True)
    ap.add_argument("--max_alt_ratio", type=float, required=True)
    ap.add_argument("--block_size", type=int, required=True)
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    S = Summary(a.sample, "01_markers")
    F = Flags(a.sample, "01_markers")

    # ---- parameters in force, so the report records what produced it --------
    for k in ("min_dp", "max_dp", "min_qual", "min_alt_ratio",
              "max_alt_ratio", "block_size"):
        S.add(f"param.{k}", getattr(a, k))

    # ---- raw het SNPs ------------------------------------------------------
    print("loading raw het SNPs...", file=sys.stderr)
    qual, dp, af = load_raw_het(a.raw_vcf)
    n_raw = len(dp)
    S.add("raw_het_snps", n_raw)

    keep = ((dp >= a.min_dp) & (dp <= a.max_dp) & (qual >= a.min_qual) &
            (af >= a.min_alt_ratio) & (af <= a.max_alt_ratio))
    S.add("raw_het_passing_all_filters", int(keep.sum()))
    S.add("retention_pct", round(100 * keep.mean(), 2))

    # independent attrition per filter (they overlap; this is the marginal view)
    rej = {
        "min_dp":  (dp < a.min_dp),
        "max_dp":  (dp > a.max_dp),
        "qual":    (qual < a.min_qual),
        "alt_ratio": (af < a.min_alt_ratio) | (af > a.max_alt_ratio),
    }
    for name, m in rej.items():
        S.add(f"rejected_by_{name}", int(m.sum()))
        S.add(f"rejected_by_{name}_pct", round(100 * m.mean(), 2))

    dp_med = float(np.median(dp))
    for q, lab in [(25, "q1"), (50, "median"), (75, "q3"), (90, "p90")]:
        S.add(f"raw_depth_{lab}", int(np.percentile(dp, q)))
    S.add("raw_depth_max", int(dp.max()))
    af_med = float(np.median(af))
    S.add("raw_alt_fraction_median", round(af_med, 4))

    # ---- final marker set --------------------------------------------------
    print("loading marker positions...", file=sys.stderr)
    pos = load_marker_positions(a.markers_vcf)
    lens = read_fai(a.fai)
    chroms = sorted(pos, key=chrom_sort_key)
    n_markers = sum(len(v) for v in pos.values())
    genome = sum(lens[c] for c in lens)
    S.add("markers_final", n_markers)
    n_pass = int(keep.sum())
    S.add("markers_removed_post_filter", n_pass - n_markers,
          "difference between filter_markers output and the final VCF = "
          "rrna_blacklist (or any other post-filter exclusion)")
    S.add("genome_size_bp", genome)
    het_pct = 100 * n_markers / genome
    S.add("implied_heterozygosity_pct", round(het_pct, 4),
          "markers / genome size; SNP-only, filtered -- expect below a k-mer estimate")

    gaps_all, dens, per_chrom = [], [], []
    for c in chroms:
        g = np.diff(pos[c])
        gaps_all.append(g)
        d = len(pos[c]) / (lens[c] / 1e6)
        dens.append(d)
        per_chrom.append((c, lens[c], len(pos[c]), d, g))
    gaps_all = np.concatenate(gaps_all) if gaps_all else np.array([0])

    S.add("marker_spacing_median_bp", int(np.median(gaps_all)))
    S.add("marker_spacing_p99_bp", int(np.percentile(gaps_all, 99)))
    S.add("marker_spacing_max_bp", int(gaps_all.max()))
    n_big = int((gaps_all > a.block_size).sum())
    S.add("gaps_exceeding_block_size", n_big,
          "each is a CO-blind region: no crossover can be called across it")
    dens = np.array(dens)
    S.add("markers_per_mb_min", round(dens.min(), 0))
    S.add("markers_per_mb_max", round(dens.max(), 0))
    S.add("markers_per_mb_ratio", round(dens.max() / dens.min(), 2))
    S.add("expected_markers_per_block", round(dens.mean() * a.block_size / 1e6, 1),
          "genome-wide density x block_size; a CELL sees far fewer")

    # ---- flags -------------------------------------------------------------
    F.check(a.max_dp < dp_med, "MAX_DP_BELOW_PEAK",
            f"max_dp={a.max_dp} is BELOW the median het-SNP depth ({dp_med:.0f}). "
            f"It is decapitating the diploid peak, not trimming collapsed "
            f"paralogs. Currently rejecting {100*rej['max_dp'].mean():.1f}% of "
            f"het sites. Set max_dp to ~1.5-2x the median.")
    F.check(abs(af_med - 0.5) > 0.02, "REFERENCE_BIAS",
            f"median ALT fraction {af_med:.3f} (expect 0.5). ALT-carrying reads "
            f"align marginally worse and are lost -- costs sensitivity via "
            f"dropout rather than creating false switches.")
    F.check(n_big > 0, "MARKER_GAPS",
            f"{n_big} marker gaps exceed block_size ({a.block_size:,} bp); "
            f"largest {gaps_all.max():,} bp. CO-blind regions.")
    F.check(het_pct < 0.1, "LOW_HET",
            f"implied heterozygosity {het_pct:.3f}% is suspiciously low. Check "
            f"max_dp placement and compare against any published estimate.")
    F.check(dens.max() / dens.min() > 3, "UNEVEN_MARKER_DENSITY",
            f"markers/Mb varies {dens.max()/dens.min():.1f}x across chromosomes "
            f"({dens.min():.0f}-{dens.max():.0f}). Detection efficiency will "
            f"vary with it.")

    # ---- plots -------------------------------------------------------------
    print("\nplotting...", file=sys.stderr)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    clip = int(min(dp.max(), np.percentile(dp, 99.9) * 2.5))
    ax.hist(dp[dp <= clip], bins=200, color=C_RAW)
    ax.axvline(a.min_dp, color="#1f77b4", ls="--", lw=1.2, label=f"min_dp={a.min_dp}")
    ax.axvline(a.max_dp, color=C_CUT, ls="--", lw=1.2, label=f"max_dp={a.max_dp}")
    ax.axvline(dp_med, color=C_REF, ls="--", lw=1.2, label=f"median={dp_med:.0f}")
    ax.set(xlabel="HiFi depth at het SNP (DP)", ylabel="het SNPs",
           title=f"{a.sample}: depth at het SNPs\n"
                 f"diploid peak + collapsed-paralog tail")
    ax.legend(fontsize=8)
    save(fig, f"{a.out_dir}/depth.png")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(af, bins=100, color=C_KEEP)
    ax.axvline(a.min_alt_ratio, color=C_CUT, ls="--", lw=1.2,
               label=f"{a.min_alt_ratio}/{a.max_alt_ratio} cut")
    ax.axvline(a.max_alt_ratio, color=C_CUT, ls="--", lw=1.2)
    ax.axvline(0.5, color=C_REF, ls=":", lw=1.2, label="expected het 0.5")
    ax.set(xlabel="ALT fraction (AD/DP)", ylabel="het SNPs",
           title=f"{a.sample}: allele ratio\nmedian={af_med:.3f}")
    ax.legend(fontsize=8)
    save(fig, f"{a.out_dir}/allele_ratio.png")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(np.clip(qual, 0, 300), bins=150, color="#9467bd")
    ax.axvline(a.min_qual, color=C_CUT, ls="--", lw=1.2,
               label=f"min_qual={a.min_qual:.0f}")
    ax.set(xlabel="QUAL (clipped at 300)", ylabel="het SNPs",
           title=f"{a.sample}: variant quality")
    ax.legend(fontsize=8)
    save(fig, f"{a.out_dir}/qual.png")

    ncol = 4
    nrow = int(np.ceil(len(chroms) / ncol))
    fig = plt.figure(figsize=(4 * ncol, 2.6 * (nrow + 1.4)))
    gs = fig.add_gridspec(nrow + 1, ncol, height_ratios=[1.5] + [1] * nrow)
    axg = fig.add_subplot(gs[0, :])
    bins = np.logspace(0, 7, 90)
    axg.hist(gaps_all, bins=bins, color="#2ca02c")
    axg.set_xscale("log"); axg.set_yscale("log")
    axg.axvline(a.block_size, color=C_RAW, ls="--",
                label=f"block_size {a.block_size/1e6:.2g} Mb")
    axg.axvline(np.median(gaps_all), color=C_CUT, ls=":",
                label=f"median {np.median(gaps_all):,.0f} bp")
    axg.set(xlabel="gap to next marker (bp)", ylabel="count",
            title=f"{a.sample}: marker spacing, genome-wide "
                  f"({n_markers:,} markers, {n_big} gaps > block_size)")
    axg.legend(fontsize=9)
    for i, (c, L, n, d, g) in enumerate(per_chrom):
        ax = fig.add_subplot(gs[1 + i // ncol, i % ncol])
        ax.hist(g, bins=np.logspace(0, 7, 50), color=C_KEEP)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.axvline(a.block_size, color=C_RAW, ls="--", lw=.8)
        ax.set_title(f"chr{c}  n={n:,}  med={np.median(g):,.0f}", fontsize=9)
        ax.tick_params(labelsize=7)
    save(fig, f"{a.out_dir}/spacing.png", dpi=110)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    labs = [f"chr{c}" for c, *_ in per_chrom]
    ax.bar(range(len(per_chrom)), dens, color=C_KEEP)
    ax.axhline(dens.mean(), color=C_CUT, ls="--", lw=1,
               label=f"mean {dens.mean():.0f}/Mb")
    ax.set_xticks(range(len(per_chrom)))
    ax.set_xticklabels(labs, rotation=45, fontsize=8)
    ax.set(ylabel="markers per Mb",
           title=f"{a.sample}: marker density by chromosome "
                 f"(max/min = {dens.max()/dens.min():.1f}x)")
    ax.legend(fontsize=8)
    save(fig, f"{a.out_dir}/per_chrom.png")

    # ---- per-chromosome table ---------------------------------------------
    with open(f"{a.out_dir}/per_chrom.tsv", "w") as f:
        f.write("chrom\tlength_bp\tmarkers\tmarkers_per_mb\t"
                "median_gap_bp\tmax_gap_bp\tgaps_over_block\n")
        for c, L, n, d, g in per_chrom:
            f.write(f"{c}\t{L}\t{n}\t{d:.1f}\t{np.median(g):.0f}\t"
                    f"{g.max()}\t{(g > a.block_size).sum()}\n")

    S.echo(); F.echo()
    S.write(f"{a.out_dir}/summary.tsv")
    F.write(f"{a.out_dir}/flags.txt")


if __name__ == "__main__":
    main()
