#!/usr/bin/env python3
"""
Convert cellsnp-lite output to per-cell TSV files in the format expected
by hapCO_identification.R.

Reads:
  cellSNP.tag.DP.mtx, cellSNP.tag.AD.mtx, cellSNP.tag.OTH.mtx
  cellSNP.samples.tsv (cell barcodes, one per row)
  cellSNP.base.vcf.gz (sites with coverage; row order matches matrix rows)

Writes:
  per-cell TSV files: chrom, pos, ref, ref_count, alt, alt_count
  chrom_map.tsv: internal_id, name, size
  conversion_summary.txt: stats

Cells with fewer than --min_markers covered sites are skipped (no file written).
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from pathlib import Path

import scipy.io
import scipy.sparse


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--cellsnp_dir", required=True,
                   help="Directory with cellSNP.tag.{DP,AD}.mtx etc")
    p.add_argument("--fai", required=True,
                   help="Reference FASTA index (chrom<TAB>size<TAB>...)")
    p.add_argument("--out_dir", required=True,
                   help="Output directory for per-cell TSVs")
    p.add_argument("--min_markers", type=int, default=50,
                   help="Skip cells with fewer than this many covered markers (default: 50)")
    return p.parse_args()


def read_sites(vcf_gz):
    """Read site table from cellsnp's base.vcf.gz. Returns list of (chrom, pos, ref, alt)."""
    sites = []
    with gzip.open(vcf_gz, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            chrom, pos, _id, ref, alt = fields[:5]
            sites.append((chrom, int(pos), ref, alt))
    return sites


def read_barcodes(samples_tsv):
    """Read cell barcodes (one per line)."""
    with open(samples_tsv) as f:
        return [line.strip() for line in f if line.strip()]


def read_chrom_sizes(fai):
    """Read chromosome sizes from FASTA index. Preserve order."""
    chrom_sizes = []  # [(name, size), ...]
    with open(fai) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                chrom_sizes.append((parts[0], int(parts[1])))
    return chrom_sizes


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cellsnp = Path(args.cellsnp_dir)

    # Load metadata
    print(f"[INFO] Reading sites from {cellsnp}/cellSNP.base.vcf.gz", file=sys.stderr)
    sites = read_sites(cellsnp / "cellSNP.base.vcf.gz")
    n_sites = len(sites)

    print(f"[INFO] Reading barcodes from {cellsnp}/cellSNP.samples.tsv", file=sys.stderr)
    barcodes = read_barcodes(cellsnp / "cellSNP.samples.tsv")
    n_cells = len(barcodes)

    print(f"[INFO] {n_sites} sites x {n_cells} cells", file=sys.stderr)

    # Load matrices (scipy returns 1-indexed COO; convert to 0-indexed CSC for column access)
    print(f"[INFO] Loading DP matrix...", file=sys.stderr)
    dp = scipy.io.mmread(str(cellsnp / "cellSNP.tag.DP.mtx")).tocsc()
    print(f"[INFO] Loading AD matrix...", file=sys.stderr)
    ad = scipy.io.mmread(str(cellsnp / "cellSNP.tag.AD.mtx")).tocsc()

    if dp.shape != (n_sites, n_cells):
        sys.exit(f"DP matrix shape {dp.shape} != ({n_sites}, {n_cells})")
    if ad.shape != (n_sites, n_cells):
        sys.exit(f"AD matrix shape {ad.shape} != ({n_sites}, {n_cells})")

    # Build chrom_map.tsv (only chromosomes that appear in sites)
    chrom_sizes_all = read_chrom_sizes(args.fai)
    site_chroms = {chrom for chrom, _, _, _ in sites}
    chrom_map = []
    for name, size in chrom_sizes_all:
        if name in site_chroms:
            chrom_map.append((len(chrom_map) + 1, name, size))

    chrom_map_path = out_dir / "chrom_map.tsv"
    with open(chrom_map_path, 'w') as f:
        f.write("internal_id\tname\tsize\n")
        for cid, name, size in chrom_map:
            f.write(f"{cid}\t{name}\t{size}\n")
    print(f"[INFO] Wrote chrom_map: {len(chrom_map)} chromosomes", file=sys.stderr)

    # Sort sites by (chrom, pos) for clean output, but maintain mapping back to matrix row indices
    site_order = sorted(range(n_sites), key=lambda i: (sites[i][0], sites[i][1]))

    # Per-cell TSVs
    n_written = 0
    n_skipped_empty = 0
    n_skipped_below_min = 0
    n_markers_per_cell = []

    for cell_idx, barcode in enumerate(barcodes):
        # Get this cell's column from DP — non-zero rows = sites with coverage
        col_dp = dp.getcol(cell_idx)
        col_ad = ad.getcol(cell_idx)

        # COO format gives (row indices, values)
        cov_rows = col_dp.nonzero()[0]  # site indices with DP > 0
        n_cov = len(cov_rows)
        n_markers_per_cell.append(n_cov)

        if n_cov == 0:
            n_skipped_empty += 1
            continue

        if n_cov < args.min_markers:
            n_skipped_below_min += 1
            continue

        # Build dense lookups for this cell (only over the covered rows)
        # col_dp and col_ad are sparse; .toarray().flatten() is fine for column selection
        dp_vals = col_dp.toarray().flatten()
        ad_vals = col_ad.toarray().flatten()

        # Sort the covered sites by (chrom, pos) using the precomputed site_order
        # Faster: filter site_order to just the covered rows for this cell
        cov_rows_set = set(cov_rows.tolist())
        cell_site_order = [i for i in site_order if i in cov_rows_set]

        out_file = out_dir / f"{barcode}.tsv"
        with open(out_file, 'w') as f:
            for i in cell_site_order:
                chrom, pos, ref, alt = sites[i]
                dp_v = int(dp_vals[i])
                ad_v = int(ad_vals[i])
                ref_count = dp_v - ad_v  # OTH ignored — fine for CO calling
                if ref_count < 0:
                    ref_count = 0  # shouldn't happen but safety
                # Format: chrom, pos, ref_base, ref_count, alt_base, alt_count
                f.write(f"{chrom}\t{pos}\t{ref}\t{ref_count}\t{alt}\t{ad_v}\n")

        n_written += 1

    # Summary
    summary_path = out_dir / "conversion_summary.txt"
    with open(summary_path, 'w') as f:
        f.write(f"Total cells in cellsnp output: {n_cells}\n")
        f.write(f"Cells with no coverage (skipped): {n_skipped_empty}\n")
        f.write(f"Cells below min_markers={args.min_markers} (skipped): {n_skipped_below_min}\n")
        f.write(f"Cells written: {n_written}\n")
        f.write(f"Total sites: {n_sites}\n")
        f.write(f"\n")
        f.write(f"Marker count distribution across all cells (n=marker count):\n")
        if n_markers_per_cell:
            sorted_counts = sorted(n_markers_per_cell)
            n = len(sorted_counts)
            f.write(f"  min:    {sorted_counts[0]}\n")
            f.write(f"  q1:     {sorted_counts[n//4]}\n")
            f.write(f"  median: {sorted_counts[n//2]}\n")
            f.write(f"  q3:     {sorted_counts[3*n//4]}\n")
            f.write(f"  max:    {sorted_counts[-1]}\n")
            f.write(f"  mean:   {sum(sorted_counts)/n:.1f}\n")

    print(f"[INFO] Done. Wrote {n_written} per-cell files. See {summary_path}.",
          file=sys.stderr)


if __name__ == "__main__":
    main()
