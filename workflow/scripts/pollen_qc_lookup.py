#!/usr/bin/env python3
"""
Pollen QC: check expression of pollen marker orthologs in STARsolo count matrix.

Inputs:
  - BLAST output (tabular, outfmt 6) from blastp of pollen markers vs proteome
  - STARsolo raw matrix dir (matrix.mtx, features.tsv, barcodes.tsv)
  - --chromosomes: comma-separated list of "real chromosome" scaffold names
                   (used to distinguish chromosome-resident hits from debris)

For each marker, reports expression stats and flags debris-only hits.
"""

import argparse
from collections import defaultdict
from pathlib import Path
import sys

import scipy.io
import scipy.sparse
import numpy as np


def parse_blast(path, evalue_max=1e-10, top_n=5):
    """Parse BLAST tabular output. Returns dict marker -> [list of (gene_id, evalue, pident, qcovs)]."""
    hits_by_marker = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            qseqid, sseqid, pident, length, evalue, bitscore, qcovs = parts[:7]
            # Strip trailing .N transcript suffix to get gene ID
            gene_id = sseqid.rsplit(".", 1)[0] if sseqid.endswith((".1", ".2", ".3", ".4", ".5")) else sseqid
            try:
                evalue_f = float(evalue)
                pident_f = float(pident)
                qcovs_f = float(qcovs)
            except ValueError:
                continue
            if evalue_f > evalue_max:
                continue
            marker = qseqid.split("|")[0]
            hits_by_marker[marker].append((gene_id, evalue_f, pident_f, qcovs_f))
    result = {}
    for m, hits in hits_by_marker.items():
        seen = set()
        unique_hits = []
        for h in sorted(hits, key=lambda x: x[1]):
            if h[0] not in seen:
                seen.add(h[0])
                unique_hits.append(h)
            if len(unique_hits) >= top_n:
                break
        result[m] = unique_hits
    return result


def gene_scaffold(gene_id):
    """Extract scaffold name from a gene ID like '_scaffold_33_000020' -> '_scaffold_33'.
    Returns the scaffold prefix or the full ID if format doesn't match."""
    parts = gene_id.split("_")
    if len(parts) >= 3 and parts[1] == "scaffold":
        return "_".join(parts[:3])  # e.g., '_scaffold_33'
    return gene_id


def load_starsolo(matrix_dir):
    md = Path(matrix_dir)
    print(f"  Loading {md / 'matrix.mtx'}...", file=sys.stderr)
    M = scipy.io.mmread(str(md / "matrix.mtx")).tocsr()
    print(f"  Loading features.tsv...", file=sys.stderr)
    gene_ids = []
    with open(md / "features.tsv") as fh:
        for line in fh:
            gene_ids.append(line.split("\t")[0])
    print(f"  Loading barcodes.tsv...", file=sys.stderr)
    with open(md / "barcodes.tsv") as fh:
        barcodes = [line.strip() for line in fh]
    print(f"  Matrix shape: {M.shape}, non-zero: {M.nnz}", file=sys.stderr)
    return M, gene_ids, barcodes


def identify_real_cells(M, n_top=10000):
    umi_per_cell = np.asarray(M.sum(axis=0)).flatten()
    if (umi_per_cell > 0).sum() < n_top:
        n_top = (umi_per_cell > 0).sum()
        print(f"  Only {n_top} barcodes have any UMIs; using all of them", file=sys.stderr)
    threshold = np.sort(umi_per_cell)[-n_top]
    mask = umi_per_cell >= threshold
    print(f"  Top {mask.sum()} barcodes have UMI >= {threshold} (used as 'real cells')",
          file=sys.stderr)
    return mask, umi_per_cell


def compute_marker_stats(M, gene_ids, marker_to_hits, cell_mask, real_chroms):
    gene_id_to_row = {g: i for i, g in enumerate(gene_ids)}
    rows = []
    for marker, hits in marker_to_hits.items():
        for rank, (gid, ev, pident, qcovs) in enumerate(hits, 1):
            scaffold = gene_scaffold(gid)
            on_chrom = (scaffold in real_chroms) if real_chroms else True
            row_idx = gene_id_to_row.get(gid)
            if row_idx is None:
                rows.append({
                    "marker": marker, "rank": rank, "gene_id": gid, "scaffold": scaffold,
                    "on_chrom": on_chrom,
                    "evalue": ev, "pident": pident, "qcovs": qcovs,
                    "in_matrix": False,
                    "n_cells_expr": 0, "total_umi_cells": 0,
                    "pct_cells_expressing": 0.0, "mean_umi_per_expr_cells": 0.0,
                })
                continue
            counts = M[row_idx, :].toarray().flatten()
            n_cells = cell_mask.sum()
            counts_in_cells = counts[cell_mask]
            mask_expr_cells = counts_in_cells > 0
            n_expr_cells = int(mask_expr_cells.sum())
            total_cells = int(counts_in_cells.sum())
            mean_cells = float(counts_in_cells[mask_expr_cells].mean()) if n_expr_cells else 0.0
            pct_cells = (n_expr_cells / n_cells * 100) if n_cells else 0.0
            rows.append({
                "marker": marker, "rank": rank, "gene_id": gid, "scaffold": scaffold,
                "on_chrom": on_chrom,
                "evalue": ev, "pident": pident, "qcovs": qcovs,
                "in_matrix": True,
                "n_cells_expr": n_expr_cells, "total_umi_cells": total_cells,
                "pct_cells_expressing": round(pct_cells, 2),
                "mean_umi_per_expr_cells": round(mean_cells, 2),
            })
    return rows


def emit_debris_warnings(rows, marker_to_hits):
    """Emit warning lines for markers where ALL hits are on debris (off-chromosome)."""
    warnings = []
    markers_seen = set()
    for marker in marker_to_hits.keys():
        marker_rows = [r for r in rows if r["marker"] == marker]
        if not marker_rows:
            continue
        on_chrom_count = sum(1 for r in marker_rows if r["on_chrom"])
        if on_chrom_count == 0:
            scaffolds = sorted(set(r["scaffold"] for r in marker_rows))
            warnings.append(f"WARNING: {marker} has no hits on real chromosomes. "
                            f"All hits are on debris contigs: {', '.join(scaffolds)}. "
                            "Expression signal for this marker is unreliable.")
    return warnings


def format_table(rows):
    cols = [
        ("marker", 8),
        ("rank", 4),
        ("scaffold", 18),
        ("chr?", 5),
        ("gene_id", 26),
        ("evalue", 11),
        ("pident", 7),
        ("qcovs", 6),
        ("n_cells_expr", 13),
        ("pct_cells", 10),
        ("mean_umi/cell", 14),
    ]
    header = "  ".join(name.ljust(width) for name, width in cols)
    sep = "  ".join("-" * width for _, width in cols)
    out = [header, sep]
    for r in rows:
        chr_flag = "YES" if r["on_chrom"] else "DEBRIS"
        cells = [
            str(r["marker"]).ljust(8),
            str(r["rank"]).ljust(4),
            str(r["scaffold"])[:18].ljust(18),
            chr_flag.ljust(5),
            str(r["gene_id"])[:26].ljust(26),
            f"{r['evalue']:.2e}".ljust(11),
            f"{r['pident']:.1f}".ljust(7),
            f"{r['qcovs']:.0f}".ljust(6),
            str(r["n_cells_expr"]).ljust(13),
            f"{r['pct_cells_expressing']:.2f}%".ljust(10),
            f"{r['mean_umi_per_expr_cells']:.2f}".ljust(14),
        ]
        out.append("  ".join(cells))
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast", required=True)
    ap.add_argument("--matrix_dir", required=True)
    ap.add_argument("--chromosomes", default="",
                    help="Comma-separated real-chromosome scaffold names. "
                         "Hits not on these are flagged as DEBRIS. "
                         "Empty = no filtering, all hits treated as chromosome.")
    ap.add_argument("--n_top_cells", type=int, default=10000)
    ap.add_argument("--top_hits", type=int, default=5)
    ap.add_argument("--out_tsv")
    args = ap.parse_args()

    real_chroms = set(args.chromosomes.split(",")) if args.chromosomes else set()
    if real_chroms:
        print(f"=== Real chromosomes: {len(real_chroms)} scaffolds ===", file=sys.stderr)
        for c in sorted(real_chroms): print(f"  {c}", file=sys.stderr)
    else:
        print(f"=== No chromosome filter set ===", file=sys.stderr)

    print(f"\n=== Parsing BLAST: {args.blast} ===", file=sys.stderr)
    marker_to_hits = parse_blast(args.blast, top_n=args.top_hits)
    print(f"  {len(marker_to_hits)} markers with hits", file=sys.stderr)

    print(f"\n=== Loading STARsolo: {args.matrix_dir} ===", file=sys.stderr)
    M, gene_ids, barcodes = load_starsolo(args.matrix_dir)

    print(f"\n=== Identifying real cells ===", file=sys.stderr)
    cell_mask, umi_per_cell = identify_real_cells(M, n_top=args.n_top_cells)

    print(f"\n=== Computing per-marker stats ===", file=sys.stderr)
    rows = compute_marker_stats(M, gene_ids, marker_to_hits, cell_mask, real_chroms)

    print(f"\n=== Pollen marker expression: {args.matrix_dir} ===")
    print(f"Real-cell set: top {cell_mask.sum()} barcodes by UMI count")
    print(f"  Total UMIs in cells: {int(umi_per_cell[cell_mask].sum()):,}")
    print(f"  Median UMIs/cell: {int(np.median(umi_per_cell[cell_mask]))}")
    print()
    print(format_table(rows))

    # Emit debris warnings
    if real_chroms:
        warnings = emit_debris_warnings(rows, marker_to_hits)
        if warnings:
            print()
            print("=" * 80)
            for w in warnings:
                print(w)
        else:
            print()
            print("All markers have at least one hit on a real chromosome.")

    if args.out_tsv:
        import csv
        with open(args.out_tsv, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"\nWrote: {args.out_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
