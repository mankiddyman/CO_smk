#!/bin/bash
# co_sweep.sh -- hapCO_identification.R over a grid of marker_num x block_size,
# for one sample's selected cells; everything else exactly as rule co_calling.
# Output: results/co_sweep/<sample>/mn<N>_bs<B>/per_cell/*_co_pred.txt, with
# PDFs kept for the first 6 cells only (examples). Rerunnable: finished grid
# points carry .done and are skipped.
# Usage: co_sweep.sh SAMPLE THREADS RSCRIPT "MARKER_NUMS" "BLOCK_SIZES"
set -uo pipefail
S=$1; P=$2; RS=$3; MNS=$4; BSS=$5
ROOT=results/co_sweep
CELLS=results/cell_qc/$S/good_cells.tsv
[ -s "$CELLS" ] || { echo "no $CELLS"; exit 1; }
for mn in $MNS; do
  for bs in $BSS; do
    O=$ROOT/$S/mn${mn}_bs${bs}
    [ -f $O/.done ] && { echo "skip $O"; continue; }
    mkdir -p $O/per_cell
    echo "$(date +%H:%M:%S)  $S  marker_num=$mn  block_size=$bs  ($(wc -l < $CELLS) cells)"
    xargs -P $P -I{} "$RS" workflow/scripts/hapCO_identification.R \
        --input results/cell_data/$S/{}.tsv --prefix {} \
        --chrom_map results/cell_data/$S/chrom_map.tsv --outpath $O/per_cell \
        --cell_markers 200 --block_size $bs --marker_num $mn \
        --baseAF 0.3 --windowAF 0.4 --genotype 0.2 < $CELLS > $O/hapco.log 2>&1
    head -6 $CELLS | sed 's/$/_co.pdf/' > $O/keep_pdfs.txt
    find $O/per_cell -name '*_co.pdf' | grep -v -F -f $O/keep_pdfs.txt | xargs -r rm -f
    echo "  $(ls $O/per_cell/*_co_pred.txt 2>/dev/null | wc -l) cells called"
    touch $O/.done
  done
done
