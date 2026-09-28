#!/bin/bash
# recall_panel.sh -- re-call ONLY the review panel's cells with other hapCO settings, then render them.
# Usage: recall_panel.sh SAMPLE SEED LABEL "hapCO args"
#   e.g. recall_panel.sh Dbinata_hap1 1 mn6_term24 "--marker_num 6 --block_size 2000000 --terminal_marker_num 24"
# Unspecified args take hapCO's defaults (block 2 Mb, marker_num 8, baseAF 0.3, windowAF 0.4,
# genotype 0.2, terminal off) -- pass everything you want to differ from those.
set -uo pipefail
S=$1; SEED=$2; LABEL=$3; ARGS=$4
RS=$PWD/.snakemake/conda/738e17a23d1043909c1a53d6c039531c_/bin/Rscript
PY=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python
PANEL=qc/review/$S/panel_seed$SEED.txt
[ -s "$PANEL" ] || { echo "no panel yet -- run: review_panel.py $S 50 $SEED"; exit 1; }
O=qc/review/$S/runs/$LABEL/per_cell; mkdir -p "$O"
echo "$ARGS" > qc/review/$S/runs/$LABEL/args.txt
xargs -P 16 -I{} "$RS" workflow/scripts/hapCO_identification.R \
    --input results/cell_data/$S/{}.tsv --prefix {} \
    --chrom_map results/cell_data/$S/chrom_map.tsv --outpath "$O" --cell_markers 200 $ARGS \
    < "$PANEL" > qc/review/$S/runs/$LABEL/hapco.log 2>&1
"$PY" workflow/scripts/review_panel.py "$S" 50 "$SEED" "$O" "$LABEL"
