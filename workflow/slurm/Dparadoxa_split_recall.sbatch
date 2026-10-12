#!/bin/bash
#SBATCH --job-name=paradoxa_split_recall
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/run_Dparadoxa_split_recall_%j.out
# Dparadoxa_split from crossover calling on, with the settings stamped 2026-10-12 (config
# co_calling.Drosera_paradoxa: marker_num 6, terminal_marker_num 4, mask_classes doubled/missing/mixed):
# hapCO per cell (new per-cell PDFs), dosage masks, aggregate, landscape. Only if that worked: the after
# picture (calling_overview.py) and every cell's hapCO page in one PDF (if pdfunite is there).
# --use-conda needs conda >= 24.7.1 first on PATH: the directory in CO_CONDA_DIR (set by the block that
# submits this), else the smk env's own. Job 1078173 failed on the base conda 22.11.1.
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin
S=Dparadoxa_split
export PATH="${CO_CONDA_DIR:-$ENV}:$ENV:$PATH"
echo "conda used: $(type -P conda) ($(conda --version 2>&1))"
rm -f results/crossovers/$S/all_cells_co.pdf
$ENV/snakemake results/crossovers/$S/co_diagnostics.pdf results/landscape/$S/01_landscape.pdf \
    --cores 16 --resources mem_mb=60000 --use-conda --rerun-triggers mtime params --rerun-incomplete --keep-going
rc=$?
echo "snakemake exit status $rc"
if [ $rc -eq 0 ]; then
    $ENV/python workflow/scripts/calling_overview.py $S > qc/calling_overview/overview_${S}_after_$(date +%F).txt 2>&1
    if command -v pdfunite > /dev/null; then
        pdfunite results/crossovers/$S/per_cell/*_co.pdf results/crossovers/$S/all_cells_co.pdf
    fi
fi
exit $rc
