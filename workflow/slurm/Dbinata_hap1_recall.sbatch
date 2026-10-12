#!/bin/bash
#SBATCH --job-name=binata_recall
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/run_Dbinata_hap1_recall_%j.out
# Dbinata_hap1 from the molecule tables on, with the settings stamped 2026-10-12 (config
# co_calling.Drosera_binata: marker_num 8, block_size 1 Mb, terminal_marker_num 4, mask_classes
# doubled/missing/mixed) and chr6_hap1:1-70000 left out (config/exclude_regions/Dbinata_hap1.csv):
# molecule tables, hapCO per cell (new per-cell PDFs), dosage masks, aggregate, landscape. Only if that
# worked: the after picture (calling_overview.py), the Python port checked against the new R calls cell by
# cell (co_sweep.py --check, which also refreshes the sweep on the new tables), and every cell's hapCO page
# in one PDF (if pdfunite is there).
# --use-conda needs conda >= 24.7.1 first on PATH: the directory in CO_CONDA_DIR (set by the block that
# submits this), else the smk env's own.
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
ENV=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin
S=Dbinata_hap1
export PATH="${CO_CONDA_DIR:-$ENV}:$ENV:$PATH"
echo "conda used: $(type -P conda) ($(conda --version 2>&1))"
rm -f results/crossovers/$S/all_cells_co.pdf
$ENV/snakemake results/crossovers/$S/co_diagnostics.pdf results/landscape/$S/01_landscape.pdf \
    --cores 16 --resources mem_mb=60000 --use-conda --rerun-triggers mtime params --rerun-incomplete --keep-going
rc=$?
echo "snakemake exit status $rc"
if [ $rc -eq 0 ]; then
    echo "molecule tables: $(tail -n 2 logs/cell_data_molecules/$S.log | tr '\n' ' ')"
    $ENV/python workflow/scripts/calling_overview.py $S > qc/calling_overview/overview_${S}_after_$(date +%F).txt 2>&1
    $ENV/python workflow/scripts/co_sweep.py $S --check --threads 16 > logs/${S}_sweep_check.log 2>&1
    echo "co_sweep.py --check exit status $? (0 = the Python port calls every cell exactly as the R did; logs/${S}_sweep_check.log)"
    if command -v pdfunite > /dev/null; then
        pdfunite results/crossovers/$S/per_cell/*_co.pdf results/crossovers/$S/all_cells_co.pdf
    fi
fi
exit $rc
