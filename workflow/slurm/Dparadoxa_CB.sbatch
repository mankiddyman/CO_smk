#!/bin/bash
#SBATCH --job-name=paradoxa_CB
#SBATCH --cpus-per-task=48
#SBATCH --mem=400G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/run_Dparadoxa_CB_%j.out
# D. paradoxa on the crossover reference joined as the pollen inherit it (config mapping_reference
# Dparadoxa_CB: chr1_hap2 + chr2_hap1 + chr3-6 hap1): the full chain from read mapping to the landscape.
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/snakemake \
    results/landscape/Dparadoxa_CB/01_landscape.pdf \
    --cores 48 --resources mem_mb=390000 --use-conda --rerun-triggers mtime params --rerun-incomplete --keep-going
