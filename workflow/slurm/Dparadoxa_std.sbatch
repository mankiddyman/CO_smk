#!/bin/bash
#SBATCH --job-name=paradoxa_std
#SBATCH --cpus-per-task=24
#SBATCH --mem=200G
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/run_Dparadoxa_std_%j.out
# D. paradoxa on the mapping reference hap1 + chr2_hap2 (config mapping_reference):
# the full chain from read mapping to the landscape, submitted 2026-09-28.
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/snakemake \
    results/landscape/Dparadoxa_std/01_landscape.pdf \
    --cores 24 --resources mem_mb=190000 --use-conda --rerun-triggers mtime params --rerun-incomplete --keep-going
