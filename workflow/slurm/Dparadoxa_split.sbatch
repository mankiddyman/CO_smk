#!/bin/bash
#SBATCH --job-name=paradoxa_split
#SBATCH --cpus-per-task=48
#SBATCH --mem=400G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/run_Dparadoxa_split_%j.out
# D. paradoxa on Dparadoxa_CB cut at its three translocation breakpoints (config mapping_reference
# Dparadoxa_split: pieces from config/split_plans/Dparadoxa_split.csv): the full chain from read mapping
# to the landscape, on a reference in which no join between translocated pieces is assumed. The
# reference itself is built before this job (rule mapping_reference), because the chromosome list
# the DAG is built from is read from its index.
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/snakemake \
    results/landscape/Dparadoxa_split/01_landscape.pdf \
    --cores 48 --resources mem_mb=390000 --use-conda --rerun-triggers mtime params --rerun-incomplete --keep-going
