#!/bin/bash
#SBATCH --job-name=CB_checks
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/Dparadoxa_CB_checks_%j.out
# After the Dparadoxa_CB pipeline job: linkage scan on C + B, then crossovers std vs CB and
# transmission ratios. Queue with: sbatch --dependency=afterok:<CB job id> <this file>
set -euo pipefail
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
PY=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python
$PY workflow/scripts/linkage_scan.py Dparadoxa_CB qc/linkage/Dparadoxa_CB
$PY workflow/scripts/paradoxa_cb_check.py qc/linkage/Dparadoxa_CB/cb_check.txt
