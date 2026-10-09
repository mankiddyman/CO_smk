#!/bin/bash
#SBATCH --job-name=split_checks
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/Dparadoxa_split_checks_%j.out
# After the Dparadoxa_split pipeline job: which pieces travel together in the pollen (linkage scan on the
# pieces), and which cells carry a piece twice or not at all (haploidness per piece: an unbalanced
# translocation gamete shows as one piece with both alleles). Both run even if the other fails. Queue with:
#   sbatch --dependency=afterok:<split job id> <this file>
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk || exit 1
PY=/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python
rc=0
$PY workflow/scripts/linkage_scan.py Dparadoxa_split qc/linkage/Dparadoxa_split || { echo "linkage_scan FAILED"; rc=1; }
$PY workflow/scripts/chrom_haploidness.py Dparadoxa_split qc/haplotypes/Dparadoxa_split/per_piece || { echo "chrom_haploidness FAILED"; rc=1; }
exit $rc
