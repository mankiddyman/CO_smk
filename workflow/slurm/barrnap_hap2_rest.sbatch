#!/bin/bash
#SBATCH --job-name=barrnap_hap2
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/barrnap_hap2_rest_%j.out
# rDNA annotation (barrnap, the same tool and settings as rule rrna_blacklist) for the five
# D. paradoxa hap2 chromosomes no sample annotates yet: chr1_hap2 (C) and chr3-6 hap2.
# Writes results/blacklist/hap2_rest/rrna.gff, read by workflow/scripts/rdna_overview.py.
cd /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk
DUAL=/netscratch/dep_mercier/grp_marques/Aaryan/reproducible_phd/results/Drosera_paradoxa/assembly_final/pass2/Drosera_paradoxa_chr.fa
O=results/blacklist/hap2_rest
mkdir -p $O
BARRNAP=$(ls .snakemake/conda/*/bin/barrnap 2>/dev/null | head -n 1)
[ -n "$BARRNAP" ] || { echo "no barrnap found under .snakemake/conda"; exit 1; }
export PATH=$(dirname $BARRNAP):$PATH
/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk/bin/python - $DUAL $O/hap2_rest.fa \
    chr1_hap2 chr3_hap2 chr4_hap2 chr5_hap2 chr6_hap2 << 'ENDPY'
import sys
fa, out, names = sys.argv[1], sys.argv[2], sys.argv[3:]
idx = {l.split("\t")[0]: [int(x) for x in l.split("\t")[1:5]] for l in open(fa + ".fai")}
with open(fa, "rb") as f, open(out, "wb") as o:
    for n in names:
        length, off, lb, lw = idx[n]
        f.seek(off)
        seq = f.read(length + ((length + lb - 1) // lb) * (lw - lb)).replace(b"\n", b"").replace(b"\r", b"")[:length]
        o.write(b">" + n.encode() + b"\n")
        for k in range(0, length, 80):
            o.write(seq[k:k + 80] + b"\n")
        print(n, length, len(seq))
ENDPY
date
barrnap --kingdom euk --threads 8 $O/hap2_rest.fa > $O/rrna.gff 2> $O/barrnap.log
echo "rRNA features: $(grep -vc '^#' $O/rrna.gff)"
date
rm -f $O/hap2_rest.fa
