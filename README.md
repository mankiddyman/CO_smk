cat > README.md <<'EOF'
# Crossover landscape inference from scRNA-seq

Snakemake pipeline to infer meiotic crossover landscapes from single-cell
RNA-seq of pollen nuclei, given a phased genome assembly + annotation.

Currently configured for *Cuscuta epithymum* (sample 6303_B) across both
haplotypes (hap1, hap2) of the HapHiC-scaffolded assembly.

## Inputs (per haplotype)

- Assembly FASTA (chromosome-scaffolded)
- Gene annotation GFF3 (Helixer)
- scRNA-seq reads (10x Genomics 5' v1, sequenced on BGI DNBSEQ-G400)
  - R1 = cell barcode (16) + UMI (10), 28 bp
  - R3 = cDNA, 100 bp
  - R2 (10 bp linker) and I1 are not used

Inputs declared in `config/samples.csv`. One row per `sample_id` =
`{library}_{haplotype}`.

## Pipeline (current state)

1. `star_index` — build STAR index per haplotype from FASTA + GFF3
2. `starsolo_align` — STARsolo alignment with chemistry-correct CB/UMI parsing,
   raw count matrix (no cell filtering — done downstream)
3. `starsolo_qc` — extract mapping/barcode metrics for comparison

Planned: emptyDrops cell calling, SNP calling per cell, crossover inference.

## Running

```bash
# Dry-run
snakemake -n --cores 16 --use-conda

# One haplotype
snakemake --cores 16 --use-conda --rerun-incomplete \
  results/starsolo/6303_B_hap1/.done

# Everything
snakemake --cores 16 --use-conda --rerun-incomplete
```

On HPC, request ≥16 CPU and ≥100 GB RAM per align job.

## Repo layout

config/
config.yaml       # paths to whitelist etc
samples.csv       # sample sheet (one row per sample × haplotype)
workflow/
Snakefile
envs/
star.yaml
scripts/
starsolo_qc.py
docs/
decisions.md      # why STARsolo not Cell Ranger, etc
results/            # gitignored
logs/               # gitignored


## Related

- Core comparative-genomics pipeline (assembly → annotation → comparative):
  https://github.com/mankiddyman/reproducible_phd
- This repo consumes assembly + annotation as inputs; it does not derive them.

## Status

Active development. Pipeline runs end-to-end as of YYYY-MM-DD; no scientific
results yet.
EOF
