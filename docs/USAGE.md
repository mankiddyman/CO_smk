# CO_smk — Usage Guide

Species-agnostic Snakemake pipeline to infer meiotic crossover (CO) landscapes from single-cell RNA-seq of pollen. Given a genome assembly, its annotation, HiFi reads, and 10x pollen scRNA-seq, it calls per-cell haplotype switches (crossovers) and produces recombination landscapes, interference statistics, and a composite karyoplot.

Repo: github.com/mankiddyman/CO_smk
Working dir on HPC: `/netscratch/dep_mercier/grp_marques/Aaryan/CO_smk`

---

## What the pipeline does (one line per stage)

```
Stage 0  starsolo_align        align scRNA to genome, get cell x gene matrix + BAM (chemistry-aware)
Stage 1  markers               HiFi -> per-chrom het SNP markers (minimap2, bcftools) -> filter
Stage 2  cell_calling          knee/EmptyDrops -> called barcodes
Stage 3  cell_snp_counting     per-cell allele counts at marker sites (cellsnp-lite)
Stage 3b cellsnp_to_per_cell   sparse matrices -> per-cell haplotype TSVs
Stage 3c switch_diagnostics    per-cell switch-rate QC
Stage 3d select_cells          apply QC thresholds -> good cells
Stage 4  co_calling            sliding-block CO detection per cell (hapCO_identification.R)
Stage 4b co_aggregate          per-cell CO calls -> landscape summary + intervals
Stage 5  recombination_landscape   landscape, Marey maps, CO classes, CoC + Gamma interference
Stage 6  composite_landscape   karyoploteR figure: ideograms + CO + gene/SNP/LTR/satellite tracks + SyRI ribbons
```

The DAG is in `.github/assets/rulegraph.svg` (regenerate with the `generate_dag` rule).

---

## Tool stack

Everything runs in per-rule conda envs under `workflow/envs/` (Snakemake `--use-conda` builds them). Key tools by env:

- **star.yaml** — STAR/STARsolo (alignment), gffread (GFF->GTF), samtools
- **markers.yaml** — minimap2, bcftools, samtools
- **cell_calling.yaml** — R + DropletUtils (EmptyDrops/knee)
- **cellsnp.yaml** — cellsnp-lite, python (scipy/numpy for matrix conversion)
- **cell_qc.yaml** — R (switch diagnostics, cell selection, CO calling, aggregation, landscape)
- **karyoploter.yaml** — R + karyoploteR (composite figure)

The Snakemake driver itself runs from the `smk` micromamba env at
`/netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk`.

> Note: `python` on this system is aliased to ipython, which breaks `argparse` with `--` flags.
> Always call helper scripts with `python3`, not `python`.

---

## Inputs: what you need per species

To run a species end-to-end you need FOUR things:

1. **Genome assembly** (`.fa`) + its `.fai` index (`samtools faidx`)
2. **Annotation** (`.gff3`) with `gene`/`mRNA`/`exon`/`CDS` features (Helixer or BRAKER output both work — must contain `gene` features)
3. **HiFi reads** (`.fasta` / `.fasta.gz`) — used to call heterozygous SNP markers that define the two haplotypes
4. **10x pollen scRNA-seq** — R1 (barcode+UMI) + cDNA read (R2 for 3' chemistry, R3 for 5' chemistry)

Optional (only needed for the Stage 6 composite figure):
5. **LTR repeat GFF** (e.g. DANTE_LTR output)
6. **Satellite repeat GFF** (e.g. TideCluster / TRASH output)
7. **SyRI output dir** (`syri.out` + `mapids.txt`) — for hap1<->hap2 synteny ribbons
8. **hap2 assembly** `.fai` — the composite rule expects `hap2.fasta.fai` next to the hap1 assembly

You can run Stages 0-5 without 5-8. The composite figure (Stage 6) needs them.

### Where inputs live

There is no fixed inputs directory inside the repo. All input paths are **absolute paths recorded in `config/samples.csv`** (one row per sample). This is deliberate — the pipeline points at data wherever it lives (genome in the Drosera/Cuscuta project trees, scRNA in `/biodata`, whitelists in the Cell Ranger install). Nothing gets copied into the repo.

Reference genome/annotation are currently pulled from each project's own results tree (e.g. `.../Drosera/pan_drosera/genome/`, `.../Cuscuta/epithymum/results/assemblies/...`).

---

## config/samples.csv — column reference

One row per sample. `sample_id` is the wildcard used throughout `results/`.
Columns (CSV, not TSV — semicolon-separate multiple FASTQs within a field if ever needed):

| Column | Meaning | Required for |
|---|---|---|
| `sample_id` | unique ID, used as `{sample}` wildcard | always |
| `species` | groups samples; used to find the marker-reference sample | always |
| `haplotype` | e.g. `hap1`, `hap2`, `primary` | bookkeeping |
| `assembly_fasta` | absolute path to genome `.fa` (must have `.fai` alongside) | always |
| `annotation_gff` | absolute path to `.gff3` | always |
| `repeats_ltr_gff` | LTR repeat GFF | composite only |
| `repeats_satellite_gff` | satellite repeat GFF | composite only |
| `syri_dir` | dir containing `syri.out` + `mapids.txt` | composite only |
| `scrna_r1` | barcode+UMI read (R1) | always |
| `scrna_r2` | cDNA read for 3' chemistry | 3' v4 samples |
| `scrna_i1` | sample index (unused by pipeline) | — |
| `scrna_r3` | cDNA read for 5' chemistry | 5' v1 samples |
| `hifi_reads` | HiFi reads for marker calling | marker reference sample |
| `centromere` | e.g. `holocentric` / `monocentric` | bookkeeping |
| `exp_ploidy` | expected ploidy | bookkeeping |
| `chr_number_2n` | 2n chromosome count | bookkeeping |
| `chemistry` | `10x_5p_v1` or `10x_3p_v4` — drives alignment branching | always |
| `cdna_read` | which read holds cDNA (`R2`/`R3`) — documentation | always |
| `notes` | free text (record library provenance, kit, date, who confirmed chemistry) | always |
| `is_marker_reference` | `yes`/`no` — is this the sample whose HiFi markers define the haplotypes for its species | always |

**Fill empty fields with nothing (blank between commas).** Blank optional fields are fine as long as you don't request an output that needs them.

### The `is_marker_reference` mechanism

Markers are called once per species from the marker-reference sample's HiFi reads, then reused for all scRNA samples of that species. `marker_sample_for(sample_id)` looks up the row where `species` matches AND `is_marker_reference == yes`. So exactly one row per species should have `is_marker_reference = yes` (usually the hap1 assembly with HiFi reads).

### Chemistry branching (Stage 0)

`starsolo_align` branches on the `chemistry` column:

| chemistry | cDNA read | barcode read | CB len | UMI len | whitelist (config key) | strand |
|---|---|---|---|---|---|---|
| `10x_5p_v1` | R3 | R1 | 16 | 10 | `whitelist_5p_v1` | Reverse |
| `10x_3p_v4` | R2 | R1 | 16 | 12 | `whitelist_3p_v4` | Reverse |

To add a new chemistry, extend `starsolo_inputs_for()` and `starsolo_params_for()` in `workflow/Snakefile` and add the whitelist path to `config/config.yaml`.

---

## config/config.yaml — the knobs

Current `config/config.yaml`:

```yaml
# Barcode whitelists (from the Cell Ranger install)
whitelist_5p_v1: /opt/share/software/packages/cellranger-9.0.1/bin/lib/python/cellranger/barcodes/737K-august-2016.txt
whitelist_3p_v4: /netscratch/dep_mercier/grp_marques/Aaryan/CO_smk/resources/whitelists/3M-3pgex-may-2023_TRU.txt

# Stage 1 — marker filtering (het SNPs from HiFi)
marker_filters:
  min_dp: 10
  max_dp: 80
  min_qual: 30
  min_alt_ratio: 0.3
  max_alt_ratio: 0.7

# Stage 3 — cell QC thresholds
cell_qc:
  ref_low: 0.2            # alt_ratio < this -> REF (genotype 0)
  alt_high: 0.8           # alt_ratio > this -> ALT (genotype 1)
  min_total_markers: 2500
  min_per_chrom_markers: 100
  max_switch_rate: 0.09

# Stage 3d — CO calling (sliding-block)
co_calling:
  cell_markers: 2500      # match upstream cell selection threshold
  block_size: 2000000     # 2 Mb min block size
  marker_num: 23          # min markers per block (23 = value used for all reported Cuscuta figures)

# Stage 4 — recombination landscape plotting
landscape:
  bin_mb: 5               # window size for sliding-window CO rate (Mb)
  coc_interval_mb: 5      # interval size for CoC binning (Mb)
  n_boot: 500             # bootstrap iterations for CIs

# Stage 5 — composite landscape figure
composite:
  win_mb: 5               # sliding window size for all density tracks (Mb)
  step_mb: 0.5            # sliding step (Mb)
```

> Stage numbering above follows the config file's own comments (landscape=Stage 4,
> composite=Stage 5). The stage list at the top of this doc splits things a bit finer
> (landscape=5, composite=6); they refer to the same rules.
>
> **`marker_num` provenance note:** the committed value is **23**, and this is the value
> that produced the reported Cuscuta figures (373 COs, 90 cells, Gamma nu=2.40, CI [1.61, 4.54]).
> Verified against `co_calling_params.txt`, `co_summary.txt`, and `gamma_interference_summary.tsv`
> from the May 4 run, which all agree. Git history shows the line only ever moved 8 -> 50 -> 30 -> 23;
> a "20" mentioned in an old config comment was never actually committed or run. **23 is canonical.**

### Parameters most likely to need per-species tuning

- **`co_calling.marker_num`** — the big one. Min markers to trust a haplotype block. Too high = miss real (esp. close double) crossovers; too low = call noise as COs. History on Cuscuta: 8 -> 50 -> 30 -> 23 (committed and used for all reported figures). Lower recovers close DCOs but can inflate CO count. Directly affects the Gamma nu / interference estimate. The reported Cuscuta numbers (373 COs, nu=2.40) are the marker_num=23 results.
- **`cell_qc.min_total_markers`** and **`co_calling.cell_markers`** — scale with heterozygosity + sequencing depth. Keep these two equal.
- **`cell_qc.max_switch_rate`** — noise ceiling per cell. 0.09 worked for Cuscuta.

Everything else (bin sizes, bootstrap, filter DP ranges) is usually fine across species.

> These are currently GLOBAL, not per-species. If two species need different `marker_num`,
> either run them in separate invocations with an edited config, or (future) add per-species
> override columns to samples.csv. Don't add that complexity until a second species actually needs it.

---

## Running the pipeline

Driver env: `micromamba activate /netscratch/dep_mercier/grp_marques/Aaryan/micromamba_envs/smk`

### Dry run first, always

```bash
snakemake -n --cores 16 --use-conda --rerun-triggers mtime <target>
```

`--rerun-triggers mtime` avoids spurious reruns from conda-env or param churn — use it routinely.

### Full pipeline for everything in `rule all`

```bash
snakemake --cores 16 --use-conda --rerun-triggers mtime
```

### One sample, one target (recommended while iterating)

```bash
# Just alignment (QC of mapping rate for a new/shallow library)
snakemake --cores 16 --use-conda \
  results/starsolo/<sample_id>/Solo.out/GeneFull/raw/matrix.mtx

# Full landscape for one sample
snakemake --cores 16 --use-conda \
  results/landscape/<sample_id>/01_landscape.pdf

# Composite figure for one sample (needs repeats + SyRI + hap2.fai)
snakemake --cores 16 --use-conda \
  results/composite/<sample_id>/composite.pdf
```

### Building conda envs ahead of time

```bash
snakemake --cores 16 --use-conda --conda-create-envs-only
```

Do this after moving the repo or wiping `.snakemake/conda` — gets the ~5-10 min env build out of the way before a real run.

---

## Adding a NEW species: checklist

The whole point of this doc. To add species X for full CO landscape:

1. **Confirm inputs exist**: genome `.fa` (+`.fai`), `.gff3` (has `gene` features), HiFi reads, 10x pollen scRNA (R1 + cDNA read).
   - `samtools faidx genome.fa` if no `.fai`.
   - Check GFF feature types: `awk '!/^#/{print $3}' X.gff3 | sort -u` — must include `gene`, `mRNA`, `exon`, `CDS`.
   - Check scaffold naming matches between genome and GFF (watch for a leading-underscore mismatch: genome `scaffold_1` vs gffread proteome `_scaffold_1`). Real chromosomes vs debris contigs: note the cutoff (e.g. Dbinata <=32, Dparadoxa <=14).

2. **Confirm chemistry**: inspect R1 length + structure. 3' v4 = R1 150bp with poly-T at ~pos 29, 16bp CB + 12bp UMI, cDNA on R2. 5' v1 = cDNA on R3, 10bp UMI. Record kit + who confirmed in `notes`.

3. **Add whitelist** to `config/config.yaml` if it's a new chemistry (decompress the Cell Ranger `.txt.gz` to `resources/whitelists/`, plain text not gzip).

4. **Add rows to `config/samples.csv`**:
   - One marker-reference row (`is_marker_reference = yes`) with `hifi_reads` populated — usually the hap1 assembly.
   - One row per scRNA library, pointing at the same `assembly_fasta`/`annotation_gff`, with `chemistry` set and `is_marker_reference = no`.

5. **QC-only first if the library is shallow** (like the Drosera QC): add the sample to `DROSERA_QC_ONLY`-style exclusion in `workflow/Snakefile` so `rule all` skips downstream, then run just the starsolo target and check `Log.final.out` mapping rate + a knee plot before committing to full analysis or ordering deep sequencing.

6. **Dry run** the landscape target, read the plan, confirm only the expected rules fire.

7. **Run markers first** (slow, ~hours; per-chrom parallel): `results/markers/<ref_sample>/markers.tsv`. Inspect `filter_summary.txt` for per-chromosome marker counts.

8. **Then cells -> CO -> landscape.** Tune `marker_num` if the CO count / interference looks off; rerun from `co_calling` onward.

9. **Composite** last, once repeats + SyRI are available.

---

## Gotchas (learned the hard way)

- **`python` is ipython** here — use `python3` for helper scripts, or IPython eats `--flags`.
- **Scaffold naming**: gffread prepends `_` to scaffold names in proteome/gene IDs (`_scaffold_1_...`) while the genome FASTA uses `scaffold_1`. Bites any join between annotation-derived IDs and genome coordinates.
- **Annotation completeness**: if a GFF was accidentally run on debris contigs only (happened with Dbinata), genes on real chromosomes won't be in the count matrix. STARsolo can only count annotated features — TBLASTN against the genome bypasses this for spot-checks.
- **Multi-haplotype reference depresses mapping rate**: aligning to hap1+hap2+debris splits allele-shared reads into multi-mapping / "too many loci". Report total alignment rate, not just uniquely-mapped %, and note the reference composition. (Cuscuta was hap1-only -> not comparable head-to-head.)
- **`--rerun-triggers mtime`** on every invocation to avoid bogus reruns.
- **cellsnp uses CB/UB tags** — the STARsolo BAM must carry `CB`/`UB` (it does via `--outSAMattributes`). Don't strip them.
- **All ggsave uses `device=cairo_pdf`** for UTF-8 safety; keep plot text ASCII where possible (nu not the greek letter, +/- not the symbol).
- **Commit workflow**: `pull --rebase --strategy-option=ours`, regen DAG, append to `docs/decisions.md`, commit, push; resolve SVG conflicts with `--strategy-option=ours`.

---

## Outputs map (per sample)

```
results/
  starsolo/<s>/         Log.final.out (mapping rate), Solo.out/GeneFull/raw/ (matrix), BAM
  markers/<ref>/        markers.vcf.gz, markers.tsv, filter_summary.txt, diagnostic_histograms.pdf
  cells/<s>/            barcodes_called.tsv, knee_plot.pdf
  snps/<s>/             cellSNP.tag.{DP,AD,OTH}.mtx
  cell_data/<s>/        per-cell haplotype TSVs, chrom_map.tsv
  cell_qc/<s>/          switches.tsv, switch_diagnostics.pdf, good_cells.tsv
  crossovers/<s>/       co_intervals.bed, co_per_cell.tsv, co_summary.txt, co_diagnostics.pdf
  landscape/<s>/        01_landscape.pdf, 02_marey_maps.pdf, 03_co_class_proportions.pdf,
                        04_coc.pdf, 05_gamma_interference.pdf, *_summary.tsv
  composite/<s>/        composite.pdf
```
