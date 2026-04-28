# Design decisions

Append-only log. Each entry: date, decision, reasoning, status.

---

## 2026-04-20 — Separate repo, not integrated with core pipeline

**Decision:** Crossover analysis lives in its own repo, coupled to the core
pipeline only via file paths (assembly FASTA, annotation GFF3 on disk).

**Reasoning:**
- Different input lifecycle: core derives assemblies de novo; crossover
  consumes any existing assembly+annotation.
- Core's `species.csv` schema is built around `run_hifiasm == yes`. Crossover
  species don't fit that contract.
- Different tool ecosystem (scRNA, SNP calling) vs (hifiasm, BUSCO).
- Likely standalone publication unit.

**Reconsider if:** ≥2 crossover analyses share substantial rule logic and
duplication starts hurting.

---

## 2026-04-20 — Two rows per library (one per haplotype)

**Decision:** `samples.csv` primary key is `sample_id = {library}_{haplotype}`.
scRNA + HiFi paths are duplicated across rows for the same library.

**Reasoning:**
- Flat schema; every rule gets all its inputs from a single row lookup.
- Premature normalisation (separate libraries / haplotypes / runs tables) is
  the wrong move while the pipeline is still one library.

**Reconsider if:** library count grows past ~5 and duplication causes
inconsistencies.

---

## 2026-04-20 — STARsolo over Cell Ranger

**Decision:** Use STARsolo as the aligner/quantifier in the pipeline.

**Reasoning:**
- Cell Ranger is monolithic, hard to wrap in Snakemake, and its cell-calling
  defaults are mammalian-tuned (called 53,925 cells on a library targeted
  at ~10,000; median UMI 278, median genes 212 — clearly including ambient
  droplets).
- STARsolo gives explicit control over CB/UMI parsing — critical for this
  library where the chemistry call is non-trivial.
- Raw matrix → emptyDrops cell calling as a separate step lets us tune
  parameters for plant scRNA without re-aligning.
- Same underlying aligner (STAR) as Cell Ranger, so mapping rates should be
  comparable; the gain is in cell calling and reproducibility.

**Reconsider if:** STARsolo mapping rate is materially worse (>5%) than
Cell Ranger's 64% on the same data.

---

## 2026-04-20 — Library is 10x 5' v1 on BGI DNBSEQ

**Established by:**
- R1=28 bp, R2=10 bp (linker, not cDNA), R3=100 bp (cDNA), I1=10 bp
- 86% barcode hit rate against `737K-august-2016.txt` (5' v1 whitelist)
- 0.2% hit rate against `3M-5pgex-jan-2023.txt.gz` (5' v2 whitelist)
- Read headers `V35017xxx` confirm BGI DNBSEQ-G400 sequencer

**Implication:** Cell Ranger's `--chemistry=auto` correctly detected
`SC5P-R2`. R1 is 2 bp longer than spec (28 vs 26); Cell Ranger trims
silently. STARsolo needs `--soloBarcodeReadLength 0` to skip the length
check.

**cDNA is on R3, not R2.** Cell Ranger's "R2-only" naming is internal
convention; the actual cDNA file fed to STAR is R3.



---

## 2026-04-21 — Convert GFF3 to GTF before STAR indexing

**Decision:** Add a `gff_to_gtf` rule (gffread) upstream of `star_index`.
`star_index` now consumes GTF, not GFF3.

**Reasoning:**
- First STARsolo run on hap1 produced 1.0% exonic gene mapping
  (`Reads Mapped to Gene: Unique Gene = 0.00998`) vs Cell Ranger's 46.6%
  on the same data.
- Genome mapping was identical (64.3%) and barcode parsing was identical
  (89.2%), so the aligner was doing the same work — the failure was
  exon→transcript→gene linkage at counting time.
- STAR's `--sjdbGTFtagExonParentTranscript Parent` flag handles splice
  junction definition from GFF3, but gene/transcript assignment for
  counting expects GTF-style `transcript_id` / `gene_id` attributes.
  Helixer GFF3 uses `ID` / `Parent` hierarchy (`exon → mRNA → gene`)
  which STAR did not resolve for counting.
- `gffread -T` converts the hierarchy to flat GTF with explicit
  `transcript_id` and `gene_id` attributes on every exon. Standard fix.

**Implementation notes:**
- `gff_to_gtf` writes to `results/annotation/{sample}.gtf` — one GTF
  per sample × haplotype combination. Cheap (<1 min) so no reason to
  cache outside the pipeline.
- `--sjdbGTFtagExonParentTranscript Parent` removed from STAR index
  call; the GTF default behavior is correct.
- Added `gffread` to `workflow/envs/star.yaml`.

**Reconsider if:** future annotation sources produce GTF directly
(skip the conversion) or produce GFF3 that fails gffread validation
(would need upstream cleaning).

---

## 2026-04-21 — Conda env hygiene: one purpose per env

**Decision:** `star.yaml` contains only alignment/indexing tools
(star, samtools, gffread). Analysis libraries (scanpy, pandas) live
in `scanalysis.yaml`. The `star_index` rule was accidentally pointed
at `scanalysis.yaml` and happened to work — fixed.

**Reasoning:** Original attempt combined star + scanpy in one env,
which was unsolvable via conda (python-igraph → glpk → libdeflate
conflict with htslib pinned by STAR). Splitting by purpose keeps
each solve small and each env's dependency surface narrow.

**Reconsider if:** two rules of different purposes need to share
state via the same env (unlikely in a reproducibility-focused pipeline).



---

## 2026-04-21 — Using GFF3 directly with GeneFull counting, accepting gap vs CR

**Decision:** Revert to GFF3 input with `--sjdbGTFtagExonParentTranscript Parent`.
Primary count matrix is `GeneFull` (exonic + intronic). Don't chase CR's 47%
exonic match.

**Reasoning:**
- GFF3 → GTF via gffread actively made counting worse (GeneFull dropped
  32% → 1.8%, Gene stayed at 1%). Added flag we removed was doing useful work.
- CR's GTF structurally identical to ours; not a GTF format issue.
- For crossover inference, we need BAM with CB tags (working) + good-enough
  count matrix for emptyDrops cell calling. 32% GeneFull is sufficient.
- Don't need to match CR's exonic rate exactly; downstream uses BAM alignments.

**Reconsider if:** cell calling produces obviously wrong cell counts, or
downstream SNP calling shows coverage too sparse, trace it to counting gaps.


---

## 2026-04-21 — Accept GeneFull ~32% vs CR's ~50%, use BAM-driven downstream

**Decision:** Use STARsolo with GFF3 + `--sjdbGTFtagExonParentTranscript Parent`.
Primary counting mode is `GeneFull` (exonic + intronic). Accept ~15pp gap
on gene counting vs Cell Ranger's output on the same data.

**Run-by-run empirical record:**

| Config | Gene % | GeneFull % | Genome map % |
|---|---|---|---|
| Cell Ranger 9.0.1 (`SC5P-R2 --include-introns`) | 46.6 | 49.7 | 64.3 |
| STARsolo + GFF3 + `--Parent` | 1.0 | 32.1 | 64.3 |
| STARsolo + gffread GTF (no `--Parent`) | 1.0 | 1.8 | 64.3 |
| STARsolo + GFF3 + `--Parent` (reverted) | 1.0 | 32.1 | 64.3 |

**Reasoning:**
- Genome mapping rate is identical across all configs (64.3%), meaning
  the aligner is doing the same work. The gap is purely in how STARsolo
  attributes reads to genes vs how Cell Ranger does.
- Gffread-converted GTF made counting worse, not better. The `Parent`
  flag is doing real work in STARsolo's internal gene model for this
  annotation style that gffread's output doesn't preserve.
- For crossover inference, the BAM (CB-tagged, 100% of first 100k
  checked) is the actual input to downstream SNP calling. The raw count
  matrix is used only for cell calling (emptyDrops), where 32% gene
  assignment is sufficient.
- Perfect match to Cell Ranger would be nice but not worth chasing. Our
  pipeline is reproducible, parameterizable, Snakemake-native, and the
  alignment is demonstrably correct.

**Reconsider if:** emptyDrops on the raw matrix produces suspiciously
few or many cells, or if downstream per-cell SNP depth is sparse.
Either could trace back to counting gaps we're accepting here.

---

## 2026-04-21 — GTF attribute order matters for STARsolo gene counting

**Decision:** `gff_to_gtf` rule pipes gffread output through awk to reorder
attributes so `gene_id` precedes `transcript_id` in column 9.

**Reasoning:**
- STARsolo 2.7.11b silently collapsed all 20,164 Helixer genes into a
  single `MissingGeneID` bucket in features.tsv when column 9 was
  `transcript_id "..."; gene_id "..."` (gffread default).
- Cell Ranger's `mkref` produces `gene_id "..."; transcript_id "...";`
  and works correctly with STAR.
- Same coordinates, same IDs, same feature counts, only attribute order
  differs. Undocumented STAR parser requirement.

**Lesson for future debugging:** inspect `geneInfo.tab` in the STAR index
and `features.tsv` in the STARsolo output after any annotation change.
Summary.csv percentages can look fine while the matrix structure is broken.

---

## 2026-04-21 — STARsolo strand for 10x 5' chemistry is Reverse

**Decision:** `--soloStrand Reverse` in `starsolo_align`.

**Reasoning:**
- 10x 5' chemistry uses template-switching at the 5' cap; cDNA read (R3)
  aligns antisense to the transcript.
- Empirical test on 1M-read subset, same index, varying `--soloStrand`:

  | Strand | Gene | GeneFull |
  |---|---|---|
  | Forward (default) | 0.9% | 1.6% |
  | Reverse | 41.0% | 45.6% |
  | Unstranded | 41.8% | 46.8% |

- Unstranded marginally higher but double-counts antisense overlaps.
  Reverse preserves real biological stranding for downstream SNP work.
- Cell Ranger handles this via internal chemistry auto-detect; STARsolo
  requires explicit config. Manual is explicit for 10x 3' but ambiguous
  for 5'.

**Validation:** GeneFull mapping matches Cell Ranger's 49.7% within ~4pp.

**Pair with:** GTF attribute-order fix above — both required simultaneously.
Either alone produced "low but plausible" numbers that hid structural failure.


---

## 2026-04-22 — Cell calling: emptyDrops with lower=300 for plant pollen

**Decision:** Use `emptyDrops(lower=300, FDR<=0.001)` as primary cell call,
yielding ~7k cells. Output multiple thresholds for downstream comparison.

**Reasoning:**
- Cuscuta epithymum pollen scRNA shows minimal transcriptional complexity:
  top 20 genes are nearly identical between ambient pool and high-UMI
  cells, meaning emptyDrops cannot distinguish cells from ambient by
  expression profile. Effectively, `lower` IS the threshold.
- Empirical cell counts at varying thresholds:

  | lower | n_cells | median UMI | comment |
  |---|---|---|---|
  | 100 | 60,413 | 268 | far too permissive |
  | 300 | 6,874 | 451 | matches expected ~10k loaded → ~5-6k recovered |
  | 500 | 2,353 | 663 | conservative |
  | 1000 | 329 | 1,382 | far too restrictive |

- Library prep targeted ~10k loaded cells; expected recovery ~5-7k.
  lower=300 brackets this expectation.
- Downstream filters (per-marker depth, per-chromosome switch rate)
  will clean residual ambient/doublet contamination. Cell calling at
  this stage is a coarse first pass.
- For pollen specifically, ambient contamination IS a haplotype
  mixture (50/50 at each marker). Stage 3 switch-rate filtering
  catches contaminated cells via excess apparent heterozygosity.

**Reconsider if:** downstream switch-rate filtering throws out >50% of
called cells, OR if marker coverage per cell is too sparse — then
lower=300 may need adjusting up or down based on what the SNP stage
needs.



---

## 2026-04-24 — Marker filter thresholds

**Decision:** Filter HiFi raw VCF with:
- Biallelic SNPs only
- Heterozygous genotype only
- DP: 10-80
- QUAL: ≥30
- ALT allele ratio: 0.3-0.7

Yields ~4.3M markers from 6.8M het biallelic SNPs (63% retention).

**Reasoning:**
- Diagnostic histograms on raw VCF showed healthy library: modal DP~40x,
  very high QUAL (median 222), ALT ratio centered near 0.5 but with
  reference-bias shoulder (median 0.39, not 0.5).
- DP>80 excludes repeat-region false hets (pileups from paralogs).
- DP<10 excludes low-coverage noise where allele ratio is unreliable.
- QUAL≥30 is a safety net; filter is not discriminating since most
  variants are QUAL>200.
- Ratio 0.3-0.7 retains the main het distribution including left-shifted
  shoulder from reference bias. Tightening to 0.4-0.6 was considered
  but would cut the left side of the real het distribution because of
  the asymmetric reference-bias shift.

**Note on marker count:**
- 4.3M markers is far more than needed for CO resolution (~60k would
  suffice for 10 kb spacing on 600 Mb genome).
- Kept the full set at this stage; downstream can subset (e.g. to
  markers in expressed genes) when Stage 2 pileup performance requires.

**Reconsider if:**
- Stage 2 per-cell pileup is too slow → subset to expressed genes or
  space markers more widely.
- CO calls look overclustered → reference bias may be misleading
  per-cell calls; tighten to 0.4-0.6 and re-examine.




---

## 2026-04-23 — Parallel HiFi variant calling

**Decision:** Split variant calling by scaffold. Three rules:
- `hifi_align` — one BAM per sample (minimap2, ~28 min for hap1)
- `hifi_variants_per_chrom` — one VCF per scaffold (bcftools mpileup + call, parallel)
- `hifi_variants_merge` — concatenate scaffolds into final VCF

**Reasoning:**
- First serial attempt ran 22h before being killed. `bcftools mpileup --threads`
  only parallelizes BCF compression, not the pileup scan itself.
- Parallel version ran in ~10h wall time, bounded by largest scaffold
  (scaffold_1 at 115 Mb ~ 10h on its own).
- Lesson: for big bcftools jobs, parallelize at the Snakemake level by region.

**Future improvement:** could split largest scaffolds into sub-regions
to break the scaffold_1 ceiling. Not implemented — 10h is acceptable
and adds complexity. Backlogged.

---

## 2026-04-24 — Marker filter thresholds

**Decision:** Filter raw HiFi VCF to biallelic het SNPs with:
- DP: 10-80
- QUAL: ≥30
- ALT allele ratio: 0.3-0.7

**Result:** 4,248,719 markers from 7,305,567 raw variants (58% retention).
Per-chrom marker counts scale roughly with chromosome length.

**Reasoning:**
- Diagnostic histograms on hap1 showed healthy library: modal DP ≈40x,
  QUAL almost always >200, ALT ratio centered near 0.5.
- DP<10 drops noisy low-coverage calls; DP>80 excludes repeat/paralog pileups.
- QUAL≥30 is a safety net — most variants pass regardless.
- Ratio 0.3-0.7 retains the heterozygous distribution including the left
  shoulder caused by reference bias (median observed at 0.39, not 0.5).
  Tightening to 0.4-0.6 was considered but would cut genuine hets due to
  the asymmetric shift.

**Reference bias context:** *C. epithymum* is highly heterozygous
(~1 het per 90bp). On a HiFi read averaging ~170 het sites, this compounds
the alignment-mismatch penalty for alt-carrying reads. Reads with the alt
allele are systematically slightly worse-aligned and undercounted.
Need to re-apply this correction at per-cell allele counting (Stage 2)
or at downstream CO calling (Stage 3/4).

**Marker count is much larger than needed:** 4.2M vs ~60k needed for
~10 kb CO resolution on a 600 Mb genome. Kept the full set; subset
later if performance requires (e.g. limit to expressed-gene markers).

---

## 2026-04-25 — Stage 2: per-cell allele counting

**Decision:** Use cellsnp-lite (not pysam, not custom bcftools post-processing,
not per-cell BAM demultiplexing) for per-cell allele counts.

**Reasoning:**
- Demultiplexing scRNA BAM into 6,874 per-cell BAMs would inflate disk
  ~3-5x, hammer filesystem with file-count, and be ~10x slower.
- cellsnp-lite walks the BAM once with parallel-by-region; built for exactly
  this case (BAM with CB tags + SNP list → sparse cell × site matrices).
- Validated on a scaffold_1:1-50Mb subset (434k markers, 25 min): UMI
  dedup confirmed (18 reads → 5 unique molecules at one spot-check site),
  AD/DP/OTH outputs match by-hand counts.

**Settings:**
- `--minMAPQ 20`, `--minLEN 30` — reject low-quality/short alignments
- `--cellTAG CB --UMItag UB` — STARsolo's corrected cell + UMI tags
- `--minCOUNT 1` — accept any (cell, site) pair with ≥1 UMI;
  per-cell DP filtering happens downstream
- `--minMAF 0` — markers are already curated, don't refilter

**Result on hap1:**
- 694,537 of 4.25M markers had any scRNA coverage (16%) — consistent with
  scRNA only seeing expressed gene regions in pollen.
- 8.26M non-zero (cell, site) entries in DP matrix.
- 3.69M non-zero entries in AD matrix (45% of DP entries had ≥1 ALT UMI).
- 28k non-zero entries in OTH matrix (0.35% — sequencing errors).
- Wall time: ~1.5h on full markers with 16 threads.

**Outputs:** `results/snps/{sample}/cellSNP.{tag.DP,tag.AD,tag.OTH}.mtx`,
`cellSNP.samples.tsv`, `cellSNP.base.vcf.gz` (sites with any coverage).

**Side fix:** `starsolo_align` rule now also outputs `Aligned.sortedByCoord.out.bam.bai`
(samtools index after STAR), since cellsnp-lite needs an indexed BAM.
