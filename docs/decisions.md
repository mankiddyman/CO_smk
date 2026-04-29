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



---

## 2026-04-28 — Stage 3a: cellsnp MTX → per-cell TSVs

**Decision:** Convert cellsnp-lite output (sparse MTX matrices) to one TSV
per cell, in the format expected by Meng's `hapCO_identification.R`:
`chrom, pos, ref, ref_count, alt, alt_count`.

**Reasoning:**
- Adopting Meng's CO calling algorithm requires per-cell tables. Bridging
  sparse-matrix → per-cell-files at this point is cleaner than rewriting
  the algorithm to consume MTX directly.
- ~6,874 files of ~1k rows each (~25 KB each) is fine for filesystem.
- Skipping cells with <50 markers — none filtered on hap1 (min cell had 463).

**Result on hap1:**
- 6,874 cells written, all retained.
- Median 1,027 markers per cell, max 6,156.
- 293 MB total disk for per-cell TSVs.
- Validated: 3 random (cell, marker) entries all match raw cellsnp output exactly.

**Side outputs:** `chrom_map.tsv` (internal_id ↔ scaffold name ↔ size).
This map is the bridge from Meng's integer-numbered chromosome convention
to our scaffold names. Per-cell files keep scaffold names directly;
chromsome-size lookups in downstream scripts go through the map.



---

## 2026-04-28 — Stage 3b/3c: cell QC and selection

**Decision:** Two-step QC for selecting cells suitable for CO calling.

**Stage 3b — switch_diagnostics:** for every cell, compute per-marker
genotype calls (alt_ratio thresholds 0.2/0.8) and count haplotype switches
per chromosome (NOT across chromosome boundaries — that was a minor flaw
in Meng's bash code). Output: switches.tsv with total/per-chrom marker
and switch counts, plus diagnostic histograms.

**Stage 3c — select_cells:** apply thresholds to switches.tsv to write
good_cells.tsv (one barcode per line) for downstream CO calling.

**Thresholds for Cuscuta epithymum hap1:**
- total markers ≥ 2,500
- per-chromosome markers ≥ 100
- switch rate ≤ 0.09

**Result:** 90 cells pass (1.3% of 6,874 called cells).
Median marker count 3,506; median switch rate 0.073.

**Reasoning:**
- Our switch rate distribution is shifted ~1.5x higher than Meng's R. breviuscula
  data: median 0.152 vs her 0.095. Likely driven by reference bias (median
  ALT ratio 0.39 in HiFi suggests ~20% of true ALT reads systematically
  undercounted; this propagates to scRNA-seq) plus dense markers amplifying
  noise per cell.
- Visual inspection of cells across switch-rate bands showed:
  - sr ≤ 0.08: cleanly callable, clear haplotype blocks after smoothing
  - sr 0.08-0.09: callable with low coverage being limiting; smoothing
    rescues most patterns
  - sr ≥ 0.10: unreliable, smoothing breaks down
- Marker count was a stronger predictor of callability than switch rate
  per se — cells with >3,000 markers were callable even at sr=0.10, while
  sparse cells (< 2,000 markers) were dicey even at sr=0.07.
- Combined filter chosen to optimize for: cells where smoothing recovers
  clean haplotype patterns AND switch rate is low enough that doublets/
  contaminated cells are excluded. Visual spot-check on calibration PDF
  with 12 cells confirmed the threshold lands where callability breaks down.

**Comparison to literature (rough sanity check on n):**
- Shi et al 2019 (pear): 12 cells published a recombination map.
- Castellani et al 2024 (R. breviuscula, our lab's prior paper): 644 cells.
- Dreissig et al 2017 (barley): 24 cells.
- Our 90 cells is comfortably in the publishable range for a single-individual
  recombination landscape.

**Optional secondary tier (not implemented):** ≥2,000 markers, sr ≤ 0.10
gives 151 cells. Could be used for sensitivity analysis if the primary
landscape looks suspiciously sparse.

**Note on what we changed from Meng's pipeline:**
- Meng used sr ≤ 0.13, ≥300 total markers; chromosome-transition switches
  not separated.
- Our thresholds are stricter (sr ≤ 0.09, ≥2,500 markers) because our
  switch-rate distribution is shifted higher.
- We count switches per chromosome separately, then sum, rather than across
  the full sorted file. Avoids spurious "switches" at chromosome boundaries.

---

## 2026-04-29 — Stage 3d: per-cell CO calling

**Decision:** Adapt Meng's `hapCO_identification.R` to call crossovers
per cell on the 90 selected good cells. Aggregate per-cell calls into
landscape diagnostics.

**Implementation:**
- `workflow/scripts/hapCO_identification.R`: Modified Meng's R script to
  read scaffold names from chrom_map.tsv (instead of integer chromosome
  IDs from genome.fai). Otherwise unchanged: same smoothing, same block
  filter, same breakpoint refinement.
- `workflow/scripts/co_aggregate.R`: New script. Reads per-cell
  `_co_pred.txt` files and produces co_summary.txt, co_intervals.bed,
  co_per_cell.tsv, and co_diagnostics.pdf.
- `co_calling` rule: parallelized via `xargs -P {threads}`. One R script
  invocation per cell, ~1-2 sec each.
- `co_aggregate` rule: single job, reads all per-cell outputs, produces
  aggregate diagnostics.

**Parameters chosen for Cuscuta epithymum hap1:**
- `cell_markers: 2500` (matches our cell selection threshold)
- `block_size: 2000000` (Meng's default for Rhynchospora)
- `marker_num: 50` (raised from Meng's default of 8)

**Reasoning for raising marker_num to 50:**

Initial run with Meng's default `-n 8` produced 505 COs (mean 5.61/cell).
Visual inspection of per-cell PDFs flagged ~4 cells with suspicious
short intercalated DCO segments. Examination of the marker_count column
in `_co_pred.txt` files showed these suspicious blocks were supported
by only 9-30 markers, while well-supported COs had blocks with 100-400+
markers.

Distribution of marker_count across all 505 calls:
  - (8,15]:    44 COs  (likely noise — sparse-block artifacts)
  - (15,30]:   41 COs  (mostly noise)
  - (30,50]:   32 COs  (border zone)
  - (50,100]:  60 COs  (mostly real)
  - (100,500]: 304 COs (solid)
  - (500,Inf]:  24 COs (very solid)

No clean bimodal cut in the distribution. Set threshold at 50 markers
as a defensible round number that:
  - Eliminates the suspicious 9-30 marker blocks
  - Retains all 100+ marker blocks (clearly real)
  - Drops 23% of calls (505 → 484... actually higher with full -n 50:
    final result is 484 COs at -n 50)
  - Per-chromosome rate: scaffold_1 dropped from 1.31 to 1.00 COs/cell
    (length-normalized: ~1.0 CO per 10Mb across all 7 chromosomes)

Choice of A vs B (calling parameter vs post-hoc filter): Implemented as
A — change `-n` in CO calling. Means per-cell PDFs reflect filtered
calls directly. Tradeoff: harder to do sensitivity analysis later,
need to re-run calling to test other thresholds.

**Final result for hap1:**
- 484 total COs across 90 cells
- Mean 5.38, median 5 COs per cell
- Per-chromosome counts roughly proportional to chromosome length
  (uniform CO density of ~1.0 per 10 Mb across all 7 scaffolds)
- Long tail: cells with 9-14 COs (n=18). Flagged for investigation.

**Unfinished work:**
- Investigate high-CO tail (cells with ≥9 COs).
- Hap2 pipeline currently broken: hap2 BAM aligned against hap1-derived
  markers VCF returns 0 sites because scaffold names don't match between
  haplotypes. Two paths: (A) skip hap2, hap1 is sufficient; (B) run
  variant calling separately on hap2. Deferred.


---

## 2026-04-29 — Stage 3d revisions: append-mode bug + threshold tuning

**Bug fixed:** The hapCO_identification.R script's append-mode write
(`append = chr_idx > 1`) didn't reset between Snakemake re-runs of the
same cell. Multiple `--forcerun co_calling` invocations triplicated
per-cell `_co_pred.txt` files. Caught when investigating high-CO tail
showed cells with exactly 3 COs on multiple chromosomes — telltale
data duplication pattern.

Fix: added `rm -f {params.out_dir}/*_co_pred.txt
{params.out_dir}/*_co_block_pred.txt {params.out_dir}/*_co.pdf` at the
start of the co_calling rule's shell block. Idempotent re-runs now
produce correct output regardless of file state.

**Threshold revised: marker_num 50 → 30.**

After fixing the duplication bug, fresh CO calls with -n 50 gave 260
total COs / mean 2.89 per cell. Visual inspection of 3 cells (zero,
mid, high) showed the algorithm was over-filtering: clear haplotype
blocks visible by eye but rejected because their marker support fell
in the 30-50 range.

Re-ran with -n 30. Distribution shows real COs are recovered:
  - Total: 350 COs (was 260)
  - Mean: 3.89 / cell (was 2.89)
  - Median: 4 (was 3)
  - 0-CO cells: 1 (was 6)
  - Visual inspection on same 3 cells confirms missed COs now recovered
  - No new spurious DCOs appeared

Reasoning: yesterday's distribution showed the artifact zone was
mostly markers <30 (suspicious DCO blocks were 9, 10, 13, 19, 19, 30
markers). The (30,50] zone was largely real-but-lower-supported COs.

**Final hap1 result:**
- 350 COs across 90 cells (89 with ≥1 CO)
- Mean 3.89 / median 4 COs per cell
- 0.55 COs / chromatid average
- Per-chromosome rates: 0.46-0.74 (length-correlated)
- scaffold_1 highest (115 Mb): 67 COs, 0.74/cell
- scaffold_7 lowest (63 Mb): 41 COs, 0.46/cell

**Hap2 deferred:** SAMPLES filter added in Snakefile to exclude hap2
from `rule all`. Hap2 needs separate variant calling on hap2 reference;
deferred until after hap1 results published.
