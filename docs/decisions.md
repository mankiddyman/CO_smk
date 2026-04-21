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
