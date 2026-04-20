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
