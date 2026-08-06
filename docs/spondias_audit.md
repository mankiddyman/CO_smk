# Spondias tuberosa CO landscape — audit

**Revised 2026-08-06.** Supersedes the 2026-08-04 version. Most items that were
open questions are now answered; what remains open is flagged as such.

**Headline result:** 153 cells, 1740 crossovers, **1137 cM over 284 Mb
(4.0 cM/Mb)**, 16 chromosomes, all above the 50 cM obligate-CO floor.

**Repo:** `/netscratch/dep_mercier/grp_marques/Aaryan/CO_smk`
**Sample:** `Spondias_hap1` · **config species key:** `Spondias_tuberosa`

---

## The three things to say out loud

If you remember nothing else from this audit:

1. **1137 cM is a LOWER BOUND.** CO detection is marker-limited and has not
   saturated. True map length plausibly 1300–1900 cM.
2. **The landscape SHAPE is sound.** Marker spacing is uniform genome-wide, so
   under-detection is uniform along chromosomes. Relative CO density is the
   reportable result.
3. **~60% of the library is rRNA** from a single 49 kb array. That is a
   library-prep issue affecting the whole 6178A batch, and it is the single
   biggest available improvement for future samples.

---

## §0 Parameter table — every value, why, and the evidence

Every threshold below was swept. Numbers are recorded inline in
`config/config.yaml` next to each parameter.

| Stage | Parameter | Value | Basis |
|---|---|---|---|
| markers | `min_dp` | 10 | rejects 0.1%, not doing work |
| markers | `max_dp` | **150** | het-SNP depth median 91; 80 sat *below* the peak and rejected 72% of real sites. Yield plateaus from ~140. Gives 0.53% het vs 0.62% published |
| markers | `min_qual` | 30 | rejects 2.4% |
| markers | `min/max_alt_ratio` | 0.3 / 0.7 | rejects 18.5% — the largest single filter. Cut sits at the edge of a broad low-AF shelf |
| cells | `default_lower` | 300 | **swept 100/300/500/1000 → result IDENTICAL.** Inconsequential |
| cell QC | `min_total_markers` | 800 | knee in the per-cell marker distribution |
| cell QC | `min_per_chrom_markers` | 30 | swept 15–40; ≤25 puts a 2nd chromosome below the floor |
| cell QC | `max_switch_rate` | **0.12** | swept 0.05–0.12; COs/cell FLAT (11.27–11.40), only 0.12 has zero floor violations |
| CO call | `block_size` | 1 Mb | ≥1.5 Mb breaks the floor; ≤250 kb falls below `marker_num` at observed density. **Report sensitivity: 1307 cM @500 kb, 981 @1.5 Mb** |
| CO call | `marker_num` | 4 | swept 3–8 → map varies 3.2%. Data-driven |
| CO call | `base_af` | 0.3 | 0.2 breaks the floor (3 chromosomes); 0.3 vs 0.4 differ 0.4% |
| CO call | `window_af` | 0.4 | **INERT** — `win_af` only takes {0,.25,.5,.75,1}, so any value in (0.25, 0.5] is identical |
| CO call | `genotype` | 0.2 | **INERT** — after smoothing, column 8 holds only 0/1. Verified 0.1/0.2/0.45 → identical |
| cellsnp | `min_mapq` | 20 | must stay equal to `filter_bam_uniq`'s `-q` |

**The pattern worth noticing:** the *obligate-CO floor* did most of the
selecting — it bounded `block_size` above, excluded `base_af 0.2`, and chose
`max_switch_rate 0.12`. One physical constraint (one CO per bivalent → 2 of 4
chromatids recombinant → ≥50 cM per chromosome) outperformed every statistical
criterion. Foreground this in the methods.

---

## §1 Provenance

- [ ] `git status` clean; working tree == HEAD
- [ ] Config values match `results/crossovers/Spondias_hap1/co_calling_params.txt`
      (written BY the run — ground truth)
- [ ] `results/markers/Spondias_hap1/filter_summary.txt` filter line matches config
- [ ] **Cuscuta untouched**: `snakemake results/landscape/6303_B_hap1/01_landscape.pdf -n`
      → nothing to be done. If it plans a rerun, per-species pinning has broken
      and Cuscuta's published figures are at risk.

**Sweep discipline — learned the hard way this week.** Two sweeps silently
produced stale results. Any sweep must carry an *independent* confirmation that
the thing actually changed (e.g. `wc -l cellSNP.samples.tsv` per iteration), and
must never hide stderr with `> /dev/null 2>&1` without checking the exit code.

**Known landmine:** `results/annotation/Spondias_hap1.gtf` is hand-staged by
`fix_spondias_gtf.py`, not produced by a rule. **Never run `--forceall`** on
this sample — `gff_to_gtf` would clobber it.

---

## §2 Inputs

### Assembly
- 16 chromosomes named `1`–`16`, 283,703,128 bp
- [ ] **Confirm published 2n for *S. tuberosa*.** `samples.csv` says 32 (n=16),
      which matches the assembly. Verify against cytology — the obligate-CO
      arithmetic depends on it.
- [ ] Collaborator-assembled. **Ask for their QC** (BUSCO, N50, Merqury).
      You are building a genetic map on their coordinates.

### Annotation
- Normalized by `fix_spondias_gtf.py`: seqnames `scaffold_N` → `N` (hap2 did
  NOT need this — collaborator was inconsistent between haplotypes), and
  `transcript_id` recovered from `ID`/`Parent` (it held a functional
  *description*, non-unique: 24,151 genes → 15,971 distinct strings)
- `gene_id` never modified → functional annotation joins back to their names
- **Only used for the STAR splice index**, not for marker or CO calling

### HiFi (→ markers)
- ~35 Gb → ~91× on 284 Mb
- [ ] **Is this the same individual as the pollen donor?** First-order
      assumption. If not, the markers don't describe the meiosis you measured.
      Verify explicitly with whoever supplied the material.

### scRNA
- lib `6178A`, mzhang, Nov 2023, 2 lanes, 935 M read pairs
- **10x 5' v1**, established three ways: library table; 70.7% barcode match to
  the 5' whitelist vs 0.3% to 3'; R1 = 28 bp with CB 1–16, UMI 17–26 (entropy
  ~2.0) and a fixed `TT` adapter at 27–28 (entropy ~1.0, 75% TT) → **UMI is
  10 bp**
- cDNA on **R2** (2-file layout), unlike Cuscuta's 5' which used R3
- [ ] **Was the pollen from one plant or pooled?** Pooling would mix meioses
      from different individuals and invalidate a single marker set.

---

## §3 Alignment — RESOLVED

| param | value | status |
|---|---|---|
| CB 1–16, UMI 17–26 (len 10) | | proven from read structure |
| `--soloStrand Reverse` | | **CONFIRMED genome-wide**: 80.9% of reads antisense to genes vs 11.2% sense (7.2:1). Unstranded would be ~50/50 |
| `--soloBarcodeReadLength 0` | | needed, R1 is 28 not 26 |
| `--sjdbOverhang` | 99 | should be read_len−1 = 89. **Harmless** (over-long junction sequences waste index memory, lose no sensitivity). Fix at next index build |
| `--alignIntronMax` | unset → ~590 kb | see below |

**`alignIntronMax` — evidence collected, fix deferred.** Measured on the
MAPQ≥20 BAM (2% sample, 539,048 N gaps) vs the Helixer annotation (121,873
introns):

| | p50 | p99 | p99.9 | max | >20 kb |
|---|---|---|---|---|---|
| annotated | 177 | 2,881 | 8,249 | 19,287 | **0** |
| observed | 201 | 3,570 | 229,194 | 766,661 | 1,600 (0.297%) |

The two agree to p99 then diverge 28×. The observed distribution is **bimodal**:
real introns to ~10 kb, a valley, then a separate artefact population from ~50 kb
to the cap. **`align_intron_max: 20000`** sits in the valley, above the largest
annotated intron. Recorded in config, **not applied** — it would invalidate the
BAM and everything downstream. Apply at next alignment.

Mapping: 25.05% unique, 60.01% multi, 14.9% too short. The multi is almost
entirely the rDNA array — see §4.

---

## §4 The rDNA array — RESOLVED, and the biggest single finding

**What it is.** 45S rDNA at **chr13:1–49,152**, five tandem units, 11,154 bp
period (barrnap 18S starts 928/12,082/23,236/34,390/45,545). Holds **1.25
billion records = 85% of the BAM** in 0.017% of the assembly, ~2.3 million×
coverage vs 41× background. 98% MAPQ 0.

**What it broke.** `cell_snp_counting` hung for **2 days 20 h**. cellsnp applies
`--minMAPQ` *after* fetching and decompressing each pileup, so every fetch near
the array decompressed a mountain of reads to keep ~none — 10.7 TB of reads
against a 65 GB file. Excluding markers by coordinate did **not** help, because
STAR's ~590 kb default `alignIntronMax` creates spurious long introns and htslib
returns any read whose *span* overlaps the query.

**The fix.** `filter_bam_uniq` applies the identical MAPQ predicate once, up
front. Output-equivalent (verified: per-region `view -c -q 20` counts identical).
**65 GB → 5.3 GB; 2d20h hung → 338 seconds.** The general lesson: when a tool
filters *after* an expensive fetch, and the filter is cheap and static, hoist it
out and materialise the filtered input once.

**Also**: `rrna_blacklist` rule (barrnap → BED, 2 kb flank, `bedtools merge -d
3000` to close the ~1,318 bp intra-array spacers) removes 915 markers (0.062%).
The chr4 5S array (305 copies) is **harmless** — Pol III, not polyadenylated,
0.9× background.

**Multimapping is NOT a genome property.** Outside the array, **98.2% of reads
map to exactly one place** (mean NH 1.05). The earlier note of "NH 2.3–3.1
elsewhere" was wrong. No WGD/paralogy problem.

### The library-prep implication
- **~60% of reads are 45S rRNA** (of 467.7 M input: 397.8 M mapped, ~278 M with
  primary in the array, ~117 M usable)
- Sequencing saturation is already **77.7%**, so resequencing this library gains
  little (~20–25% more unique molecules)
- [ ] **Raise rRNA depletion with mzhang for batch 6178A.** Same budget → ~2.5×
      usable reads. This is the highest-value action available.

---

## §5 Markers — RESOLVED

- 1,471,619 markers; **median spacing 33 bp**, p99 2,373 bp, **only 2 gaps
  > 1 Mb genome-wide**. Marker density is not limiting anywhere.
- Filter attrition: AF window rejects **18.5%** (largest), QUAL 2.4%,
  `max_dp` 2.2%, `min_dp` 0.1%
- Markers/Mb ranges 3,742 (chr12) → 6,923 (chr11), ~1.85× spread

**`max_dp` rationale.** The filter exists to reject **collapsed paralogs** —
duplicated regions merged into one assembly copy, where paralogous differences
mimic a perfect 50/50 het but do **not segregate in meiosis**. MAPQ cannot
detect this (the second copy isn't in the assembly to compete), so depth is the
only handle. Het-SNP depth: q1 79, **median 91**, q3 104, p90 117, max 2,280.
`max_dp: 80` sat *below the median* and rejected **72%** of real het sites,
giving 0.13% heterozygosity — irreconcilable with the collaborators' 0.62%.
Yield plateaus from ~140 (140→250 gains 1.3%), indicating **few collapsed
2-copy paralogs** in this assembly. 150 chosen.

- [ ] **Read the collaborators' manuscript** to confirm the 0.62% figure and its
      method. Your 0.53% is 86% of it; the gap is expected (SNPs only, no
      indels; ratio filter; k-mer vs mapping-based; different denominator).

---

## §6 Cells and genotypes — RESOLVED

- 8,497 barcodes called; knee 674, inflection 442; median 1,127 UMI, 863 genes
- **Knee is soft** (no plateau) — expected for dehydrated, transcriptionally
  quiescent mature pollen. **Swept and shown inconsequential** (§0).
- **6,001,202 cell×marker observations, mean DP 1.20, 89.5% single-read.** This
  is the fundamental fragility: most genotype calls rest on one read.

### Contamination — QUANTIFIED at 0.57%
Three independent estimates agree and all are low:
- **0.57%** minor-allele rate at well-behaved markers (DP≥5)
- **~3.5%** ambient per droplet (mean ambient UMI 38.9 / median cell UMI 1127)
- **0.032** median switch rate in selected cells

The initial alarm (an apparent 26% contamination) was an artefact of using a
*fixed* 0.10/0.90 window regardless of DP — at DP=10 a single contaminating read
gives exactly 0.10. The fix was to classify *markers* by their behaviour across
all cells rather than classify individual observations.

**Unresolved sub-finding (footnote-level).** 652 of 1,172 markers with ≥5
high-depth observations show a persistent minor allele near **0.30** with a hard
ceiling at 0.42. They sit in very highly expressed genes (median 5,277× vs
1,419× for clean markers). Four hypotheses tested and **all excluded**:
collapsed paralogy (HiFi depth identical, 95 vs 98), RNA editing (substitution
spectra identical), multimapping (bad has *less*: 2.6% vs 7.4%), ambient RNA
(bad is *less* ambient-enriched: 0.864 vs 1.276). Cause unknown. Scale: 652 of
1.47 M markers (0.04%), at DP 1.20 contributing one read each, in a pipeline
stable to 3.2% across a 2.7× `marker_num` change. **Report honestly, do not
chase further.**

### Reference bias
ALT fraction 0.473 (HiFi raw), 0.4729 (scRNA global), 0.4704 (single-read) —
a consistent ~2.7% ALT deficit. Mechanism is **ALT dropout**: reads carrying the
ALT allele have one extra mismatch, align marginally worse, and are lost. They
are *not* miscalled as REF. So it costs **sensitivity**, not false switches.
Worth one methods sentence.

---

## §7 Detection limitation — THE key finding

**CO detection is marker-limited and has not saturated.** Confirmed three ways:

1. **Between-cell**: mean COs/chromosome rises 0.53 → 1.03 across marker bins
   (30–50 up to 500+), still climbing at the top
2. **Within-cell** (controls for cell quality/doublets — splits each cell's own
   chromosomes into low/high marker-density halves): 0.0376 → 0.0410 CO/Mb,
   +9%, 95% CI [+0.0001, +0.0067] excluding zero
3. **Excess zero-CO chromosomes**: 38% observed vs 25% expected under
   binomial(n=2, p=0.5)

**Efficiency is uncertain.** A saturating fit gives asymptote 0.854 COs/chr
(→1367 cM, 83% efficiency), but the top marker bin *observes* 1.03 — above the
fitted ceiling — so the exponential form is wrong and 83% is optimistic. True
efficiency plausibly 60–87%; **true map plausibly 1300–1900 cM**.

Median cell carries **132 markers/chromosome**; ~148 would give ~90% detection.

**Correction to an earlier error:** zero-CO chromosomes are *expected*, not
impossible. One CO involves 2 of 4 chromatids, so ~50% of gametes carry zero for
that chromosome even with the obligate CO satisfied. An earlier draft called
this "meiotically impossible" — that was wrong.

**The fix is depth per cell, not more cells.** More cells improves precision of
the mean; more markers per cell improves accuracy. Since the library is 77.7%
saturated, that means rRNA depletion, not resequencing.

---

## §8 Biological sanity — PASSED

- **All 16 chromosomes above 50 cM** (obligate-CO floor)
- **1.42 COs per bivalent** (11.37 COs/cell ÷ 16 × 2) — normal plant range 1–3
- **4.0 cM/Mb** — high end of the typical 1–5 range
- **chr1 is NOT elevated.** Once matched on marker *density* and expressed per
  Mb it sits **4th of 16** with overlapping CIs (0.0518 vs mean 0.0420). An
  earlier matched-*total-markers* test suggested 1.8× elevation, but that
  didn't control for length — chr1's 200 markers span 26.4 Mb vs ~18 Mb for
  others. Its higher total is length.
- Across the other 15, r(length, CO/Mb) = −0.22 — weakly negative, the textbook
  direction (obligate CO as a fixed cost)

---

## §9 Figures to review

- [ ] `results/cells/Spondias_hap1/knee_plot.pdf` — soft knee, shown inconsequential
- [ ] `results/cell_qc/Spondias_hap1/switch_diagnostics.pdf` — **bimodal**: noise
      mode ~0.26, clean haploid shoulder 0.01–0.10. This justifies the QC.
- [ ] `results/crossovers/Spondias_hap1/co_diagnostics.pdf` — COs/cell unimodal,
      no long right tail; CO positions per chromosome
- [ ] `results/landscape/Spondias_hap1/01_landscape.pdf`
- [ ] **`02_marey_maps.pdf`** — most diagnostic. Should be monotonic sigmoid.
      chr1's 2 Mb-bin profile shows distal peaks (2–10, 20–24 Mb) and a trough
      at 14–18 Mb = centromere. Titles now read the sample (were hardcoded
      "Cuscuta epithymum hap1").
- [ ] `diag/*.png` — saturation, marker filters, spacing, introns, contamination

---

## §10 Limitations for the paper

- **1137 cM is a lower bound**; detection ~60–87% efficient, not saturated
- **Block-size sensitivity**: 1307 cM @500 kb, 1137 @1 Mb, 981 @1.5 Mb. Report
  the range. The 280 extra COs at 500 kb are *indistinguishable* from the shared
  1,719 (support 7948 vs 8054, width 347 vs 365 kb) — likely real; 1 Mb kept as
  conservative on density (~10 markers/block vs ~5)
- **~12% of the genome is CO-blind** (terminal 1 Mb of each arm, by construction)
- **~60% of the library is rRNA** → ~117 M usable reads → mean DP 1.20
- **153 cells** — small; per-chromosome precision limited
- **Interference (ν) is NOT reportable.** 1 Mb blocks → ±1 Mb CO localisation;
  ν depends on inter-CO distances. Precedent: Cuscuta's ν moved 3.70 → 2.40 on
  one threshold change *with far better resolution*.
- Markers from one HiFi individual; assembly and annotation collaborator-produced

---

## §11 External comparison — STILL OPEN

- [ ] Published genetic maps for Anacardiaceae (mango, cashew, pistachio) —
      is 1137 cM in range?
- [ ] Any prior linkage map for *S. tuberosa*? Even low-density would validate
- [ ] **Cytological chiasma counts per bivalent** — a direct independent check
      on 1.42 COs/bivalent. Worth searching for.
- [ ] Is 4.0 cM/Mb defensible for a woody perennial?

---

## §12 Technical debt

- [x] ~~`baseAF`/`windowAF`/`genotype` not in config~~ — exposed and swept
- [x] ~~`default_lower` hardcoded~~ — exposed and swept
- [x] ~~plot titles hardcoded to Cuscuta~~ — now from `--sample_id`
- [ ] `--alignIntronMax` unset — value chosen (20 kb), apply at next alignment
- [ ] `--sjdbOverhang 99` should be 89 — apply at next index build
- [ ] `ref_low`/`alt_high` in config but the rule never passes them (harmless:
      config values coincide with script defaults)
- [ ] cellsnp filters live in the shell string, not `params:` → no rerun trigger
- [ ] `hap2_chroms` hardcoded to Cuscuta scaffolds in `composite_landscape.R`
- [ ] `composite_landscape` NaN-crashes for samples with empty repeats/SyRI
      columns → blocks `rule all` and bare `snakemake --unlock`
- [ ] `Spondias_hap1.gtf` hand-staged, not rule-produced
- [ ] `marker_diagnostics.R:118` duplicates the config allele-ratio filter
- [ ] `window_af` and `genotype` are inert at this depth — consider removing or
      documenting in the script itself

---

## §13 The decision

- [ ] Do I believe **1137 cM (lower bound), 4.0 cM/Mb, 153 cells**?
- [ ] Weakest link: single-read genotyping at DP 1.20, and the fact that the
      obligate-CO floor was used to *select* parameters and is then cited as
      evidence the map is sound. Is that circular? (Defence: the floor is a hard
      physical constraint, not a soft expectation — but say so explicitly.)
- [ ] "How do you know the COs are real?" → uniform marker spacing, stable
      across `marker_num` 3–8, contamination 0.57%, all chromosomes clear the
      obligate-CO floor, COs/bivalent 1.42 within the normal plant range
- [ ] Unresolved items I will state up front:
      1. Map length is a detection-limited lower bound
      2. The ~0.30 minor-allele marker subset (cause unknown)
      3. No external validation yet (§11)
