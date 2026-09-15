# Dbinata_hap1 marker audit

Audited 2026-09-15. Companion to `docs/spondias_audit.md`.
Config block: `marker_filters: Drosera_binata:` in `config/config.yaml`.

## Verdict

Markers are sound. 10,843,674 markers at 15,662/Mb, median spacing 13 bp, 0 gaps
exceeding `block_size`, 1.2x max/min density across chromosomes. 84.4% of markers
in alignable sequence are independently confirmed by an assembly-based
hap1-vs-hap2 alignment -- an orthogonal method sharing no failure mode with
read mapping.

Thresholds are unchanged from `default:`. The point of the block is that they are
now binata's own, with evidence.

| param | value | basis |
|---|---|---|
| min_dp | 10 | below the diploid peak |
| max_dp | 80 | 1.6x het-SNP median 51; above the peak |
| min_qual | 30 | near-inert; QUAL saturates ~220 |
| min_alt_ratio | 0.3 | below it the real:artefact trade goes worse than 1:1 |
| max_alt_ratio | 0.7 | matches Cuscuta and Spondias; consistency over 2.4% yield |

## Provenance of the HiFi

`hifi_reads` is two files labelled with two accessions (5546D, 6718B), and the
BAM contains two PacBio movies: `m84191_241130_003922_s3` (Revio, Nov 2024, ~74%
of depth) and `m64074_230718_093436` (Sequel II, Jul 2023, ~26%).

Tested by asking which movie supplied the ALT reads at 655 sites:

| class | movie | alt share | depth share |
|---|---|---|---|
| clean | m64074 | 25.5% | 26.0% |
| clean | m84191 | 74.5% | 74.0% |
| shelf | m64074 | 25.9% | 26.7% |
| shelf | m84191 | 74.1% | 73.3% |

ALT share tracks depth share in both classes. One individual sequenced twice,
not two plants pooled. A pool would have shelf sites drawing ALT reads
overwhelmingly from one movie.

## Assembly integrity

minimap2 asm20, both directions. asm5 under-aligns badly at this divergence
(31.9% vs 59.2% query coverage); its numbers should not be used.

- every `chrN_hap1` has `chrN_hap2` as dominant partner at 87.5-93.8%; no
  splitting, no chimeras
- reciprocal coverage symmetric: 59.2% hap1->hap2, 59.7% hap2->hap1. Asymmetry
  would indicate unresolved duplication in one haplotype.
- hap1 691.8 Mb, hap2 672.2 Mb -- 2.8% apart

## Heterozygosity

6,683,993 substitutions over 242.4 Mb alignable at de<0.05 = 2.76%.

Consistent with D. paradoxa GenomeScope 2.85% (same genus). GenomeScope FAILED
on binata -- the model breaks at this divergence -- so the remembered "~2%" prior
was never a measurement. 2.76% is the number to cite.

Genome splits three ways:

| fraction of hap1 | what it is | markers |
|---|---|---|
| 35% alignable at de<0.05 | orthologous pairs; 86% of all exons | 4.08 M, validatable |
| ~24% alignable at de>=0.05 | more divergent, still paired | mixed |
| 41% not alignable | hemizygous TEs, divergent repeat families | 6.76 M, unvalidatable |

## The low-AF shelf

Raw ALT fraction is bimodal: clean peak at 0.5 plus a wide asymmetric shelf below
0.3. Three hypotheses tested, two rejected.

- Two plants pooled -- REJECTED, see above.
- Collapsed repeats / paralogs -- REJECTED. That mechanism gives high depth AND
  low AF together. The observed shelf ANTI-correlates with depth: 47.2% of sites
  below AF 0.3 at DP 10-25, falling to 22.8% above DP 250.
- Reference bias localised to divergent regions -- ACCEPTED. The reference IS
  hap1. Where hap2 is too divergent to map, hap2 reads are lost, and those are
  precisely the ALT-carrying reads. Depth falls and ALT fraction falls together.
  The failure can only push AF down, never up, which is why the shelf is
  one-sided.

Everything fits: shelf is low-side only; shrinks with depth; 10.3% exonic vs
32.9% non-exonic (genes are conserved between haplotypes, TEs are not); 25.2%
alignable vs 37.2% non-alignable; both movies contribute proportionally;
GenomeScope failed.

The global REFERENCE_BIAS flag at median ALT fraction 0.400 is a FALSE POSITIVE
-- a median taken over a bimodal mixture. The exonic subset sits at 0.485 with
58.7% inside the 0.4-0.6 window. Do not report 0.400 as bias.

## Where the cut belongs

Raw calls in alignable sequence, confirmed by an assembly difference at the same
position and alleles:

| AF bin | n | confirmed | est. real | real | artefact |
|---|---|---|---|---|---|
| <0.10 | 62,315 | 9.1% | 10% | 6k | 56k |
| 0.10-0.20 | 1,030,823 | 14.4% | 16% | 170k | 861k |
| 0.20-0.30 | 1,055,814 | 39.1% | 45% | 471k | 585k |
| 0.30-0.40 | 1,385,057 | 70.6% | 81% | 1,116k | 269k |
| 0.40-0.60 | 3,867,113 | 87.5% | ceiling | -- | -- |
| 0.60-0.70 | 820,620 | 87.6% | ceiling | -- | -- |
| 0.70-0.80 | 309,868 | 84.1% | 96% | 298k | 12k |
| >0.80 | 99,329 | 77.6% | 89% | 88k | 11k |

The 0.40-0.70 plateau at 87.6% is the assembly alignment's own sensitivity
ceiling, not a defect in those markers; est. real divides observed confirmation
by it. The middle bins are therefore calibration, not independent measurement.

0.3 is where the trade flips. Admitting 0.20-0.30 buys ~471k real markers for
~585k artefacts. 0.30-0.40 is ~1,116k real for ~269k.

Artefacts here are worse than their count suggests: they cluster spatially in
divergent sequence, so they produce the SAME wrong genotype at the SAME loci
across many cells. Correlated error reads as a crossover; random error reads as
noise. More cells does not average it away.

## Limitation for the methods section

41% of hap1 has no hap2 alignment at de<0.05. Marker under-detection there is NOT
uniform along chromosomes -- it concentrates in divergent, repeat-rich sequence,
which in a monocentric genome means pericentromeric.

Marker-blindness and genuine crossover suppression therefore predict the same
signal and cannot be separated from these data. The Spondias principle --
"uniform marker spacing means uniform under-detection, so landscape shape is
robust" -- DOES NOT HOLD for binata. Pericentromeric suppression must be reported
with that caveat attached.

Mitigating consideration, not a substitute for the caveat: mismatch repair aborts
strand invasion between divergent sequences (heteroduplex rejection), so the
non-alignable fraction is expected to be genuinely CO-poor. Markers not recovered
there are worth less than their count implies.

## Open items

- chr3 carries ~46,893 reads/Mb against a ~20,000 median across the other 15
  chromosomes -- 2.3x. Possible collapsed satellite or rDNA array. Unexamined.
  Check before producing a landscape figure.
- MANIFEST.txt records repo_commit 07f75d1 (DIRTY WORKING TREE). The reference is
  not reproducible from a commit; republish from a clean tree.
- Accession karyotype unconfirmed. exp_ploidy=2, chr_number_2n=32 and 16
  chromosomes per haplotype are mutually consistent, and 1:1 chromosome pairing
  plus paradoxa's successful diploid GenomeScope fit both support a diploid.
  Worth confirming against accession records for 5546D / 6718B.

## Reproduction

    minimap2 -cx asm20 --cs -t 24 -I8G hap1.fa hap2.fa > h2_to_h1.paf
    # alignable: primary alignments with de<0.05, merged -> 242.4 Mb
    # asm SNPs:  cs-tag substitutions in hap1 coordinates -> 6,683,993
    # join raw VCF to asm SNPs, restrict to alignable, stratify by AD/DP
