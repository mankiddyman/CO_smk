# Dparadoxa_hap1 marker audit

Audited 2026-09-23. Companion to `docs/binata_marker_audit.md`.
Config block: `marker_filters: Drosera_paradoxa:` in `config/config.yaml`.

## Verdict

Thresholds are NOT binata's. Two differ, both because paradoxa's artefact
population is larger and reaches further up the ALT ratio.

| param | paradoxa | binata | basis |
|---|---|---|---|
| min_dp | 10 | 10 | inert; DP<10 is 0.3% of het SNPs |
| max_dp | 100 | 80 | 1.75x het-SNP median 57, at p90, below the 2x paralog line |
| min_qual | 30 | 30 | costs 0.01%; 99.9% of QUAL<30 already removed by depth or ratio |
| min_alt_ratio | **0.4** | 0.3 | where the real:artefact trade flips HERE |
| max_alt_ratio | 0.7 | 0.7 | identical ascertainment across species |

## Source calls

`results/markers/Dparadoxa_hap1/hifi_raw.vcf.gz`: 38,343,168 records, 98.0%
SNPs, 2.0% indels; 37,314,339 biallelic het SNPs.

Depth at het SNPs: p1 20, p5 30, p25 45, p50 57, p75 75, p90 98, p95 115,
p99 158. QUAL saturates at 222 as it does for binata.

## The shelf is two populations, not one

| depth band | n | ratio <0.3 | ratio 0.4-0.6 |
|---|---|---|---|
| <=45 | 9,050,840 | 53.7% | 24.2% |
| 45-57 | 8,942,113 | 21.0% | 55.7% |
| 57-75 | 9,795,294 | 33.5% | 45.9% |
| 75+ | 9,526,092 | 53.5% | 16.5% |

Binata's shelf anti-correlated with depth throughout: reference bias, where
divergent hap2 reads fail to map and only ALT reads are lost. That population
is here too (the lowest band). The top band is the opposite sign and is
collapsed paralogs: minority alleles from a second locus piled onto one. They
also raise the ratio, which is why median depth is 111 in the 0.80-0.85 bin
against 55 at the 0.5 peak.

## Confirmation against the assembly (chr5)

Method as in the binata audit: raw calls inside alignable sequence, confirmed
by an assembly difference at the same position.

    sort -k6,6 -k8,8n results/hap_align/Dparadoxa_hap1/by_chrom/chr5_hap2.paf > chr5_sorted.paf
    k8 paftools.js call -f refs/Drosera_paradoxa_hap1/genome.fa chr5_sorted.paf > chr5_asm.vcf
    # alignable: primary alignments, de<0.05, blocks >=50 kb (paftools -L) -> 111.4 Mb

757,782 assembly SNVs on chr5. The 0.40-0.60 plateau reaches 83.4%, which is
the alignment's own sensitivity ceiling rather than a defect in those markers;
estimated real divides the observed rate through by it.

| ratio | n (genome) | in alignable (chr5) | confirmed | est. real |
|---|---|---|---|---|
| 0.10-0.15 | 2,999,279 | 167,371 | 1.1% | 1.3% |
| 0.15-0.20 | 4,406,302 | 239,543 | 1.9% | 2.3% |
| 0.20-0.25 | 4,032,297 | 219,341 | 4.0% | 4.8% |
| 0.25-0.30 | 3,677,590 | 194,223 | 9.2% | 11.0% |
| 0.30-0.35 | 3,287,750 | 156,397 | 20.4% | 24.5% |
| 0.35-0.40 | 3,131,628 | 134,862 | 41.2% | 49.4% |
| 0.40-0.45 | 3,683,745 | 155,538 | 64.8% | 77.7% |
| 0.45-0.50 | 3,553,236 | 156,263 | 77.4% | 92.8% |
| 0.50-0.60 | 5,887,665 | 259,914 | 83.2% | ceiling |
| 0.60-0.65 | 1,292,990 | 54,069 | 78.7% | 94.3% |
| 0.65-0.70 | 681,148 | 25,596 | 72.2% | 86.5% |
| 0.70-0.85 | 653,477 | 22,403 | 67.0% | 80.5% |

The trade per bin: 0.30-0.35 buys ~806k real for ~2.48M artefact (1:3),
0.35-0.40 ~1.55M for ~1.58M (1:1), 0.40-0.45 ~2.86M for ~821k (3.5:1). Hence
0.40.

## Open

- chr5 only. The other five chromosomes exceeded minimap2's 32-bit chaining
  limit at asm20 and are unaligned; chr5 is the least repetitive of the six,
  so these artefact fractions are a lower bound. Re-run when they land.
- `gap_flag_bp` is provisional until `block_size` is chosen at the co_calling
  gate.
- Marker density lands near 6,700-7,000/Mb against binata's 15,662/Mb, which
  reflects the tighter ratio window and paradoxa's larger genome, not a
  shortage: at a 2 Mb block that is still >13,000 markers per block.
