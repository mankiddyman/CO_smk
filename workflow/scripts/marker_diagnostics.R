#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(VariantAnnotation)
    library(optparse)
})

option_list <- list(
    make_option("--vcf", type="character"),
    make_option("--out_stats", type="character"),
    make_option("--out_pdf", type="character"),
    make_option("--out_summary", type="character")
)
opt <- parse_args(OptionParser(option_list=option_list))

message("Reading VCF: ", opt$vcf)
vcf <- readVcf(opt$vcf)
n_total <- nrow(vcf)
message(sprintf("Total variants: %d", n_total))

# Pull what we need into a data frame
gt <- geno(vcf)$GT
ad <- geno(vcf)$AD  # list-column: c(ref_count, alt_count) per variant
dp <- geno(vcf)$DP

is_biallelic <- lengths(alt(vcf)) == 1
is_snp <- isSNV(vcf)

# Compute allele depths and ratios
# Compute allele depths and ratios (vectorized — works on 10M+ variants)
ad_list <- ad[,1]
ref_count <- vapply(ad_list, function(x) if (length(x) >= 1) x[1] else NA_integer_, integer(1))
alt_count <- vapply(ad_list, function(x) if (length(x) >= 2) x[2] else NA_integer_, integer(1))
total_dp <- ref_count + alt_count
alt_ratio <- alt_count / total_dp


stats <- data.frame(
    chrom = as.character(seqnames(vcf)),
    pos = start(vcf),
    qual = qual(vcf),
    dp = total_dp,
    ref_count = ref_count,
    alt_count = alt_count,
    alt_ratio = alt_ratio,
    gt = as.character(gt[,1]),
    is_biallelic_snp = is_biallelic & is_snp
)

# Subset to biallelic SNPs for downstream analysis
snp_stats <- stats[stats$is_biallelic_snp & !is.na(stats$dp), ]
message(sprintf("Biallelic SNPs: %d", nrow(snp_stats)))

# Hets only (the marker candidates)
het_gts <- c("0/1", "1/0", "0|1", "1|0")
hets <- snp_stats[snp_stats$gt %in% het_gts, ]
message(sprintf("Heterozygous biallelic SNPs: %d", nrow(hets)))

write.table(snp_stats, opt$out_stats, sep="\t", quote=FALSE, row.names=FALSE)

# --- Histograms ---
pdf(opt$out_pdf, width=10, height=8)
par(mfrow=c(2,2))

# Depth histogram (cap at 99th percentile for visibility)
dp_cap <- quantile(snp_stats$dp, 0.99, na.rm=TRUE)
hist(pmin(snp_stats$dp, dp_cap), breaks=100,
     main="Depth at biallelic SNPs",
     xlab=sprintf("DP (capped at 99th pctile = %g)", dp_cap),
     col="lightblue")
abline(v=c(10, 20), col="red", lty=2)
legend("topright", c("DP=10","DP=20"), col="red", lty=2, cex=0.8)

# Quality histogram
qual_cap <- quantile(snp_stats$qual, 0.99, na.rm=TRUE)
hist(pmin(snp_stats$qual, qual_cap), breaks=100,
     main="Variant quality (QUAL)",
     xlab=sprintf("QUAL (capped at %g)", qual_cap),
     col="lightgreen")
abline(v=c(20, 30), col="red", lty=2)

# Allele ratio (hets only)
hist(hets$alt_ratio, breaks=50,
     main=sprintf("ALT allele ratio (heterozygous SNPs, n=%d)", nrow(hets)),
     xlab="ALT / (REF+ALT)",
     col="orange")
abline(v=c(0.3, 0.5, 0.7), col="red", lty=c(2,1,2))

# Depth at hets specifically
hist(pmin(hets$dp, dp_cap), breaks=100,
     main=sprintf("Depth at heterozygous SNPs (n=%d)", nrow(hets)),
     xlab=sprintf("DP (capped)"),
     col="pink")
abline(v=c(10, 20), col="red", lty=2)

dev.off()

# --- Summary text ---
sink(opt$out_summary)
cat("=== Marker definition diagnostics ===\n\n")
cat(sprintf("Total variants in raw VCF: %d\n", n_total))
cat(sprintf("Biallelic SNPs: %d\n", nrow(snp_stats)))
cat(sprintf("Heterozygous biallelic SNPs: %d\n", nrow(hets)))
cat("\n--- Depth distribution at SNPs ---\n")
print(summary(snp_stats$dp))
cat("\n--- Depth distribution at heterozygous SNPs ---\n")
print(summary(hets$dp))
cat("\n--- ALT allele ratio at heterozygous SNPs ---\n")
print(summary(hets$alt_ratio))
cat("\n--- Quality score distribution ---\n")
print(summary(snp_stats$qual))
cat("\n--- Per-chromosome het SNP counts ---\n")
print(table(hets$chrom))
cat("\n--- Hets passing typical thresholds ---\n")
for (dp_thr in c(10, 15, 20)) {
    for (qual_thr in c(20, 30)) {
        n_pass <- sum(hets$dp >= dp_thr & hets$qual >= qual_thr &
                      hets$alt_ratio >= 0.3 & hets$alt_ratio <= 0.7,
                      na.rm=TRUE)
        cat(sprintf("DP>=%d, QUAL>=%d, 0.3<=ratio<=0.7:  %d markers\n",
                    dp_thr, qual_thr, n_pass))
    }
}
sink()

message("Done. Inspect ", opt$out_pdf, " and ", opt$out_summary)
