#!/usr/bin/env Rscript
# Compute per-cell switch rates from per-cell TSVs.
# Produces switches.tsv + switch_diagnostics.pdf.
#
# Faithful to Meng's CO_calling.sh logic, with one improvement:
# switches counted per-chromosome then summed, not naively across the whole cell file.
# This avoids spurious "switches" at chromosome transitions.

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("-d", "--cell_data_dir"), type = "character", default = NULL,
              help = "Directory containing per-cell TSVs and chrom_map.tsv"),
  make_option(c("-o", "--out_dir"), type = "character", default = NULL,
              help = "Output directory for switches.tsv and switch_diagnostics.pdf"),
  make_option("--ref_low", type = "double", default = 0.2,
              help = "alt_ratio < this -> genotype 0 (REF) [default: %default]"),
  make_option("--alt_high", type = "double", default = 0.8,
              help = "alt_ratio > this -> genotype 1 (ALT) [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$cell_data_dir) || is.null(opt$out_dir)) {
  stop("Both --cell_data_dir and --out_dir required.")
}

dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# Read chrom_map
chrom_map <- read.table(file.path(opt$cell_data_dir, "chrom_map.tsv"),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
cat("Chromosomes:", nrow(chrom_map), "\n")

# Find per-cell TSVs (excluding chrom_map and conversion_summary)
all_files <- list.files(opt$cell_data_dir, pattern = "\\.tsv$", full.names = TRUE)
exclude <- c("chrom_map.tsv")
cell_files <- all_files[!basename(all_files) %in% exclude]
cat("Cells to process:", length(cell_files), "\n")

# Process each cell
process_cell <- function(filepath) {
  bc <- sub("\\.tsv$", "", basename(filepath))
  df <- read.table(filepath, header = FALSE, sep = "\t",
                   col.names = c("chrom", "pos", "ref", "ref_count", "alt", "alt_count"),
                   stringsAsFactors = FALSE)

  # Keep only informative markers (at least one read)
  df <- df[df$ref_count > 0 | df$alt_count > 0, ]
  total_markers <- nrow(df)
  if (total_markers == 0) {
    return(data.frame(barcode = bc, total_markers = 0, total_switches = 0,
                      switch_rate = NA, stringsAsFactors = FALSE))
  }

  # Per-marker genotype: 0=REF, 1=ALT, 0.5=ambiguous (matches Meng's encoding)
  alt_ratio <- df$alt_count / (df$ref_count + df$alt_count)
  geno <- ifelse(alt_ratio < opt$ref_low, 0,
          ifelse(alt_ratio > opt$alt_high, 1, 0.5))

  # Per-chromosome switch counting — switches do NOT cross chromosome boundaries
  per_chrom_switches <- integer(nrow(chrom_map))
  per_chrom_markers <- integer(nrow(chrom_map))
  names(per_chrom_switches) <- chrom_map$name
  names(per_chrom_markers) <- chrom_map$name

  for (i in seq_len(nrow(chrom_map))) {
    chrom_name <- chrom_map$name[i]
    chrom_geno <- geno[df$chrom == chrom_name]
    per_chrom_markers[i] <- length(chrom_geno)
    if (length(chrom_geno) >= 2) {
      # A "switch" = adjacent markers have different genotype values
      per_chrom_switches[i] <- sum(chrom_geno[-1] != chrom_geno[-length(chrom_geno)])
    }
  }
  total_switches <- sum(per_chrom_switches)

  # Build output row
  result <- data.frame(
    barcode = bc,
    total_markers = total_markers,
    total_switches = total_switches,
    switch_rate = total_switches / total_markers,
    stringsAsFactors = FALSE
  )
  for (chr in chrom_map$name) {
    result[[paste0("markers_", chr)]] <- per_chrom_markers[chr]
    result[[paste0("switches_", chr)]] <- per_chrom_switches[chr]
  }
  return(result)
}
cat("Processing cells...\n")
results <- do.call(rbind, lapply(seq_along(cell_files), function(i) {
  if (i %% 500 == 0) cat("  ", i, "/", length(cell_files), "\n")
  process_cell(cell_files[i])
}))

# Save switches table
out_table <- file.path(opt$out_dir, "switches.tsv")
write.table(results, out_table, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
cat("Wrote", out_table, "\n")

# Diagnostic plots
out_pdf <- file.path(opt$out_dir, "switch_diagnostics.pdf")
pdf(out_pdf, width = 8, height = 11, family = "Helvetica")
par(mfrow = c(3, 1), mar = c(5, 5, 3, 2))

# Plot 1: Marker count distribution
mk <- results$total_markers
hist(mk[mk < quantile(mk, 0.99)], breaks = 50,
     main = "Marker count distribution per cell (top 1% trimmed)",
     xlab = "Markers per cell",
     ylab = "Cells",
     col = "gray40", border = "white")
abline(v = c(300, 500), col = c("orange", "red"), lty = "dashed", lwd = 2)
legend("topright", c("300 (loose)", "500 (strict)"),
       lty = "dashed", col = c("orange", "red"), bty = "n")
mtext(sprintf("Median: %d, mean: %.0f, n cells: %d",
              median(mk), mean(mk), length(mk)), side = 3, line = 0.5, cex = 0.8)

# Plot 2: Switch rate distribution (cells with >=300 markers)
sub <- results[results$total_markers >= 300, ]
sr <- sub$switch_rate
hist(sr, breaks = 50,
     main = sprintf("Switch rate (cells with >=300 markers, n=%d)", nrow(sub)),
     xlab = "Switch rate (switches / markers)",
     ylab = "Cells",
     col = "gray40", border = "white")
abline(v = c(0.10, 0.13), col = c("red", "orange"), lty = "dashed", lwd = 2)
legend("topright", c("0.10 (strict)", "0.13 (loose)"),
       lty = "dashed", col = c("red", "orange"), bty = "n")
mtext(sprintf("Median: %.3f, IQR: %.3f-%.3f",
              median(sr, na.rm = TRUE),
              quantile(sr, 0.25, na.rm = TRUE),
              quantile(sr, 0.75, na.rm = TRUE)),
      side = 3, line = 0.5, cex = 0.8)

# Plot 3: Markers vs switches scatter
plot(sub$total_markers, sub$total_switches,
     pch = 16, cex = 0.5, col = adjustcolor("black", alpha.f = 0.3),
     xlab = "Markers per cell",
     ylab = "Switches per cell",
     main = "Markers vs switches per cell")
abline(a = 0, b = 0.13, col = "orange", lwd = 2, lty = "dashed")
abline(a = 0, b = 0.10, col = "red", lwd = 2, lty = "dashed")
legend("topright",
       c("y=0.13x (Meng's loose)", "y=0.10x (Meng's strict)"),
       lty = "dashed", col = c("orange", "red"), bty = "n")

dev.off()
cat("Wrote", out_pdf, "\n")

# Console summary
cat("\n=== Summary ===\n")
cat(sprintf("Total cells processed: %d\n", nrow(results)))
cat(sprintf("Cells with >=300 markers: %d\n",
            sum(results$total_markers >= 300)))
cat(sprintf("Cells with >=500 markers: %d\n",
            sum(results$total_markers >= 500)))

cat(sprintf("Switch rate quartiles (cells >=300 markers): %.3f / %.3f / %.3f / %.3f / %.3f\n",
            quantile(sub$switch_rate, 0, na.rm = TRUE),
            quantile(sub$switch_rate, 0.25, na.rm = TRUE),
            quantile(sub$switch_rate, 0.5, na.rm = TRUE),
            quantile(sub$switch_rate, 0.75, na.rm = TRUE),
            quantile(sub$switch_rate, 1, na.rm = TRUE)))

cat(sprintf("\nApplying Meng's loose filter (>=300 markers, switch_rate <=0.13): %d cells pass\n",
            sum(results$total_markers >= 300 & results$switch_rate <= 0.13, na.rm = TRUE)))
cat(sprintf("Applying Meng's strict filter (>=500 markers, switch_rate <=0.10): %d cells pass\n",
            sum(results$total_markers >= 500 & results$switch_rate <= 0.10, na.rm = TRUE)))
