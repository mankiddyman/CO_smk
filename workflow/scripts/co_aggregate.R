#!/usr/bin/env Rscript
# Aggregate per-cell CO predictions into a landscape summary.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--co_dir", type = "character", help = "Per-cell CO output directory"),
  make_option("--good_cells", type = "character", help = "good_cells.tsv (one barcode per line)"),
  make_option("--chrom_map", type = "character", help = "chrom_map.tsv"),
  make_option("--out_dir", type = "character", help = "Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

barcodes <- readLines(opt$good_cells)
chrom_map <- read.table(opt$chrom_map, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Read all per-cell co_pred files
all_cos <- do.call(rbind, lapply(barcodes, function(bc) {
  f <- file.path(opt$co_dir, paste0(bc, "_co_pred.txt"))
  if (!file.exists(f) || file.info(f)$size == 0) return(NULL)
  df <- read.table(f, header = FALSE, sep = " ",
                   col.names = c("chrom", "bp_start", "bp_end",
                                 "geno_pre", "geno_post", "marker_count"),
                   stringsAsFactors = FALSE, fill = TRUE)
  df$barcode <- bc
  df
}))

if (is.null(all_cos) || nrow(all_cos) == 0) {
  cat("No COs found across all cells!\n")
  writeLines(c("No COs called."), file.path(opt$out_dir, "co_summary.txt"))
  quit(status = 0)
}

# Per-cell counts
co_per_cell <- table(factor(all_cos$barcode, levels = barcodes))
co_counts <- as.integer(co_per_cell)
co_counts_full <- c(co_counts, rep(0, length(barcodes) - length(co_counts)))

# Per-chromosome counts
co_per_chrom <- table(factor(all_cos$chrom, levels = chrom_map$name))

# Per-cell-per-chrom counts
co_per_cell_chrom <- table(
  factor(all_cos$barcode, levels = barcodes),
  factor(all_cos$chrom, levels = chrom_map$name)
)

# Write summary
sink(file.path(opt$out_dir, "co_summary.txt"))
cat("=== CO calling summary ===\n")
cat(sprintf("Cells processed: %d\n", length(barcodes)))
cat(sprintf("Total COs called: %d\n", nrow(all_cos)))
cat(sprintf("Cells with at least one CO: %d\n", length(unique(all_cos$barcode))))
cat(sprintf("Mean COs per cell: %.2f\n", mean(co_counts_full)))
cat(sprintf("Median COs per cell: %.0f\n", median(co_counts_full)))
cat("\n=== COs per cell distribution ===\n")
print(table(co_counts_full))
cat("\n=== COs per chromosome (total across cells) ===\n")
print(co_per_chrom)
cat("\n=== Mean COs per cell, per chromosome ===\n")
print(round(colMeans(co_per_cell_chrom), 2))
sink()

# Write CO BED
all_cos_bed <- all_cos[, c("chrom", "bp_start", "bp_end", "barcode")]
write.table(all_cos_bed, file.path(opt$out_dir, "co_intervals.bed"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Write per-cell counts
write.table(data.frame(barcode = barcodes, n_cos = co_counts_full),
            file.path(opt$out_dir, "co_per_cell.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# Diagnostic plots
pdf(file.path(opt$out_dir, "co_diagnostics.pdf"), width = 10, height = 12)
par(mfrow = c(3, 1), mar = c(5, 5, 3, 2))

# Plot 1: COs per cell
hist(co_counts_full, breaks = 0:max(20, max(co_counts_full) + 1) - 0.5,
     main = sprintf("COs per cell (n=%d cells)", length(barcodes)),
     xlab = "Number of COs", ylab = "Cells",
     col = "steelblue", border = "white")
mtext(sprintf("Mean: %.1f, median: %.0f", mean(co_counts_full), median(co_counts_full)),
      side = 3, line = 0.3, cex = 0.9)

# Plot 2: COs per chromosome (total)
barplot(co_per_chrom, las = 2,
        main = "Total COs called per chromosome (across all cells)",
        ylab = "CO count", col = "steelblue", border = "white")

# Plot 3: CO position scatter per chromosome
plot(NA, xlim = c(0, max(chrom_map$size)),
     ylim = c(0, nrow(chrom_map) + 0.5),
     xaxt = "n", yaxt = "n",
     xlab = "Position (Mb)", ylab = "Chromosome",
     main = "CO breakpoint positions across chromosomes")
axis(1, at = pretty(c(0, max(chrom_map$size))),
     labels = pretty(c(0, max(chrom_map$size))) / 1e6)
axis(2, at = 1:nrow(chrom_map), labels = chrom_map$name, las = 1)
for (i in 1:nrow(chrom_map)) {
  chrom <- chrom_map$name[i]
  cos_here <- all_cos[all_cos$chrom == chrom, ]
  if (nrow(cos_here) > 0) {
    points(cos_here$bp_start, rep(i, nrow(cos_here)),
           pch = "|", col = "red", cex = 0.8)
  }
  segments(0, i, chrom_map$size[i], i, col = "gray", lwd = 1)
}

dev.off()
cat("Done. Outputs in", opt$out_dir, "\n")
