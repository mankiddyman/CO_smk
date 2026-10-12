#!/usr/bin/env Rscript
# Aggregate per-cell CO predictions into a landscape summary.
#
# --masks (optional; results/crossovers/{sample}/masks.tsv from co_masks.py): cell x chromosome (or piece)
# pairs whose dosage class is masked (config co_calling.<species>.mask_classes: a nucleus that carries
# the sequence twice, not at all, or mixed). Their calls are noise, not crossovers: they are moved to
# co_intervals_masked.bed, and that cell does not count for that sequence (cells_per_chrom.tsv: the
# denominator the landscape scripts use per sequence). Without --masks every cell counts everywhere.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--co_dir", type = "character", help = "Per-cell CO output directory"),
  make_option("--good_cells", type = "character", help = "good_cells.tsv (one barcode per line)"),
  make_option("--chrom_map", type = "character", help = "chrom_map.tsv"),
  make_option("--out_dir", type = "character", help = "Output directory"),
  make_option("--masks", type = "character", default = NULL,
              help = "masks.tsv (barcode, chrom, ..., class, masked yes/no) [optional]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

barcodes <- readLines(opt$good_cells)
barcodes <- barcodes[nchar(barcodes) > 0]
chrom_map <- read.table(opt$chrom_map, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Masked cell x chromosome pairs (only cells in good_cells count)
masked <- data.frame(barcode = character(), chrom = character(), class = character(),
                     stringsAsFactors = FALSE)
if (!is.null(opt$masks)) {
  m <- read.table(opt$masks, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                  colClasses = "character")
  masked <- m[m$masked == "yes" & m$barcode %in% barcodes, c("barcode", "chrom", "class")]
}
mask_key <- paste(masked$barcode, masked$chrom, sep = "\t")
n_masked <- sapply(chrom_map$name, function(c) sum(masked$chrom == c))
cells_per_chrom <- data.frame(
  chrom    = chrom_map$name,
  n_cells  = length(barcodes) - as.integer(n_masked),
  n_masked = as.integer(n_masked),
  masked_classes = sapply(chrom_map$name, function(c) {
    x <- table(masked$class[masked$chrom == c])
    if (length(x) == 0) "" else paste(names(x), as.integer(x), collapse = ", ")
  }),
  stringsAsFactors = FALSE
)
write.table(cells_per_chrom, file.path(opt$out_dir, "cells_per_chrom.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

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

dropped <- NULL
if (!is.null(all_cos) && nrow(all_cos) > 0) {
  is_masked <- paste(all_cos$barcode, all_cos$chrom, sep = "\t") %in% mask_key
  dropped <- all_cos[is_masked, ]
  all_cos <- all_cos[!is_masked, ]
}
mfile <- file.path(opt$out_dir, "co_intervals_masked.bed")
if (is.null(dropped) || nrow(dropped) == 0) {
  invisible(file.create(mfile))
} else {
  write.table(dropped[, c("chrom", "bp_start", "bp_end", "barcode")], mfile,
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

if (is.null(all_cos) || nrow(all_cos) == 0) {
  cat("No COs found across all cells!\n")
  writeLines(c("No COs called."), file.path(opt$out_dir, "co_summary.txt"))
  invisible(file.create(file.path(opt$out_dir, "co_intervals.bed")))
  write.table(data.frame(barcode = barcodes, n_cos = rep(0L, length(barcodes))),
              file.path(opt$out_dir, "co_per_cell.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  pdf(file.path(opt$out_dir, "co_diagnostics.pdf")); plot.new(); title("No COs called"); dev.off()
  quit(status = 0)
}

# Per-cell counts
co_per_cell <- table(factor(all_cos$barcode, levels = barcodes))
co_counts_full <- as.integer(co_per_cell)

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
cat("\n=== Masked cell x chromosome (dosage classes; calls there dropped) ===\n")
if (nrow(masked) == 0) {
  cat("none\n")
} else {
  cat(sprintf("%d cell x chromosome masked, %d calls dropped (co_intervals_masked.bed)\n",
              nrow(masked), if (is.null(dropped)) 0L else nrow(dropped)))
  for (i in seq_len(nrow(cells_per_chrom))) {
    if (cells_per_chrom$n_masked[i] > 0)
      cat(sprintf("  %-20s %3d masked (%s), %d cells count; %d calls dropped\n",
                  cells_per_chrom$chrom[i], cells_per_chrom$n_masked[i], cells_per_chrom$masked_classes[i],
                  cells_per_chrom$n_cells[i],
                  if (is.null(dropped)) 0L else sum(dropped$chrom == cells_per_chrom$chrom[i])))
  }
}
cat("\n=== COs per cell distribution ===\n")
print(table(co_counts_full))
cat("\n=== COs per chromosome (total across cells) ===\n")
print(co_per_chrom)
cat("\n=== Mean COs per cell, per chromosome (cells that count for it) ===\n")
per_chrom_mean <- as.numeric(co_per_chrom) / cells_per_chrom$n_cells
names(per_chrom_mean) <- chrom_map$name
print(round(per_chrom_mean, 2))
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
mtext(sprintf("Mean: %.1f, median: %.0f%s", mean(co_counts_full), median(co_counts_full),
              if (nrow(masked) > 0) sprintf("; %d cell x chromosome masked", nrow(masked)) else ""),
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
