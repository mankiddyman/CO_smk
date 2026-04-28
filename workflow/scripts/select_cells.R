#!/usr/bin/env Rscript
# Apply cell QC thresholds to switches.tsv to produce good_cells.tsv.

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option("--switches", type = "character", help = "Path to switches.tsv"),
  make_option("--out", type = "character", help = "Path for good_cells.tsv"),
  make_option("--out_summary", type = "character", help = "Path for selection_summary.txt"),
  make_option("--min_total_markers", type = "integer", default = 2500),
  make_option("--min_per_chrom_markers", type = "integer", default = 100),
  make_option("--max_switch_rate", type = "double", default = 0.09)
)
opt <- parse_args(OptionParser(option_list = option_list))

sw <- read.table(opt$switches, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chrom_cols <- grep("^markers_", names(sw), value = TRUE)
sw$min_chrom_markers <- apply(sw[, chrom_cols], 1, min)

# Apply filters
pass_total <- sw$total_markers >= opt$min_total_markers
pass_chrom <- sw$min_chrom_markers >= opt$min_per_chrom_markers
pass_sr    <- sw$switch_rate <= opt$max_switch_rate & !is.na(sw$switch_rate)
pass_all   <- pass_total & pass_chrom & pass_sr

good <- sw[pass_all, ]

# Write good_cells.tsv (just barcodes, one per line — used by downstream rules)
write.table(good$barcode, opt$out,
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Selection summary
sink(opt$out_summary)
cat("=== Cell selection summary ===\n")
cat(sprintf("Input cells: %d\n", nrow(sw)))
cat(sprintf("Filters applied:\n"))
cat(sprintf("  total markers >= %d\n", opt$min_total_markers))
cat(sprintf("  min per-chromosome markers >= %d\n", opt$min_per_chrom_markers))
cat(sprintf("  switch rate <= %.3f\n", opt$max_switch_rate))
cat("\n")
cat(sprintf("Cells passing total marker filter: %d\n", sum(pass_total)))
cat(sprintf("Cells passing per-chrom marker filter: %d\n", sum(pass_chrom)))
cat(sprintf("Cells passing switch rate filter: %d\n", sum(pass_sr, na.rm = TRUE)))
cat(sprintf("Cells passing ALL filters: %d (%.1f%%)\n",
            sum(pass_all), 100 * sum(pass_all) / nrow(sw)))
cat("\n")
cat("Distribution of selected cells' switch rates:\n")
print(summary(good$switch_rate))
cat("\nDistribution of selected cells' marker counts:\n")
print(summary(good$total_markers))
sink()

cat(sprintf("Wrote %d cells to %s\n", sum(pass_all), opt$out))
cat(sprintf("Summary in %s\n", opt$out_summary))
