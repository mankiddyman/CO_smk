#!/usr/bin/env Rscript
# select_haploid_cells.R -- select_cells' haploidness mode.
# Keeps barcodes whose informative markers follow ONE haplotype:
# haploidness >= MIN_HAPLOIDNESS with >= MIN_WINDOWS full 15-molecule windows
# (table from rule cell_haploidness). Writes one barcode per line, no header:
# the format co_calling reads.
# Usage: select_haploid_cells.R TABLE GOOD_CELLS SUMMARY MIN_HAPLOIDNESS MIN_WINDOWS
a <- commandArgs(trailingOnly = TRUE)
d <- read.delim(gzfile(a[1]), stringsAsFactors = FALSE)
h <- as.numeric(a[4]); w <- as.integer(a[5])
s <- d[!is.na(d$haploidness), ]
k <- s[s$haploidness >= h & s$windows >= w, ]
writeLines(k$barcode, a[2])
out <- c("select_on: haploidness",
         sprintf("min_haploidness: %s", h), sprintf("min_windows: %d", w),
         sprintf("barcodes counted by cellsnp: %d", nrow(d)),
         sprintf("scored (>= 5 full windows): %d", nrow(s)),
         sprintf("selected: %d (%.1f%% of counted)", nrow(k), 100 * nrow(k) / nrow(d)),
         sprintf("selected, median molecules: %.0f", median(k$molecules)),
         sprintf("selected, median haploidness: %.3f", median(k$haploidness)))
writeLines(out, a[3])
cat(out, sep = "\n")
