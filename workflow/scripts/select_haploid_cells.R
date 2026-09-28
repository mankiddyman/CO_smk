#!/usr/bin/env Rscript
# select_haploid_cells.R -- select_cells' haploidness mode.
# Keeps barcodes whose informative markers follow ONE haplotype
#   haploidness >= MIN_HAPLOIDNESS with >= MIN_WINDOWS full 15-molecule windows
# AND that carry enough evidence to call crossovers on EVERY chromosome
#   molecules >= MIN_MOLECULES, weakest main chromosome >= MIN_CHROM_MOLECULES
# (table from rule cell_haploidness). Writes one barcode per line, no header:
# the format co_calling reads.
# Usage: select_haploid_cells.R TABLE GOOD_CELLS SUMMARY MIN_HAPLOIDNESS MIN_WINDOWS [MIN_MOLECULES MIN_CHROM_MOLECULES]
a <- commandArgs(trailingOnly = TRUE)
d <- read.delim(gzfile(a[1]), stringsAsFactors = FALSE)
h <- as.numeric(a[4]); w <- as.integer(a[5])
mm <- if (length(a) >= 6 && nzchar(a[6])) as.numeric(a[6]) else 0
mc <- if (length(a) >= 7 && nzchar(a[7])) as.numeric(a[7]) else 0
if (mc > 0 && !("min_chrom_molecules" %in% names(d)))
  stop("table has no min_chrom_molecules column -- rerun cell_haploidness")
s <- d[!is.na(d$haploidness), ]
k <- s[s$haploidness >= h & s$windows >= w, ]
k2 <- k[k$molecules >= mm & (if (mc > 0) k$min_chrom_molecules >= mc else TRUE), ]
writeLines(k2$barcode, a[2])
out <- c("select_on: haploidness",
         sprintf("min_haploidness: %s", h), sprintf("min_windows: %d", w),
         sprintf("min_molecules: %s", mm), sprintf("min_chrom_molecules: %s", mc),
         sprintf("barcodes counted by cellsnp: %d", nrow(d)),
         sprintf("scored (>= 5 full windows): %d", nrow(s)),
         sprintf("pass haploidness and windows: %d", nrow(k)),
         sprintf("selected, also passing the coverage floors: %d (%.1f%% of counted)",
                 nrow(k2), 100 * nrow(k2) / nrow(d)),
         sprintf("selected, median molecules: %.0f", median(k2$molecules)),
         if ("min_chrom_molecules" %in% names(k2))
           sprintf("selected, median weakest-chromosome molecules: %.0f", median(k2$min_chrom_molecules)) else NULL,
         sprintf("selected, median haploidness: %.3f", median(k2$haploidness)))
writeLines(out, a[3])
cat(out, sep = "\n")
