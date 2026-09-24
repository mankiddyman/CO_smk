#!/usr/bin/env Rscript
# Switch-rate diagnostics for one sample: the histogram cell_qc's threshold is
# read from, rate against cell depth, and per-chromosome rates of kept cells.
# Usage: Rscript plot_switch_diagnostics.R <sample> [threshold]
args <- commandArgs(trailingOnly = TRUE)
s <- args[1]; thr <- if (length(args) > 1) as.numeric(args[2]) else 0.10
d <- read.delim(file.path("results/cell_qc", s, "switches.tsv"))
st <- read.delim(file.path("results/cells", s, "barcode_stats.tsv"))
d$umi <- st$total_umi[match(d$barcode, st$barcode)]
d <- d[d$total_markers > 0, ]
out <- file.path("qc/cell_qc", s); dir.create(out, recursive = TRUE, showWarnings = FALSE)
pdf(file.path(out, "switch_diagnostics.pdf"), width = 11, height = 8.5)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

h <- hist(d$switch_rate, breaks = 120, plot = FALSE)
plot(h, col = ifelse(h$mids <= thr, "#2E7D32", "#D85A30"), border = NA,
     main = sprintf("%s switch rate: green = kept at %.2f", s, thr),
     xlab = "switches / adjacent marker pairs")
abline(v = thr, lty = 2)
legend("topright", bty = "n", cex = .85,
       legend = sprintf("<=%.2f: %s\n>%.2f: %s", thr, format(sum(d$switch_rate <= thr), big.mark = ","),
                        thr, format(sum(d$switch_rate > thr), big.mark = ",")))

plot(d$total_markers, d$switch_rate, log = "x", pch = 16, cex = .25, col = "#88878040",
     xlab = "markers in the cell", ylab = "switch rate", main = "rate against cell depth")
q <- cut(log10(d$total_markers), 25)
lines(10 ^ tapply(log10(d$total_markers), q, median), tapply(d$switch_rate, q, median),
      col = "#D85A30", lwd = 2.5)
abline(h = thr, lty = 3)

keep <- d[d$switch_rate <= thr, ]
if (nrow(keep) > 1) {
  hist(log10(keep$total_markers), breaks = 40, col = "#2E7D32", border = NA,
       main = sprintf("markers per kept cell (n = %s)", format(nrow(keep), big.mark = ",")),
       xlab = "log10(markers)")
  abline(v = log10(c(800, 2500)), lty = 2, col = c("grey30", "#534AB7"))
  legend("topright", bty = "n", cex = .8, lty = 2, col = c("grey30", "#534AB7"),
         legend = c(sprintf(">=800: %s", format(sum(keep$total_markers >= 800), big.mark = ",")),
                    sprintf(">=2500: %s", format(sum(keep$total_markers >= 2500), big.mark = ","))))
  mk <- grep("^markers_", names(d), value = TRUE); sww <- sub("^markers_", "switches_", mk)
  rate <- as.matrix(keep[, sww]) / pmax(1, as.matrix(keep[, mk]) - 1)
  boxplot(rate, names = sub("^markers_chr", "", mk), las = 2, outline = FALSE, col = "#2E7D32",
          border = "grey30", main = "per-chromosome rate, kept cells", xlab = "chromosome",
          ylab = "switch rate")
}
invisible(dev.off())
cat(sprintf("%s: %s cells, %s at <= %.2f (%.1f%%); wrote %s\n", s, format(nrow(d), big.mark = ","),
            format(nrow(keep), big.mark = ","), thr, 100 * nrow(keep) / nrow(d),
            file.path(out, "switch_diagnostics.pdf")))
