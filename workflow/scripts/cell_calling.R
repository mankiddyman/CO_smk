#!/usr/bin/env Rscript

# Cell calling for STARsolo raw output using DropletUtils.
# Outputs multiple thresholds; downstream picks via the symlink.

suppressPackageStartupMessages({
    library(DropletUtils)
    library(Matrix)
    library(BiocParallel)
    library(optparse)
})

option_list <- list(
    make_option("--matrix_dir", type="character"),
    make_option("--out_dir", type="character"),
    make_option("--sample", type="character"),
    make_option("--default_lower", type="integer", default=300,
                help="Default emptyDrops lower threshold for canonical output"),
    make_option("--fdr", type="numeric", default=0.001),
    make_option("--niters", type="integer", default=10000),
    make_option("--seed", type="integer", default=42)
)
opt <- parse_args(OptionParser(option_list=option_list))

set.seed(opt$seed)
dir.create(opt$out_dir, showWarnings=FALSE, recursive=TRUE)

# --- Load ---
message("Loading raw matrix from ", opt$matrix_dir)
sce <- read10xCounts(opt$matrix_dir, col.names=TRUE)
m <- counts(sce)
message(sprintf("Matrix: %d genes x %d barcodes, %d non-zero",
                nrow(m), ncol(m), length(m@x)))

totals <- colSums(m)

# --- barcodeRanks ---
message("Computing barcode ranks...")
br <- barcodeRanks(m, lower=100)
knee_val <- metadata(br)$knee
inflection_val <- metadata(br)$inflection
message(sprintf("Knee: %g UMI, Inflection: %g UMI", knee_val, inflection_val))

knee_called <- names(totals)[totals >= knee_val]
inflection_called <- names(totals)[totals >= inflection_val]

# --- emptyDrops at multiple thresholds ---
ed_thresholds <- c(100, 300, 500, 1000)
ed_results <- list()
for (lo in ed_thresholds) {
    set.seed(opt$seed)
    message(sprintf("Running emptyDrops with lower=%d...", lo))
    ed <- emptyDrops(m, lower=lo, niters=opt$niters, BPPARAM=SerialParam())
    is_cell <- ed$FDR <= opt$fdr & !is.na(ed$FDR)
    ed_results[[as.character(lo)]] <- list(
        ed = ed,
        called = colnames(m)[is_cell]
    )
    message(sprintf("  -> %d cells", sum(is_cell)))
}

# --- Write barcode lists ---
writeLines(knee_called, file.path(opt$out_dir, "barcodes_knee.tsv"))
writeLines(inflection_called, file.path(opt$out_dir, "barcodes_inflection.tsv"))
for (lo in ed_thresholds) {
    writeLines(ed_results[[as.character(lo)]]$called,
               file.path(opt$out_dir, sprintf("barcodes_emptydrops_lower%d.tsv", lo)))
}

# --- Canonical "use this one" output ---
default_path <- file.path(opt$out_dir,
                          sprintf("barcodes_emptydrops_lower%d.tsv", opt$default_lower))
canonical_path <- file.path(opt$out_dir, "barcodes_called.tsv")
if (file.exists(canonical_path)) file.remove(canonical_path)
file.copy(default_path, canonical_path)
message(sprintf("Canonical barcodes: %s -> barcodes_called.tsv", default_path))

# --- Threshold comparison table ---
build_row <- function(method, called) {
    if (length(called) == 0) {
        return(data.frame(method=method, n_cells=0,
                          median_umi=NA, median_genes=NA,
                          min_umi=NA, max_umi=NA))
    }
    cells_m <- m[, called, drop=FALSE]
    data.frame(
        method = method,
        n_cells = length(called),
        median_umi = median(colSums(cells_m)),
        median_genes = median(colSums(cells_m > 0)),
        min_umi = min(colSums(cells_m)),
        max_umi = max(colSums(cells_m))
    )
}

comparison <- rbind(
    build_row(sprintf("knee_%g", knee_val), knee_called),
    build_row(sprintf("inflection_%g", inflection_val), inflection_called),
    do.call(rbind, lapply(ed_thresholds, function(lo) {
        build_row(sprintf("ed_lower%d", lo),
                  ed_results[[as.character(lo)]]$called)
    }))
)
comparison$is_default <- comparison$method == sprintf("ed_lower%d", opt$default_lower)
write.table(comparison, file.path(opt$out_dir, "threshold_comparison.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
print(comparison)

# --- Per-barcode stats (for the default threshold) ---
ed_default <- ed_results[[as.character(opt$default_lower)]]$ed
stats_df <- data.frame(
    barcode = colnames(m),
    total_umi = totals,
    n_genes = colSums(m > 0),
    ed_pvalue = ed_default$PValue,
    ed_fdr = ed_default$FDR,
    is_cell_default = ed_default$FDR <= opt$fdr & !is.na(ed_default$FDR),
    is_cell_knee = totals >= knee_val,
    is_cell_inflection = totals >= inflection_val
)
write.table(stats_df, file.path(opt$out_dir, "barcode_stats.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

# --- Knee plot (PDF — no font issues) ---
pdf(file.path(opt$out_dir, "knee_plot.pdf"), width=9, height=7)
totals_sorted <- sort(totals, decreasing=TRUE)
totals_pos <- totals_sorted[totals_sorted > 0]

plot(seq_along(totals_pos), totals_pos, log="xy",
     type="l", lwd=2,
     xlab="Barcode rank (log)", ylab="UMI count (log)",
     main=sprintf("Barcode rank: %s", opt$sample))
abline(h=knee_val, col="red", lty=2)
abline(h=inflection_val, col="orange", lty=2)
for (lo in ed_thresholds) {
    abline(h=lo, col="blue", lty=3)
}
n_default <- length(ed_results[[as.character(opt$default_lower)]]$called)
abline(v=n_default, col="darkgreen", lty=1, lwd=2)
legend("topright",
       legend=c(sprintf("knee = %g", knee_val),
                sprintf("inflection = %g", inflection_val),
                sprintf("ed thresholds (%s)", paste(ed_thresholds, collapse=",")),
                sprintf("default (lower=%d): %d cells",
                        opt$default_lower, n_default)),
       col=c("red","orange","blue","darkgreen"),
       lty=c(2,2,3,1), lwd=c(1,1,1,2),
       bty="n", cex=0.85)
dev.off()

# --- One-line QC summary ---
default_called <- ed_results[[as.character(opt$default_lower)]]$called
default_m <- m[, default_called, drop=FALSE]
qc <- data.frame(
    sample = opt$sample,
    n_barcodes_total = ncol(m),
    knee_umi = knee_val,
    inflection_umi = inflection_val,
    default_threshold = sprintf("ed_lower%d", opt$default_lower),
    n_cells_default = length(default_called),
    median_umi_default = median(colSums(default_m)),
    median_genes_default = median(colSums(default_m > 0))
)
write.table(qc, file.path(opt$out_dir, "qc_summary.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)

message("=== DONE ===")
message(sprintf("Default cells (lower=%d): %d", opt$default_lower, length(default_called)))
