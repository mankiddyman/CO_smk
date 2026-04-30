#!/usr/bin/env Rscript
# =============================================================================
# Composite recombination landscape figure with annotation tracks
# =============================================================================
# Reads a single row from config/samples.csv to find all input paths.
# Layout (per chromosome pair, top to bottom):
#   - SyRI synteny ribbons (hap1 <-> hap2)
#   - hap1 ideogram
#   - 5 heatmap tracks below hap1: CO rate, gene density, SNP density,
#     satellite density, LTR density
#   - hap2 ideogram (no tracks)
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(optparse)
  library(karyoploteR)
  library(GenomicRanges)
  library(data.table)
})

# =============================================================================
# CLI
# =============================================================================

option_list <- list(
  make_option("--samples_csv", type = "character",
              default = "config/samples.csv"),
  make_option("--sample_id",   type = "character"),
  make_option("--co_intervals", type = "character", default = NULL),
  make_option("--co_per_cell",  type = "character", default = NULL),
  make_option("--chrom_map",    type = "character", default = NULL),
  make_option("--markers_vcf",  type = "character", default = NULL,
              help = "markers.vcf.gz [default: results/markers/{sample}/markers.vcf.gz]"),
  make_option("--out_pdf",      type = "character", default = NULL),
  make_option("--win_mb",       type = "double", default = 5),
  make_option("--step_mb",      type = "double", default = 0.5),
  make_option("--hap2_fai",     type = "character", default = NULL),
  make_option("--hap2_chroms",  type = "character",
              default = "scaffold_2,scaffold_4,scaffold_6,scaffold_8,scaffold_10,scaffold_12,scaffold_14")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$sample_id)) stop("--sample_id is required")

# =============================================================================
# READ SAMPLES.CSV ROW + DERIVE PATHS
# =============================================================================

cat("Reading samples.csv...\n")
samples <- read.table(opt$samples_csv, sep = ",", header = TRUE,
                      stringsAsFactors = FALSE, quote = "\"")
row <- samples[samples$sample_id == opt$sample_id, ]
if (nrow(row) != 1)
  stop(sprintf("Expected 1 row for sample_id '%s', got %d",
               opt$sample_id, nrow(row)))

sid <- opt$sample_id
if (is.null(opt$co_intervals)) opt$co_intervals <- sprintf("results/crossovers/%s/co_intervals.bed", sid)
if (is.null(opt$co_per_cell))  opt$co_per_cell  <- sprintf("results/crossovers/%s/co_per_cell.tsv", sid)
if (is.null(opt$chrom_map))    opt$chrom_map    <- sprintf("results/cell_data/%s/chrom_map.tsv", sid)
if (is.null(opt$markers_vcf))  opt$markers_vcf  <- sprintf("results/markers/%s/markers.vcf.gz", sid)
if (is.null(opt$out_pdf))      opt$out_pdf      <- sprintf("results/composite/%s/composite.pdf", sid)

genes_gff <- row$annotation_gff
ltr_gff   <- row$repeats_ltr_gff
sat_gff   <- row$repeats_satellite_gff
syri_dir  <- row$syri_dir

if (is.null(opt$hap2_fai)) {
  opt$hap2_fai <- file.path(dirname(row$assembly_fasta), "hap2.fasta.fai")
}
hap2_chrom_names <- strsplit(opt$hap2_chroms, ",")[[1]]

cat(sprintf("Sample: %s (%s, %s)\n", sid, row$species, row$haplotype))
cat("  Genes:      ", genes_gff, "\n")
cat("  LTRs:       ", ltr_gff, "\n")
cat("  Satellites: ", sat_gff, "\n")
cat("  SyRI:       ", syri_dir, "\n")
cat("  Hap2 fai:   ", opt$hap2_fai, "\n")
cat("  Markers:    ", opt$markers_vcf, "\n")

dir.create(dirname(opt$out_pdf), showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# DATA LOADING
# =============================================================================

load_data <- function(opt) {
  co <- read.table(opt$co_intervals, sep = "\t",
                   col.names = c("chrom", "start", "end", "barcode"))
  co$chrom <- as.character(co$chrom)
  co$start <- as.numeric(co$start)
  co$end   <- as.numeric(co$end)
  co$mid   <- (co$start + co$end) / 2

  cells <- read.table(opt$co_per_cell, sep = "\t", header = TRUE)

  chrom_map <- read.table(opt$chrom_map, sep = "\t", header = TRUE)
  chrom_map$name <- as.character(chrom_map$name)
  chrom_map$size <- as.numeric(chrom_map$size)

  list(co = co, cells = cells, chrom_map = chrom_map,
       n_cells = length(unique(cells$barcode)))
}

# =============================================================================
# GENERIC HELPERS: SLIDING WINDOW + GFF + VCF
# =============================================================================

compute_sliding_rate <- function(co_positions, chrom_size, n_cells,
                                 win_mb, step_mb) {
  win_size  <- win_mb * 1e6
  step_size <- step_mb * 1e6
  centres   <- seq(win_size / 2, chrom_size - win_size / 2, by = step_size)
  rates     <- numeric(length(centres))
  for (i in seq_along(centres)) {
    lo <- centres[i] - win_size / 2
    hi <- centres[i] + win_size / 2
    n  <- sum(co_positions >= lo & co_positions < hi)
    rates[i] <- (n * 100 / n_cells) / win_mb
  }
  list(centres = centres, rates = rates)
}

compute_feature_density <- function(features_df, chrom_size, win_mb, step_mb,
                                    metric = c("count", "bp_covered")) {
  metric <- match.arg(metric)
  win_size  <- win_mb * 1e6
  step_size <- step_mb * 1e6
  centres   <- seq(win_size / 2, chrom_size - win_size / 2, by = step_size)
  vals      <- numeric(length(centres))

  for (i in seq_along(centres)) {
    lo <- centres[i] - win_size / 2
    hi <- centres[i] + win_size / 2
    in_win <- features_df$start < hi & features_df$end > lo
    if (metric == "count") {
      vals[i] <- sum(in_win)
    } else {
      sub <- features_df[in_win, ]
      if (nrow(sub) == 0) {
        vals[i] <- 0
      } else {
        s <- pmax(sub$start, lo)
        e <- pmin(sub$end, hi)
        vals[i] <- sum(e - s) / win_size
      }
    }
  }
  list(centres = centres, vals = vals)
}

compute_point_density <- function(positions, chrom_size, win_mb, step_mb) {
  win_size  <- win_mb * 1e6
  step_size <- step_mb * 1e6
  centres   <- seq(win_size / 2, chrom_size - win_size / 2, by = step_size)
  vals      <- numeric(length(centres))
  for (i in seq_along(centres)) {
    lo <- centres[i] - win_size / 2
    hi <- centres[i] + win_size / 2
    vals[i] <- sum(positions >= lo & positions < hi)
  }
  list(centres = centres, vals = vals)
}

load_gff <- function(gff_path, feature_types, valid_chroms) {
  gff <- read.table(gff_path, sep = "\t", comment.char = "#",
                    quote = "", stringsAsFactors = FALSE,
                    col.names = c("chrom", "source", "type", "start", "end",
                                  "score", "strand", "phase", "attributes"))
  gff$chrom <- as.character(gff$chrom)
  gff <- gff[gff$chrom %in% valid_chroms, ]
  gff <- gff[gff$type %in% feature_types, ]
  gff$start <- as.numeric(gff$start)
  gff$end   <- as.numeric(gff$end)
  gff
}

load_vcf_positions <- function(vcf_path, valid_chroms) {
  cat(sprintf("  Loading VCF (%s)...\n", basename(vcf_path)))
  # data.table fread can decompress + skip ## header lines via the cmd= interface
  vcf <- fread(cmd = sprintf("zcat %s | grep -v '^##'", vcf_path),
               sep = "\t", header = TRUE, select = c(1, 2),
               col.names = c("chrom", "pos"),
               showProgress = FALSE)
  vcf$chrom <- as.character(vcf$chrom)
  vcf <- vcf[vcf$chrom %in% valid_chroms]
  cat(sprintf("  %d variants on hap1 chromosomes\n", nrow(vcf)))
  vcf
}

# =============================================================================
# HAP2 + SYRI LOADERS
# =============================================================================

load_hap2_chrom_map <- function(hap2_fai_path, hap2_chroms) {
  fai <- read.table(hap2_fai_path, sep = "\t",
                    col.names = c("name", "size", "offset",
                                  "linebases", "linewidth"))
  fai$name <- as.character(fai$name)
  fai$size <- as.numeric(fai$size)
  fai <- fai[fai$name %in% hap2_chroms, ]
  fai[, c("name", "size")]
}

load_syri_ribbons <- function(syri_path, types_keep, hap1_chroms, hap2_chroms,
                              min_size_bp = 100000) {
  cat("  Loading SyRI...\n")
  syri <- read.table(syri_path, sep = "\t", quote = "", comment.char = "",
                     fill = TRUE, stringsAsFactors = FALSE,
                     col.names = c("chrom_ref","start_ref","end_ref",
                                   "ref_seq","alt_seq",
                                   "chrom_qry","start_qry","end_qry",
                                   "id","parent_id","type","status"))
  syri <- syri[syri$type %in% types_keep, ]
  syri <- syri[syri$chrom_ref %in% hap1_chroms, ]
  syri <- syri[syri$chrom_qry %in% hap2_chroms, ]
  syri$start_ref <- as.numeric(syri$start_ref)
  syri$end_ref   <- as.numeric(syri$end_ref)
  syri$start_qry <- as.numeric(syri$start_qry)
  syri$end_qry   <- as.numeric(syri$end_qry)
  syri <- syri[(syri$end_ref - syri$start_ref) >= min_size_bp |
               (syri$end_qry - syri$start_qry) >= min_size_bp, ]
  cat(sprintf("  %d ribbon records (>= %d bp, types: %s)\n",
              nrow(syri), min_size_bp, paste(types_keep, collapse=",")))
  syri
}

# =============================================================================
# TRACK DATA BUILDERS
# =============================================================================

build_co_track_data <- function(d, win_mb, step_mb) {
  out_list <- lapply(seq_len(nrow(d$chrom_map)), function(i) {
    chrom      <- d$chrom_map$name[i]
    chrom_size <- d$chrom_map$size[i]
    cos_chrom  <- d$co[d$co$chrom == chrom, ]
    sr <- compute_sliding_rate(cos_chrom$mid, chrom_size, d$n_cells,
                               win_mb, step_mb)
    data.frame(chrom = chrom, pos = sr$centres, val = sr$rates)
  })
  do.call(rbind, out_list)
}

build_gff_density_track <- function(d, gff_path, feature_types, metric,
                                    win_mb, step_mb, label) {
  cat(sprintf("  Loading %s gff...\n", label))
  feats <- load_gff(gff_path, feature_types, d$chrom_map$name)
  cat(sprintf("  %d %s records on %d chromosomes\n",
              nrow(feats), label, length(unique(feats$chrom))))
  out_list <- lapply(seq_len(nrow(d$chrom_map)), function(i) {
    chrom      <- d$chrom_map$name[i]
    chrom_size <- d$chrom_map$size[i]
    f_chrom    <- feats[feats$chrom == chrom, ]
    fd <- compute_feature_density(f_chrom, chrom_size, win_mb, step_mb, metric)
    data.frame(chrom = chrom, pos = fd$centres, val = fd$vals)
  })
  do.call(rbind, out_list)
}

build_snp_density_track <- function(d, vcf_path, win_mb, step_mb) {
  vcf <- load_vcf_positions(vcf_path, d$chrom_map$name)
  out_list <- lapply(seq_len(nrow(d$chrom_map)), function(i) {
    chrom      <- d$chrom_map$name[i]
    chrom_size <- d$chrom_map$size[i]
    pos_chrom  <- vcf$pos[vcf$chrom == chrom]
    pd <- compute_point_density(pos_chrom, chrom_size, win_mb, step_mb)
    data.frame(chrom = chrom, pos = pd$centres, val = pd$vals)
  })
  do.call(rbind, out_list)
}

# =============================================================================
# HEATMAP TRACK PLOTTER
# =============================================================================

plot_heatmap_track <- function(kp, df, r0, r1, palette_fn, label,
                               step_mb, ymax = NULL,
                               percentile_cap = 0.99) {
  # Cap at chosen percentile so a few outliers don't wash out the rest
  if (is.null(ymax)) {
    ymax <- quantile(df$val[df$val > 0], percentile_cap, na.rm = TRUE)
    if (is.na(ymax) || ymax == 0) ymax <- max(df$val, na.rm = TRUE)
    if (is.na(ymax) || ymax == 0) ymax <- 1
  }

  step_size <- step_mb * 1e6
  pal <- palette_fn(100)

  for (chrom in unique(df$chrom)) {
    sub <- df[df$chrom == chrom, ]
    if (nrow(sub) == 0) next

    norm_vals <- pmin(pmax(sub$val / ymax, 0), 1)
    cols <- pal[ pmax(1, ceiling(norm_vals * 100)) ]

    kpRect(kp, chr = chrom,
           x0 = sub$pos - step_size / 2,
           x1 = sub$pos + step_size / 2,
           y0 = 0, y1 = 1,
           r0 = r0, r1 = r1,
           col = cols, border = NA)

    kpAddLabels(kp, chr = chrom, labels = label,
                r0 = r0, r1 = r1, side = "left",
                cex = 0.5, label.margin = 0.025)
  }

  invisible(ymax)
}

# =============================================================================
# SYRI RIBBON PLOTTER
# =============================================================================

plot_syri_ribbons <- function(kp, syri_df, type_colors, alpha = 0.4) {
  src <- toGRanges(data.frame(
    chr   = syri_df$chrom_ref,
    start = syri_df$start_ref,
    end   = syri_df$end_ref))
  dst <- toGRanges(data.frame(
    chr   = syri_df$chrom_qry,
    start = syri_df$start_qry,
    end   = syri_df$end_qry))
  cols <- type_colors[syri_df$type]
  kpPlotLinks(kp, data = src, data2 = dst,
              col = adjustcolor(cols, alpha.f = alpha),
              border = NA, y = 0)
}

# =============================================================================
# COLOR-BAR LEGEND HELPER
# =============================================================================

draw_colorbar <- function(x_left, x_right, y_bot, y_top,
                          palette_fn, label, max_val, n = 100) {
  cols <- palette_fn(n)
  step <- (y_top - y_bot) / n
  for (i in 1:n) {
    rect(x_left, y_bot + (i - 1) * step,
         x_right, y_bot + i * step,
         col = cols[i], border = NA, xpd = NA)
  }
  rect(x_left, y_bot, x_right, y_top, col = NA, border = "gray50", xpd = NA)
  text(x_right + 0.005, y_top, sprintf("%.2g", max_val),
       cex = 0.55, adj = c(0, 0.5), xpd = NA)
  text(x_right + 0.005, y_bot, "0",
       cex = 0.55, adj = c(0, 0.5), xpd = NA)
  text((x_left + x_right) / 2, y_top + 0.012, label,
       cex = 0.6, adj = c(0.5, 0), font = 2, xpd = NA)
}

# =============================================================================
# MAIN
# =============================================================================

cat("Loading data...\n")
d <- load_data(opt)
cat(sprintf("  %d cells, %d COs, %d hap1 chromosomes\n",
            d$n_cells, nrow(d$co), nrow(d$chrom_map)))

cat("Loading hap2 chrom sizes...\n")
hap2_map <- load_hap2_chrom_map(opt$hap2_fai, hap2_chrom_names)
cat(sprintf("  %d hap2 chromosomes\n", nrow(hap2_map)))

# Interleaved hap1/hap2 chromosome order
mapids <- read.table(file.path(syri_dir, "mapids.txt"), sep = "\t",
                     col.names = c("hap1", "hap2"), stringsAsFactors = FALSE)
mapids <- mapids[match(d$chrom_map$name, mapids$hap1), ]
interleaved <- c(rbind(mapids$hap1, mapids$hap2))

all_sizes <- c(setNames(d$chrom_map$size, d$chrom_map$name),
               setNames(hap2_map$size, hap2_map$name))
custom_genome <- toGRanges(data.frame(
  chr   = interleaved,
  start = 1,
  end   = all_sizes[interleaved]))

cat("Building track data...\n")
co_track_df <- build_co_track_data(d, opt$win_mb, opt$step_mb)
gene_df <- build_gff_density_track(d, genes_gff,
                                   feature_types = c("gene"),
                                   metric = "count",
                                   opt$win_mb, opt$step_mb, "gene")
sat_df  <- build_gff_density_track(d, sat_gff,
                                   feature_types = c("tandem_repeat"),
                                   metric = "bp_covered",
                                   opt$win_mb, opt$step_mb, "satellite")
ltr_df  <- build_gff_density_track(d, ltr_gff,
                                   feature_types = c("transposable_element"),
                                   metric = "bp_covered",
                                   opt$win_mb, opt$step_mb, "LTR")
snp_df  <- build_snp_density_track(d, opt$markers_vcf,
                                   opt$win_mb, opt$step_mb)

cat("Loading SyRI ribbons...\n")
syri_df <- load_syri_ribbons(file.path(syri_dir, "syri.out"),
                             types_keep  = c("SYN", "INV"),
                             hap1_chroms = d$chrom_map$name,
                             hap2_chroms = hap2_map$name,
                             min_size_bp = 100000)

# =============================================================================
# RENDER
# =============================================================================

cat("Rendering plot...\n")

# A4 portrait: 8.27 x 11.69 inches; we use 7.5 x 10.5 to leave margin
pdf(opt$out_pdf, width = 7.5, height = 10.5)

pp <- getDefaultPlotParams(plot.type = 2)
pp$leftmargin     <- 0.13
pp$rightmargin    <- 0.13
pp$topmargin      <- 30
pp$bottommargin   <- 35
pp$ideogramheight <- 9
pp$data1height    <- 95   # bumped slightly for 5 tracks
pp$data1inmargin  <- 2
pp$data2height    <- 50
pp$data2inmargin  <- 2

kp <- plotKaryotype(genome = custom_genome, plot.type = 2,
                    plot.params = pp, chromosomes = "all",
                    cex = 0.55, main = "")
kpAddBaseNumbers(kp, tick.dist = 20e6, tick.len = 4, cex = 0.4,
                 minor.tick.dist = 5e6, minor.tick.len = 2)

# Palettes (light = low, saturated = high)
co_pal   <- colorRampPalette(c("#ffffff", "#fde0dd", "#fa9fb5",
                               "#c51b8a", "#7a0177"))
gene_pal <- colorRampPalette(c("#ffffff", "#e5f5e0", "#a1d99b",
                               "#41ab5d", "#005a32"))
snp_pal  <- colorRampPalette(c("#ffffff", "#deebf7", "#9ecae1",
                               "#3182bd", "#08306b"))
sat_pal  <- colorRampPalette(c("#ffffff", "#efedf5", "#bcbddc",
                               "#807dba", "#3f007d"))
ltr_pal  <- colorRampPalette(c("#ffffff", "#fee6ce", "#fdae6b",
                               "#e6550d", "#7f2704"))

# Stack 5 heatmap bands below ideogram with small gaps between bands
# Band height = 0.16, gap = 0.03, top of topmost band = 0.97
band_h <- 0.16
gap    <- 0.03
top    <- 0.97

co_max  <- plot_heatmap_track(kp, co_track_df,
                              r0 = top - 1*band_h - 0*gap,
                              r1 = top - 0*band_h - 0*gap,
                              palette_fn = co_pal, label = "CO",
                              step_mb = opt$step_mb)
gene_max <- plot_heatmap_track(kp, gene_df,
                               r0 = top - 2*band_h - 1*gap,
                               r1 = top - 1*band_h - 1*gap,
                               palette_fn = gene_pal, label = "Gene",
                               step_mb = opt$step_mb)
snp_max  <- plot_heatmap_track(kp, snp_df,
                               r0 = top - 3*band_h - 2*gap,
                               r1 = top - 2*band_h - 2*gap,
                               palette_fn = snp_pal, label = "SNP",
                               step_mb = opt$step_mb)
sat_max  <- plot_heatmap_track(kp, sat_df,
                               r0 = top - 4*band_h - 3*gap,
                               r1 = top - 3*band_h - 3*gap,
                               palette_fn = sat_pal, label = "Sat",
                               step_mb = opt$step_mb)
ltr_max  <- plot_heatmap_track(kp, ltr_df,
                               r0 = top - 5*band_h - 4*gap,
                               r1 = top - 4*band_h - 4*gap,
                               palette_fn = ltr_pal, label = "LTR",
                               step_mb = opt$step_mb)

# Synteny ribbons
type_colors <- c(SYN = "gray60", INV = "tomato")
plot_syri_ribbons(kp, syri_df, type_colors)


# Color-bar legends in the right margin
op <- par(no.readonly = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0),
    new = TRUE)
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")

# Title in the top margin (NDC overlay so karyoploteR cannot push it around)
text(0.5, 0.985,
     sprintf("Composite landscape - %s   (n=%d cells, %d COs)",
             sid, d$n_cells, nrow(d$co)),
     cex = 1.0, font = 2, adj = c(0.5, 0.5), xpd = NA)

# Stack 5 colorbars + ribbon legend vertically along the right side
bar_x_left  <- 0.905
bar_x_right <- 0.935
bar_h       <- 0.085
bar_gap     <- 0.045
top_y       <- 0.88

draw_colorbar_v <- function(y_top, palette_fn, label, max_val) {
  draw_colorbar(bar_x_left, bar_x_right,
                y_top - bar_h, y_top,
                palette_fn, label, max_val)
}

draw_colorbar_v(top_y - 0 * (bar_h + bar_gap), co_pal,
                "CO\ncM/Mb", co_max)
draw_colorbar_v(top_y - 1 * (bar_h + bar_gap), gene_pal,
                sprintf("Genes\n/%g Mb", opt$win_mb), gene_max)
draw_colorbar_v(top_y - 2 * (bar_h + bar_gap), snp_pal,
                sprintf("SNPs\n/%g Mb", opt$win_mb), snp_max)
draw_colorbar_v(top_y - 3 * (bar_h + bar_gap), sat_pal,
                "Sat\nbp frac", sat_max)
draw_colorbar_v(top_y - 4 * (bar_h + bar_gap), ltr_pal,
                "LTR\nbp frac", ltr_max)

# Ribbon legend at the bottom of the right column
ribbon_y <- top_y - 5 * (bar_h + bar_gap) + 0.01
text((bar_x_left + bar_x_right) / 2, ribbon_y + 0.04,
     "Ribbons", cex = 0.6, font = 2, adj = c(0.5, 0), xpd = NA)
rect(bar_x_left, ribbon_y + 0.020, bar_x_right, ribbon_y + 0.030,
     col = adjustcolor("gray60", alpha.f = 0.6), border = NA, xpd = NA)
text(bar_x_right + 0.005, ribbon_y + 0.025,
     "SYN", cex = 0.55, adj = c(0, 0.5), xpd = NA)
rect(bar_x_left, ribbon_y + 0.005, bar_x_right, ribbon_y + 0.015,
     col = adjustcolor("tomato", alpha.f = 0.6), border = NA, xpd = NA)
text(bar_x_right + 0.005, ribbon_y + 0.010,
     "INV", cex = 0.55, adj = c(0, 0.5), xpd = NA)

par(op)
dev.off()
cat("Wrote:", opt$out_pdf, "\n")
