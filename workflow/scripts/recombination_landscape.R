#!/usr/bin/env Rscript
# =============================================================================
# Recombination Landscape Analysis
# =============================================================================
# Produces 7 PDF plots and 3 TSV summaries from CO calls.
#
# Inputs:
#   --co_intervals   results/crossovers/{sample}/co_intervals.bed
#   --co_per_cell    results/crossovers/{sample}/co_per_cell.tsv
#   --chrom_map      results/cell_data/{sample}/chrom_map.tsv
#   --out_dir        results/landscape/{sample}/
#
# Outputs:
#   01_landscape.pdf
#   02_marey_maps.pdf
#   03_co_class_proportions.pdf
#   04_coc_full.pdf
#   05_coc_zoom.pdf
#   06_gamma_interference.pdf
#   07_inter_co_distances.pdf
#   landscape_summary.tsv
#   coc_table.tsv
#   gamma_interference_summary.tsv
#
# Each plot lives in its own function. To change aesthetics of plot N:
# replace plot_N() function only.
# =============================================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(optparse))

# =============================================================================
# CLI
# =============================================================================

option_list <- list(
  make_option("--co_intervals", type = "character"),
  make_option("--co_per_cell",  type = "character"),
  make_option("--chrom_map",    type = "character"),
  make_option("--out_dir",      type = "character"),
  make_option("--bin_mb",       type = "double",  default = 5,
              help = "Bin size for landscape histogram in Mb [default: %default]"),
  make_option("--coc_interval_mb", type = "double", default = 5,
              help = "Target interval size for CoC binning [default: %default]"),
  make_option("--n_boot",       type = "integer", default = 500,
              help = "Bootstrap iterations [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# DATA LOADING
# =============================================================================

load_data <- function(opt) {
  co <- read.table(opt$co_intervals, sep = "\t",
                   col.names = c("chrom", "start", "end", "barcode"))
  co$chrom   <- as.character(co$chrom)
  co$start   <- as.numeric(co$start)
  co$end     <- as.numeric(co$end)
  co$barcode <- as.character(co$barcode)
  co$mid     <- (co$start + co$end) / 2

  cells <- read.table(opt$co_per_cell, sep = "\t", header = TRUE)
  cells$barcode <- as.character(cells$barcode)
  cells$n_cos   <- as.numeric(cells$n_cos)

  chrom_map <- read.table(opt$chrom_map, sep = "\t", header = TRUE)
  chrom_map$name <- as.character(chrom_map$name)
  chrom_map$size <- as.numeric(chrom_map$size)

  list(co = co, cells = cells, chrom_map = chrom_map,
       n_cells = length(unique(cells$barcode)))
}

# =============================================================================
# PLOT 1: PER-CHROMOSOME RECOMBINATION LANDSCAPE (cM/Mb per bin)
# =============================================================================

# =============================================================================
# PLOT 1: PER-CHROMOSOME RECOMBINATION LANDSCAPE
# Sliding-window CO rate (cM/Mb) with bootstrap ±1 SD ribbon.
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
    rates[i] <- (n * 100 / n_cells) / win_mb  # cM/Mb
  }
  list(centres = centres, rates = rates)
}

plot_landscape <- function(d, out_path, bin_mb,
                           win_mb = 5, step_mb = 0.5, n_boot = 200) {
  suppressPackageStartupMessages(library(ggplot2))

  total_cM      <- nrow(d$co) * 100 / d$n_cells
  total_size_Mb <- sum(d$chrom_map$size) / 1e6
  genome_rate   <- total_cM / total_size_Mb

  # Build long-format data: per chromosome, per window
  build_chrom_df <- function(chrom, chrom_size, cos_chrom) {
    real <- compute_sliding_rate(cos_chrom$mid, chrom_size, d$n_cells,
                                 win_mb, step_mb)
    boot_rates <- matrix(NA, nrow = n_boot, ncol = length(real$centres))
    for (b in seq_len(n_boot)) {
      boot_bcs <- sample(d$cells$barcode, replace = TRUE)
      boot_co  <- cos_chrom[cos_chrom$barcode %in% boot_bcs, ]
      br <- compute_sliding_rate(boot_co$mid, chrom_size, d$n_cells,
                                 win_mb, step_mb)
      boot_rates[b, ] <- br$rates
    }
    sd_rate <- apply(boot_rates, 2, sd, na.rm = TRUE)
    chr_avg <- nrow(cos_chrom) * 100 / d$n_cells / (chrom_size / 1e6)

    data.frame(
      chrom       = chrom,
      pos_Mb      = real$centres / 1e6,
      rate        = real$rates,
      ymin        = pmax(real$rates - sd_rate, 0),
      ymax        = real$rates + sd_rate,
      chrom_size  = chrom_size / 1e6,
      chr_avg     = chr_avg
    )
  }

  per_chrom <- do.call(rbind, lapply(seq_len(nrow(d$chrom_map)), function(i) {
    chrom      <- d$chrom_map$name[i]
    chrom_size <- d$chrom_map$size[i]
    build_chrom_df(chrom, chrom_size, d$co[d$co$chrom == chrom, ])
  }))

  # Preserve chromosome order
  per_chrom$chrom <- factor(per_chrom$chrom, levels = d$chrom_map$name)

  # Per-chromosome means and ends as a separate frame for h/v lines
  ref_lines <- data.frame(
    chrom      = factor(d$chrom_map$name, levels = d$chrom_map$name),
    chr_avg    = sapply(d$chrom_map$name, function(c) {
      n <- sum(d$co$chrom == c)
      sz <- d$chrom_map$size[d$chrom_map$name == c]
      n * 100 / d$n_cells / (sz / 1e6)
    }),
    chrom_end  = d$chrom_map$size / 1e6
  )

  # Legend handled via colour aesthetic mapping with manual scale
  p <- ggplot(per_chrom, aes(x = pos_Mb)) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = "+/-1 SD"),
                alpha = 0.5) +
    geom_line(aes(y = rate, colour = "CO rate"), linewidth = 0.6) +
    geom_hline(aes(yintercept = chr_avg, colour = "Chromosome mean"),
               data = ref_lines, linewidth = 0.4) +
    geom_hline(aes(yintercept = genome_rate, colour = "Genome-wide mean"),
               linewidth = 0.4) +
    geom_vline(aes(xintercept = chrom_end, colour = "Chromosome end"),
               data = ref_lines, linewidth = 0.4) +
    facet_wrap(~ chrom, ncol = 1, strip.position = "top",
               scales = "free_x") +
    scale_colour_manual(
      name = NULL,
      values = c(
        "CO rate"          = "black",
        "Chromosome mean"  = "darkgreen",
        "Genome-wide mean" = "magenta",
        "Chromosome end"   = "blue"
      ),
      breaks = c("CO rate", "Chromosome mean",
                 "Genome-wide mean", "Chromosome end")
    ) +
    scale_fill_manual(name = NULL, values = c("+/-1 SD" = "gray70")) +
    scale_x_continuous(expand = c(0, 0),
                       limits = c(0, max(d$chrom_map$size) / 1e6)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
   labs(
      x        = "Chromosome length (Mb)",
      y        = expression(paste("CO rate (cM Mb"^-1, ")")),
      title    = sprintf("Recombination landscape - Cuscuta epithymum hap1"),
      subtitle = sprintf(
        paste0(
          "n = %d pollen cells, %d crossovers (mean %.2f / cell)\n",
          "Expected from CO assurance: ~%.1f / cell (1 CO per bivalent x %d bivalents x 0.5 chromatid)\n",
          "Genome-wide mean rate: %.2f cM/Mb (%.0f cM total over %.0f Mb assembled)"
        ),
        d$n_cells,
        nrow(d$co),
        nrow(d$co) / d$n_cells,
        0.5 * nrow(d$chrom_map),
        nrow(d$chrom_map),
        genome_rate,
        nrow(d$co) * 100 / d$n_cells,
        sum(d$chrom_map$size) / 1e6
      )
    ) + theme_classic(base_size = 9) +
    theme(
      strip.background = element_blank(),
      strip.text.x     = element_text(angle = 0, hjust = 0, size = 8),
      panel.spacing.y  = unit(0.3, "lines"),
      legend.position  = "top",
      legend.box       = "horizontal",
      legend.margin    = margin(0, 0, 0, 0),
      legend.key.size  = unit(0.8, "lines"),
      plot.margin      = margin(4, 4, 4, 4),
       plot.title       = element_text(size = 10, face = "bold"),
      plot.subtitle    = element_text(size = 8, lineheight = 1.1,
                                      colour = "gray30")
    )

  # Sizing: A4 width ~7 in usable, ~1 in per chromosome panel
  n_panels <- nrow(d$chrom_map)
  total_h  <- 0.9 * n_panels + 1.6   # +1 inch for legend + axis
  ggsave(out_path, p, width = 7, height = total_h, units = "in",
         device = cairo_pdf)
}


# =============================================================================
# PLOT 2: MAREY MAPS (cumulative cM vs physical Mb)
# =============================================================================

plot_marey <- function(d, out_path) {
  suppressPackageStartupMessages(library(ggplot2))

  # Long-format data: cumulative cM per chromosome
  marey_df <- do.call(rbind, lapply(seq_len(nrow(d$chrom_map)), function(i) {
    chrom      <- d$chrom_map$name[i]
    chrom_size <- d$chrom_map$size[i]
    cos_chrom  <- d$co[d$co$chrom == chrom, ]
    cos_chrom  <- cos_chrom[order(cos_chrom$mid), ]

    if (nrow(cos_chrom) == 0) {
      return(data.frame(chrom = chrom,
                        pos_Mb = c(0, chrom_size / 1e6),
                        cum_cM = c(0, 0)))
    }

    data.frame(
      chrom  = chrom,
      pos_Mb = c(0, cos_chrom$mid / 1e6, chrom_size / 1e6),
      cum_cM = c(0, seq_along(cos_chrom$mid) * 100 / d$n_cells,
                 nrow(cos_chrom) * 100 / d$n_cells)
    )
  }))
  marey_df$chrom <- factor(marey_df$chrom, levels = d$chrom_map$name)

  # Final point per chromosome for label placement
  ref_df <- do.call(rbind, lapply(seq_len(nrow(d$chrom_map)), function(i) {
    chrom      <- d$chrom_map$name[i]
    chrom_size <- d$chrom_map$size[i]
    n_cos      <- sum(d$co$chrom == chrom)
    total_cM   <- n_cos * 100 / d$n_cells
    data.frame(
      chrom        = chrom,
      chrom_end_Mb = chrom_size / 1e6,
      total_cM     = total_cM,
      cM_per_Mb    = total_cM / (chrom_size / 1e6)
    )
  }))
  ref_df$chrom <- factor(ref_df$chrom, levels = d$chrom_map$name)

  total_cM_genome <- nrow(d$co) * 100 / d$n_cells
  total_size_Mb   <- sum(d$chrom_map$size) / 1e6
  genome_rate     <- total_cM_genome / total_size_Mb

  # Distinct, colour-blind-friendly palette
  n_chrom <- nrow(d$chrom_map)
  palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[seq_len(n_chrom)]
  names(palette) <- d$chrom_map$name

  p <- ggplot(marey_df, aes(x = pos_Mb, y = cum_cM, colour = chrom)) +
    geom_step(linewidth = 0.7) +
    geom_point(aes(x = chrom_end_Mb, y = total_cM),
               data = ref_df, shape = 16, size = 2) +
    scale_colour_manual(name = NULL, values = palette) +
    scale_x_continuous(expand = expansion(mult = c(0.005, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0.005, 0.05))) +
    labs(
      x        = "Physical position (Mb)",
      y        = "Cumulative genetic distance (cM)",
      title    = "Marey map - Cuscuta epithymum hap1",
      subtitle = sprintf(
        paste0(
          "n = %d pollen cells, %d crossovers\n",
          "Total genetic length: %.0f cM over %.0f Mb (mean %.2f cM/Mb)\n",
          "Each step is one observed CO; slope = local recombination rate"
        ),
        d$n_cells, nrow(d$co),
        total_cM_genome, total_size_Mb, genome_rate
      )
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position    = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background  = element_rect(fill = "white", colour = "gray80"),
      legend.key.size    = unit(0.9, "lines"),
      legend.text        = element_text(size = 8),
      plot.title         = element_text(size = 11, face = "bold"),
      plot.subtitle      = element_text(size = 8.5, lineheight = 1.1,
                                        colour = "gray30"),
      plot.margin        = margin(6, 6, 6, 6)
    )

  # Single square-ish panel suitable as figure subpanel
  ggsave(out_path, p, width = 6, height = 5, units = "in", device = cairo_pdf)
}

# =============================================================================
# PLOT 3: CO CLASS PROPORTIONS PER CHROMOSOME
# =============================================================================

plot_co_classes <- function(d, out_path) {
  suppressPackageStartupMessages(library(ggplot2))

  chrom_names <- d$chrom_map$name
  all_bcs     <- d$cells$barcode

  co_matrix <- matrix(0, nrow = length(all_bcs), ncol = length(chrom_names),
                      dimnames = list(all_bcs, chrom_names))
  for (i in seq_len(nrow(d$co))) {
    co_matrix[d$co$barcode[i], d$co$chrom[i]] <-
      co_matrix[d$co$barcode[i], d$co$chrom[i]] + 1
  }

  max_class    <- 4
  class_levels <- c("0", "1", "2", "3", "4+")

  long_df <- do.call(rbind, lapply(chrom_names, function(c) {
    counts      <- co_matrix[, c]
    capped      <- pmin(counts, max_class)
    tab         <- tabulate(capped + 1, nbins = max_class + 1)
    props       <- tab / sum(tab)
    chr_size_Mb <- d$chrom_map$size[d$chrom_map$name == c] / 1e6
    # Each segment's height = proportion × chromosome size
    segment_Mb  <- props * chr_size_Mb
    data.frame(chrom         = c,
               chrom_size_Mb = chr_size_Mb,
               co_class      = class_levels,
               proportion    = props,
               segment_Mb    = segment_Mb)
  }))

  # Order by natural chromosome number (1, 3, 5, 7, 9, 11, 13)
  long_df$chrom    <- factor(long_df$chrom,    levels = chrom_names)
  long_df$co_class <- factor(long_df$co_class, levels = class_levels)

  # Per-chromosome annotation above each bar: size + mean COs/cell
  ann_df <- do.call(rbind, lapply(chrom_names, function(c) {
    chr_size_Mb <- d$chrom_map$size[d$chrom_map$name == c] / 1e6
    n           <- sum(d$co$chrom == c)
    chr_mean    <- n / d$n_cells
    data.frame(chrom         = c,
               chrom_size_Mb = chr_size_Mb,
               label         = sprintf("%.0f Mb\nmean %.2f", chr_size_Mb, chr_mean))
  }))
  ann_df$chrom <- factor(ann_df$chrom, levels = chrom_names)

  fill_palette <- c("0"  = "#F0F0F0",
                    "1"  = "#9ECAE1",
                    "2"  = "#4292C6",
                    "3"  = "#08519C",
                    "4+" = "#08306B")

  max_height <- max(d$chrom_map$size) / 1e6

  p <- ggplot(long_df, aes(x = chrom, y = segment_Mb, fill = co_class)) +
    geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
    # Percent labels inside segments tall enough to fit them
    geom_text(aes(label = ifelse(segment_Mb >= 5,
                                 sprintf("%.0f%%", proportion * 100), "")),
              position = position_stack(vjust = 0.5),
              size = 2.7, colour = "white", fontface = "bold") +
    # Per-chromosome label above the bar
    geom_text(aes(x = chrom, y = chrom_size_Mb + max_height * 0.04,
                  label = label),
              data = ann_df, inherit.aes = FALSE,
              size = 2.6, colour = "gray30", lineheight = 0.9) +
    scale_fill_manual(name = "COs per cell", values = fill_palette) +
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, max_height * 1.18),
                       breaks = pretty(c(0, max_height))) +
    labs(
      x        = NULL,
      y        = "Chromosome size (Mb)",
      title    = "Per-chromosome CO class distribution - Cuscuta epithymum hap1",
      subtitle = sprintf(
        paste0(
          "n = %d pollen cells, %d total crossovers (mean %.2f / cell)\n",
          "Bar height = chromosome size; segment height = proportion x size"
        ),
        d$n_cells, nrow(d$co), nrow(d$co) / d$n_cells
      )
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x       = element_text(size = 9),
      legend.position   = "right",
      legend.title      = element_text(size = 9),
      legend.key.size   = unit(0.9, "lines"),
      plot.title        = element_text(size = 11, face = "bold"),
      plot.subtitle     = element_text(size = 8.5, lineheight = 1.1,
                                       colour = "gray30"),
      plot.margin       = margin(6, 6, 6, 6)
    )

  ggsave(out_path, p, width = 8, height = 5.5, units = "in", device = cairo_pdf)
}

# =============================================================================
# COC COMPUTATION (used by plots 4 and 5)
# =============================================================================

compute_coc <- function(d, interval_mb, n_boot) {
  N <- ceiling((max(d$chrom_map$size) / 1e6) / interval_mb)
  cell_ids   <- d$cells$barcode
  cell_index <- seq_along(cell_ids); names(cell_index) <- cell_ids
  n_cells    <- d$n_cells

  collect <- function(co_subset) {
    ds <- numeric(0); os <- numeric(0); es <- numeric(0)
    for (chr in d$chrom_map$name) {
      chr_size <- d$chrom_map$size[d$chrom_map$name == chr]
      chr_co   <- co_subset[co_subset$chrom == chr, ]
      if (nrow(chr_co) == 0) next
      breaks <- seq(0, chr_size, length.out = N + 1)
      mat <- matrix(0, nrow = n_cells, ncol = N)
      for (i in seq_len(nrow(chr_co))) {
        bc <- chr_co$barcode[i]
        if (!bc %in% cell_ids) next
        iv <- findInterval(chr_co$mid[i], breaks)
        if (iv >= 1 && iv <= N) mat[cell_index[[bc]], iv] <- 1
      }
      P <- colSums(mat) / n_cells
      iv_size_Mb <- (chr_size / N) / 1e6
      for (a in 1:(N - 1)) {
        for (b in (a + 1):N) {
          obs <- sum(mat[, a] * mat[, b])
          exp <- P[a] * P[b] * n_cells
          if (exp > 0) {
            ds <- c(ds, (b - a) * iv_size_Mb)
            os <- c(os, obs)
            es <- c(es, exp)
          }
        }
      }
    }
    list(d = ds, o = os, e = es)
  }

  bin_agg <- function(p) {
    if (length(p$d) == 0) return(NULL)
    bins <- round(p$d / interval_mb) * interval_mb
    ub <- sort(unique(bins))
    coc <- numeric(length(ub)); n_pairs <- numeric(length(ub))
    tot_e <- numeric(length(ub))
    for (i in seq_along(ub)) {
      idx <- which(bins == ub[i])
      coc[i]     <- sum(p$o[idx]) / sum(p$e[idx])
      n_pairs[i] <- length(idx)
      tot_e[i]   <- sum(p$e[idx])
    }
    list(d_bin = ub, coc = coc, n_pairs = n_pairs, total_exp = tot_e)
  }

  real <- bin_agg(collect(d$co))

  boot_mat <- matrix(NA, nrow = n_boot, ncol = length(real$d_bin))
  for (b in seq_len(n_boot)) {
    boot_bcs <- sample(cell_ids, replace = TRUE)
    boot_co  <- d$co[d$co$barcode %in% boot_bcs, ]
    bp <- bin_agg(collect(boot_co))
    if (is.null(bp)) next
    for (i in seq_along(real$d_bin)) {
      idx <- which(bp$d_bin == real$d_bin[i])
      if (length(idx) > 0) boot_mat[b, i] <- bp$coc[idx]
    }
  }

  ci_lo <- apply(boot_mat, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  ci_hi <- apply(boot_mat, 2, function(x) quantile(x, 0.975, na.rm = TRUE))

  data.frame(d_bin = real$d_bin, coc = real$coc,
             ci_lo = as.numeric(ci_lo), ci_hi = as.numeric(ci_hi),
             n_pairs = real$n_pairs, total_exp = real$total_exp)
}

# =============================================================================
# PLOT 4: CoC FULL RANGE
# =============================================================================
# =============================================================================
# PLOT 4: CoC PLOT WITH ADAPTIVE OPACITY
# =============================================================================
# Single panel: CoC vs distance. Opacity scales with statistical reliability
# (total expected count per bin) - well-supported bins are solid, low-support
# bins fade gracefully without being hidden.

plot_coc <- function(coc_df, out_path, n_cells, interval_mb,
                     gamma_summary = NULL) {
  suppressPackageStartupMessages(library(ggplot2))

  # Map total_exp to alpha: 1 at total_exp=10+, scales down to 0.2 below
  coc_df$alpha_val <- pmin(1, pmax(0.2, coc_df$total_exp / 10))

  subtitle_text <- sprintf(
    paste0(
      "n = %d cells, interval size = %g Mb (MADpattern style)\n",
      "CoC < 1: positive interference; opacity proportional to ",
      "statistical power per bin"
    ),
    n_cells, interval_mb
  )
  if (!is.null(gamma_summary)) {
    subtitle_text <- paste0(
      subtitle_text,
      sprintf("\nGamma model: nu = %.2f, 95%% CI [%.2f, %.2f]",
              gamma_summary$nu, gamma_summary$ci_lo, gamma_summary$ci_hi)
    )
  }

  ymax <- min(max(c(coc_df$ci_hi, coc_df$coc), na.rm = TRUE) * 1.1, 4)

  p <- ggplot(coc_df, aes(x = d_bin, y = coc)) +
    geom_hline(yintercept = 1, colour = "red", linetype = "dashed") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, alpha = alpha_val),
                  width = interval_mb * 0.3, colour = "gray50") +
    geom_line(aes(alpha = alpha_val), colour = "darkblue", linewidth = 0.6,
              group = 1) +
    geom_point(aes(alpha = alpha_val), colour = "darkblue", size = 2.2) +
    annotate("text", x = max(coc_df$d_bin) * 0.97, y = 1.05,
             label = "no interference (CoC = 1)", hjust = 1,
             colour = "red", size = 2.7) +
    scale_alpha_identity() +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                       limits = c(0, ymax)) +
    labs(
      x        = "Inter-interval distance (Mb)",
      y        = "Coefficient of coincidence (CoC)",
      title    = "Crossover interference (CoC) - Cuscuta epithymum hap1",
      subtitle = subtitle_text
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 8.5, lineheight = 1.1,
                                    colour = "gray30"),
      plot.margin   = margin(6, 6, 6, 6)
    )

  ggsave(out_path, p, width = 7, height = 4.5, units = "in", device = cairo_pdf)
}

# =============================================================================
# GAMMA MODEL (used by plots 6 and 7)
# =============================================================================

get_inter_co_dists <- function(co_subset) {
  out <- numeric(0)
  for (bc in unique(co_subset$barcode)) {
    pos <- sort(co_subset$mid[co_subset$barcode == bc])
    if (length(pos) >= 2) out <- c(out, diff(pos) / 1e6)
  }
  out
}

fit_gamma <- function(dists) {
  if (length(dists) < 3) return(list(nu = NA, theta = NA, n = length(dists)))
  nll <- function(p) {
    if (p[1] <= 0 || p[2] <= 0) return(Inf)
    -sum(dgamma(dists, shape = p[1], scale = p[2], log = TRUE))
  }
  m <- mean(dists); v <- var(dists)
  init <- if (v > 0) c(m^2 / v, v / m) else c(2, m / 2)
  res <- tryCatch(
    optim(init, nll, method = "L-BFGS-B",
          lower = c(0.05, 0.05), upper = c(50, 1000)),
    error = function(e) NULL
  )
  if (is.null(res)) return(list(nu = NA, theta = NA, n = length(dists)))
  list(nu = res$par[1], theta = res$par[2], n = length(dists))
}

compute_gamma <- function(d, n_boot) {
  per_chrom <- data.frame(chrom = character(), nu = numeric(),
                          theta = numeric(), n = integer())
  pooled_dists <- numeric(0)
  for (chrom in d$chrom_map$name) {
    dists <- get_inter_co_dists(d$co[d$co$chrom == chrom, ])
    pooled_dists <- c(pooled_dists, dists)
    fit <- fit_gamma(dists)
    per_chrom <- rbind(per_chrom,
                       data.frame(chrom = chrom, nu = fit$nu,
                                  theta = fit$theta, n = fit$n))
  }
  pooled_fit <- fit_gamma(pooled_dists)
  per_chrom <- rbind(per_chrom,
                     data.frame(chrom = "POOLED", nu = pooled_fit$nu,
                                theta = pooled_fit$theta, n = pooled_fit$n))

  boot_nus <- replicate(n_boot, {
    s <- sample(pooled_dists, replace = TRUE)
    fit_gamma(s)$nu
  })
  boot_nus <- boot_nus[!is.na(boot_nus)]
  ci <- quantile(boot_nus, c(0.025, 0.975))

  list(per_chrom    = per_chrom,
       pooled_fit   = pooled_fit,
       pooled_dists = pooled_dists,
       ci_pooled    = ci)
}

# =============================================================================
# PLOT 6: GAMMA NU PER CHROMOSOME
# =============================================================================
# =============================================================================
# PLOT 5: GAMMA INTERFERENCE - combined histogram + per-chromosome nu
# =============================================================================

plot_gamma_combined <- function(g, out_path, n_cells) {
  suppressPackageStartupMessages(library(ggplot2))

  hist_df <- data.frame(distance_Mb = g$pooled_dists)
  xs <- seq(0.1, max(g$pooled_dists), length.out = 200)
  curve_df <- data.frame(
    x       = xs,
    gamma_y = dgamma(xs, shape = g$pooled_fit$nu, scale = g$pooled_fit$theta),
    exp_y   = dexp(xs, rate = 1 / mean(g$pooled_dists))
  )

  p <- ggplot(hist_df, aes(x = distance_Mb)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 25, fill = "lightblue", colour = "white") +
    geom_line(aes(x = x, y = exp_y, colour = "Exponential null (nu = 1)"),
              data = curve_df, linewidth = 0.8, linetype = "dashed") +
    geom_line(aes(x = x, y = gamma_y, colour = "Gamma fit"),
              data = curve_df, linewidth = 1.0) +
    scale_colour_manual(name = NULL,
                        values = c("Gamma fit" = "darkblue",
                                   "Exponential null (nu = 1)" = "red")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      x        = "Inter-CO distance (Mb)",
      y        = "Density",
      title    = "Crossover interference - Cuscuta epithymum hap1",
      subtitle = sprintf(
        paste0(
          "Gamma model fit to %d inter-CO physical distances (pooled across chromosomes)\n",
          "Pooled nu = %.2f, 95%% CI [%.2f, %.2f] - values >1 indicate positive interference\n",
          "Per-chromosome fits omitted (n_distances per chromosome too small for reliable fit)"
        ),
        length(g$pooled_dists),
        g$pooled_fit$nu, g$ci_pooled[1], g$ci_pooled[2]
      )
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position    = c(0.98, 0.95),
      legend.justification = c(1, 1),
      legend.background  = element_rect(fill = "white", colour = "gray80"),
      legend.key.size    = unit(0.9, "lines"),
      legend.text        = element_text(size = 8.5),
      plot.title         = element_text(size = 11, face = "bold"),
      plot.subtitle      = element_text(size = 8.5, lineheight = 1.1,
                                        colour = "gray30"),
      plot.margin        = margin(6, 6, 6, 6)
    )

  ggsave(out_path, p, width = 7, height = 4.5, units = "in",
         device = cairo_pdf)
}


# =============================================================================
# SUMMARY TABLES
# =============================================================================

write_landscape_summary <- function(d, out_path) {
  out <- data.frame(
    chrom              = d$chrom_map$name,
    size_Mb            = d$chrom_map$size / 1e6,
    n_cos              = sapply(d$chrom_map$name,
                                function(c) sum(d$co$chrom == c)),
    total_genetic_cM   = sapply(d$chrom_map$name,
                                function(c) sum(d$co$chrom == c) * 100 /
                                            d$n_cells),
    mean_cos_per_cell  = sapply(d$chrom_map$name,
                                function(c) sum(d$co$chrom == c) / d$n_cells)
  )
  out$cM_per_Mb <- out$total_genetic_cM / out$size_Mb
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_coc_table <- function(coc_df, out_path) {
  write.table(coc_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_gamma_summary <- function(g, out_path) {
  out <- g$per_chrom
  out$ci_lo <- NA; out$ci_hi <- NA
  out$ci_lo[out$chrom == "POOLED"] <- g$ci_pooled[1]
  out$ci_hi[out$chrom == "POOLED"] <- g$ci_pooled[2]
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# =============================================================================
# MAIN
# =============================================================================

cat("Loading data...\n")
d <- load_data(opt)
cat(sprintf("  %d cells, %d COs, %d chromosomes\n",
            d$n_cells, nrow(d$co), nrow(d$chrom_map)))

cat("Plot 1: landscape...\n")
plot_landscape(d, file.path(opt$out_dir, "01_landscape.pdf"), opt$bin_mb)

cat("Plot 2: Marey maps...\n")
plot_marey(d, file.path(opt$out_dir, "02_marey_maps.pdf"))

cat("Plot 3: CO class proportions...\n")
plot_co_classes(d, file.path(opt$out_dir, "03_co_class_proportions.pdf"))

cat(sprintf("Computing CoC (n_boot=%d)...\n", opt$n_boot))
coc_df <- compute_coc(d, opt$coc_interval_mb, opt$n_boot)

cat(sprintf("Computing Gamma (n_boot=%d)...\n", opt$n_boot * 2))
g <- compute_gamma(d, opt$n_boot * 2)

# Pull out pooled gamma summary for use in CoC plot subtitle
gamma_summary <- list(
  nu    = g$pooled_fit$nu,
  ci_lo = as.numeric(g$ci_pooled[1]),
  ci_hi = as.numeric(g$ci_pooled[2])
)

cat("Plot 4: CoC...\n")
plot_coc(coc_df, file.path(opt$out_dir, "04_coc.pdf"),
         d$n_cells, opt$coc_interval_mb, gamma_summary)

cat("Plot 5: gamma combined...\n")
plot_gamma_combined(g, file.path(opt$out_dir, "05_gamma_interference.pdf"),
                    d$n_cells)


cat("Writing summary tables...\n")
write_landscape_summary(d, file.path(opt$out_dir, "landscape_summary.tsv"))
write_coc_table(coc_df, file.path(opt$out_dir, "coc_table.tsv"))
write_gamma_summary(g, file.path(opt$out_dir, "gamma_interference_summary.tsv"))

cat(sprintf("\nDone. Pooled gamma nu = %.2f, 95%% CI [%.2f, %.2f]\n",
            g$pooled_fit$nu, g$ci_pooled[1], g$ci_pooled[2]))
cat("All outputs in:", opt$out_dir, "\n")
