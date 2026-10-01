#!/usr/bin/env Rscript
# =============================================================================
# talk_landscape.R -- recombination landscape restyled for a talk slide
# =============================================================================
# Plotting only. Crossovers and window rates are recombination_landscape.R's own:
# load_data() and compute_sliding_rate() are read out of that script and run
# unchanged, with plot_landscape()'s windows (5 Mb, 0.5 Mb step;
# cM/Mb = 100 x crossovers / cells / Mb). Run from the pipeline root.
#
#   Rscript workflow/scripts/talk_landscape.R --species Dbinata_hap1 \
#       [--species Dparadoxa_std --in-progress Dparadoxa_std] \
#       [--highlight CHR] [--genes GFF] [--ymax N] [--out results/talk]
#
#   --species      sample_id, repeatable; one panel each, left to right. A sample
#                  with no crossover results yet gets a dashed placeholder panel.
#   --in-progress  sample_id: data drawn, panel flagged "in progress"
#   --highlight    chromosome drawn at 2 pt, in every panel that has it
#                  (SAMPLE=CHR for one panel only)
#   --genes        gene GFF for the summary: one per --species in the same
#                  order, or SAMPLE=GFF
#   --ymax         pin the shared y maximum (default: from the data, printed)
#   --out          output directory [results/talk]
#
# Panels keep one fixed size: REF_SLOTS panels fill FIG_A, fewer make a
# narrower figure with the same panel geometry (slides can build on each other).
#
# Writes to --out:
#   recomb_metaplot.png/.svg          A  one panel per species, 0-100 % of length
#   recomb_stacked_<sample>.png/.svg  B  per-chromosome panels, same y-limits
#   recomb_summary.csv                C  one row per species
# =============================================================================

# ---- parameters -------------------------------------------------------------
LANDSCAPE_R <- "workflow/scripts/recombination_landscape.R"  # functions reused
SAMPLES_CSV <- "config/samples.csv"                            # species, centromere
INPUT_FMT   <- c(co_intervals = "results/crossovers/%s/co_intervals.bed",
                 co_per_cell  = "results/crossovers/%s/co_per_cell.tsv",
                 chrom_map    = "results/cell_data/%s/chrom_map.tsv")
CHECK_FMT   <- "results/landscape/%s/landscape_summary.tsv"    # cross-check only

WIN_MB  <- 5     # plot_landscape() defaults; keep, so rates match 01_landscape.pdf
STEP_MB <- 0.5

N_BINS        <- 100         # relative bins per chromosome
SMOOTH_KERNEL <- c(1, 1, 1)  # 3-bin running mean, the same for every species
YMAX_Q        <- 0.99        # y max: 99th percentile of window rates per species,
YMAX_STEP     <- 1           #   max over species, rounded up to this (cM/Mb)
COLD_FRAC     <- 0.25        # cold zone: mean curve below 25 % of genome mean
COLD_MIN_BINS <- 5           # shorter runs of cold bins are not a zone
OUTER_FRAC    <- 0.10        # chromosome ends, for the end / genome-mean ratio

COL <- c(monocentric = "#D95F02", holocentric = "#1B9E77",
         text = "#334155", text2 = "#64748B",
         zone = "#E9EDF2", rule = "#94A3B8", pale = "#CBD5E1")
TXT <- c(x = "position along chromosome (%)", y = "crossover rate (cM/Mb)",
         cold = "recombination-cold zone", genome = "genome average",
         hot = "hot", pending = "in progress",
         foot = "n = %s pollen cells \u00b7 %s crossovers \u00b7 %d chromosomes")
PT  <- c(axis_title = 15, tick = 12.5, label = 14, cold = 15, header = 16,
         tag = 12, foot = 9, strip = 12)               # font sizes, pt
LW  <- c(chrom = 1, mean = 4, highlight = 2, genome = 1.5, axis = 0.8,
         stacked = 1.5)                                # line widths, pt
FIG_A     <- c(w = 10, h = 4.6)                # inches, for REF_SLOTS panels
REF_SLOTS <- 2
SPACING_PT <- 30                               # gap between panels
FIG_B     <- c(w = 10, row_h = 0.58, ncol = 2) # inches per row of panels
DPI       <- 300

# ---- setup ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
lw  <- function(pt) pt * 96 / 72 / ggplot2::.pt   # stroke in pt -> ggplot linewidth
fsz <- function(pt) pt / ggplot2::.pt             # font in pt -> geom_text size
CAIRO <- if (capabilities("cairo")) "cairo" else getOption("bitmapType")

args  <- commandArgs(trailingOnly = TRUE)
FLAGS <- c("--species", "--in-progress", "--highlight", "--genes", "--ymax", "--out")
if (!length(args) || length(args) %% 2 || !all(args[c(TRUE, FALSE)] %in% FLAGS))
  stop("usage: talk_landscape.R --species SAMPLE [--species SAMPLE] [--in-progress SAMPLE] ",
       "[--highlight CHR] [--genes GFF] [--ymax N] [--out DIR]", call. = FALSE)
arg <- function(f) args[which(args == f & seq_along(args) %% 2 == 1) + 1]

SAMPLES  <- unique(unlist(strsplit(arg("--species"), ",")))
FLAGGED  <- unique(unlist(strsplit(arg("--in-progress"), ",")))
OUT      <- c(arg("--out"), "results/talk")[1]
YMAX_PIN <- as.numeric(c(arg("--ymax"), NA)[1])
if (!length(SAMPLES)) stop("give at least one --species", call. = FALSE)
if (length(setdiff(FLAGGED, SAMPLES)))
  stop("--in-progress ", toString(setdiff(FLAGGED, SAMPLES)), ": not a --species", call. = FALSE)

# SAMPLE=VALUE binds to one panel; a bare value binds by position (genes)
# or to every panel (highlight).
bind <- function(vals, what, by_position) {
  out <- setNames(rep(NA_character_, length(SAMPLES)), SAMPLES)
  for (k in seq_along(vals)) {
    named <- grepl("^[^=/]+=", vals[k])
    s <- if (named) sub("=.*", "", vals[k]) else if (by_position) SAMPLES[k] else "*"
    v <- if (named) sub("^[^=]*=", "", vals[k]) else vals[k]
    if (is.na(s)) stop(what, ": more values than --species", call. = FALSE)
    if (s == "*") { out[] <- v; next }
    if (!s %in% SAMPLES) stop(what, " ", vals[k], ": not a --species", call. = FALSE)
    out[s] <- v
  }
  out
}
HIGHLIGHT <- bind(arg("--highlight"), "--highlight", by_position = FALSE)
GENES     <- bind(arg("--genes"),     "--genes",     by_position = TRUE)

# The pipeline's own functions, parsed out of its script (its main body never runs).
reuse <- function(path, fns) {
  env <- new.env()
  for (e in parse(path, keep.source = FALSE))
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) &&
        as.character(e[[2]]) %in% fns) eval(e, env)
  miss <- setdiff(fns, ls(env))
  if (length(miss)) stop("not found in ", path, ": ", toString(miss), call. = FALSE)
  env
}
RL <- reuse(LANDSCAPE_R, c("load_data", "compute_sliding_rate"))

sheet <- read.csv(SAMPLES_CSV, colClasses = "character", check.names = FALSE)

# ---- per-species data -------------------------------------------------------
resample_rel <- function(rel, rate, n = N_BINS, sub = 20) {
  # average the window curve over each relative bin; before the first and after
  # the last window centre the nearest window (which covers that stretch) holds
  if (length(rel) < 2) return(rep(rate[1], n))
  x <- (seq_len(n * sub) - 0.5) / (n * sub)
  colMeans(matrix(approx(rel, rate, xout = x, rule = 2)$y, nrow = sub))
}
smooth_bins <- function(y, k = SMOOTH_KERNEL) {
  h <- (length(k) - 1) %/% 2
  vapply(seq_along(y), function(i) {
    j <- (i - h):(i + h); ok <- j >= 1 & j <= length(y)
    sum(k[ok] * y[j[ok]]) / sum(k[ok])
  }, numeric(1))
}
longest_run <- function(flag) {        # first and last bin of the longest TRUE run
  r <- rle(flag); end <- cumsum(r$lengths); start <- end - r$lengths + 1
  k <- which(r$values & r$lengths >= COLD_MIN_BINS)
  if (!length(k)) return(NULL)
  k <- k[which.max(r$lengths[k])]
  c(start[k], end[k])
}

load_species <- function(s) {
  r   <- sheet[sheet$sample_id == s, , drop = FALSE]
  cen <- if (nrow(r) && r$centromere[1] %in% c("monocentric", "holocentric")) r$centromere[1] else NA
  sp  <- list(sample = s, flagged = s %in% FLAGGED, centromere = cen,
              species = sub("^([A-Za-z])[a-z]*_", "\\1. ", if (nrow(r)) r$species[1] else s),
              colour = if (is.na(cen)) COL[["text2"]] else COL[[cen]])
  f <- setNames(sprintf(INPUT_FMT, s), names(INPUT_FMT))
  for (k in names(f)) cat(sprintf("%-16s %-12s %s  %s\n", s, k, f[[k]],
      if (file.exists(f[[k]])) format(file.mtime(f[[k]]), "%Y-%m-%d %H:%M") else "MISSING"))
  if (!all(file.exists(f)) || any(file.size(f) == 0)) {
    cat(sprintf("%-16s no crossover results -> placeholder panel\n", s))
    return(c(sp, ok = FALSE))
  }
  d  <- RL$load_data(as.list(f))
  cm <- d$chrom_map
  short <- cm$size < WIN_MB * 1e6
  if (any(short)) cat(sprintf("%-16s not drawn (< %g Mb): %s\n", s, WIN_MB, toString(cm$name[short])))

  win <- do.call(rbind, lapply(which(!short), function(i) {
    w <- RL$compute_sliding_rate(d$co$mid[d$co$chrom == cm$name[i]], cm$size[i],
                                 d$n_cells, WIN_MB, STEP_MB)
    data.frame(chrom = cm$name[i], size = cm$size[i], centre = w$centres, rate = w$rates)
  }))
  curves <- do.call(rbind, lapply(split(win, factor(win$chrom, unique(win$chrom))), function(w)
    data.frame(chrom = w$chrom[1], x = (seq_len(N_BINS) - 0.5) * 100 / N_BINS,
               y = smooth_bins(resample_rel(w$centre / w$size[1], w$rate)))))
  mean_y <- as.numeric(tapply(curves$y, curves$x, mean))
  genome <- nrow(d$co) * 100 / d$n_cells / (sum(cm$size) / 1e6)   # as plot_landscape()
  run    <- longest_run(mean_y < COLD_FRAC * genome)
  zone   <- if (is.null(run)) NULL else c(lo = run[1] - 1, hi = run[2]) * 100 / N_BINS
  cat(sprintf("%-16s %d cells, %d crossovers, %d chromosomes; genome mean %.2f cM/Mb; cold zone %s\n",
              s, d$n_cells, nrow(d$co), nrow(cm), genome,
              if (is.null(zone)) "none" else sprintf("%g-%g %% (mean < %.2f cM/Mb)",
                                                     zone[["lo"]], zone[["hi"]], COLD_FRAC * genome)))

  # summary numbers (crossovers by interval midpoint, as the window rates)
  size  <- setNames(cm$size, cm$name)
  rel   <- d$co$mid / size[d$co$chrom]
  outer <- rel < OUTER_FRAC | rel >= 1 - OUTER_FRAC
  outer_rate <- sum(outer, na.rm = TRUE) * 100 / d$n_cells / (2 * OUTER_FRAC * sum(cm$size) / 1e6)
  in_zone <- function(p) if (is.null(zone)) NA else
    100 * mean(p * 100 >= zone[["lo"]] & p * 100 < zone[["hi"]], na.rm = TRUE)
  genes_pct <- NA; n_genes <- NA
  if (!is.na(GENES[[s]])) {
    if (!file.exists(GENES[[s]])) {
      cat(sprintf("%-16s --genes %s: file not found, gene column left empty\n", s, GENES[[s]]))
    } else {
      g <- read.delim(GENES[[s]], header = FALSE, comment.char = "#", quote = "", fill = TRUE,
                      colClasses = c("character", "NULL", "character", "numeric", "numeric",
                                     rep("NULL", 4)))
      g <- g[g$V3 %in% "gene" & g$V1 %in% cm$name, ]
      n_genes <- nrow(g)
      if (n_genes) genes_pct <- in_zone(((g$V4 + g$V5) / 2) / size[g$V1])
      else cat(sprintf("%-16s no 'gene' rows on this sample's chromosomes in %s\n", s, GENES[[s]]))
    }
  }
  summary <- data.frame(
    sample = s, species = sp$species, centromere = cen,
    status = if (sp$flagged) "in progress" else "final",
    n_cells = d$n_cells, n_crossovers = nrow(d$co), n_chromosomes = nrow(cm),
    genome_mean_cM_per_Mb = genome,
    cold_zone_from_pct = if (is.null(zone)) NA else zone[["lo"]],
    cold_zone_to_pct   = if (is.null(zone)) NA else zone[["hi"]],
    pct_length_in_cold_zone     = if (is.null(zone)) 0 else zone[["hi"]] - zone[["lo"]],
    pct_crossovers_in_cold_zone = in_zone(rel),
    outer10_cM_per_Mb = outer_rate, outer10_over_genome_mean = outer_rate / genome,
    pct_genes_in_cold_zone = genes_pct, n_genes = n_genes)

  # cross-check against the pipeline's own landscape table, if present
  chk <- sprintf(CHECK_FMT, s)
  if (file.exists(chk)) {
    ls <- read.delim(chk)
    same <- identical(as.character(ls$chrom), cm$name) &&
            all(ls$n_cos == as.numeric(table(factor(d$co$chrom, cm$name))))
    cat(sprintf("%-16s crossovers per chromosome %s %s\n", s,
                if (same) "match" else "DIFFER FROM", chk))
  }
  c(sp, list(ok = TRUE, win = win, curves = curves, genome = genome, zone = zone,
             mean = data.frame(x = unique(curves$x), y = mean_y), summary = summary,
             foot = sprintf(TXT[["foot"]], format(d$n_cells, big.mark = ","),
                            format(nrow(d$co), big.mark = ","), nrow(cm))))
}

SP <- lapply(SAMPLES, load_species); names(SP) <- SAMPLES
OK <- SP[vapply(SP, `[[`, logical(1), "ok")]
if (!length(OK) && is.na(YMAX_PIN)) stop("no crossover results for any --species", call. = FALSE)

q99  <- vapply(OK, function(sp) quantile(sp$win$rate, YMAX_Q, names = FALSE), numeric(1))
YMAX <- if (is.na(YMAX_PIN)) ceiling(max(q99) / YMAX_STEP) * YMAX_STEP else YMAX_PIN
cat(sprintf("y max = %g cM/Mb  [%s 99th pct of %g-Mb window rates: %s]\n", YMAX,
            if (is.na(YMAX_PIN)) "rounded up from" else "pinned by --ymax; data", WIN_MB,
            paste(sprintf("%s %.2f", names(q99), q99), collapse = ", ")))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
stack_df  <- function(f) do.call(rbind, lapply(OK, f))
add_panel <- function(df) if (is.null(df) || !nrow(df)) NULL else
  transform(df, panel = factor(sample, levels = SAMPLES))

clip_marks <- function(df, xcol = "x") {   # one marker per run above YMAX, at its peak
  if (is.null(df) || !nrow(df)) return(NULL)
  do.call(rbind, lapply(split(df, list(df$sample, df$chrom), drop = TRUE), function(d) {
    d <- d[order(d[[xcol]]), ]
    over <- d$y > YMAX
    if (!any(over)) return(NULL)
    id <- cumsum(c(TRUE, diff(over) != 0))
    do.call(rbind, lapply(unique(id[over]), function(k) {
      dd <- d[id == k, ]; dd[which.max(dd$y), ] }))
  }))
}

save_both <- function(g, base, w, h) {
  png(paste0(base, ".png"), width = w, height = h, units = "in", res = DPI, bg = "white",
      type = CAIRO)
  grid.newpage(); grid.draw(g); invisible(dev.off())
  if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite(paste0(base, ".svg"), w, h, bg = "white")
  else svg(paste0(base, ".svg"), w, h, bg = "white")   # cairo svg: text as outlines
  grid.newpage(); grid.draw(g); invisible(dev.off())
  cat(sprintf("wrote %s.png + .svg  (%.2f x %.2f in)\n", normalizePath(base, mustWork = FALSE), w, h))
}

theme_talk <- function() {
  theme_classic(base_size = PT[["tick"]], base_family = "sans") +
    theme(axis.line  = element_line(colour = COL[["rule"]], linewidth = lw(LW[["axis"]])),
          axis.ticks = element_line(colour = COL[["rule"]], linewidth = lw(LW[["axis"]])),
          axis.ticks.length = unit(4, "pt"),
          axis.text  = element_text(colour = COL[["text2"]], size = PT[["tick"]]),
          axis.title = element_text(colour = COL[["text"]], size = PT[["axis_title"]]),
          axis.title.x = element_text(margin = margin(t = 6)),
          axis.title.y = element_text(margin = margin(r = 8)),
          strip.background = element_blank(),
          panel.background = element_rect(fill = "white", colour = NA),
          plot.background  = element_rect(fill = "white", colour = NA))
}

# ---- A: metaplot, one panel per species -------------------------------------
curves <- stack_df(function(sp) cbind(sample = sp$sample, colour = sp$colour, sp$curves))
means  <- stack_df(function(sp) cbind(sample = sp$sample, colour = sp$colour, chrom = "mean", sp$mean))
hl     <- if (is.null(curves)) NULL else curves[(curves$chrom == HIGHLIGHT[curves$sample]) %in% TRUE, ]
for (s in names(OK)) if (!is.na(HIGHLIGHT[[s]]) && !any(hl$sample == s))
  cat(sprintf("%-16s --highlight %s: no such chromosome here\n", s, HIGHLIGHT[[s]]))
zones  <- stack_df(function(sp) if (!is.null(sp$zone))
  data.frame(sample = sp$sample, lo = sp$zone[["lo"]], hi = sp$zone[["hi"]], top = YMAX))
genome <- stack_df(function(sp) data.frame(sample = sp$sample, rate = sp$genome))
hot    <- stack_df(function(sp) if (!is.null(sp$zone)) {
  ends <- list(sp$mean$x >= 4 & sp$mean$x <= sp$zone[["lo"]] - 2,      # left hot end
               sp$mean$x <= 96 & sp$mean$x >= sp$zone[["hi"]] + 2)     # right hot end
  do.call(rbind, lapply(ends, function(e) if (any(e)) {
    m  <- sp$mean[e, ]; m <- m[which.max(m$y), ]                       # its peak
    yc <- min(m$y, YMAX); txt <- min(yc + 0.22 * YMAX, 0.93 * YMAX)
    data.frame(sample = sp$sample, colour = sp$colour, x = m$x, tip = yc + 0.05 * YMAX,
               txt = txt, tail = if (txt - 0.04 * YMAX - (yc + 0.05 * YMAX) > 0.05 * YMAX)
                                   txt - 0.04 * YMAX else NA)
  }))
})
holder  <- data.frame(sample = setdiff(SAMPLES, names(OK)))
gl <- if (length(OK)) genome[genome$sample == names(OK)[1], ] else NULL   # labelled once
if (!is.null(gl)) gl$x <- if (!is.null(OK[[1]]$zone)) OK[[1]]$zone[["hi"]] - 1.5 else 98.5
clipA <- rbind(clip_marks(curves), clip_marks(means))
frame <- function(df) geom_rect(data = add_panel(transform(df, top = YMAX)),
                                aes(xmin = 0, xmax = 100, ymin = 0, ymax = top), fill = NA,
                                colour = COL[["pale"]], linetype = "22", linewidth = lw(1.2))

pA <- ggplot() +
  { if (!is.null(zones)) geom_rect(data = add_panel(zones), aes(xmin = lo, xmax = hi, ymin = 0, ymax = top),
                                   fill = COL[["zone"]]) } +
  { if (nrow(holder)) list(frame(holder),
      geom_text(data = add_panel(transform(holder, y = YMAX / 2)), aes(x = 50, y = y),
                label = TXT[["pending"]], colour = COL[["text2"]], size = fsz(PT[["cold"]]),
                family = "sans")) } +
  { if (!is.null(genome)) geom_hline(data = add_panel(genome), aes(yintercept = rate),
                                     colour = COL[["rule"]], linetype = "22",
                                     linewidth = lw(LW[["genome"]])) } +
  { if (!is.null(curves)) geom_line(data = add_panel(curves),
                                    aes(x, pmin(y, YMAX), group = chrom, colour = colour),
                                    linewidth = lw(LW[["chrom"]]), alpha = 0.3) } +
  { if (!is.null(hl) && nrow(hl)) geom_line(data = add_panel(hl),
                                            aes(x, pmin(y, YMAX), group = chrom, colour = colour),
                                            linewidth = lw(LW[["highlight"]])) } +
  { if (!is.null(means)) geom_line(data = add_panel(means), aes(x, pmin(y, YMAX), colour = colour),
                                   linewidth = lw(LW[["mean"]]), lineend = "round") } +
  { if (!is.null(clipA)) geom_point(data = add_panel(clipA), aes(x, YMAX, colour = colour),
                                    shape = 17, size = 1.8) } +
  { if (!is.null(zones)) geom_text(data = add_panel(zones),
                                   aes(x = (lo + hi) / 2, y = 0.62 * YMAX,
                                       label = ifelse(hi - lo < 60, sub("-", "-\n", TXT[["cold"]]), TXT[["cold"]])),
                                   colour = COL[["text"]], fontface = "bold", lineheight = 0.9,
                                   size = fsz(PT[["cold"]]), family = "sans") } +
  { if (!is.null(hot)) list(
      if (any(!is.na(hot$tail))) geom_segment(data = add_panel(hot[!is.na(hot$tail), ]),
                   aes(x = x, xend = x, y = tail, yend = tip, colour = colour),
                   linewidth = lw(1.6), arrow = arrow(length = unit(6, "pt"), type = "closed")),
      geom_text(data = add_panel(hot), aes(x = x, y = txt, label = TXT[["hot"]], colour = colour),
                vjust = 0, fontface = "bold", size = fsz(PT[["label"]]), family = "sans")) } +
  { if (!is.null(gl)) geom_text(data = add_panel(gl), aes(x = x, y = rate + 0.03 * YMAX),
                                label = TXT[["genome"]], hjust = 1, vjust = 0, colour = COL[["text2"]],
                                fontface = "bold", size = fsz(PT[["label"]]), family = "sans") } +
  scale_colour_identity() +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, YMAX), breaks = Filter(function(b) b <= YMAX, pretty(c(0, YMAX), n = 4)),
                     expand = expansion(mult = c(0.015, 0))) +
  coord_cartesian(clip = "off") +
  facet_wrap(~panel, nrow = 1, drop = FALSE) +
  labs(x = TXT[["x"]], y = TXT[["y"]]) +
  theme_talk() +
  theme(strip.text = element_blank(), panel.spacing.x = unit(SPACING_PT, "pt"),
        plot.margin = margin(4, 16, 2, 4))

# headers above and footnotes below each panel, aligned to its left edge
gA  <- ggplotGrob(pA)
pan <- gA$layout[grepl("^panel", gA$layout$name), ]
pan <- pan[order(pan$l), ]
gA  <- gtable::gtable_add_rows(gA, unit(PT[["header"]] * 2.1, "pt"), pos = 0)
gA  <- gtable::gtable_add_rows(gA, unit(PT[["foot"]] * 2.4, "pt"), pos = -1)
for (k in seq_along(SAMPLES)) {
  sp  <- SP[[SAMPLES[k]]]
  a   <- textGrob(sp$species, x = 0, y = 0.3, hjust = 0, vjust = 0,
                  gp = gpar(fontsize = PT[["header"]], fontface = "italic", col = COL[["text"]]))
  b   <- textGrob(if (is.na(sp$centromere)) "" else sp$centromere,
                  x = unit(0, "npc") + grobWidth(a) + unit(10, "pt"), y = 0.3, hjust = 0, vjust = 0,
                  gp = gpar(fontsize = PT[["header"]] - 1, fontface = "bold", col = sp$colour))
  hdr <- list(a, b)
  if (isTRUE(sp$flagged) && isTRUE(sp$ok)) {           # "in progress" badge
    x0  <- unit(0, "npc") + grobWidth(a) + grobWidth(b) + unit(24, "pt")
    tag <- textGrob(TXT[["pending"]], x = x0 + unit(7, "pt"), y = unit(0.3, "npc") + unit(1.5, "pt"),
                    hjust = 0, vjust = 0, gp = gpar(fontsize = PT[["tag"]], fontface = "bold",
                                                     col = COL[["text2"]]))
    pill <- roundrectGrob(x = x0, y = unit(0.3, "npc") - unit(2.5, "pt"), just = c(0, 0),
                          width = grobWidth(tag) + unit(14, "pt"),
                          height = unit(PT[["tag"]] + 6, "pt"), r = unit(0.5, "snpc"),
                          gp = gpar(fill = COL[["zone"]], col = COL[["pale"]], lwd = 1))
    hdr <- c(hdr, list(pill, tag))
  }
  gA <- gtable::gtable_add_grob(gA, do.call(grobTree, hdr), t = 1, l = pan$l[k], r = pan$r[k],
                                clip = "off", name = paste0("header-", k))
  if (isTRUE(sp$ok))
    gA <- gtable::gtable_add_grob(gA, textGrob(sp$foot, x = 0, y = 0.35, hjust = 0,
                                               gp = gpar(fontsize = PT[["foot"]], col = COL[["text2"]])),
                                  t = nrow(gA), l = pan$l[k], r = pan$r[k], clip = "off",
                                  name = paste0("foot-", k))
}

# same panel size whatever the number of panels: REF_SLOTS panels fill FIG_A
png(tmp <- tempfile(fileext = ".png"), width = FIG_A[["w"]], height = FIG_A[["h"]],
    units = "in", res = DPI, type = CAIRO)
fixed <- sum(vapply(seq_along(gA$widths), function(i)
  if (unitType(gA$widths[i]) == "null") 0 else convertWidth(gA$widths[i], "in", TRUE), numeric(1)))
invisible(dev.off()); unlink(tmp)
n_pan   <- nrow(pan)
spacing <- SPACING_PT / 72
panel_w <- (FIG_A[["w"]] - (fixed - (n_pan - 1) * spacing) - (REF_SLOTS - 1) * spacing) / REF_SLOTS
save_both(gA, file.path(OUT, "recomb_metaplot"), fixed + n_pan * panel_w, FIG_A[["h"]])

# ---- B: stacked per-chromosome backup, same y-limits and colour -------------
for (sp in OK) {
  w <- transform(sp$win, sample = sp$sample, pos = centre / 1e6, y = rate)
  short <- sub("_hap[0-9]+$", "", w$chrom)                   # chr1_hap1 -> chr1
  if (length(unique(short)) != length(unique(w$chrom))) short <- w$chrom
  w$lab <- factor(short, unique(short))
  clipB  <- clip_marks(w, xcol = "pos")
  nrow_b <- ceiling(nlevels(w$lab) / FIG_B[["ncol"]])
  pB <- ggplot(w, aes(pos, pmin(y, YMAX))) +
    geom_hline(yintercept = sp$genome, colour = COL[["rule"]], linetype = "22", linewidth = lw(1)) +
    geom_line(colour = sp$colour, linewidth = lw(LW[["stacked"]])) +
    { if (!is.null(clipB)) geom_point(data = clipB, aes(pos, YMAX), colour = sp$colour,
                                      shape = 17, size = 1.4) } +
    facet_wrap(~lab, ncol = FIG_B[["ncol"]], dir = "v", strip.position = "left") +
    scale_x_continuous(limits = c(0, max(w$size) / 1e6), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, YMAX), breaks = c(0, YMAX), expand = expansion(mult = c(0.03, 0))) +
    coord_cartesian(clip = "off") +
    labs(x = "position (Mb)", y = TXT[["y"]],
         caption = paste0(if (sp$flagged) paste0(TXT[["pending"]], " \u00b7 ") else "",
                          sp$species, " \u00b7 ", sp$foot)) +
    theme_talk() +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, colour = COL[["text"]],
                                           size = PT[["strip"]], margin = margin(r = 4)),
          axis.text = element_text(size = PT[["strip"]] - 1),
          panel.spacing.y = unit(14, "pt"), panel.spacing.x = unit(26, "pt"),
          plot.caption = element_text(hjust = 0, size = PT[["foot"]], colour = COL[["text2"]]),
          plot.caption.position = "plot", plot.margin = margin(6, 16, 4, 4))
  save_both(ggplotGrob(pB), file.path(OUT, paste0("recomb_stacked_", sp$sample)),
            FIG_B[["w"]], FIG_B[["row_h"]] * nrow_b + 1.3)
}

# ---- C: numbers for the notes -----------------------------------------------
if (length(OK)) {
  tab <- do.call(rbind, lapply(OK, `[[`, "summary"))
  tab$y_max_cM_per_Mb <- YMAX
  num <- vapply(tab, is.double, logical(1))
  tab[num] <- lapply(tab[num], round, 3)
  write.csv(tab, file.path(OUT, "recomb_summary.csv"), row.names = FALSE)
  cat(sprintf("wrote %s\n", normalizePath(file.path(OUT, "recomb_summary.csv"))))
  print(t(tab), quote = FALSE)
}
