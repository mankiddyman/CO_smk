#!/usr/bin/env Rscript
# CO identification in gamete nuclei
# Original by Hequan Sun (sunhequan@gmail.com), MPIPZ
# Modified by Meng Zhang (mzhang@mpipz.mpg.de) for Rhynchospora scRNA-seq pollen
# Adapted by Aaryan Bhatia for Cuscuta epithymum:
#   - Reads chrom_map.tsv (internal_id, name, size) instead of genome.fai
#   - Filters rows by scaffold name string (not integer chromosome ID)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "Per-cell TSV with chrom, pos, ref, ref_count, alt, alt_count"),
  make_option(c("-p", "--prefix"), type = "character", default = "out",
              help = "Output filename prefix [default: %default]"),
  make_option(c("-g", "--chrom_map"), type = "character", default = NULL,
              help = "chrom_map.tsv (internal_id<TAB>name<TAB>size)"),
  make_option(c("-o", "--outpath"), type = "character", default = getwd(),
              help = "Output directory [default: cwd]"),
  make_option(c("-c", "--cell_markers"), type = "integer", default = 200,
              help = "Min informative markers required [default: %default]"),
  make_option(c("-s", "--block_size"), type = "integer", default = 2000000,
              help = "Min block size in bp [default: %default]"),
  make_option(c("-n", "--marker_num"), type = "integer", default = 8,
              help = "Min markers per block [default: %default]"),
  make_option("--baseAF", type = "double", default = 0.3,
              help = "Smoothing pass-1 threshold [default: %default]"),
  make_option("--windowAF", type = "double", default = 0.4,
              help = "Smoothing pass-2 threshold [default: %default]"),
  make_option("--genotype", type = "double", default = 0.2,
              help = "Block genotype threshold [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$chrom_map)) {
  stop("Need --input and --chrom_map")
}

dir.create(opt$outpath, showWarnings = FALSE, recursive = TRUE)

### functions: smooth, blocks, breakpoints (unchanged from Meng) ###############

allele_count_smoother <- function(acnt_chr) {
  this_dim <- dim(acnt_chr)
  if (this_dim[2] != 6) stop("Wrong dim of acnt_chr")
  acnt_plus_af <- cbind(acnt_chr,
                        acnt_chr[, 4] / (acnt_chr[, 4] + acnt_chr[, 6]),
                        rep(0, this_dim[1]))
  colnames(acnt_plus_af) <- c(paste("V", 1:6, sep = ""), "rawAF", "smtAF")
  adj <- 2
  for (row in 1:this_dim[1]) {
    surrounding <- max(row - adj, 1):min(row + adj, this_dim[1])
    if (row == 1) {
      right <- (row + 1):min(row + adj, this_dim[1])
      win_af <- sum(acnt_plus_af[right, 7] > 0) / (length(surrounding) - 1)
    } else if (row == this_dim[1]) {
      left <- max(row - adj, 1):(row - 1)
      win_af <- sum(acnt_plus_af[left, 7] > 0) / (length(surrounding) - 1)
    } else {
      left <- max(row - adj, 1):(row - 1)
      right <- (row + 1):min(row + adj, this_dim[1])
      win_af <- (sum(acnt_plus_af[left, 7]) + sum(acnt_plus_af[right, 7])) /
                (length(surrounding) - 1)
    }
    if (win_af > 1 - opt$baseAF)        acnt_plus_af[row, 8] <- 1
    else if (win_af < opt$baseAF)       acnt_plus_af[row, 8] <- 0
    else                                 acnt_plus_af[row, 8] <- acnt_plus_af[row, 7]
  }
  for (row in 1:this_dim[1]) {
    surrounding <- max(row - adj, 1):min(row + adj, this_dim[1])
    if (row == 1) {
      right <- (row + 1):min(row + adj, this_dim[1])
      win_af <- sum(acnt_plus_af[right, 8] > 0) / (length(surrounding) - 1)
    } else if (row == this_dim[1]) {
      left <- max(row - adj, 1):(row - 1)
      win_af <- sum(acnt_plus_af[left, 8] > 0) / (length(surrounding) - 1)
    } else {
      left <- max(row - adj, 1):(row - 1)
      right <- (row + 1):min(row + adj, this_dim[1])
      win_af <- (sum(acnt_plus_af[left, 8] > 0) + sum(acnt_plus_af[right, 8] > 0)) /
                (length(surrounding) - 1)
    }
    if (opt$windowAF <= win_af && win_af <= 1 - opt$windowAF) {
      acnt_plus_af[row, 8] <- acnt_plus_af[row, 7]
    }
  }
  acnt_plus_af
}

get_genotype_block <- function(acnt_smoothed) {
  this_dim <- dim(acnt_smoothed)
  this_block <- matrix(NA, 1, 4)
  this_block[1] <- acnt_smoothed[1, 2]
  this_block[2] <- acnt_smoothed[1, 2]
  this_block[3] <- -1
  this_block[4] <- 1
  blocks <- matrix(NA, 0, 4)
  if (acnt_smoothed[1, 8] >= 1 - opt$genotype)       this_block[3] <- 1
  else if (acnt_smoothed[1, 8] <= opt$genotype)      this_block[3] <- 0
  else                                               this_block[3] <- -1
  if (this_dim[1] > 1) {
    for (row in 2:this_dim[1]) {
      if (acnt_smoothed[row, 8] != acnt_smoothed[row - 1, 8]) {
        blocks <- rbind(blocks, this_block)
        this_block[1] <- acnt_smoothed[row, 2]
        this_block[2] <- acnt_smoothed[row, 2]
        this_block[3] <- -1
        this_block[4] <- 1
        if (acnt_smoothed[row, 8] >= 1 - opt$genotype)       this_block[3] <- 1
        else if (acnt_smoothed[row, 8] <= opt$genotype)      this_block[3] <- 0
        else                                                 this_block[3] <- -1
      } else {
        this_block[2] <- acnt_smoothed[row, 2]
        this_block[4] <- this_block[4] + 1
      }
    }
  }
  blocks <- rbind(blocks, this_block)
  colnames(blocks) <- c("sta", "end", "genotype", "markernum")
  blocks
}

filter_blocks <- function(blocks, min_block_size, min_marker) {
  new_blocks <- matrix(NA, 0, 4)
  colnames(new_blocks) <- c("sta", "end", "genotype", "markernum")
  for (row in 1:dim(blocks)[1]) {
    if (blocks[row, 2] - blocks[row, 1] >= min_block_size &&
        blocks[row, 4] >= min_marker) {
      new_blocks <- rbind(new_blocks, blocks[row, ])
    }
  }
  new_blocks
}

get_breakpoints <- function(blocks, chr_id) {
  breakpoints <- matrix(NA, 0, 6)
  if (dim(blocks)[1] > 1) {
    marker_acc <- blocks[1, 4]
    for (row in 2:dim(blocks)[1]) {
      if (blocks[row, 3] != blocks[row - 1, 3]) {
        tmp <- matrix(NA, 0, 6)
        tmp[1] <- chr_id
        tmp[2] <- blocks[row - 1, 2]
        tmp[3] <- blocks[row, 1]
        tmp[4] <- blocks[row - 1, 3]
        tmp[5] <- blocks[row, 3]
        tmp[6] <- marker_acc
        breakpoints <- rbind(breakpoints, tmp)
        marker_acc <- blocks[row, 4]
      } else {
        marker_acc <- marker_acc + blocks[row, 4]
      }
    }
  }
  breakpoints
}

fine_breakpoints <- function(breakpoints, acnt_smoothed) {
  bp_updated <- breakpoints
  for (bp in 1:dim(breakpoints)[1]) {
    interval <- acnt_smoothed[acnt_smoothed[, 2] >= breakpoints[bp, 2] &
                              acnt_smoothed[, 2] <= breakpoints[bp, 3], ]
    seq <- interval[, 8]
    seq_len <- length(seq)
    max_score <- 0
    max_pos <- 1
    for (pos in 1:(seq_len - 1)) {
      l0 <- sum(seq[1:pos] == 0); l1 <- sum(seq[1:pos] == 1)
      r0 <- sum(seq[(pos + 1):seq_len] == 0); r1 <- sum(seq[(pos + 1):seq_len] == 1)
      af0L <- l0 / pos; af1L <- l1 / pos
      af0R <- r0 / (seq_len - pos); af1R <- r1 / (seq_len - pos)
      score <- max(af0L * af1R, af1L * af0R)
      if (score > max_score) { max_score <- score; max_pos <- pos }
    }
    bp_updated[bp, 2] <- interval[max_pos, 2]
    bp_updated[bp, 3] <- interval[max_pos + 1, 2]
  }
  bp_updated
}

fine_blocks <- function(bp_updated, final_blocks) {
  for (bl in 2:dim(final_blocks)[1]) {
    final_blocks[bl - 1, 2] <- bp_updated[bl - 1, 2]
    final_blocks[bl, 1]     <- bp_updated[bl - 1, 3]
  }
  final_blocks
}

make_final_blocks <- function(filtered, chr_size) {
  final_blocks <- matrix(NA, 0, 4)
  this_block <- matrix(NA, 1, 4)
  this_block[1] <- 1; this_block[2] <- 1
  acc <- 0
  row <- 1
  if (dim(filtered)[1] > 1) {
    for (row in 2:dim(filtered)[1]) {
      if (filtered[row, 3] != filtered[row - 1, 3]) {
        this_block[2] <- filtered[row - 1, 2]
        this_block[3] <- filtered[row - 1, 3]
        this_block[4] <- acc
        final_blocks <- rbind(final_blocks, this_block)
        this_block[1] <- filtered[row, 1]
        this_block[2] <- filtered[row, 2]
      } else {
        acc <- acc + filtered[row, 4]
      }
    }
  }
  if (dim(filtered)[1] > 0) {
    this_block[2] <- filtered[row, 2]
    this_block[3] <- filtered[row, 3]
    this_block[4] <- acc
    final_blocks <- rbind(final_blocks, this_block)
    if (final_blocks[length(final_blocks[, 1]), 2] < chr_size) {
      final_blocks[length(final_blocks[, 1]), 2] <- chr_size
    }
  }
  final_blocks
}

### main ######################################################################

chrom_map <- read.table(opt$chrom_map, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
chr_num <- nrow(chrom_map)
chrsize <- chrom_map$size
min_chr_idx <- which.min(chrsize)
max_chr_idx <- which.max(chrsize)

barcode <- opt$prefix
cat("Info: barcode", barcode, "\n")

pdf(file.path(opt$outpath, paste0(barcode, "_co.pdf")),
    family = "Helvetica", height = 11.7, width = 8.3)
par(mai = c(0.4, 1, 0.1, 0.7))
m <- cbind(1:chr_num)
layout(m, heights = rep(1, chr_num))

acnt <- read.table(opt$input, header = FALSE, sep = "\t",
                   stringsAsFactors = FALSE)
informative <- acnt$V4 > 0 | acnt$V6 > 0

if (sum(informative) >= opt$cell_markers) {
  acnt_valid <- acnt[informative, ]
  total_smoothed <- 0

  for (chr_idx in 1:chr_num) {
    chrom_name <- chrom_map$name[chr_idx]
    this_size  <- chrsize[chr_idx]

    acnt_chr <- acnt_valid[acnt_valid$V1 == chrom_name, ]
    if (nrow(acnt_chr) < 2) {
      # too few markers — empty plot
      plot.new()
      title(main = paste(chrom_name, "(too few markers)"), cex.main = 0.8)
      next
    }

    acnt_chr_smoothed <- allele_count_smoother(acnt_chr)
    smt <- acnt_chr_smoothed[acnt_chr_smoothed[, 7] != acnt_chr_smoothed[, 8], ]
    total_smoothed <- total_smoothed + length(smt$V1)

    blocks <- get_genotype_block(acnt_chr_smoothed)
    filtered <- filter_blocks(blocks, opt$block_size, opt$marker_num)

    if (nrow(filtered) == 0) {
      # no blocks pass filter — empty plot
      plot.new()
      title(main = paste(chrom_name, "(no blocks pass filter)"), cex.main = 0.8)
      next
    }

    breakpoints <- get_breakpoints(filtered, chr_idx)
    if (nrow(breakpoints) > 0) {
      bp_updated <- fine_breakpoints(breakpoints, acnt_chr_smoothed)
    } else {
      bp_updated <- breakpoints
    }
    final_blocks2 <- make_final_blocks(filtered, this_size)
    if (nrow(breakpoints) > 0) {
      final_blocks <- fine_blocks(bp_updated, final_blocks2)
    } else {
      final_blocks <- final_blocks2
    }

    ## visualize
    this_ylim <- 20
    plot(acnt_chr[acnt_chr$V4 > 0, 2],
         acnt_chr[acnt_chr$V4 > 0, 4],
         col = "red", type = "h", cex.lab = 1.2,
         ylim = c(-this_ylim, this_ylim),
         xlim = c(1, max(chrsize)),
         axes = FALSE, xlab = "", ylab = "Allele count")
    axis(1, cex.axis = 1.0,
         at = c(seq(0, this_size, 1e7), this_size),
         labels = round(c(seq(0, this_size, 1e7), this_size) / 1e6))
    axis(2, cex.axis = 1.0,
         at = c(-this_ylim, 0, this_ylim),
         labels = c(this_ylim, 0, this_ylim))
    points(acnt_chr[acnt_chr$V6 > 0, 2],
           acnt_chr[acnt_chr$V6 > 0, 6] * -1,
           col = "blue", type = "h")
    lines(c(0, this_size), c(0, 0), col = "black", lty = 2)
    mtext(paste0("(Mb ", chrom_name, ")\n"),
          at = this_size + 6.5e6, padj = 0.1, side = 1, line = 2.0, cex = 0.7)
    rect(1, this_ylim - 2, this_size, this_ylim + 3, col = "gray", border = NA)

    if (nrow(final_blocks) > 0) {
      rect(as.vector(final_blocks[, 1]), this_ylim - 2,
           as.vector(final_blocks[, 2]), this_ylim + 3,
           col = ifelse(final_blocks[, 3] == 1, "red",
                  ifelse(final_blocks[, 3] == 0, "blue", "purple")),
           border = NA)
    }

    if (chr_idx == min_chr_idx) {
      legend(this_size + 0.4 * (chrsize[max_chr_idx] - this_size), 10,
             pch = c(15, 15), col = c("red", "blue"),
             legend = c("genotype A", "genotype B"),
             horiz = FALSE, border = "gray", cex = 1.0)
    }

    ## write outputs
    bp_out <- file.path(opt$outpath, paste0(barcode, "_co_pred.txt"))
    bk_out <- file.path(opt$outpath, paste0(barcode, "_co_block_pred.txt"))
    append_flag <- chr_idx > 1
    # Replace integer chr_idx with scaffold name in breakpoints output
    if (nrow(bp_updated) > 0) {
      bp_named <- bp_updated
      bp_named[, 1] <- chrom_name  # write scaffold name, not internal_id
      write.table(bp_named, file = bp_out, quote = FALSE,
                  row.names = FALSE, col.names = FALSE, append = append_flag)
    }
    blocks_named <- cbind(rep(chrom_name, nrow(final_blocks)), final_blocks)
    write.table(blocks_named, file = bk_out, quote = FALSE,
                row.names = FALSE, col.names = FALSE, append = append_flag)
  }

  cat("Info: total smoothed sites:", total_smoothed,
      "out of", length(acnt_valid$V1),
      "for barcode", barcode, "\n")
} else {
  par(mai = c(1, 1, 0.1, 0.7))
  plot(1, 1, col = "red", type = "h", axes = FALSE,
       xlab = paste0("Markers too few (", sum(informative),
                     " < ", opt$cell_markers, ") -- no CO calls"),
       ylab = "")
  cat("Info: skipping", barcode, "-- only", sum(informative), "markers\n")
}
dev.off()
