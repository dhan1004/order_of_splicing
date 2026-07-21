#!/usr/bin/env Rscript
# =============================================================================
# plot_gene_splicing_lib.R
#
# Refactored core of plot_gene_splicing_order.R: the same intron reconstruction,
# arc-building, colours and layer order, but wrapped as a callable function that
# RETURNS a ggplot (no ggsave, no commandArgs, no file reading).
#
# Plotting logic is preserved verbatim from the original Sections 4-7 so output
# is identical to the standalone script. Differences: (a) title takes an
# optional suffix so paired panels can be labelled "GENE (human)" / "(mouse)";
# (b) the per-gene TSV pair table and console printing are dropped (not wanted
# for a scan of thousands); (c) legend/axis kept but title sizing eased so
# panels read at contact-sheet scale.
#
# source() this from a driver; see render_scan_contactsheets.R.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggforce)
})

COL_DOWNSTREAM <- "#E05C5C"
COL_UPSTREAM   <- "#4B8BBE"
COL_EXON       <- "#7a5195"
EXON_STUB      <- 200

# ── normalise a raw pairs data.frame (Section 2 logic, no file IO) ────────────
normalise_pairs <- function(raw) {
  rename_col <- function(df, old, new) {
    if (old %in% names(df) && !new %in% names(df))
      names(df)[names(df) == old] <- new
    df
  }
  raw <- rename_col(raw, "upstream",   "upstream_count")
  raw <- rename_col(raw, "downstream", "downstream_count")
  raw <- rename_col(raw, "total",      "total_reads")

  if (!"fraction_downstream" %in% names(raw))
    raw$fraction_downstream <- raw$downstream_count / raw$total_reads
  if (!"intron1_length" %in% names(raw))
    raw$intron1_length <- raw$intron1_end - raw$intron1_start
  if (!"intron2_length" %in% names(raw))
    raw$intron2_length <- raw$intron2_end - raw$intron2_start
  raw
}

# ── build arcs (Section 5 logic, unchanged) ──────────────────────────────────
.build_arcs <- function(pairs_df, introns) {
  if (nrow(pairs_df) == 0) return(NULL)
  idx_map <- introns %>% select(start, end, number)

  pairs_indexed <- pairs_df %>%
    left_join(idx_map %>% rename(i1_num = number),
              by = c("intron1_start" = "start", "intron1_end" = "end")) %>%
    left_join(idx_map %>% rename(i2_num = number),
              by = c("intron2_start" = "start", "intron2_end" = "end")) %>%
    filter(!is.na(i1_num), !is.na(i2_num)) %>%
    mutate(span = abs(i2_num - i1_num)) %>%
    arrange(span, i1_num)

  if (nrow(pairs_indexed) == 0) return(NULL)

  arc_list <- vector("list", nrow(pairs_indexed))
  for (k in seq_len(nrow(pairs_indexed))) {
    row   <- pairs_indexed[k, ]
    i_mid <- introns$mid[introns$number == row$i1_num]
    j_mid <- introns$mid[introns$number == row$i2_num]
    if (length(i_mid) == 0 || length(j_mid) == 0) next

    above  <- (k %% 2 == 1)
    dir_y  <- if (above) 1 else -1
    height <- dir_y * (0.40 + 0.18 * (row$span - 1))

    x_from <- if (row$direction == "downstream_biased") i_mid else j_mid
    x_to   <- if (row$direction == "downstream_biased") j_mid else i_mid

    padj_val <- if ("binom_padj" %in% names(row) && !is.na(row$binom_padj))
      row$binom_padj else NA_real_

    arc_list[[k]] <- tibble(
      x = c(x_from, (x_from + x_to) / 2, x_to),
      y = c(0, height, 0),
      group = k, direction = row$direction, span = row$span, above = above,
      x_tip = x_to, padj = padj_val, fd = row$fraction_downstream,
      total_reads = row$total_reads, constitutive = row$span == 1
    )
  }
  bind_rows(arc_list)
}

# ── main callable: returns a ggplot for one gene, or NULL if absent ──────────
# raw        : already-normalised pairs data.frame (whole file, many genes)
# gene_name  : gene_symbol to plot
# title_suffix: e.g. "human" / "mouse" -> title "GENE (human)"
plot_gene_splicing <- function(raw, gene_name, title_suffix = NULL) {

  df_gene_sig <- raw %>%
    filter(gene_symbol == gene_name,
           direction %in% c("downstream_biased", "upstream_biased"))
  df_gene_all <- raw %>% filter(gene_symbol == gene_name)
  if (nrow(df_gene_all) == 0) return(NULL)

  # ── Section 4: reconstruct introns + exons (unchanged) ──────────────────────
  introns_tbl <- df_gene_all %>%
    select(intron1_start, intron1_end, intron2_start, intron2_end) %>%
    pivot_longer(everything(),
                 names_to = c("which", ".value"),
                 names_pattern = "(intron[12])_(start|end)") %>%
    distinct(start, end) %>% arrange(start) %>%
    mutate(number = row_number(), mid = (start + end) / 2, length = end - start)

  exons_internal <- if (nrow(introns_tbl) >= 2) {
    tibble(start = introns_tbl$end[-nrow(introns_tbl)],
           end   = introns_tbl$start[-1],
           number = seq_len(nrow(introns_tbl) - 1) + 1L) %>%
      filter(end > start)
  } else tibble(start = numeric(), end = numeric(), number = integer())

  exons_all <- bind_rows(
    tibble(start = min(introns_tbl$start) - EXON_STUB,
           end = min(introns_tbl$start), number = 1L),
    exons_internal,
    tibble(start = max(introns_tbl$end),
           end = max(introns_tbl$end) + EXON_STUB,
           number = nrow(introns_tbl) + 1L)
  ) %>% arrange(start) %>%
    mutate(mid = (start + end) / 2,
           constitutive = number > 1 & number < max(number))

  arcs <- .build_arcs(df_gene_sig, introns_tbl)

  # ── Section 7: plot (layer order preserved) ─────────────────────────────────
  x_range  <- range(c(exons_all$start, exons_all$end))
  x_pad    <- diff(x_range) * 0.03
  x_limits <- c(x_range[1] - x_pad, x_range[2] + x_pad)
  arc_extent <- if (!is.null(arcs) && nrow(arcs) > 0)
    max(abs(arcs$y)) + 0.12 else 0.55
  y_limits <- c(-arc_extent - 0.10, arc_extent + 0.10)

  p <- ggplot() +
    annotate("segment", x = min(exons_all$start), xend = max(exons_all$end),
             y = 0, yend = 0, colour = "black", linewidth = 0.4) +
    annotate("segment",
             x = max(exons_all$end) - diff(x_range) * 0.015,
             xend = max(exons_all$end), y = 0, yend = 0,
             colour = "grey40", linewidth = 0.9,
             arrow = arrow(length = unit(5, "pt"), type = "closed"))

  if (!is.null(arcs) && nrow(arcs) > 0) {
    max_reads <- max(arcs$total_reads, na.rm = TRUE)
    for (g in unique(arcs$group)) {
      ag <- filter(arcs, group == g)
      col_g <- if (ag$direction[1] == "downstream_biased")
        COL_DOWNSTREAM else COL_UPSTREAM
      lwd_g <- pmax(0.3, pmin(1.5,
        0.3 + 1.2 * sqrt(ag$total_reads[1]) / sqrt(max_reads)))
      alpha_g <- 0.4
      p <- p +
        geom_bezier(data = ag, aes(x = x, y = y, group = group),
                    colour = col_g, linewidth = lwd_g, alpha = alpha_g,
                    lineend = "round", show.legend = FALSE) +
        annotate("point", x = ag$x_tip[1], y = 0,
                 colour = col_g, size = 2.0, shape = 17, alpha = alpha_g)
    }
  }

  p <- p +
    geom_rect(data = exons_all,
              aes(xmin = start, xmax = end, ymin = -0.28, ymax = 0.28),
              fill = COL_EXON, colour = COL_EXON, linewidth = 0.3) +
    geom_text(data = exons_all, aes(x = mid, y = -0.5, label = number),
              size = 2.6, colour = COL_EXON, fontface = "bold", vjust = 1) +
    geom_text(data = introns_tbl, aes(x = mid, y = 0.5, label = number),
              size = 2.4, colour = "black", fontface = "bold", vjust = 0)

  leg_df <- tibble(x = c(NA_real_, NA_real_), y = c(NA_real_, NA_real_),
    bias = factor(c("Downstream-biased", "Upstream-biased")))

  title_txt <- if (is.null(title_suffix)) gene_name
               else paste0(gene_name, " (", title_suffix, ")")

  p +
    geom_line(data = leg_df, aes(x = x, y = y, colour = bias),
              linewidth = 1.3, na.rm = TRUE) +
    scale_colour_manual(
      values = c("Downstream-biased" = COL_DOWNSTREAM,
                 "Upstream-biased" = COL_UPSTREAM),
      name = NULL, drop = FALSE) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
    scale_x_continuous(labels = function(x) paste0(round(x / 1e3, 1), " kb"),
                       expand = c(0, 0)) +
    labs(title = title_txt, x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(
      axis.line.y = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(colour = "grey70", linewidth = 0.4),
      axis.text.x = element_text(colour = "grey40", size = 6),
      axis.ticks.x = element_line(colour = "grey70"),
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      legend.position = "none",
      plot.margin = margin(4, 8, 2, 8))
}