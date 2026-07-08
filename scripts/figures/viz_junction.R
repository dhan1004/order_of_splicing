#!/usr/bin/env Rscript
# =============================================================================
# figures_style.R  —  Fairbrother Lab
# Produces two  figure types:
#  Binned mean ± SEM line plot (Figure 1-style)
#  Overlaid histogram by splicing order group (Figure 2-style)
#
# Usage: source("viz_junction.R")  or  Rscript viz_junction.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(scales)
})

# ── Config ────────────────────────────────────────────────────────────────────
INPUT_TSV  <- "merged/sig_structure_features_final.tsv"   # ← adjust path
OUT_DIR    <- "figures/0408b"
MIN_READS  <- 10
MIN_INF    <- 5
N_BINS     <- 20   # number of equal-width bins for the line plot

dir.create(OUT_DIR, showWarnings = FALSE)

save_fig <- function(p, name, w = 7, h = 6) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h,
         dpi = 200, bg = "white")
  message("  saved: ", name)
}

# ── Theme matching original figures ──────────────────────────────────────────
theme_thesis <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
  theme(
    panel.grid.major   = element_line(color = "grey88", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "grey70", fill = NA),
    axis.text          = element_text(color = "black", size = base_size - 1),
    axis.title         = element_text(color = "black", size = base_size,
                                      face = "bold"),
    plot.title         = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0),
    plot.subtitle      = element_text(size = base_size - 2, color = "grey40"),
    legend.background  = element_rect(fill = "white", color = "grey80",
                                      linewidth = 0.3),
    legend.key.size    = unit(0.45, "cm"),
    legend.title       = element_text(size = base_size - 1, face = "bold"),
    legend.text        = element_text(size = base_size - 2)
  )
}

# =============================================================================
# FIGURE A: Binned mean ± SEM line plot
#   x = diff column
#   y = mean fraction_downstream ± SEM
#   Matches the style in Image 1
# =============================================================================
plot_binned_mean_sem <- function(df,
                                  diff_col,
                                  xlab        = "Boundary H-bond Difference (binned)\n(Intron1 \u2212 Intron2)",
                                  title       = NULL,
                                  n_bins      = N_BINS,
                                  negate_x    = FALSE,   # flip to intron1-intron2
                                  point_color = "#5b9bd5",
                                  line_color  = "#5b9bd5",
                                  err_color   = "grey45") {

  if (!diff_col %in% names(df)) {
    message("Column not found: ", diff_col); return(invisible(NULL))
  }

  d <- df %>%
    select(x_raw = all_of(diff_col), y = fraction_downstream) %>%
    drop_na() %>%
    mutate(x = if (negate_x) -x_raw else x_raw)  # flip sign if needed

  # Equal-width bins on the (possibly negated) x
  breaks  <- seq(min(d$x), max(d$x), length.out = n_bins + 1)
  d <- d %>%
    mutate(bin = cut(x, breaks = breaks, include.lowest = TRUE))

  bin_sum <- d %>%
    group_by(bin) %>%
    summarise(
      mean_y  = mean(y),
      sem     = sd(y) / sqrt(n()),
      n       = n(),
      x_mid   = mean(x),   # midpoint of bin
      .groups = "drop"
    ) %>%
    filter(n >= 5)  # drop near-empty tail bins

  ggplot(bin_sum, aes(x_mid, mean_y)) +
    geom_errorbar(aes(ymin = mean_y - sem, ymax = mean_y + sem),
                  width = diff(range(bin_sum$x_mid)) / (n_bins * 2),
                  color = err_color, linewidth = 0.7) +
    geom_line(color = line_color, linewidth = 1.1) +
    geom_point(shape = 21, size = 3.5, fill = "white",
               color = point_color, stroke = 1.8) +
    scale_y_continuous(
      breaks = pretty_breaks(5),
      labels = label_number(accuracy = 0.01)
    ) +
    scale_x_continuous(breaks = pretty_breaks(7)) +
    labs(
      x     = xlab,
      y     = "Mean Fraction Downstream \u00b1 SEM",
      title = title
    ) +
    theme_thesis()
}

# =============================================================================
# FIGURE B: Overlaid histogram by splicing order group
#   Upstream spliced first  = fraction_downstream < 0.5
#   Downstream spliced first = fraction_downstream > 0.5
#   Matches the style in Image 2
# =============================================================================
plot_order_histogram <- function(df,
                                  diff_col,
                                  xlab     = "Boundary H-bond Difference",
                                  title    = "Distribution by Splicing Order",
                                  n_breaks = 60,
                                  negate_x = FALSE,
                                  col_up   = "#5b7fcb",   # blue
                                  col_dn   = "#f5c842") { # yellow-gold

  if (!diff_col %in% names(df)) {
    message("Column not found: ", diff_col); return(invisible(NULL))
  }

  d <- df %>%
    select(x_raw = all_of(diff_col), fd = fraction_downstream) %>%
    drop_na() %>%
    filter(fd != 0.5) %>% 
    mutate(
      x     = if (negate_x) -x_raw else x_raw,
      group = if_else(fd < 0.5, "Upstream spliced first",
                                "Downstream spliced first"),
      group = factor(group, levels = c("Upstream spliced first",
                                       "Downstream spliced first"))
    )

  # Shared bin breaks
  bin_breaks <- seq(
    floor(quantile(d$x, 0.001) / 10) * 10,
    ceiling(quantile(d$x, 0.999) / 10) * 10,
    length.out = n_breaks + 1
  )

  ggplot(d, aes(x, fill = group, color = group)) +
    geom_histogram(breaks = bin_breaks,
                   alpha = 0.55, position = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "#cc2222", linewidth = 0.9) +
    scale_fill_manual(
      values = c("Upstream spliced first"   = col_up,
                 "Downstream spliced first" = col_dn),
      name = NULL
    ) +
    scale_color_manual(
      values = c("Upstream spliced first"   = col_up,
                 "Downstream spliced first" = col_dn),
      name = NULL
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
    labs(
      x     = xlab,
      y     = "Count",
      title = title
    ) +
    theme_thesis() +
    theme(
      legend.position   = c(0.02, 0.98),
      legend.justification = c(0, 1),
      panel.grid.major.x = element_blank()
    )
}

# =============================================================================
# LOAD & FILTER
# =============================================================================
message("Loading: ", INPUT_TSV)
raw <- read_tsv(INPUT_TSV, show_col_types = FALSE)

df <- raw %>%
  filter(!is.na(fraction_downstream),
         total_reads >= MIN_READS,
         (upstream_count + downstream_count) >= MIN_INF)

message(sprintf("  Filtered to %s pairs", format(nrow(df), big.mark = ",")))

# =============================================================================
# GENERATE FIGURES
# diff_junction_hbonds_per_nt is the primary target (normalized);
# also generate for raw boundary diff and per individual SS.
# =============================================================================

targets <- list(
  list(
    col    = "diff_junction_hbonds_per_nt",
    xlab_a = "Boundary H-bond Difference per nt (binned)\n(Downstream \u2212 Upstream Intron)",
    xlab_b = "Boundary H-bond Difference per nt\n(Downstream \u2212 Upstream Intron)",
    title  = "Boundary H-bond Difference",
    stem   = "junction_norm"
  ),
  list(
    col    = "diff_junction_hbonds_raw",
    xlab_a = "Boundary H-bond Difference (binned)\n(Downstream \u2212 Upstream Intron)",
    xlab_b = "Boundary H-bond Difference\n(Downstream \u2212 Upstream Intron)",
    title  = "Boundary H-bond Difference (raw)",
    stem   = "junction_raw"
  ),
  list(
    col    = "diff_5ss_hbonds_per_nt",
    xlab_a = "5\u2019SS H-bond Difference per nt (binned)\n(Downstream \u2212 Upstream Intron)",
    xlab_b = "5\u2019SS H-bond Difference per nt\n(Downstream \u2212 Upstream Intron)",
    title  = "5\u2019 Splice Site H-bond Difference",
    stem   = "5ss_norm"
  ),
  list(
    col    = "diff_3ss_hbonds_per_nt",
    xlab_a = "3\u2019SS H-bond Difference per nt (binned)\n(Downstream \u2212 Upstream Intron)",
    xlab_b = "3\u2019SS H-bond Difference per nt\n(Downstream \u2212 Upstream Intron)",
    title  = "3\u2019 Splice Site H-bond Difference",
    stem   = "3ss_norm"
  ),
  list(
    col    = "diff_hbonds_per_nt",
    xlab_a = "Internal H-bond Difference per nt (binned)\n(Intron2 \u2212 Intron1)",
    xlab_b = "Internal H-bond Difference per nt\n(Intron2 \u2212 Intron1)",
    title  = "Internal H-bond Difference",
    stem   = "internal_norm"
  )
)

for (t in targets) {
  if (!t$col %in% names(df)) {
    message("Skipping (column absent): ", t$col); next
  }

  # Figure A: binned mean ± SEM
  pa <- plot_binned_mean_sem(df,
    diff_col = t$col,
    xlab     = t$xlab_a,
    title    = paste("Binned Mean ±SEM —", t$title)
  )
  if (!is.null(pa)) save_fig(pa, paste0("A_", t$stem, "_binned_sem"))

  # Figure B: histogram by splicing order
  pb <- plot_order_histogram(df,
    diff_col = t$col,
    xlab     = t$xlab_b,
    title    = paste("Distribution by Splicing Order —", t$title)
  )
  if (!is.null(pb)) save_fig(pb, paste0("B_", t$stem, "_histogram"))
}

message("\nDone. Figures in: ", OUT_DIR, "/")