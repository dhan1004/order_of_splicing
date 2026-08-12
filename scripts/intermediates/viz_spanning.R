#!/usr/bin/env Rscript
# viz_spanning.R
# Visualize the pooled frozen-intermediate signal from pool_spanning.py.
# Histogram facets for the bounded [0,1] frozen_score (per your preference),
# plus a length/depth diagnostic panel to eyeball the confounds up front.
#
# Usage:
#   Rscript viz_spanning.R spanning_analysis_set.tsv spanning_plots.png

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
infile  <- ifelse(length(args) >= 1, args[1], "spanning_analysis_set.tsv")
outfile <- ifelse(length(args) >= 2, args[2], "spanning_plots.png")

theme_thesis <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        plot.title = element_text(face = "bold"))

d <- read.delim(infile, stringsAsFactors = FALSE)
d$frozen_score <- as.numeric(d$frozen_score)
d <- d[!is.na(d$frozen_score), ]

# ---- P1: overall frozen_score distribution -------------------------------
p1 <- ggplot(d, aes(frozen_score)) +
  geom_histogram(binwidth = 0.05, boundary = 0,
                 fill = "steelblue", colour = "white") +
  geom_vline(xintercept = 0.5, linetype = 2, colour = "grey40") +
  labs(title = "Frozen signature (pooled)",
       subtitle = "3'ss span / (3'ss + 5'ss span) per retained intron",
       x = "3'ss span / (3'ss + 5'ss span)",
       y = "introns") +
  theme_thesis

# ---- P2: split by direction (upstream/downstream retained) ---------------
if ("direction" %in% names(d)) {
  p2 <- ggplot(d, aes(frozen_score)) +
    geom_histogram(binwidth = 0.05, boundary = 0,
                   fill = "darkorange", colour = "white") +
    facet_wrap(~direction, ncol = 1, scales = "free_y") +
    labs(title = "By pair direction", x = "3'ss span / (3'ss + 5'ss span)", y = "introns") +
    theme_thesis
} else {
  p2 <- ggplot() + theme_void()
}

# ---- P3: confound check — length vs frozen_score -------------------------
d$len_kb <- d$retained_intron_length / 1000
p3 <- ggplot(d, aes(len_kb + 0.01, frozen_score)) +
  geom_bin2d(bins = 40) +
  scale_x_log10(labels = label_number()) +
  scale_fill_viridis_c(trans = "log10") +
  labs(title = "Confound: intron length",
       x = "retained intron length (kb, log)", y = "3'ss span / (3'ss + 5'ss span)") +
  theme_thesis

# ---- P4: confound check — coverage depth vs frozen_score -----------------
p4 <- ggplot(d, aes(total_span, frozen_score)) +
  geom_bin2d(bins = 40) +
  scale_x_log10() +
  scale_fill_viridis_c(trans = "log10") +
  labs(title = "Confound: coverage depth",
       x = "total span reads (log)", y = "3'ss span / (3'ss + 5'ss span)") +
  theme_thesis

combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Pairwise splice-site spanning — frozen intermediate signal",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

ggsave(outfile, combined, width = 12, height = 9, dpi = 200)
cat("Wrote", outfile, "\n")