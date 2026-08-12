#!/usr/bin/env Rscript
# viz_lariat_overlap.R
# Visualize the two outputs of lariat_order_overlap.py.
#
#   1. gene-ordinal : lariat rate vs order_rank (with Wilson CIs) + length control
#   2. pairwise     : lariat rate for first- vs second-spliced introns
#
# Usage:
  # Rscript scripts/figures/viz_lariat_overlap.R \
  #   results/lariat_overlap/gene_ordinal_intron_lariat.tsv \
  #   results/lariat_overlap/pairwise_intron_lariat.tsv \
  #   results/figures/lariat_overlap

suppressMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(patchwork)
})
 
args <- commandArgs(trailingOnly = TRUE)
gene_path <- args[1]
pair_path <- args[2]
outdir    <- ifelse(length(args) >= 3, args[3], ".")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
 
# palette from your thesis convention
col_down <- "#E05C5C"   # downstream / later
col_up   <- "#4B8BBE"   # upstream / earlier
 
theme_thesis <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )
 
# Wilson score interval for a proportion (better than normal approx at low rates)
wilson <- function(k, n, z = 1.96) {
  p <- k / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  data.frame(lo = pmax(0, centre - half), hi = pmin(1, centre + half))
}
 
# ------------------------------------------------------------------
# 1. GENE-ORDINAL: lariat rate by order_rank
# ------------------------------------------------------------------
g <- read.delim(gene_path)
g$lariat_present <- as.logical(g$lariat_present)  # pandas writes "True"/"False"
g <- subset(g, n_introns_in_gene > 1)
 
# collapse the noisy long tail: cap rank for display
RANK_CAP <- 8
g$rank_disp <- pmin(g$order_rank, RANK_CAP)
g$rank_lab  <- ifelse(g$order_rank >= RANK_CAP, paste0(RANK_CAP, "+"),
                      as.character(g$order_rank))
 
rate_by_rank <- g %>%
  group_by(rank_disp) %>%
  summarise(k = sum(lariat_present), n = n(), .groups = "drop") %>%
  mutate(rate = k / n) %>%
  bind_cols(wilson(.$k, .$n)) %>%
  mutate(rank_lab = ifelse(rank_disp >= RANK_CAP, paste0(RANK_CAP, "+"),
                           as.character(rank_disp)))
 
p1 <- ggplot(rate_by_rank, aes(x = factor(rank_disp), y = rate)) +
  geom_col(fill = col_down, alpha = 0.85) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.25, color = "grey30") +
  scale_x_discrete(labels = rate_by_rank$rank_lab) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Lariat rate rises with later splicing order",
       subtitle = "Rank 1 = earliest-spliced intron (within tested pairs); Wilson 95% CI",
       x = "Order rank within gene", y = "Fraction with lariat support") +
  theme_thesis
 
# length control: is the rank effect just intron length?
# bin length, facet the rank trend within length bins
g$len_bin <- cut(g$intron_length,
                 breaks = quantile(g$intron_length, c(0, .25, .5, .75, 1), na.rm = TRUE),
                 labels = c("Q1 (short)", "Q2", "Q3", "Q4 (long)"),
                 include.lowest = TRUE)
 
rate_rank_len <- g %>%
  filter(!is.na(len_bin)) %>%
  group_by(len_bin, rank_disp) %>%
  summarise(k = sum(lariat_present), n = n(), .groups = "drop") %>%
  filter(n >= 30) %>%
  mutate(rate = k / n)
 
p2 <- ggplot(rate_rank_len, aes(x = rank_disp, y = rate, color = len_bin)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_x_continuous(breaks = 1:RANK_CAP) +
  scale_color_viridis_d(option = "D", end = 0.9, name = "Intron length quartile") +
  labs(title = "Order effect within intron-length quartiles",
       subtitle = "If lines still slope up, position matters beyond length",
       x = "Order rank within gene", y = "Fraction with lariat support") +
  theme_thesis
 
# ------------------------------------------------------------------
# 2. PAIRWISE: first- vs second-spliced
# ------------------------------------------------------------------
p <- read.delim(pair_path)
p$lariat_present <- as.logical(p$lariat_present)  # pandas writes "True"/"False"
pair_rate <- p %>%
  group_by(first_spliced) %>%
  summarise(k = sum(lariat_present), n = n(), .groups = "drop") %>%
  mutate(rate = k / n,
         role = ifelse(first_spliced, "First-spliced", "Second-spliced")) %>%
  bind_cols(wilson(.$k, .$n))
 
p3 <- ggplot(pair_rate, aes(x = role, y = rate, fill = role)) +
  geom_col(alpha = 0.85, width = 0.6) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, color = "grey30") +
  scale_fill_manual(values = c("First-spliced" = col_up,
                               "Second-spliced" = col_down), guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1),
                     expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Lariats favor the later-spliced intron",
       subtitle = "Per intron-in-pair membership (rows not independent)",
       x = NULL, y = "Fraction with lariat support") +
  theme_thesis
 
# ------------------------------------------------------------------
# save
# ------------------------------------------------------------------
ggsave(file.path(outdir, "gene_ordinal_rate_by_rank.png"), p1,
       width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "gene_ordinal_rank_by_length.png"), p2,
       width = 7, height = 5, dpi = 300)
ggsave(file.path(outdir, "pairwise_first_vs_second.png"), p3,
       width = 5, height = 5, dpi = 300)
 
combined <- (p1 | p3) / p2 + plot_annotation(tag_levels = "A")
ggsave(file.path(outdir, "lariat_overlap_combined.png"), combined,
       width = 12, height = 9, dpi = 300)
 
cat("wrote figures to", outdir, "\n")