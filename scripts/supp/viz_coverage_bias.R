#!/usr/bin/env Rscript
# =============================================================================
# plot_coverage_bias.R
#
# Tests whether the downstream-first splicing bias correlates with intron pair
# position within the transcript — which would suggest a 3' coverage artifact.
#
# Input:  TSV produced by coverage_bias_check.py
#         (pooled splicing order TSV with normalized_position column added)
#
# Usage:
#   Rscript plot_coverage_bias.R \
#     --input  /users/dhan30/scratch/splicing_order/coverage_bias/pooled_with_position.tsv \
#     --outdir /users/dhan30/scratch/splicing_order/coverage_bias/figures
#
# Outputs:
#   01_position_vs_fraction_downstream_binned.pdf/.png
#   02_position_distribution_by_direction.pdf/.png
#   03_coverage_by_position.pdf/.png
#   coverage_bias_summary.txt
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# --- CLI args ----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_file <- NULL
out_dir    <- "coverage_bias_figures"

for (i in seq_along(args)) {
  if (args[i] == "--input"  && i < length(args)) input_file <- args[i + 1]
  if (args[i] == "--outdir" && i < length(args)) out_dir    <- args[i + 1]
}

if (is.null(input_file)) {
  # Fallback: edit this path directly
  input_file <- "/users/dhan30/scratch/splicing_order/coverage_bias/pooled_with_position.tsv"
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading:", input_file, "\n")
df <- read.delim(input_file, stringsAsFactors = FALSE, check.names = FALSE)
cat("  Rows:", nrow(df), "\n\n")

# --- Theme -------------------------------------------------------------------
theme_thesis <- function() {
  theme_classic(base_size = 12) +
    theme(
      axis.text       = element_text(size = 10),
      axis.title      = element_text(size = 11),
      plot.title      = element_text(size = 12, face = "bold"),
      plot.subtitle   = element_text(size = 10, color = "grey40"),
      legend.position = "right"
    )
}

save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300)
  cat("  Saved:", name, "\n")
}

# =============================================================================
# FIGURE 1: Binned mean fraction_downstream vs normalized position
# This is the key test — if flat, no coverage bias
# =============================================================================
cat("[1/3] Binned mean fraction_downstream vs position...\n")

n_bins <- 20

bin_summary <- df %>%
  mutate(pos_bin = cut(normalized_position,
                       breaks = seq(0, 1, length.out = n_bins + 1),
                       include.lowest = TRUE)) %>%
  group_by(pos_bin) %>%
  summarise(
    mean_fd  = mean(fraction_downstream, na.rm = TRUE),
    sem      = sd(fraction_downstream,   na.rm = TRUE) / sqrt(n()),
    n        = n(),
    pos_mid  = mean(normalized_position, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  filter(n >= 100)   # drop sparse edge bins

# Pearson correlation
ct <- cor.test(df$normalized_position, df$fraction_downstream, method = "pearson")
r_val <- round(ct$estimate, 4)
p_val <- signif(ct$p.value, 3)
subtitle_text <- paste0("Pearson r = ", r_val, "   p = ", p_val,
                        "   n = ", scales::comma(nrow(df)), " pairs")

p1 <- ggplot(bin_summary, aes(x = pos_mid, y = mean_fd)) +
  geom_hline(yintercept = 0.5,   linetype = "dashed", color = "grey60", linewidth = 0.8) +
  geom_hline(yintercept = mean(df$fraction_downstream),
             linetype = "dotted", color = "#E05C5C", linewidth = 0.8) +
  geom_errorbar(aes(ymin = mean_fd - sem, ymax = mean_fd + sem),
                width = 0.02, color = "grey50", linewidth = 0.6) +
  geom_line(color = "#4B8BBE", linewidth = 1) +
  geom_point(shape = 21, size = 2.8, fill = "white",
             color = "#4B8BBE", stroke = 1.5) +
  annotate("text", x = 0.75, y = mean(df$fraction_downstream) + 0.001,
           label = paste0("Global mean = ", round(mean(df$fraction_downstream), 4)),
           color = "#E05C5C", size = 3.2, hjust = 0) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     labels = c("5' end\n(0)", "0.2", "0.4", "0.6", "0.8",
                                "3' end\n(1.0)")) +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
    title    = "Coverage Bias Check: Splicing Order vs. Transcript Position",
    subtitle = subtitle_text,
    x        = "Normalized position within transcript (5' → 3')",
    y        = "Mean fraction downstream ± SEM"
  ) +
  theme_thesis()

save_fig(p1, "01_position_vs_fraction_downstream_binned", w = 8, h = 5)

# =============================================================================
# FIGURE 2: Distribution of normalized positions by splicing direction
# If 3' bias: downstream-biased pairs should cluster toward position 1.0
# =============================================================================
cat("[2/3] Position distribution by splicing direction...\n")

df_dir <- df %>%
  mutate(direction = case_when(
    fraction_downstream >= 0.75 ~ "Downstream-biased\n(fd ≥ 0.75)",
    fraction_downstream <= 0.25 ~ "Upstream-biased\n(fd ≤ 0.25)",
    TRUE ~ "No strong bias\n(0.25 < fd < 0.75)"
  ))

dir_colors <- c(
  "Downstream-biased\n(fd ≥ 0.75)" = "#E05C5C",
  "No strong bias\n(0.25 < fd < 0.75)" = "#AAAAAA",
  "Upstream-biased\n(fd ≤ 0.25)"   = "#4B8BBE"
)

p2 <- ggplot(df_dir, aes(x = normalized_position, fill = direction,
                          color = direction)) +
  geom_density(alpha = 0.35, linewidth = 0.8) +
  scale_fill_manual(values  = dir_colors) +
  scale_color_manual(values = dir_colors) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     labels = c("5'", "0.2", "0.4", "0.6", "0.8", "3'")) +
  labs(
    title    = "Intron Pair Position Within Transcript by Splicing Direction",
    subtitle = "If 3' coverage bias drives downstream signal, red should shift right",
    x        = "Normalized position within transcript (5' → 3')",
    y        = "Density",
    fill     = NULL, color = NULL
  ) +
  theme_thesis()

save_fig(p2, "02_position_distribution_by_direction", w = 8, h = 5)

# =============================================================================
# FIGURE 3: Read coverage (total informative reads) vs transcript position
# Tests whether coverage itself is 3'-biased
# =============================================================================
cat("[3/3] Read coverage vs transcript position...\n")

cov_bins <- df %>%
  mutate(pos_bin = cut(normalized_position,
                       breaks = seq(0, 1, length.out = n_bins + 1),
                       include.lowest = TRUE)) %>%
  group_by(pos_bin) %>%
  summarise(
    median_cov = median(total, na.rm = TRUE),
    mean_cov   = mean(total,   na.rm = TRUE),
    n          = n(),
    pos_mid    = mean(normalized_position, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  filter(n >= 100)

p3 <- ggplot(cov_bins, aes(x = pos_mid, y = median_cov)) +
  geom_line(color = "#4B8BBE", linewidth = 1) +
  geom_point(shape = 21, size = 2.8, fill = "white",
             color = "#4B8BBE", stroke = 1.5) +
  scale_x_continuous(breaks = seq(0, 1, 0.2),
                     labels = c("5'", "0.2", "0.4", "0.6", "0.8", "3'")) +
  labs(
    title    = "Read Coverage vs. Transcript Position",
    subtitle = "Flat line = uniform coverage; rising toward 3' = coverage bias",
    x        = "Normalized position within transcript (5' → 3')",
    y        = "Median informative read pairs per intron pair"
  ) +
  theme_thesis()

save_fig(p3, "03_coverage_by_position", w = 8, h = 5)

# =============================================================================
# Summary text file
# =============================================================================
sink(file.path(out_dir, "coverage_bias_summary.txt"))
cat("Coverage Bias Check Summary\n")
cat("===========================\n\n")
cat("Input file:", input_file, "\n")
cat("Total pairs analyzed:", scales::comma(nrow(df)), "\n\n")
cat("Pearson correlation (normalized_position vs fraction_downstream):\n")
cat("  r =", round(ct$estimate, 5), "\n")
cat("  p =", signif(ct$p.value, 4), "\n")
cat("  95% CI: [", round(ct$conf.int[1], 5), ",",
    round(ct$conf.int[2], 5), "]\n\n")
cat("Global mean fraction_downstream:", round(mean(df$fraction_downstream), 5), "\n\n")
cat("Interpretation:\n")
if (abs(ct$estimate) < 0.02) {
  cat("  No positional trend detected. Coverage bias is unlikely to explain\n")
  cat("  the observed global downstream-first bias.\n")
} else if (abs(ct$estimate) < 0.05) {
  cat("  Weak positional trend. Unlikely to fully explain the downstream bias\n")
  cat("  but worth noting in a limitations section.\n")
} else {
  cat("  Non-trivial positional trend detected. Coverage bias may be a\n")
  cat("  contributing factor. Inspect figures carefully.\n")
}
sink()

cat("\nDone. Figures and summary written to:", out_dir, "\n")

