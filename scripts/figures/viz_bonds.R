#!/usr/bin/env Rscript
# =============================================================================
# structure_junction_analysis.R
#
# Analyzes RNA secondary structure data with focus on:
#   - Intron-exon junction hydrogen bonds (per-nt normalized)
#   - Internal intron pairing interactions (pairing_fraction, hbonds_per_nt)
#   - 5'SS vs 3'SS breakdown
#   - Upstream vs downstream comparisons for significant pairs
#   - Correlations with fraction_downstream
#
# Input: structure_features_final.tsv (from combine_structure_chunks.py)
# Output: figures as PDF + PNG in ./figures/
#
# Column conventions (from analyze_transcript.py):
#   intron1_* = upstream intron    intron2_* = downstream intron
#   diff_*    = downstream - upstream (intron2 - intron1)
#   *_per_nt  = raw count / intron length
#   diff_*_accessibility: positive → downstream more accessible → +r expected
#   diff_*_hbonds_per_nt: positive → downstream more occluded  → -r expected

# how to run: Rscript scripts/figures/viz_bonds.R
# =============================================================================


suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(stringr); library(purrr); library(forcats)
  library(scales); library(patchwork); library(ggridges)
})

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_TSV   <- "merged/sig_structure_features_final.tsv"   # change path as needed
OUT_DIR     <- "figures/0408"
MIN_READS   <- 10     # minimum total_reads for a pair to be "significant", should be redundant
MIN_INF     <- 5      # minimum upstream + downstream informative count
N_QUANTILES <- 5      # for binned accessibility/hbond plots

# Palette: upstream = teal, downstream = coral, diff = purple
COL_UP   <- "#2D9CDB"
COL_DN   <- "#F2994A"
COL_DIFF <- "#9B51E0"
COL_5SS  <- "#27AE60"
COL_3SS  <- "#EB5757"

dir.create(OUT_DIR, showWarnings = FALSE)

# Helper: save figure
save_fig <- function(p, name, w = 8, h = 6) {
  base <- file.path(OUT_DIR, name)
  ggsave(paste0(base, ".pdf"), p, width = w, height = h)
  ggsave(paste0(base, ".png"), p, width = w, height = h, dpi = 200)
  message("  saved: ", name)
}

# Helper: annotated scatter with Pearson r
cor_scatter <- function(df, x_col, y_col = "fraction_downstream",
                        xlab = x_col, title = NULL,
                        color = COL_DIFF, alpha = 0.08) {
  d <- df %>% select(x = all_of(x_col), y = all_of(y_col)) %>% drop_na()
  ct <- cor.test(d$x, d$y, method = "pearson")
  r  <- ct$estimate; p <- ct$p.value; n <- nrow(d)
  lab <- sprintf("r = %.3f\np = %s\nn = %s", r,
                 ifelse(p < 1e-300, "< 1e-300", formatC(p, format = "e", digits = 2)),
                 format(n, big.mark = ","))

  # Downsample for plotting speed (keep all for stats)
  d_plot <- if (n > 20000) slice_sample(d, n = 20000) else d

  ggplot(d_plot, aes(x, y)) +
    geom_hex(bins = 60, show.legend = TRUE) +
    scale_fill_gradient(low = "grey92", high = color, trans = "log10",
                        name = "count") +
    geom_smooth(data = d, aes(x, y), method = "lm", se = TRUE,
                color = "black", linewidth = 0.8, linetype = "dashed") +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
             label = lab, size = 3.2, family = "mono") +
    labs(x = xlab, y = "Fraction downstream spliced first",
         title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(size = 11, face = "bold"))
}

# Helper: paired violin / boxplot upstream vs downstream
paired_violin <- function(df, col1, col2, lab1 = "Upstream (intron1)",
                           lab2 = "Downstream (intron2)", title = NULL,
                           ylab = "Value") {
  d <- bind_rows(
    tibble(value = df[[col1]], intron = lab1),
    tibble(value = df[[col2]], intron = lab2)
  ) %>% drop_na()

  wt <- wilcox.test(df[[col1]], df[[col2]], paired = FALSE, exact = FALSE)
  med_diff <- median(df[[col2]], na.rm = TRUE) - median(df[[col1]], na.rm = TRUE)
  lab <- sprintf("Wilcoxon p = %s\nΔmedian = %.4f",
                 formatC(wt$p.value, format = "e", digits = 2), med_diff)

  ggplot(d, aes(intron, value, fill = intron)) +
    geom_violin(alpha = 0.7, trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
                 linewidth = 0.6) +
    scale_fill_manual(values = setNames(c(COL_UP, COL_DN), c(lab1, lab2)),
                      guide = "none") +
    annotate("text", x = 1.5, y = Inf, vjust = 1.4, size = 3,
             label = lab, family = "mono") +
    labs(x = NULL, y = ylab, title = title) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(size = 11, face = "bold"))
}

# =============================================================================
# LOAD DATA
# =============================================================================
message("Loading data: ", INPUT_TSV)
raw <- read_tsv(INPUT_TSV, show_col_types = FALSE)
message(sprintf("  Raw rows: %s  Columns: %s", format(nrow(raw), big.mark=","), ncol(raw)))

# Filter to significant (well-covered) pairs
df <- raw %>%
  filter(!is.na(fraction_downstream),
         total_reads >= MIN_READS,
         (upstream_count + downstream_count) >= MIN_INF)

message(sprintf("  After filtering (total_reads >= %d, inf >= %d): %s pairs",
                MIN_READS, MIN_INF, format(nrow(df), big.mark=",")))

# Confirm key columns present
junction_cols <- c(
  "intron1_junction_hbonds_per_nt", "intron2_junction_hbonds_per_nt",
  "intron1_5ss_hbonds_per_nt",      "intron2_5ss_hbonds_per_nt",
  "intron1_3ss_hbonds_per_nt",      "intron2_3ss_hbonds_per_nt",
  "diff_junction_hbonds_per_nt",    "diff_5ss_hbonds_per_nt",
  "diff_3ss_hbonds_per_nt",
  "intron1_pairing_fraction",       "intron2_pairing_fraction",
  "intron1_hbonds_per_nt",          "intron2_hbonds_per_nt",
  "diff_hbonds_per_nt"
)
missing <- setdiff(junction_cols, names(df))
if (length(missing) > 0) {
  warning("Missing columns (skipping related plots): ",
          paste(missing, collapse = ", "))
}

# =============================================================================
# SECTION 1: CORRELATION SUMMARY TABLE
# =============================================================================
message("\n-- Section 1: Correlation table --")

diff_cols    <- grep("^diff_",          names(df), value = TRUE)
per_nt_cols  <- grep("_per_nt$",        names(df), value = TRUE)
acc_cols     <- grep("ss[35]_accessibility", names(df), value = TRUE)
all_predictor_cols <- unique(c(diff_cols, per_nt_cols, acc_cols))

cor_rows <- map_dfr(all_predictor_cols, function(col) {
  d <- df %>% select(x = all_of(col), y = fraction_downstream) %>% drop_na()
  if (nrow(d) < 30) return(NULL)
  ct_p <- cor.test(d$x, d$y, method = "pearson")
  ct_s <- cor.test(d$x, d$y, method = "spearman", exact = FALSE)
  tibble(
    column     = col,
    r_pearson  = round(ct_p$estimate, 4),
    r_spearman = round(ct_s$estimate, 4),
    p_pearson  = ct_p$p.value,
    n          = nrow(d)
  )
}) %>% arrange(desc(abs(r_pearson)))

print(cor_rows, n = 40)
write_tsv(cor_rows, file.path(OUT_DIR, "correlation_table.tsv"))

# Forest plot
cor_plot_df <- cor_rows %>%
  filter(str_detect(column,
    "junction|5ss|3ss|hbonds_per_nt|pairing_fraction|accessibility")) %>%
  mutate(
    category = case_when(
      str_detect(column, "junction")              ~ "Junction (5+3SS)",
      str_detect(column, "5ss")                   ~ "5' Splice Site",
      str_detect(column, "3ss")                   ~ "3' Splice Site",
      str_detect(column, "accessibility")         ~ "Accessibility",
      str_detect(column, "hbonds_per_nt|pairing") ~ "Internal H-bonds",
      TRUE ~ "Other"
    ),
    col_short = str_remove_all(column,
      "_hbonds_per_nt|_per_nt|_accessibility_|diff_|intron[12]_"),
    col_label = ifelse(str_starts(column, "diff_"),
                       paste0("Δ ", col_short), col_short),
    col_label = paste0(col_label, " (n=", format(n, big.mark=","), ")")
  )

p_forest <- ggplot(cor_plot_df,
    aes(r_pearson, reorder(col_label, r_pearson), color = category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 3) +
  geom_errorbarh(aes(
    xmin = r_pearson - 1.96 / sqrt(n - 3),
    xmax = r_pearson + 1.96 / sqrt(n - 3)
  ), height = 0.3) +
  scale_color_brewer(palette = "Dark2", name = "Feature type") +
  labs(x = "Pearson r (with fraction_downstream)",
       y = NULL,
       title = "Structural Predictors of Splicing Order",
       subtitle = "Points = r; bars = 95% CI. Δ = downstream - upstream") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold"))
save_fig(p_forest, "01_correlation_forest", w = 9, h = 7)

# =============================================================================
# SECTION 2: JUNCTION H-BONDS SCATTER PLOTS
# =============================================================================
message("\n-- Section 2: Junction bond scatter plots --")

if ("diff_junction_hbonds_per_nt" %in% names(df)) {
  p_junc <- cor_scatter(df, "diff_junction_hbonds_per_nt",
    xlab = "Δ Junction H-bonds per nt (downstream - upstream)",
    title = "Junction H-bond Difference vs Splicing Order",
    color = COL_DIFF)
  save_fig(p_junc, "02a_junction_hbonds_diff_scatter")
}

if ("diff_5ss_hbonds_per_nt" %in% names(df)) {
  p5 <- cor_scatter(df, "diff_5ss_hbonds_per_nt",
    xlab = "Δ 5'SS H-bonds per nt",
    title = "5' Splice Site H-bond Difference vs Splicing Order",
    color = COL_5SS)
  save_fig(p5, "02b_5ss_hbonds_diff_scatter")
}

if ("diff_3ss_hbonds_per_nt" %in% names(df)) {
  p3 <- cor_scatter(df, "diff_3ss_hbonds_per_nt",
    xlab = "Δ 3'SS H-bonds per nt",
    title = "3' Splice Site H-bond Difference vs Splicing Order",
    color = COL_3SS)
  save_fig(p3, "02c_3ss_hbonds_diff_scatter")
}

if (all(c("diff_5ss_hbonds_per_nt", "diff_3ss_hbonds_per_nt") %in% names(df))) {
  p_comb <- p5 + p3 + plot_annotation(
    title = "Splice Site-Specific H-bond Differences",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
  save_fig(p_comb, "02d_5ss_vs_3ss_scatter", w = 14, h = 6)
}

# =============================================================================
# SECTION 3: UPSTREAM vs DOWNSTREAM VIOLINS
# =============================================================================
message("\n-- Section 3: Upstream vs downstream violin plots --")

if (all(c("intron1_junction_hbonds_per_nt", "intron2_junction_hbonds_per_nt") %in% names(df))) {
  p_viol_junc <- paired_violin(df,
    "intron1_junction_hbonds_per_nt", "intron2_junction_hbonds_per_nt",
    title = "Junction H-bonds per nt: Upstream vs Downstream",
    ylab = "Junction H-bonds per nt")
  save_fig(p_viol_junc, "03a_junction_violin")
}

if (all(c("intron1_5ss_hbonds_per_nt", "intron2_5ss_hbonds_per_nt") %in% names(df))) {
  p_viol_5ss <- paired_violin(df,
    "intron1_5ss_hbonds_per_nt", "intron2_5ss_hbonds_per_nt",
    title = "5'SS H-bonds per nt: Upstream vs Downstream",
    ylab = "5'SS H-bonds per nt")
  save_fig(p_viol_5ss, "03b_5ss_violin")
}

if (all(c("intron1_3ss_hbonds_per_nt", "intron2_3ss_hbonds_per_nt") %in% names(df))) {
  p_viol_3ss <- paired_violin(df,
    "intron1_3ss_hbonds_per_nt", "intron2_3ss_hbonds_per_nt",
    title = "3'SS H-bonds per nt: Upstream vs Downstream",
    ylab = "3'SS H-bonds per nt")
  save_fig(p_viol_3ss, "03c_3ss_violin")
}

# =============================================================================
# SECTION 4: PAIRING FRACTION
# =============================================================================
message("\n-- Section 4: Pairing fraction --")

if (all(c("intron1_pairing_fraction", "intron2_pairing_fraction") %in% names(df))) {
  p_pf <- paired_violin(df,
    "intron1_pairing_fraction", "intron2_pairing_fraction",
    title = "Fraction of Nucleotides Base-Paired: Upstream vs Downstream",
    ylab = "Pairing fraction")
  save_fig(p_pf, "04a_pairing_fraction_violin")

  df <- df %>%
    mutate(diff_pairing_fraction = intron2_pairing_fraction - intron1_pairing_fraction)

  p_pf_sc <- cor_scatter(df, "diff_pairing_fraction",
    xlab = "Δ Pairing fraction (downstream - upstream)",
    title = "Pairing Fraction Difference vs Splicing Order",
    color = COL_DIFF)
  save_fig(p_pf_sc, "04b_pairing_fraction_scatter")
}

# =============================================================================
# SECTION 5: BINNED JUNCTION BOND ANALYSIS
# =============================================================================
message("\n-- Section 5: Binned junction bond analysis --")

bin_and_plot <- function(df, diff_col, title, col = COL_DIFF) {
  if (!diff_col %in% names(df)) return(invisible(NULL))
  d <- df %>%
    select(x = all_of(diff_col), y = fraction_downstream) %>%
    drop_na() %>%
    mutate(quintile = ntile(x, N_QUANTILES))

  bin_sum <- d %>%
    group_by(quintile) %>%
    summarise(
      med   = median(y),
      lo    = quantile(y, 0.025),
      hi    = quantile(y, 0.975),
      n     = n(),
      x_mid = median(x),
      .groups = "drop"
    )

  ggplot(bin_sum, aes(x_mid, med)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = col, alpha = 0.15) +
    geom_line(color = col, linewidth = 1) +
    geom_point(aes(size = n), color = col, shape = 21, fill = "white", stroke = 1.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    labs(x = paste(diff_col, "(quintile midpoints)"),
         y = "Median fraction downstream",
         title = title,
         subtitle = "Shading = 2.5-97.5 percentile; point size ~ n pairs") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 11))
}

p_bin_junc <- bin_and_plot(df, "diff_junction_hbonds_per_nt",
  "Splicing Order by Junction H-bond Quintile", COL_DIFF)
if (!is.null(p_bin_junc)) save_fig(p_bin_junc, "05a_junction_binned")

p_bin_5ss <- bin_and_plot(df, "diff_5ss_hbonds_per_nt",
  "Splicing Order by 5'SS H-bond Quintile", COL_5SS)
if (!is.null(p_bin_5ss)) save_fig(p_bin_5ss, "05b_5ss_binned")

p_bin_3ss <- bin_and_plot(df, "diff_3ss_hbonds_per_nt",
  "Splicing Order by 3'SS H-bond Quintile", COL_3SS)
if (!is.null(p_bin_3ss)) save_fig(p_bin_3ss, "05c_3ss_binned")

# =============================================================================
# SECTION 6: JUNCTION BONDS x ACCESSIBILITY INTERACTION
# =============================================================================
message("\n-- Section 6: Junction bonds x accessibility interaction --")

acc_col  <- "diff_ss5_accessibility_25nt"
junc_col <- "diff_junction_hbonds_per_nt"

if (all(c(acc_col, junc_col) %in% names(df))) {
  df_int <- df %>%
    select(fd = fraction_downstream,
           acc  = all_of(acc_col),
           junc = all_of(junc_col)) %>%
    drop_na() %>%
    mutate(
      acc_bin  = cut(acc,
        breaks = quantile(acc, probs = 0:3/3, na.rm = TRUE),
        labels = c("Low\n(similar acc.)", "Mid\n(moderate Δacc.)", "High\n(large Δacc.)"),
        include.lowest = TRUE),
      junc_label = factor(ntile(junc, 3),
        labels = c("Low Δjunc", "Mid Δjunc", "High Δjunc"))
    )

  bin_int <- df_int %>%
    group_by(acc_bin, junc_label) %>%
    summarise(med_fd = median(fd), n = n(), .groups = "drop")

  p_interact <- ggplot(bin_int,
      aes(acc_bin, med_fd, color = junc_label, group = junc_label)) +
    geom_line(linewidth = 1) +
    geom_point(aes(size = n)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = c("#2D9CDB","#F2994A","#EB5757"),
                       name = "Junction Δ H-bonds") +
    scale_size_continuous(range = c(3, 7), guide = "none") +
    labs(x = "Accessibility difference tertile (Δ 5'SS acc, 25nt)",
         y = "Median fraction downstream",
         title = "Junction H-bond Effect Conditional on Accessibility",
         subtitle = "H-bond effects most visible at intermediate accessibility differences") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  save_fig(p_interact, "06_junction_x_accessibility_interaction", w = 8, h = 5)
}

# =============================================================================
# SECTION 7: DISTRIBUTION OF DIFF METRICS
# =============================================================================
message("\n-- Section 7: Distributions --")

dist_cols <- intersect(
  c("diff_junction_hbonds_per_nt", "diff_5ss_hbonds_per_nt",
    "diff_3ss_hbonds_per_nt", "diff_hbonds_per_nt"),
  names(df)
)

if (length(dist_cols) > 0) {
  long_dist <- df %>%
    select(all_of(dist_cols)) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value") %>%
    drop_na() %>%
    mutate(metric = str_remove_all(metric, "diff_|_per_nt") %>%
             str_replace_all("_", " ") %>% str_to_title())

  p_dists <- ggplot(long_dist, aes(value, fill = metric, color = metric)) +
    geom_density(alpha = 0.4, linewidth = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~metric, scales = "free", ncol = 2) +
    scale_fill_brewer(palette = "Set1", guide = "none") +
    scale_color_brewer(palette = "Set1", guide = "none") +
    labs(x = "Δ value (downstream - upstream)",
         y = "Density",
         title = "Distributions of Structural Difference Metrics",
         subtitle = "Centered near 0 = no systematic upstream/downstream bias") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
  save_fig(p_dists, "07_diff_metric_distributions", w = 10, h = 6)
}

# =============================================================================
# SECTION 8: 5'SS vs 3'SS BAR CHART BY INTRON POSITION
# =============================================================================
message("\n-- Section 8: 5'SS vs 3'SS breakdown --")

ss_cols <- c("intron1_5ss_hbonds_per_nt", "intron1_3ss_hbonds_per_nt",
             "intron2_5ss_hbonds_per_nt", "intron2_3ss_hbonds_per_nt")

if (all(ss_cols %in% names(df))) {
  ss_long <- df %>%
    select(fraction_downstream, all_of(ss_cols)) %>%
    drop_na() %>%
    pivot_longer(all_of(ss_cols),
                 names_to  = c("intron", "boundary"),
                 names_pattern = "(intron[12])_(5ss|3ss)_hbonds_per_nt",
                 values_to = "hbonds_per_nt") %>%
    mutate(
      intron   = if_else(intron == "intron1", "Upstream (intron1)", "Downstream (intron2)"),
      boundary = if_else(boundary == "5ss", "5' Splice Site", "3' Splice Site")
    )

  ss_sum <- ss_long %>%
    group_by(intron, boundary) %>%
    summarise(
      med = median(hbonds_per_nt, na.rm = TRUE),
      q25 = quantile(hbonds_per_nt, .25, na.rm = TRUE),
      q75 = quantile(hbonds_per_nt, .75, na.rm = TRUE),
      .groups = "drop"
    )

  p_ss_bar <- ggplot(ss_sum, aes(boundary, med, fill = intron)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = q25, ymax = q75),
                  position = position_dodge(0.8), width = 0.25) +
    scale_fill_manual(values = c("Upstream (intron1)" = COL_UP,
                                 "Downstream (intron2)" = COL_DN),
                      name = "Intron position") +
    labs(x = "Splice site boundary",
         y = "H-bonds per nt (median +/- IQR)",
         title = "H-bond Load by Splice Site: Upstream vs Downstream",
         subtitle = "Which boundary contributes more to junction occlusion?") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  save_fig(p_ss_bar, "08_5ss_vs_3ss_by_intron_position", w = 7, h = 5)
}

# =============================================================================
# SECTION 9: CROSS-INTRON PAIRING (intron1-intron2 H-bonds)
# =============================================================================
message("\n-- Section 9: Cross-intron pairing --")

if ("intron1_intron2_hbonds" %in% names(df)) {
  df <- df %>%
    mutate(
      total_intron_len       = intron1_length + intron2_length,
      intron12_hbonds_per_nt = intron1_intron2_hbonds / total_intron_len
    )

  p_i12_hist <- ggplot(df %>% drop_na(intron12_hbonds_per_nt),
      aes(intron12_hbonds_per_nt)) +
    geom_histogram(bins = 60, fill = COL_DIFF, color = "white", linewidth = 0.2) +
    labs(x = "Intron1-Intron2 H-bonds per nt (combined length)",
         y = "Count",
         title = "Cross-Intron Base Pairing Distribution",
         subtitle = "Bonds bridging upstream and downstream introns directly") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  save_fig(p_i12_hist, "09a_cross_intron_distribution")

  p_i12_sc <- cor_scatter(df, "intron12_hbonds_per_nt",
    xlab = "Intron1-Intron2 H-bonds per nt",
    title = "Cross-Intron Pairing vs Splicing Order",
    color = COL_DIFF)
  save_fig(p_i12_sc, "09b_cross_intron_scatter")
}

# =============================================================================
# SECTION 10: STRONGLY ORDERED PAIRS — feature comparison
# =============================================================================
message("\n-- Section 10: Strongly ordered pairs --")

df_sig <- df %>%
  filter(fraction_downstream < 0.2 | fraction_downstream > 0.8) %>%
  mutate(splicing_order_bias = if_else(fraction_downstream > 0.5,
                                       "downstream_first", "upstream_first"))

message(sprintf("  Strongly ordered pairs (fd < 0.2 or > 0.8): %s",
                format(nrow(df_sig), big.mark = ",")))

feature_cols <- intersect(
  c("diff_junction_hbonds_per_nt", "diff_5ss_hbonds_per_nt",
    "diff_3ss_hbonds_per_nt", "diff_hbonds_per_nt",
    "diff_ss5_accessibility_25nt", "diff_ss3_accessibility_25nt",
    "diff_pairing_fraction"),
  names(df_sig)
)

if (nrow(df_sig) > 100 && length(feature_cols) > 0) {
  sig_long <- df_sig %>%
    select(splicing_order_bias, all_of(feature_cols)) %>%
    pivot_longer(all_of(feature_cols),
                 names_to = "feature", values_to = "value") %>%
    drop_na() %>%
    mutate(feature = str_remove_all(feature, "diff_|_per_nt") %>%
             str_replace_all("_", " ") %>% str_to_title())

  p_sig <- ggplot(sig_long,
      aes(splicing_order_bias, value, fill = splicing_order_bias)) +
    geom_violin(alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
    facet_wrap(~feature, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("downstream_first" = COL_DN,
                                 "upstream_first"   = COL_UP),
                      guide = "none") +
    scale_x_discrete(labels = c("downstream_first" = "Down first",
                                 "upstream_first"   = "Up first")) +
    labs(x = "Splicing order bias",
         y = "Delta value (downstream - upstream)",
         title = "Structural Differences in Strongly Ordered Pairs",
         subtitle = "Pairs with fraction_downstream < 0.2 or > 0.8") +
    theme_classic(base_size = 10) +
    theme(plot.title = element_text(face = "bold"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 9))
  save_fig(p_sig, "10_strongly_ordered_pairs_features", w = 12, h = 8)
}

# =============================================================================
# DONE
# =============================================================================
message("\n================================================")
message("All figures written to: ", OUT_DIR, "/")
message("Correlation table: ", file.path(OUT_DIR, "correlation_table.tsv"))
message("================================================\n")