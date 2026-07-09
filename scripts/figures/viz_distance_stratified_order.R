#!/usr/bin/env Rscript
# =============================================================================
# distance_stratified_order.R
#
# Plots fraction_downstream vs. inter-intron distance for the lariat paper,
# measured two ways:
#   (A) discrete separation  — adjacent (0), 1 apart, 2 apart, 3+ apart
#   (B) continuous genomic distance (intron2_start - intron1_end), binned
#
# Input:  significant_pairs.tsv  (pairs NOT necessarily adjacent)
# Requires annotation BED to compute how many introns apart each pair is.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(data.table); library(scales)
})

# ── Config ──────────────────────────────────────────────────────────────────
PAIRS_TSV  <- "/users/dhan30/splicing_order/data/significant_pairs.tsv"
INTRON_BED <- "/users/dhan30/reference/hg38.gencode.basic.v43.introns.bed.gz"
OUT_DIR    <- "figures/distance"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

theme_thesis <- function(base_size = 13) {
  theme_bw(base_size = base_size) +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "none",
          plot.title       = element_text(face = "bold"))
}
save_fig <- function(p, name, w = 7, h = 5) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT_DIR, paste0(name, ".png")), p, width = w, height = h,
         dpi = 200, bg = "white")
}

# ── Load pairs ──────────────────────────────────────────────────────────────
df <- read_tsv(PAIRS_TSV, show_col_types = FALSE) %>%
  mutate(pair_distance = intron2_start - intron1_end)   # genomic gap (bp)

# ── Compute separation (# annotated introns between the pair) ───────────────
# BED name field: ENST00000506640.3_intron_8_0_chr1_746819_r
# Parse transcript with the _intron_ convention (NOT the full name).
ann <- fread(INTRON_BED, header = FALSE,
             col.names = c("bchr","bstart","bend","name","score","strand"))
ann[, transcript := sub("_intron_.*", "", name)]

# Per-transcript sorted intron starts, keyed for fast lookup
ann_by_tx <- ann[, .(starts = list(sort(bstart))), by = transcript]
setkey(ann_by_tx, transcript)

df <- df %>% mutate(transcript = sub("_intron_.*", "", gene_id))

count_between <- function(tx, lo, hi) {
  s <- ann_by_tx[.(tx), starts][[1]]
  if (is.null(s)) return(NA_integer_)
  sum(s > lo & s < hi)                 # strictly between = introns skipped over
}

df <- df %>%
  rowwise() %>%
  mutate(introns_between = count_between(transcript, intron1_end, intron2_start)) %>%
  ungroup() %>%
  mutate(separation = introns_between)   # 0 = adjacent

n_unmapped <- sum(is.na(df$separation))
if (n_unmapped > 0)
  message(sprintf("Warning: %d pairs (%.1f%%) had no transcript match in BED",
                  n_unmapped, 100 * n_unmapped / nrow(df)))

# =============================================================================
# PLOT A — fraction_downstream by discrete separation (violin + box)
# =============================================================================
df_sep <- df %>%
  filter(!is.na(separation)) %>%
  mutate(sep_bin = case_when(
           separation == 0 ~ "Adjacent",
           separation == 1 ~ "1 apart",
           separation == 2 ~ "2 apart",
           separation == 3 ~ "3 apart",
           separation == 4 ~ "4 apart",
           separation == 5 ~ "5 apart",
           separation >= 6 ~ "6+ apart"),
         sep_bin = factor(sep_bin,
           levels = c("Adjacent","1 apart","2 apart","3 apart",
                      "4 apart","5 apart","6+ apart")))

# Per-bin n for axis labels
sep_n <- df_sep %>% count(sep_bin) %>%
  mutate(lab = paste0(sep_bin, "\n(n=", comma(n), ")"))

pA <- ggplot(df_sep, aes(sep_bin, fraction_downstream, fill = sep_bin)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey30") +
  scale_x_discrete(labels = setNames(sep_n$lab, sep_n$sep_bin)) +
  scale_fill_viridis_d(option = "mako", begin = 0.15, end = 0.9) +
  labs(x = "Intron separation within pair",
       y = "Fraction downstream spliced first",
       title = "Splicing order by intron separation") +
  theme_thesis()

save_fig(pA, "fd_by_separation_violin")

# =============================================================================
# PLOT B — fraction_downstream vs. continuous genomic distance
# Option B1: violin over distance quantile bins (recommended — compares medians)
# =============================================================================
df_dist <- df %>% filter(pair_distance > 0)

# Quintile bins of genomic distance, labelled with their bp ranges
df_dist <- df_dist %>%
  mutate(dist_bin = cut(pair_distance,
                        breaks = quantile(pair_distance, probs = seq(0, 1, 0.2),
                                          na.rm = TRUE),
                        include.lowest = TRUE, dig.lab = 6))

dist_n <- df_dist %>% count(dist_bin) %>%
  mutate(lab = paste0(dist_bin, "\n(n=", comma(n), ")"))

pB1 <- ggplot(df_dist, aes(dist_bin, fraction_downstream, fill = dist_bin)) +
  geom_violin(alpha = 0.5, scale = "width", trim = TRUE, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey30") +
  scale_x_discrete(labels = setNames(dist_n$lab, dist_n$dist_bin)) +
  scale_fill_viridis_d(option = "mako", begin = 0.15, end = 0.9) +
  labs(x = "Genomic distance between introns (bp, quintile bins)",
       y = "Fraction downstream spliced first",
       title = "Splicing order by genomic distance") +
  theme_thesis() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_fig(pB1, "fd_by_genomic_distance_violin", w = 8)

# =============================================================================
# Option B2: faceted histograms (if you prefer the faceted look you mentioned)
# =============================================================================
pB2 <- ggplot(df_dist, aes(fraction_downstream)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60,
                 fill = "#2E8B57", alpha = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ dist_bin, ncol = 1, strip.position = "right") +
  labs(x = "Fraction downstream spliced first", y = "Density",
       title = "Splicing order across genomic distance bins") +
  theme_thesis()

save_fig(pB2, "fd_by_genomic_distance_facet", w = 7, h = 8)

# =============================================================================
# Summary table for the paper
# =============================================================================
summ <- df_sep %>%
  group_by(sep_bin) %>%
  summarise(n = n(),
            mean_fd     = mean(fraction_downstream, na.rm = TRUE),
            median_fd   = median(fraction_downstream, na.rm = TRUE),
            median_dist = median(pair_distance, na.rm = TRUE),
            pct_downstream_biased = mean(direction == "downstream_biased",
                                         na.rm = TRUE),
            .groups = "drop")
print(summ)
write_tsv(summ, file.path(OUT_DIR, "separation_summary.tsv"))

message("Done. Figures + summary written to ", OUT_DIR)