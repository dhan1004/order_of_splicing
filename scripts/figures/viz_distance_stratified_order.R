#!/usr/bin/env Rscript
# =============================================================================
# viz_distance_stratified_order.R
#
# Splicing order stratified by inter-intron distance, measured two ways:
#   (A) discrete separation  — rank2 - rank1 (1 = adjacent)
#   (B) continuous genomic distance (intron2_start - intron1_end), quintiles
#
# Input:  all_pairs_tested.tsv  (UNFILTERED — see note below)
# Requires annotation BED to assign true intron ranks.
#
# NOTE ON INPUT: this script uses all_pairs_tested.tsv, not significant_pairs.tsv.
# The latter is filtered on |fd - 0.5| >= EFFECT_THRESH, which carves an
# artificial hole around 0.5 and manufactures bimodality. Any distribution
# plotted from it shows the filter, not the biology.
#
# NOTE ON RANKS: BED intron_num is 1-based and in TRANSCRIPTION order
# (verified: on a minus-strand transcript, intron_num decreases as start
# increases). So rank = intron_num, and on minus-strand genes rank1 > rank2
# because intron1 is the leftmost by coordinate.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(readr)
  library(stringr); library(data.table); library(scales); library(patchwork)
})

# ── Config ──────────────────────────────────────────────────────────────────
PAIRS_TSV  <- "/users/dhan30/splicing_order/data/significant_pairs.tsv"
INTRON_BED <- "/users/dhan30/reference/hg38.gencode.basic.v43.introns.bed.gz"
OUT_DIR    <- "figures/distance/significant_pairs"
MIN_READS  <- 10       # drop low-coverage pairs; their fd is 0/1 by noise
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

# ── Annotation: coordinate -> (rank, strand) lookup ─────────────────────────
message("Loading annotation ...")
ann <- fread(INTRON_BED, header = FALSE,
             col.names = c("chr","start","end","name","score","strand"))

# name: ENST00000506640.3_intron_8_0_chr1_746819_r
# Split on _intron_ (critical — using the full name silently drops everything)
ann[, tx         := sub("_intron_.*", "", name)]
ann[, intron_num := as.integer(str_extract(sub(".*_intron_", "", name), "^[0-9]+"))]
ann <- ann[!is.na(intron_num)]
ann[, rank := intron_num]   # already 1-based, transcription order

# One row per coordinate. An intron shared across isoforms can carry different
# ranks; keep the modal rank so the join stays 1:1.
lookup <- ann[, .N, by = .(chr, start, end, rank, strand)][
  order(-N), .SD[1], by = .(chr, start, end)][, .(chr, start, end, rank, strand)]

message("  ", comma(nrow(lookup)), " unique introns")

# ── Load pairs, attach ranks ────────────────────────────────────────────────
df <- read_tsv(PAIRS_TSV, show_col_types = FALSE) %>%
  left_join(lookup %>% as_tibble() %>%
              rename(intron1_start = start, intron1_end = end,
                     rank1 = rank, strand1 = strand),
            by = c("chr", "intron1_start", "intron1_end")) %>%
  left_join(lookup %>% as_tibble() %>%
              rename(intron2_start = start, intron2_end = end,
                     rank2 = rank, strand2 = strand),
            by = c("chr", "intron2_start", "intron2_end")) %>%
  mutate(
    pair_distance = intron2_start - intron1_end,
    separation    = abs(rank2 - rank1),   # abs(): minus strand has rank1 > rank2
    is_first_pair = pmin(rank1, rank2) == 1L
  )

n_matched <- sum(!is.na(df$separation))
message(sprintf("  ranks matched: %s / %s (%.1f%%)",
                comma(n_matched), comma(nrow(df)), 100 * n_matched / nrow(df)))
if (n_matched / nrow(df) < 0.8)
  warning("Low rank match rate. Exact endpoint join may be failing — ",
          "consider bedtools intersect -f 1.0 -r instead.")

# ── SANITY CHECK: did the pipeline strand-correct fraction_downstream? ──────
# If it did, plus and minus strand means should agree. If they mirror around
# 0.5, half the data has inverted fd and everything downstream is wrong.
strand_check <- df %>%
  filter(!is.na(strand1)) %>%
  group_by(strand1) %>%
  summarise(mean_fd = mean(fraction_downstream, na.rm = TRUE), n = n(),
            .groups = "drop")
message("\nStrand check (means should AGREE, not mirror around 0.5):")
print(strand_check)
if (nrow(strand_check) == 2 &&
    abs(sum(strand_check$mean_fd) - 1) < 0.01 &&
    abs(diff(strand_check$mean_fd)) > 0.01) {
  warning("Plus/minus strand fd means mirror around 0.5. The pipeline likely ",
          "did NOT strand-correct. Fix upstream before trusting these figures.")
}

# ── Filter ──────────────────────────────────────────────────────────────────
df <- df %>%
  filter(!is.na(separation), total >= MIN_READS, pair_distance > 0)
message("\nPairs after filtering: ", comma(nrow(df)))

# =============================================================================
# PLOT A — fraction_downstream by discrete separation
# Faceted histograms, not violins: fd is bounded [0,1] and U-shaped, so KDE
# smears density past the bounds and invents shape that isn't there.
# =============================================================================
df_sep <- df %>%
  mutate(sep_bin = case_when(
           separation == 1 ~ "Adjacent",
           separation == 2 ~ "1 between",
           separation == 3 ~ "2 between",
           separation == 4 ~ "3 between",
           separation >= 5 ~ "4+ between"),
         sep_bin = factor(sep_bin,
           levels = c("Adjacent","1 between","2 between",
                      "3 between","4+ between")))

sep_n <- df_sep %>% count(sep_bin) %>%
  mutate(lab = paste0(sep_bin))

df_sep <- bind_rows(
df_sep %>% mutate(pair_type = "All pairs"),
df_sep %>% mutate(pair_type = ifelse(is_first_pair, "First intron", "Internal"))
) %>%
mutate(sep_bin   = factor(sep_bin, levels = c("Adjacent","1 between","2 between",
                                              "3 between","4+ between")),
        pair_type = factor(pair_type, levels = c("All pairs","First intron","Internal")))

pA <- ggplot(df_sep, aes(fraction_downstream)) +
geom_histogram(aes(y = after_stat(density)), bins = 50,
                fill = "#2E8B57", alpha = 0.85) +
geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey30") +
facet_grid(sep_bin ~ pair_type, switch = "y") +
labs(x = "Fraction downstream spliced first", y = "Density",
      title = "Splicing order by intron separation") +
theme_thesis() +
theme(strip.text.y.left = element_text(angle = 0))

save_fig(pA, "fd_by_separation_facet", w = 9, h = 8)

pA <- ggplot(df_sep, aes(fraction_downstream)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#2E8B57", alpha = 0.85) +
  facet_wrap(~ sep_bin, ncol = 1, strip.position = "right",
             labeller = labeller(sep_bin = setNames(sep_n$lab, sep_n$sep_bin))) +
  labs(x = "Fraction downstream spliced first", y = "Density",
       title = "Splicing order by intron separation") +
  theme_thesis()

save_fig(pA, "fd_by_separation_facet_all", w = 7, h = 8)

df <- df %>%
  mutate(dist_bin = cut(pair_distance,
                        breaks = quantile(pair_distance, probs = seq(0, 1, 0.2),
                                          na.rm = TRUE),
                        include.lowest = TRUE, dig.lab = 6))

summ_dist <- df %>%
  group_by(dist_bin, is_first_pair) %>%
  summarise(n       = n(),
            mean_fd = mean(fraction_downstream, na.rm = TRUE),
            se      = sd(fraction_downstream, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(lo = mean_fd - 1.96 * se,
         hi = mean_fd + 1.96 * se,
         pair_type = ifelse(is_first_pair, "First intron pairs",
                            "Internal pairs"))

pB_all <- ggplot(df, aes(fraction_downstream)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#2E8B57", alpha = 0.85) +
  facet_wrap(~ dist_bin, ncol = 1, strip.position = "right") +
  labs(x = "Fraction downstream spliced first", y = "Density",
       title = "Splicing order across genomic distance bins",
       subtitle = "All pairs pooled") +
  theme_thesis()

save_fig(pB_all, "fd_by_genomic_distance_facet_pooled", w = 7, h = 8)

df_dist <- bind_rows(
  df %>% mutate(pair_type = "All pairs"),
  df %>% mutate(pair_type = ifelse(is_first_pair, "First intron", "Internal"))
) %>%
  mutate(pair_type = factor(pair_type, levels = c("All pairs","First intron","Internal")))

pC <- ggplot(df_dist, aes(fraction_downstream)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "#2E8B57", alpha = 0.85) +
  facet_grid(dist_bin ~ pair_type, switch = "y") +
  labs(x = "Fraction downstream spliced first", y = "Density",
       title = "Splicing order across genomic distance bins") +
  theme_thesis() +
  theme(strip.text.y.left = element_text(angle = 0))

save_fig(pC, "fd_by_genomic_distance_facet_all", w = 9, h = 9)

# =============================================================================
# Summary tables
# =============================================================================
summ_sep <- df_sep %>%
  group_by(sep_bin, is_first_pair) %>%
  summarise(n           = n(),
            mean_fd     = mean(fraction_downstream, na.rm = TRUE),
            median_fd   = median(fraction_downstream, na.rm = TRUE),
            median_dist = median(pair_distance, na.rm = TRUE),
            # fraction of pairs above 0.5 — more interpretable than the mean
            # for a bimodal bounded variable
            pct_above_half = mean(fraction_downstream > 0.5, na.rm = TRUE),
            .groups = "drop")

print(summ_sep)
write_tsv(summ_sep,   file.path(OUT_DIR, "separation_summary.tsv"))
write_tsv(summ_dist,  file.path(OUT_DIR, "distance_summary.tsv"))
write_tsv(strand_check, file.path(OUT_DIR, "strand_check.tsv"))

message("\nDone. Figures + summaries written to ", OUT_DIR)