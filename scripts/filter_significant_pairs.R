#!/usr/bin/env Rscript
# =============================================================================
# filter_significant_pairs.R
#
# Runs per-pair binomial significance testing on the pooled splicing order TSV
# produced by merge_and_split.py, applies FDR correction and an effect size
# filter, and writes:
#   - significant_pairs.tsv   filtered to sig pairs only (input for structure array)
#   - gene_summary.tsv        per-gene counts of sig/total pairs
#   - all_pairs_tested.tsv    full table with binom_p, binom_padj, effect, direction
#
# Usage:
#   Rscript filter_significant_pairs.R \
#       --input  merged/splicing_order_pooled.tsv \
#       --outdir merged/
#
# Options:
#   --input            Path to pooled TSV from merge_and_split.py  [required]
#   --outdir           Output directory                             [required]
#   --min-reads        Min total informative reads per pair         [default: 10]
#   --fdr              BH-adjusted p-value threshold                [default: 0.05]
#   --effect-threshold Min |fraction_downstream - 0.5|             [default: 0.25]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(scales)
})

option_list <- list(
  make_option("--input", type = "character", default = NULL,
              help = "Path to splicing_order_pooled.tsv [required]"),
  make_option("--outdir", type = "character", default = NULL,
              help = "Output directory [required]"),
  make_option("--min-reads", type = "integer", default = 10,
              help = "Min total informative reads per pair [default: 10]"),
  make_option("--fdr", type = "double", default = 0.05,
              help = "BH FDR threshold [default: 0.05]"),
  make_option("--effect-threshold", type = "double", default = 0.25,
              help = "Min |fraction_downstream - 0.5| [default: 0.25]")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input) || is.null(opt$outdir)) {
  stop("--input and --outdir are required. Run with --help for usage.")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Load

message("Loading: ", opt$input)
df <- read.delim(opt$input, stringsAsFactors = FALSE, check.names = FALSE)
message("  Rows loaded: ", scales::comma(nrow(df)))
message("  Columns: ", paste(colnames(df), collapse = ", "))

# Normalise column names — support both per-sample and pooled naming
rename_if_exists <- function(df, old, new) {
  if (old %in% names(df) && !new %in% names(df)) names(df)[names(df) == old] <- new
  df
}
df <- rename_if_exists(df, "upstream",   "upstream_count")
df <- rename_if_exists(df, "downstream", "downstream_count")
df <- rename_if_exists(df, "total",      "total_reads")

# Recompute fraction_downstream robustly
if (!"fraction_downstream" %in% names(df)) {
  df$fraction_downstream <- df$downstream_count / df$total_reads
}
df <- df[is.finite(df$fraction_downstream), ]

# Intron lengths
if (!"intron1_length" %in% names(df))
  df$intron1_length <- df$intron1_end - df$intron1_start
if (!"intron2_length" %in% names(df))
  df$intron2_length <- df$intron2_end - df$intron2_start

# Minimum reads filter
before <- nrow(df)
df <- filter(df, total_reads >= opt$`min-reads`)
message("\nAfter min-reads filter (>=", opt$`min-reads`, "): ",
        scales::comma(nrow(df)), " pairs (dropped ", scales::comma(before - nrow(df)), ")")

# ---------------------------------------------------------------------------
# Binomial test

message("\nRunning per-pair binomial tests (H0: fraction_downstream = 0.5)...")
message("  (this may take a minute for large datasets)")

df <- df |>
  rowwise() |>
  mutate(
    binom_p    = binom.test(downstream_count, total_reads, p = 0.5)$p.value,
    effect     = fraction_downstream - 0.5,   # signed: +ve = downstream-biased
    abs_effect = abs(effect)
  ) |>
  ungroup()

# BH correction across all pairs
df$binom_padj <- p.adjust(df$binom_p, method = "BH")

# Classify
df <- df |>
  mutate(
    sig = binom_padj < opt$fdr & abs_effect >= opt$`effect-threshold`,
    direction = case_when(
      sig & effect > 0 ~ "downstream_biased",
      sig & effect < 0 ~ "upstream_biased",
      TRUE             ~ "not_significant"
    )
  )

n_sig        <- sum(df$sig)
n_downstream <- sum(df$direction == "downstream_biased")
n_upstream   <- sum(df$direction == "upstream_biased")

message("  Significant pairs (FDR<", opt$fdr,
        ", |effect|>=", opt$`effect-threshold`, "): ",
        scales::comma(n_sig),
        "  (", scales::comma(n_downstream), " downstream-biased, ",
        scales::comma(n_upstream), " upstream-biased)")

# ---------------------------------------------------------------------------
# Gene-level summary

gene_summary <- df |>
  group_by(gene_id) |>
  summarise(
    n_pairs              = n(),
    n_sig                = sum(sig),
    n_downstream_biased  = sum(direction == "downstream_biased"),
    n_upstream_biased    = sum(direction == "upstream_biased"),
    mean_fraction        = mean(fraction_downstream),
    median_fraction      = median(fraction_downstream),
    mean_total_reads     = mean(total_reads),
    total_reads          = sum(total_reads),
    .groups = "drop"
  ) |>
  mutate(pct_sig = n_sig / n_pairs) |>
  arrange(desc(n_sig))

# Outputs

# 1. Significant pairs only
sig_pairs <- df |>
  filter(sig) |>
  arrange(desc(abs_effect), binom_padj)

out_sig <- file.path(opt$outdir, "significant_pairs.tsv")
write.table(sig_pairs, out_sig, sep = "\t", row.names = FALSE, quote = FALSE)
message("\nWritten: ", out_sig, "  (", scales::comma(nrow(sig_pairs)), " rows)")

# 2. Gene summary
out_gene <- file.path(opt$outdir, "gene_summary.tsv")
write.table(gene_summary, out_gene, sep = "\t", row.names = FALSE, quote = FALSE)
message("Written: ", out_gene, "  (", scales::comma(nrow(gene_summary)), " genes)")

# 3. Full tested table
out_all <- file.path(opt$outdir, "all_pairs_tested.tsv")
write.table(df, out_all, sep = "\t", row.names = FALSE, quote = FALSE)
message("Written: ", out_all, "  (", scales::comma(nrow(df)), " rows)")

# ---------------------------------------------------------------------------
# Summary

cat("\n", strrep("=", 60), "\n", sep = "")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n", sep = "")
cat(sprintf("  Total pairs tested:              %s\n",  scales::comma(nrow(df))))
cat(sprintf("  Significant (FDR<%.2f, |e|>=%.2f): %s  (%.1f%%)\n",
            opt$fdr, opt$`effect-threshold`,
            scales::comma(n_sig), 100 * n_sig / nrow(df)))
cat(sprintf("    Downstream-biased:             %s\n",  scales::comma(n_downstream)))
cat(sprintf("    Upstream-biased:               %s\n",  scales::comma(n_upstream)))
cat(sprintf("  Unique genes (all pairs):        %s\n",
            scales::comma(length(unique(df$gene_id)))))
cat(sprintf("  Unique genes (sig pairs):        %s\n",
            scales::comma(length(unique(sig_pairs$gene_id)))))
cat(sprintf("\n  Outputs written to: %s/\n", opt$outdir))
cat(strrep("=", 60), "\n", sep = "")