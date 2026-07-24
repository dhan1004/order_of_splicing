#!/usr/bin/env Rscript
# =============================================================================
# viz_order_enrichment.R  (combined v3)
#
# Merged from viz_order_enrichment.R and viz_order_enrichment_new.R.
#
# Significance testing (binomial + BH) is performed UPSTREAM of this script.
# This script expects padj / direction / abs_effect columns already present.
# If they are absent it will compute them inline as a fallback.
#
# Figures produced:
#   01  fraction_downstream_histogram
#   02  first_vs_internal_pairs
#   03  effect_size_volcano
#   04  top_downstream_pairs          (with Jeffreys CI error bars)
#   05  top_upstream_pairs            (with Jeffreys CI error bars)
#   05a gene_sig_pair_counts          (stacked bar, ranked by n_sig)
#   05b gene_sig_pair_counts_norm     (pct_sig, normalized for gene size)
#   05c gene_sig_pairs_scatter        (n_pairs vs n_sig, coloured by direction)
#   05d gene_effect_size_ranked       (mean effect, downstream + upstream panels)
#   05e gene_effect_ranked_norm       (pct_sig ranked by |mean effect|)
#   05f gene_sig_pair_counts_effect   (n_sig ranked by mean |effect|)
#   06a intron_length_density         (all vs sig, log10 x)
#   06b intron_length_violin          (all vs sig, violin + boxplot)
#   06c paired_length_hexbin          (hexbin + patchwork marginals)
#   06d paired_length_scatter_sig     (sig pairs only, coloured by frac_down)
#   06e length_ratio_by_direction     (log2 ratio density by direction)
#   07  read_coverage                 (histogram, all vs sig)
#   08  GO_enrichment_downstream      (GeneRatio x-axis)
#   09  GO_enrichment_upstream
#   10  GO_enrichment_combined        (faceted)
#   significant_pairs.tsv
#   gene_summary.tsv
#
# Usage:
# Rscript scripts/figures/viz_order_enrichment.R \
#     --input  results/mouse/splicing_order_pooled.tsv \
#     --outdir figures/mouse \
#     --gene-bed /users/dhan30/reference/mm39.gencode.basic.vM36.genes.bed \
#     --species mouse
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(stringr)
  library(ggrepel)
  library(patchwork)
})

select <- dplyr::select

# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

option_list <- list(
  make_option("--input",  type = "character", default = NULL,
              help = "Path to splicing order TSV (pooled)"),
  make_option("--outdir", type = "character", default = "enrichment_figures"),
  make_option("--min-reads",        type = "integer", default = 10),
  make_option("--fdr",              type = "double",  default = 0.05),
  make_option("--effect-threshold", type = "double",  default = 0.25),
  make_option("--top-n",            type = "integer", default = 30),
  make_option("--gene-bed",         type = "character", default = NULL,
              help = "6-col BED for coordinate-based gene symbol mapping"),
  make_option("--species", default = "human", help = "human or mouse")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$input)) stop("--input is required.")
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

orgdb_name <- if (opt$species == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
if (!requireNamespace(orgdb_name, quietly = TRUE))
  stop("Install ", orgdb_name)
suppressPackageStartupMessages(library(orgdb_name, character.only = TRUE))
ORGDB <- get(orgdb_name)

# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

save_fig <- function(p, stem, w = 9, h = 7) {
  ggsave(file.path(opt$outdir, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(opt$outdir, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
  message("  Saved: ", stem)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

theme_thesis <- function() {
  theme_bw(base_size = 12) +
    theme(panel.grid.minor  = element_blank(),
          strip.background  = element_rect(fill = "grey92"),
          legend.position   = "right")
}

# Cohesive palette used throughout
col_downstream <- "#E05C5C"
col_upstream   <- "#4B8BBE"
col_nonsig     <- "#AAAAAA"

direction_colors <- c(downstream_biased = col_downstream,
                      upstream_biased   = col_upstream,
                      not_significant   = col_nonsig)
direction_labels <- c(downstream_biased = "Downstream-biased",
                      upstream_biased   = "Upstream-biased",
                      not_significant   = "Not significant")

MIN_PAIRS <- 10   # min intron pairs per gene for normalized figures

# --------------------------------------------------------------------------- #
# Load & prepare data
# --------------------------------------------------------------------------- #

message("Loading ", opt$input, "...")
df <- read.table(opt$input, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
message("  Rows loaded: ", nrow(df))

if (!"intron1_length" %in% names(df))
  df$intron1_length <- df$intron1_end - df$intron1_start
if (!"intron2_length" %in% names(df))
  df$intron2_length <- df$intron2_end - df$intron2_start

df <- df |> mutate(fraction_downstream = downstream / total)
df$gene_id <- sub("\\.\\d+$", "", df$gene_id)

df <- df |> filter(total >= opt$`min-reads`)
message("  After min_reads>=", opt$`min-reads`, ": ", nrow(df), " pairs")

# --------------------------------------------------------------------------- #
# Significance columns
# Use pre-computed columns if present; compute inline as fallback.
# --------------------------------------------------------------------------- #

if (!"padj" %in% names(df)) {
  message("  'padj' not found — computing binomial tests inline...")
  df <- df |>
    rowwise() |>
    mutate(pval = binom.test(downstream, total, p = 0.5,
                             alternative = "two.sided")$p.value) |>
    ungroup() |>
    mutate(padj = p.adjust(pval, method = "BH"))
}

if (!"abs_effect" %in% names(df))
  df <- df |> mutate(abs_effect = abs(fraction_downstream - 0.5))

if (!"effect" %in% names(df))
  df <- df |> mutate(effect = fraction_downstream - 0.5)

if (!"direction" %in% names(df)) {
  df <- df |> mutate(
    direction = case_when(
      padj < opt$fdr & fraction_downstream > 0.5 ~ "downstream_biased",
      padj < opt$fdr & fraction_downstream < 0.5 ~ "upstream_biased",
      TRUE ~ "not_significant"
    ))
}

if (!"sig" %in% names(df))
  df <- df |> mutate(sig = direction != "not_significant")

df$direction <- factor(df$direction,
                       levels = c("downstream_biased", "upstream_biased", "not_significant"))

n_sig  <- sum(df$sig)
n_down <- sum(df$direction == "downstream_biased")
n_up   <- sum(df$direction == "upstream_biased")
message("  Sig pairs: ", n_sig, "  (downstream: ", n_down, ", upstream: ", n_up, ")")

# --------------------------------------------------------------------------- #
# Gene symbol mapping: BED > org.Hs.eg.db > biomaRt > gene_id fallback
# --------------------------------------------------------------------------- #

message("\nMapping gene symbols...")

do_coordinate_mapping <- function(df, bed_path) {
  message("  Method 1: coordinate overlap from ", bed_path)
  first_line   <- readLines(bed_path, n = 1)
  first_fields <- strsplit(first_line, "\t")[[1]]
  has_header   <- is.na(suppressWarnings(as.integer(first_fields[2])))
  genes <- tryCatch({
    raw <- read.table(bed_path, header = has_header, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "#")
    if (ncol(raw) < 4) return(NULL)
    colnames(raw)[1:4] <- c("chr", "start", "end", "gene_symbol")
    raw$start <- as.integer(raw$start); raw$end <- as.integer(raw$end)
    raw[!is.na(raw$start) & !is.na(raw$end), ]
  }, error = function(e) { message("  BED load failed: ", e$message); NULL })
  if (is.null(genes) || nrow(genes) == 0) return(NULL)
  message("  Loaded ", scales::comma(nrow(genes)), " gene records")
  genes_by_chr <- split(genes, genes$chr)
  mapply(function(ch, ps, pe) {
    g <- genes_by_chr[[ch]]
    if (is.null(g)) return(NA_character_)
    hits <- g[g$start < pe & g$end > ps, "gene_symbol", drop = TRUE]
    if (length(hits) > 0) hits[1] else NA_character_
  }, df$chr, df$intron1_start, df$intron2_end, SIMPLIFY = TRUE, USE.NAMES = FALSE)
}

do_orgdb_mapping <- function(df) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) return(NULL)
  message("  Method 2: org.Hs.eg.db mapping")
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  select <- dplyr::select
  all_ids   <- unique(df$gene_id)
  sample_id <- all_ids[1]
  id_type <- dplyr::case_when(
    grepl("^ENST",    sample_id) ~ "ENSEMBLTRANS",
    grepl("^ENSG",    sample_id) ~ "ENSEMBL",
    grepl("^[0-9]+$", sample_id) ~ "ENTREZID",
    TRUE                         ~ "SYMBOL"
  )
  if (id_type == "SYMBOL") return(setNames(all_ids, all_ids))
  if (id_type == "ENSEMBLTRANS") {
    ensg <- tryCatch(suppressMessages(
      AnnotationDbi::mapIds(ORGDB, keys = all_ids,
                            column = "ENSEMBL", keytype = "ENSEMBLTRANS",
                            multiVals = "first")),
      error = function(e) NULL)
    if (is.null(ensg)) return(NULL)
    sym <- tryCatch(suppressMessages(
      AnnotationDbi::mapIds(ORGDB, keys = na.omit(unique(as.character(ensg))),
                            column = "SYMBOL", keytype = "ENSEMBL",
                            multiVals = "first")),
      error = function(e) NULL)
    if (is.null(sym)) return(NULL)
    setNames(as.character(sym[as.character(ensg[all_ids])]), all_ids)
  } else {
    tryCatch(suppressMessages(
      AnnotationDbi::mapIds(ORGDB, keys = all_ids,
                            column = "SYMBOL", keytype = id_type,
                            multiVals = "first")),
      error = function(e) { message("  org.Hs.eg.db failed: ", e$message); NULL })
  }
}

do_biomart_mapping <- function(df) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) return(NULL)
  message("  Method 3: biomaRt coordinate lookup (requires internet)...")
  tryCatch({
    mart    <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    regions <- df |> dplyr::select(chr, intron1_start, intron2_end) |>
      dplyr::distinct() |>
      dplyr::mutate(region = paste0(chr, ":", intron1_start, ":", intron2_end))
    res <- biomaRt::getBM(
      attributes = c("chromosome_name", "start_position", "end_position",
                     "external_gene_name"),
      filters    = "chromosomal_region",
      values     = regions$region, mart = mart)
    if (nrow(res) == 0) return(NULL)
    result <- character(nrow(df))
    for (i in seq_len(nrow(df))) {
      ch   <- sub("^chr", "", df$chr[i])
      hits <- res[res$chromosome_name == ch &
                  res$start_position  <  df$intron2_end[i] &
                  res$end_position    >  df$intron1_start[i], ]
      result[i] <- if (nrow(hits) > 0) hits$external_gene_name[1] else NA_character_
    }
    result
  }, error = function(e) { message("  biomaRt failed: ", e$message); NULL })
}

gene_symbols <- NULL
if (!is.null(opt$`gene-bed`) && nchar(opt$`gene-bed`) > 0)
  gene_symbols <- do_coordinate_mapping(df, opt$`gene-bed`)

if (is.null(gene_symbols) || mean(!is.na(gene_symbols)) < 0.5) {
  lookup <- do_orgdb_mapping(df)
  if (!is.null(lookup)) {
    candidate <- as.character(lookup[df$gene_id])
    if (sum(!is.na(candidate)) > sum(!is.na(gene_symbols %||% rep(NA, nrow(df)))))
      gene_symbols <- candidate
  }
}

if (is.null(gene_symbols) || mean(!is.na(gene_symbols)) < 0.5)
  gene_symbols <- do_biomart_mapping(df)

if (!is.null(gene_symbols)) {
  df$gene_symbol <- ifelse(is.na(gene_symbols) | gene_symbols == "",
                           df$gene_id, gene_symbols)
} else {
  message("  All mapping methods failed — using gene_id as gene_symbol")
  df$gene_symbol <- df$gene_id
}

n_mapped <- sum(df$gene_symbol != df$gene_id)
message("  Mapped ", scales::comma(n_mapped), "/", scales::comma(nrow(df)),
        " rows (", round(100 * n_mapped / nrow(df), 1), "%)")

mapping_succeeded <- n_mapped / nrow(df) > 0.5

# First intron pair flag
df <- df |>
  group_by(gene_symbol) |>
  mutate(is_first_pair = intron1_start == min(intron1_start)) |>
  ungroup()

# --------------------------------------------------------------------------- #
# Gene-level summary (computed once, used in many figures)
# --------------------------------------------------------------------------- #

gene_summary <- df |>
  group_by(gene_symbol) |>
  summarise(
    n_pairs             = n(),
    n_sig               = sum(sig),
    n_downstream_biased = sum(direction == "downstream_biased"),
    n_upstream_biased   = sum(direction == "upstream_biased"),
    mean_fraction       = mean(fraction_downstream, na.rm = TRUE),
    median_fraction     = median(fraction_downstream, na.rm = TRUE),
    mean_total_reads    = mean(total, na.rm = TRUE),
    total_reads         = sum(total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    pct_sig = n_sig / n_pairs,
    dominant_direction = case_when(
      n_downstream_biased > n_upstream_biased ~ "downstream",
      n_upstream_biased   > n_downstream_biased ~ "upstream",
      TRUE ~ "mixed"
    )
  ) |>
  arrange(desc(n_sig))

write.table(gene_summary,
            file.path(opt$outdir, "gene_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("  gene_summary.tsv written")

# --------------------------------------------------------------------------- #
# FIG 01: Fraction downstream histogram
# --------------------------------------------------------------------------- #

message("\n--- Fig 01: Fraction downstream histogram ---")
{
  mu   <- mean(df$fraction_downstream, na.rm = TRUE)
  tval <- t.test(df$fraction_downstream, mu = 0.5)
  p <- ggplot(df, aes(x = fraction_downstream)) +
    geom_histogram(bins = 60, fill = "steelblue", color = "black",
                   linewidth = 0.2, alpha = 0.85) +
    geom_vline(xintercept = 0.5, color = "red",     linewidth = 1.2, linetype = "dashed") +
    geom_vline(xintercept = mu,  color = "darkblue", linewidth = 1.2) +
    annotate("text", x = mu + 0.01, y = Inf, vjust = 1.5, hjust = 0,
             label = sprintf("Mean = %.3f\nt = %.1f, p = %.2e",
                             mu, tval$statistic, tval$p.value), size = 3.5) +
    labs(x = "Fraction Downstream Spliced First", y = "Number of Intron Pairs",
         title = "Distribution of Splicing Order",
         subtitle = sprintf("n = %s intron pairs", scales::comma(nrow(df)))) +
    theme_thesis()
  save_fig(p, "01_fraction_downstream_histogram")
}

# --------------------------------------------------------------------------- #
# FIG 02: First intron pairs vs internal pairs
# --------------------------------------------------------------------------- #

message("--- Fig 02: First vs internal intron pairs ---")
{
  mu_first    <- mean(df$fraction_downstream[df$is_first_pair],  na.rm = TRUE)
  mu_internal <- mean(df$fraction_downstream[!df$is_first_pair], na.rm = TRUE)
  p <- ggplot(df, aes(x = fraction_downstream,
                      fill = is_first_pair, color = is_first_pair)) +
    geom_histogram(aes(y = after_stat(density)), bins = 60,
                   position = "identity", alpha = 0.55, linewidth = 0.2) +
    geom_vline(xintercept = 0.5,        color = "black",        linewidth = 1,
               linetype = "dashed", alpha = 0.6) +
    geom_vline(xintercept = mu_first,    color = col_downstream, linewidth = 1.2) +
    geom_vline(xintercept = mu_internal, color = col_upstream,   linewidth = 1.2) +
    scale_fill_manual(values = c(`TRUE` = col_downstream, `FALSE` = col_upstream),
                      labels = c(`TRUE` = "First intron pair", `FALSE` = "All other pairs")) +
    scale_color_manual(values = c(`TRUE` = col_downstream, `FALSE` = col_upstream),
                       labels = c(`TRUE` = "First intron pair", `FALSE` = "All other pairs")) +
    annotate("text", x = 0.03, y = Inf, vjust = 1.5, hjust = 0,
             label = sprintf("First pairs mean = %.3f\nInternal pairs mean = %.3f",
                             mu_first, mu_internal), size = 3.5) +
    labs(x = "Fraction Downstream Spliced First", y = "Density",
         fill = NULL, color = NULL,
         title = "First Intron Pairs Drive Global Downstream Bias",
         subtitle = sprintf("First: n = %s  |  Internal: n = %s",
                            scales::comma(sum(df$is_first_pair)),
                            scales::comma(sum(!df$is_first_pair)))) +
    theme_thesis()
  save_fig(p, "02_first_vs_internal_pairs")
}

# --------------------------------------------------------------------------- #
# FIG 03: Volcano
# --------------------------------------------------------------------------- #

message("--- Fig 03: Volcano ---")
{
  plot_df <- df |>
    mutate(neg_logp = pmin(-log10(padj + 1e-300), 50))

  top_labels <- bind_rows(
    plot_df |> filter(direction == "downstream_biased") |>
      slice_max(abs_effect, n = opt$`top-n` %/% 2),
    plot_df |> filter(direction == "upstream_biased") |>
      slice_max(abs_effect, n = opt$`top-n` %/% 2)
  )

  p <- ggplot(plot_df, aes(x = effect, y = neg_logp,
                            color = direction, size = log10(total + 1))) +
    geom_point(data = filter(plot_df, direction == "not_significant"),
               alpha = 0.2, size = 0.8) +
    geom_point(data = filter(plot_df, direction != "not_significant"),
               alpha = 0.6) +
    geom_hline(yintercept = -log10(opt$fdr), linetype = "dashed",
               color = "grey40", linewidth = 0.8) +
    geom_vline(xintercept = c(-opt$`effect-threshold`, opt$`effect-threshold`),
               linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_text_repel(data = top_labels, aes(label = gene_symbol),
                    size = 2.5, max.overlaps = 20, min.segment.length = 0,
                    color = "black") +
    scale_color_manual(values = direction_colors, labels = direction_labels) +
    scale_size_continuous(name = "log\u2081\u2080(reads)", range = c(0.5, 4)) +
    labs(x = "Effect size (fraction_downstream \u2212 0.5)",
         y = expression(-log[10](p[adj])),
         color = NULL,
         title = "Volcano: Splicing Order Effect Size vs Significance",
         subtitle = sprintf("Downstream: %s  |  Upstream: %s  |  FDR < %.2f",
                            scales::comma(n_down), scales::comma(n_up), opt$fdr)) +
    theme_thesis()
  save_fig(p, "03_effect_size_volcano", w = 11, h = 8)
}

# --------------------------------------------------------------------------- #
# FIG 04 & 05: Top biased pairs with Jeffreys CI error bars
# --------------------------------------------------------------------------- #

message("--- Fig 04/05: Top biased pairs ---")
make_top_pairs_bar <- function(direction_filter, fill_col, stem) {
  top <- df |>
    filter(direction == direction_filter) |>
    arrange(desc(abs_effect)) |>
    slice_head(n = min(opt$`top-n`, 40)) |>
    mutate(
      pair_label = paste0(gene_symbol, "\n",
                          chr, ":", format(intron1_start, big.mark = ","),
                          "-", format(intron2_end, big.mark = ",")),
      pair_label = factor(pair_label, levels = rev(pair_label))
    )
  if (nrow(top) == 0) return(NULL)
  p <- ggplot(top, aes(x = pair_label, y = fraction_downstream)) +
    geom_col(fill = fill_col, alpha = 0.85) +
    geom_errorbar(aes(
      ymin = qbeta(0.025, downstream + 0.5, upstream + 0.5),
      ymax = qbeta(0.975, downstream + 0.5, upstream + 0.5)
    ), width = 0.3, color = "grey30") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1), labels = percent_format()) +
    labs(x = NULL, y = "Fraction downstream spliced first",
         title = paste0("Top ", nrow(top), " ", gsub("_", "-", direction_filter), " pairs"),
         subtitle = "Error bars: 95% Jeffreys credible interval | Sorted by |effect size|") +
    theme_thesis() +
    theme(axis.text.y = element_text(size = 8))
  save_fig(p, stem, w = 10, h = max(5, 0.35 * nrow(top) + 2))
}

make_top_pairs_bar("downstream_biased", col_downstream, "04_top_downstream_pairs")
make_top_pairs_bar("upstream_biased",   col_upstream,   "05_top_upstream_pairs")

# --------------------------------------------------------------------------- #
# FIG 05a: Gene sig pair counts — stacked bar ranked by n_sig
# --------------------------------------------------------------------------- #

message("--- Fig 05a: Gene sig pair counts (n_sig) ---")
{
  top_genes <- gene_summary |>
    filter(n_pairs >= 2) |>
    arrange(desc(n_sig)) |>
    slice_head(n = 40) |>
    mutate(gene_symbol = factor(gene_symbol, levels = rev(gene_symbol))) |>
    pivot_longer(c(n_downstream_biased, n_upstream_biased),
                 names_to = "direction", values_to = "count") |>
    mutate(direction = recode(direction,
                              n_downstream_biased = "Downstream-biased",
                              n_upstream_biased   = "Upstream-biased"))
  p <- ggplot(top_genes, aes(x = gene_symbol, y = count, fill = direction)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("Downstream-biased" = col_downstream,
                                 "Upstream-biased"   = col_upstream)) +
    labs(x = NULL, y = "Number of significant intron pairs", fill = NULL,
         title = "Genes with most significant intron pairs") +
    theme_thesis() + theme(axis.text.y = element_text(size = 8))
  save_fig(p, "05a_gene_sig_pair_counts", w = 10, h = 9)
}

# --------------------------------------------------------------------------- #
# FIG 05b: Normalized — pct_sig, ranked by pct_sig
# --------------------------------------------------------------------------- #

message("--- Fig 05b: Gene sig pair counts normalized ---")
{
  top_norm <- gene_summary |>
    filter(n_pairs >= MIN_PAIRS) |>
    arrange(desc(pct_sig)) |>
    slice_head(n = 40) |>
    mutate(gene_symbol    = factor(gene_symbol, levels = rev(gene_symbol)),
           pct_downstream = n_downstream_biased / n_pairs,
           pct_upstream   = n_upstream_biased   / n_pairs)
  long <- top_norm |>
    pivot_longer(c(pct_downstream, pct_upstream),
                 names_to = "direction", values_to = "pct") |>
    mutate(direction = recode(direction,
                              pct_downstream = "Downstream-biased",
                              pct_upstream   = "Upstream-biased"))
  label_df <- top_norm |>
    mutate(label = paste0(round(pct_sig * 100, 1), "%  (n=", n_sig, ")"))
  p <- ggplot(long, aes(x = gene_symbol, y = pct, fill = direction)) +
    geom_col() +
    geom_text(data = label_df,
              aes(x = gene_symbol, y = pct_sig, label = label),
              hjust = -0.08, size = 2.6, inherit.aes = FALSE) +
    coord_flip(clip = "off") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.22))) +
    scale_fill_manual(values = c("Downstream-biased" = col_downstream,
                                 "Upstream-biased"   = col_upstream)) +
    labs(x = NULL, fill = NULL,
         y = "Fraction of intron pairs with significant splicing order",
         title = "Genes with highest fraction of strongly ordered intron pairs",
         subtitle = paste0("Genes with \u2265", MIN_PAIRS,
                           " testable pairs; ranked by % significant")) +
    theme_thesis() + theme(axis.text.y = element_text(size = 8))
  save_fig(p, "05b_gene_sig_pair_counts_norm", w = 11, h = 9)
}

# --------------------------------------------------------------------------- #
# FIG 05c: Scatter — n_pairs vs n_sig
# --------------------------------------------------------------------------- #

message("--- Fig 05c: n_pairs vs n_sig scatter ---")
{
  null_rate  <- sum(gene_summary$n_sig) / sum(gene_summary$n_pairs)
  scatter_df <- gene_summary |>
    filter(n_pairs >= MIN_PAIRS) |>
    mutate(
      label_gene = case_when(
        pct_sig >= quantile(pct_sig, 0.92) ~ gene_symbol,
        n_sig   >= quantile(n_sig,   0.92) ~ gene_symbol,
        TRUE ~ NA_character_)
    )
  p <- ggplot(scatter_df, aes(x = n_pairs, y = n_sig,
                               color = dominant_direction, size = pct_sig)) +
    geom_abline(slope = null_rate, intercept = 0,
                linetype = "dashed", color = "grey50", linewidth = 0.7) +
    geom_point(alpha = 0.6) +
    geom_text_repel(aes(label = label_gene), size = 2.8, max.overlaps = 20,
                    segment.color = "grey60", segment.size = 0.3,
                    show.legend = FALSE) +
    scale_color_manual(
      values = c(downstream = col_downstream, upstream = col_upstream,
                 mixed = col_nonsig),
      labels = c(downstream = "Downstream-biased", upstream = "Upstream-biased",
                 mixed = "Mixed / equal"),
      name = "Dominant direction") +
    scale_size_continuous(name = "Fraction significant",
                          range = c(1, 6), labels = percent_format(accuracy = 1)) +
    scale_x_continuous(labels = comma_format()) +
    scale_y_continuous(labels = comma_format()) +
    annotate("text", x = max(scatter_df$n_pairs) * 0.72,
             y = max(scatter_df$n_pairs) * null_rate * 0.85,
             label = paste0("Null: ", round(null_rate * 100, 1), "% sig rate"),
             color = "grey40", size = 3, hjust = 0) +
    labs(x = "Total adjacent intron pairs tested",
         y = "Significant intron pairs",
         title = "Ordered splicing scales with intron count but varies in enrichment",
         subtitle = paste0("Each point = one gene (\u2265", MIN_PAIRS,
                           " pairs); dashed = genome-wide sig rate")) +
    theme_thesis()
  save_fig(p, "05c_gene_sig_pairs_scatter", w = 11, h = 8)
}

# --------------------------------------------------------------------------- #
# FIG 05d: Mean effect size per gene — downstream + upstream panels
# --------------------------------------------------------------------------- #

message("--- Fig 05d: Gene effect size ranked ---")
{
  gene_effect <- df |>
    filter(sig) |>
    group_by(gene_symbol) |>
    summarise(
      mean_effect         = mean(fraction_downstream - 0.5),
      n_sig               = n(),
      n_downstream_biased = sum(direction == "downstream_biased"),
      n_upstream_biased   = sum(direction == "upstream_biased"),
      .groups = "drop"
    ) |>
    left_join(gene_summary |> select(gene_symbol, n_pairs), by = "gene_symbol") |>
    mutate(pct_sig = n_sig / n_pairs) |>
    filter(n_sig >= MIN_PAIRS)

  effect_lim <- max(abs(gene_effect$mean_effect), na.rm = TRUE) * 1.15

  make_effect_panel <- function(dat, fill_col, nudge_dir = 1) {
    label_dat <- dat |>
      mutate(label  = paste0(round(pct_sig * 100, 0), "% sig  (n=", n_sig, ")"),
             label_x = mean_effect + nudge_dir * 0.003)
    ggplot(dat, aes(x = gene_symbol, y = mean_effect)) +
      geom_col(fill = fill_col, alpha = 0.85, width = 0.7) +
      geom_text(data = label_dat,
                aes(x = gene_symbol, y = label_x, label = label),
                hjust = ifelse(nudge_dir > 0, -0.05, 1.05),
                size = 2.5, inherit.aes = FALSE) +
      coord_flip(clip = "off") +
      scale_y_continuous(
        limits = c(ifelse(nudge_dir > 0, 0, -effect_lim),
                   ifelse(nudge_dir > 0,  effect_lim, 0)),
        labels = function(x) sprintf("%+.2f", x),
        expand = expansion(mult = c(0.02, 0.28))) +
      labs(x = NULL, y = "Mean effect size (fraction downstream \u2212 0.5)") +
      theme_thesis() +
      theme(axis.text.y = element_text(size = 8.5))
  }

  top_down <- gene_effect |> arrange(desc(mean_effect)) |> slice_head(n = 20) |>
    mutate(gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)))
  top_up   <- gene_effect |> arrange(mean_effect)        |> slice_head(n = 20) |>
    mutate(gene_symbol = factor(gene_symbol, levels = gene_symbol))

  p5d <- (make_effect_panel(top_down, col_downstream,  1) |
          make_effect_panel(top_up,   col_upstream,   -1)) +
    plot_annotation(
      title    = "Genes with strongest directional splicing order bias",
      subtitle = paste0("Top 20 downstream-biased (left) and upstream-biased (right)",
                        " | Mean effect across significant pairs | \u2265", MIN_PAIRS, " sig pairs"),
      theme    = theme_thesis() + theme(plot.title    = element_text(size = 13, face = "bold"),
                                        plot.subtitle = element_text(size = 9, color = "grey30")))
  save_fig(p5d, "05d_gene_effect_size_ranked", w = 16, h = 8)
}

# --------------------------------------------------------------------------- #
# FIG 05e: pct_sig ranked by |mean_effect|
# --------------------------------------------------------------------------- #

message("--- Fig 05e: Gene effect ranked normalized ---")
{
  gene_effect2 <- df |> filter(sig) |>
    group_by(gene_symbol) |>
    summarise(mean_effect = mean(fraction_downstream - 0.5),
              n_sig = n(),
              n_downstream_biased = sum(direction == "downstream_biased"),
              n_upstream_biased   = sum(direction == "upstream_biased"),
              .groups = "drop") |>
    left_join(gene_summary |> select(gene_symbol, n_pairs), by = "gene_symbol") |>
    mutate(pct_sig = n_sig / n_pairs) |>
    filter(n_pairs >= MIN_PAIRS)

  top_eff <- gene_effect2 |>
    arrange(desc(abs(mean_effect))) |> slice_head(n = 40) |>
    arrange(mean_effect) |>
    mutate(gene_symbol    = factor(gene_symbol, levels = gene_symbol),
           pct_downstream = n_downstream_biased / n_pairs,
           pct_upstream   = n_upstream_biased   / n_pairs) |>
    pivot_longer(c(pct_downstream, pct_upstream),
                 names_to = "direction", values_to = "pct") |>
    mutate(direction = recode(direction, pct_downstream = "Downstream-biased",
                              pct_upstream = "Upstream-biased"))

  label_eff <- gene_effect2 |> arrange(desc(abs(mean_effect))) |> slice_head(n = 40) |>
    arrange(mean_effect) |>
    mutate(gene_symbol = factor(gene_symbol, levels = gene_symbol),
           pct_sig2 = n_sig / n_pairs,
           label    = paste0(round(pct_sig2 * 100, 1), "%  (n=", n_sig, "/", n_pairs, ")"))

  p <- ggplot(top_eff, aes(x = gene_symbol, y = pct, fill = direction)) +
    geom_col() +
    geom_text(data = label_eff,
              aes(x = gene_symbol, y = pct_sig2, label = label),
              hjust = -0.08, size = 2.6, inherit.aes = FALSE) +
    coord_flip(clip = "off") +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0.25))) +
    scale_fill_manual(values = c("Downstream-biased" = col_downstream,
                                 "Upstream-biased"   = col_upstream)) +
    labs(x = NULL, fill = NULL,
         y = "Fraction of intron pairs with significant splicing order",
         title = "Genes with strongest bias, normalized for intron count",
         subtitle = paste0("Top 40 by |mean effect size| | \u2265", MIN_PAIRS, " testable pairs")) +
    theme_thesis() + theme(axis.text.y = element_text(size = 8))
  save_fig(p, "05e_gene_effect_ranked_norm", w = 11, h = 9)
}

# --------------------------------------------------------------------------- #
# FIG 05f: n_sig ranked by mean |effect|
# --------------------------------------------------------------------------- #

message("--- Fig 05f: Gene sig pairs by effect ---")
{
  gene_eff_abs <- df |> filter(sig) |>
    group_by(gene_symbol) |>
    summarise(mean_effect_abs = mean(abs(fraction_downstream - 0.5)), .groups = "drop")

  top_effect_long <- gene_summary |>
    filter(n_pairs >= 2) |>
    left_join(gene_eff_abs, by = "gene_symbol") |>
    filter(!is.na(mean_effect_abs)) |>
    arrange(desc(mean_effect_abs)) |>
    slice_head(n = 40) |>
    mutate(gene_symbol = factor(gene_symbol, levels = rev(gene_symbol))) |>
    pivot_longer(c(n_downstream_biased, n_upstream_biased),
                 names_to = "direction", values_to = "count") |>
    mutate(direction = recode(direction,
                              n_downstream_biased = "Downstream-biased",
                              n_upstream_biased   = "Upstream-biased"))
  p <- ggplot(top_effect_long, aes(x = gene_symbol, y = count, fill = direction)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c("Downstream-biased" = col_downstream,
                                 "Upstream-biased"   = col_upstream)) +
    labs(x = NULL, y = "Number of significant intron pairs", fill = NULL,
         title = "Genes with highest mean effect size among significant pairs",
         subtitle = "Top 40 ranked by mean |fraction_downstream - 0.5|") +
    theme_thesis() + theme(axis.text.y = element_text(size = 8))
  save_fig(p, "05f_gene_sig_pairs_by_effect", w = 10, h = 9)
}

# --------------------------------------------------------------------------- #
# FIG 06a: Intron length density — all vs significant
# --------------------------------------------------------------------------- #

message("--- Fig 06a: Intron length density ---")
{
  len_df <- bind_rows(
    df |> mutate(intron = "upstream (intron1)",   length = intron1_length, subset = "All pairs"),
    df |> mutate(intron = "downstream (intron2)", length = intron2_length, subset = "All pairs"),
    df |> filter(sig) |>
      mutate(intron = "upstream (intron1)",   length = intron1_length, subset = "Significant pairs"),
    df |> filter(sig) |>
      mutate(intron = "downstream (intron2)", length = intron2_length, subset = "Significant pairs")
  ) |> filter(!is.na(length), length > 0)

  len_medians <- len_df |>
    group_by(subset, intron) |>
    summarise(med = median(length), .groups = "drop")

  p <- ggplot(len_df, aes(x = length, fill = subset, color = subset)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(data = len_medians, aes(xintercept = med, color = subset),
               linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~intron, ncol = 1, scales = "free_y") +
    scale_x_log10(labels = label_comma()) +
    scale_fill_manual(values  = c("All pairs" = "steelblue",
                                  "Significant pairs" = col_downstream)) +
    scale_color_manual(values = c("All pairs" = "steelblue",
                                  "Significant pairs" = col_downstream)) +
    labs(x = "Intron length (bp, log\u2081\u2080 scale)", y = "Density",
         fill = NULL, color = NULL,
         title = "Intron length: all pairs vs significant pairs") +
    theme_thesis()
  save_fig(p, "06a_intron_length_density", w = 10, h = 8)
}

# --------------------------------------------------------------------------- #
# FIG 06b: Intron length violin — all vs significant
# --------------------------------------------------------------------------- #

message("--- Fig 06b: Intron length violin ---")
{
  len_box <- bind_rows(
    df |> pivot_longer(c(intron1_length, intron2_length),
                       names_to = "intron", values_to = "length") |>
      mutate(subset = "All pairs"),
    df |> filter(sig) |>
      pivot_longer(c(intron1_length, intron2_length),
                   names_to = "intron", values_to = "length") |>
      mutate(subset = "Significant pairs")
  ) |>
    mutate(intron = recode(intron,
                           intron1_length = "Upstream intron",
                           intron2_length = "Downstream intron")) |>
    filter(!is.na(length), length > 0)

  p <- ggplot(len_box, aes(x = subset, y = length, fill = subset)) +
    geom_violin(trim = TRUE, alpha = 0.5) +
    geom_boxplot(width = 0.15, outlier.size = 0.5, outlier.alpha = 0.3,
                 fill = "white", linewidth = 0.5) +
    facet_wrap(~intron) +
    scale_y_log10(labels = label_comma()) +
    scale_fill_manual(values = c("All pairs" = "steelblue",
                                 "Significant pairs" = col_downstream)) +
    labs(x = NULL, y = "Intron length (bp)", fill = NULL,
         title = "Intron length: all vs significant pairs") +
    theme_thesis() + theme(legend.position = "none")
  save_fig(p, "06b_intron_length_violin", w = 8, h = 6)
}

# --------------------------------------------------------------------------- #
# FIG 06c: Paired length hexbin + patchwork marginals
# --------------------------------------------------------------------------- #

message("--- Fig 06c: Paired length hexbin ---")
{
  pair_len <- df |>
    filter(!is.na(intron1_length), !is.na(intron2_length),
           intron1_length > 0, intron2_length > 0)
  len_breaks <- 10^(1:5)
  x_lim <- range(pair_len$intron1_length) * c(0.9, 1.1)
  y_lim <- range(pair_len$intron2_length) * c(0.9, 1.1)

  p_top <- ggplot(pair_len, aes(x = intron1_length)) +
    geom_histogram(bins = 60, fill = col_upstream, color = "white", linewidth = 0.15) +
    scale_x_log10(limits = x_lim, breaks = len_breaks) +
    theme_thesis() +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank(), plot.margin = margin(5, 5, 0, 5))

  p_right <- ggplot(pair_len, aes(x = intron2_length)) +
    geom_histogram(bins = 60, fill = col_downstream, color = "white", linewidth = 0.15) +
    scale_x_log10(limits = y_lim, breaks = len_breaks) +
    coord_flip() +
    theme_thesis() +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), plot.margin = margin(5, 5, 5, 0))

  p_main <- ggplot(pair_len, aes(x = intron1_length, y = intron2_length)) +
    geom_hex(bins = 60, aes(fill = after_stat(log10(count)))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "white", linewidth = 0.8) +
    scale_x_log10(labels = scales::comma, breaks = len_breaks, limits = x_lim,
                  name = "Upstream intron length (bp)") +
    scale_y_log10(labels = scales::comma, breaks = len_breaks, limits = y_lim,
                  name = "Downstream intron length (bp)") +
    scale_fill_viridis_c(name = expression(log[10](pairs)), option = "plasma") +
    labs(title = "Upstream vs downstream intron length",
         subtitle = "Dashed = equal length | Hex bins, log\u2081\u2080 axes") +
    theme_thesis() + theme(plot.margin = margin(0, 5, 5, 5))

  p_hex <- (p_top + plot_spacer() + plot_layout(widths = c(4, 1))) /
           (p_main + p_right  + plot_layout(widths = c(4, 1))) +
           plot_layout(heights = c(1, 4))

  pdf(file.path(opt$outdir, "06c_paired_length_hexbin.pdf"), width = 9, height = 8)
  print(p_hex); dev.off()
  png(file.path(opt$outdir, "06c_paired_length_hexbin.png"),
      width = 9, height = 8, units = "in", res = 300)
  print(p_hex); dev.off()
  message("  Saved: 06c_paired_length_hexbin")
}

# --------------------------------------------------------------------------- #
# FIG 06d: Paired length scatter — significant pairs coloured by fraction_downstream
# --------------------------------------------------------------------------- #

message("--- Fig 06d: Paired length scatter (sig) ---")
{
  p <- ggplot(df |> filter(sig, intron1_length > 0, intron2_length > 0),
              aes(x = intron1_length, y = intron2_length,
                  color = fraction_downstream, size = log10(total + 1))) +
    geom_point(alpha = 0.6) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey40", linewidth = 0.7) +
    scale_x_log10(labels = label_comma(), name = "Upstream intron length (bp)") +
    scale_y_log10(labels = label_comma(), name = "Downstream intron length (bp)") +
    scale_color_gradient2(low = col_upstream, mid = "grey85", high = col_downstream,
                          midpoint = 0.5,
                          name = "Fraction\ndownstream\nspliced first") +
    scale_size_continuous(name = "log\u2081\u2080(reads)", range = c(1, 4)) +
    labs(title = "Paired intron lengths — significant pairs",
         subtitle = paste0("n = ", scales::comma(n_sig),
                           " | Dashed = equal length")) +
    theme_thesis()
  save_fig(p, "06d_paired_length_scatter_sig", w = 9, h = 7)
}

# --------------------------------------------------------------------------- #
# FIG 06e: log2 length ratio by direction
# --------------------------------------------------------------------------- #

message("--- Fig 06e: Length ratio by direction ---")
{
  pair_len2 <- df |>
    filter(intron1_length > 0, intron2_length > 0) |>
    mutate(log2_length_ratio = log2(intron2_length / intron1_length))

  dir_labels_n <- c(
    downstream_biased = paste0("Downstream-biased (n=", scales::comma(n_down), ")"),
    not_significant   = paste0("Not significant (n=",
                               scales::comma(sum(df$direction == "not_significant")), ")"),
    upstream_biased   = paste0("Upstream-biased (n=", scales::comma(n_up), ")")
  )
  ratio_medians <- pair_len2 |>
    group_by(direction) |>
    summarise(med = median(log2_length_ratio), .groups = "drop")

  p <- ggplot(pair_len2, aes(x = log2_length_ratio,
                              fill = direction, color = direction)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
    geom_vline(data = ratio_medians, aes(xintercept = med, color = direction),
               linetype = "dashed", linewidth = 0.9) +
    scale_fill_manual(values = direction_colors,  labels = dir_labels_n) +
    scale_color_manual(values = direction_colors, labels = dir_labels_n) +
    labs(x = "log\u2082(intron2 length / intron1 length)", y = "Density",
         fill = NULL, color = NULL,
         title = "log\u2082(downstream / upstream) intron length ratio by direction",
         subtitle = "Positive = downstream intron longer | 0 = equal length") +
    theme_thesis()
  save_fig(p, "06e_length_ratio_by_direction", w = 10, h = 6)
}

# --------------------------------------------------------------------------- #
# FIG 07: Read coverage histogram
# --------------------------------------------------------------------------- #

message("--- Fig 07: Read coverage ---")
{
  coverage_breaks <- 10^seq(0, 7, by = 0.15)
  cov_df <- bind_rows(df |> mutate(subset = "All pairs"),
                      df |> filter(sig) |> mutate(subset = "Significant pairs"))
  cov_medians <- cov_df |>
    group_by(subset) |>
    summarise(med = median(total), n = n(), .groups = "drop")
  p <- ggplot(cov_df, aes(x = total, fill = subset, color = subset)) +
    geom_histogram(breaks = coverage_breaks, alpha = 0.55,
                   position = "identity", linewidth = 0.15) +
    geom_vline(data = cov_medians, aes(xintercept = med, color = subset),
               linetype = "dashed", linewidth = 0.8) +
    geom_text(data = cov_medians,
              aes(x = med, y = Inf,
                  label = paste0("n=", scales::comma(n), "\nmed=", round(med)),
                  color = subset),
              angle = 90, vjust = 1.5, hjust = 1.1, size = 3.2, inherit.aes = FALSE) +
    scale_x_log10(labels = scales::comma,
                  name = "Total informative reads (log\u2081\u2080 scale)") +
    scale_y_continuous(labels = scales::comma) +
    scale_fill_manual(values  = c("All pairs" = col_upstream,
                                  "Significant pairs" = col_downstream), name = NULL) +
    scale_color_manual(values = c("All pairs" = col_upstream,
                                  "Significant pairs" = col_downstream), name = NULL) +
    labs(y = "Intron pairs",
         title = "Read coverage: all pairs vs significant pairs",
         subtitle = "Higher coverage enriched in significant pairs") +
    theme_thesis()
  save_fig(p, "07_read_coverage", w = 10, h = 6)
}

# --------------------------------------------------------------------------- #
# FIG 08-10: GO enrichment
# --------------------------------------------------------------------------- #

message("--- Fig 08-10: GO enrichment ---")

go_ok <- requireNamespace("clusterProfiler", quietly = TRUE) &&
         requireNamespace("org.Hs.eg.db",    quietly = TRUE)

if (!go_ok) {
  message("  Skipping GO — install with BiocManager::install(",
          "c('clusterProfiler','org.Hs.eg.db'))")
} else {
  suppressPackageStartupMessages({
    library(clusterProfiler); library(org.Hs.eg.db)
  })
  select <- dplyr::select

  all_genes <- unique(df$gene_symbol)

  map_to_entrez <- function(keys) {
    id_type <- dplyr::case_when(
      grepl("^ENSG",    keys[1]) ~ "ENSEMBL",
      grepl("^[0-9]+$", keys[1]) ~ "ENTREZID",
      TRUE ~ "SYMBOL"
    )
    if (id_type == "ENTREZID") return(keys)
    tryCatch(suppressMessages(as.character(na.omit(unique(
      AnnotationDbi::mapIds(ORGDB, keys = keys,
                            column = "ENTREZID", keytype = id_type,
                            multiVals = "first"))))),
      error = function(e) { message("  mapIds failed: ", e$message); character(0) })
  }

  run_go_df <- function(target_genes, background, direction_label) {
    if (length(target_genes) < 5) return(NULL)
    use_genes <- map_to_entrez(target_genes)
    use_bg    <- map_to_entrez(background)
    if (length(use_genes) < 5) return(NULL)
    message("  GO: ", direction_label, " | ", length(use_genes), " genes")
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = use_genes, universe = use_bg, OrgDb = ORGDB,
        keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
        pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) { message("  enrichGO failed: ", e$message); NULL })
    if (is.null(ego)) return(NULL)
    ego_df <- as.data.frame(ego)
    if (nrow(ego_df) == 0) return(NULL)
    ego_df |>
      dplyr::mutate(
        gene_ratio_num = vapply(strsplit(GeneRatio, "/"),
                                function(x) as.numeric(x[1]) / as.numeric(x[2]),
                                numeric(1)),
        neg_log10_padj = -log10(p.adjust),
        direction      = direction_label,
        Description    = stringr::str_trunc(Description, 45)
      )
  }

  plot_go_single <- function(ego_df, direction_label, fill_color, stem) {
    if (is.null(ego_df) || nrow(ego_df) == 0) return(NULL)
    top <- ego_df |> dplyr::arrange(p.adjust) |> dplyr::slice_head(n = 25) |>
      dplyr::mutate(Description = factor(Description, levels = rev(Description)))
    top <- top |>
      dplyr::distinct(Description, .keep_all = TRUE) |>
      dplyr::mutate(Description = factor(Description, levels = rev(Description)))
    p <- ggplot(top, aes(x = gene_ratio_num, y = Description,
                          size = Count, color = neg_log10_padj)) +
      geom_point() +
      scale_color_gradient(low = "#7B9EC9", high = fill_color,
                           name = expression(-log[10](p[adj]))) +
      scale_size_continuous(name = "Gene count", range = c(2, 8)) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(title = paste0("GO Biological Process\n(", direction_label, " pairs)"),
           subtitle = paste0("Target: ", length(target_genes),
                             " genes | Background: ", length(all_genes)),
           x = "Gene ratio", y = NULL) +
      theme_thesis() + theme(axis.text.y = element_text(size = 9))
    save_fig(p, stem, w = 11, h = 9)
    write.table(ego_df,
                file.path(opt$outdir, paste0(stem, "_table.tsv")),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }

  downstream_genes <- df |>
    filter(direction == "downstream_biased", abs_effect >= opt$`effect-threshold`) |>
    pull(gene_symbol) |> unique()
  upstream_genes <- df |>
    filter(direction == "upstream_biased", abs_effect >= opt$`effect-threshold`) |>
    pull(gene_symbol) |> unique()

  target_genes <- downstream_genes
  go_down <- run_go_df(downstream_genes, all_genes, "Downstream-biased")
  plot_go_single(go_down, "Downstream-biased", col_downstream, "08_GO_enrichment_downstream")

  target_genes <- upstream_genes
  go_up <- run_go_df(upstream_genes, all_genes, "Upstream-biased")
  plot_go_single(go_up, "Upstream-biased", col_upstream, "09_GO_enrichment_upstream")

  # Combined faceted figure
  if (!is.null(go_down) || !is.null(go_up)) {
    go_combined <- bind_rows(
      if (!is.null(go_down)) go_down |> dplyr::arrange(p.adjust) |>
        dplyr::slice_head(n = 20) |> dplyr::mutate(direction_label = "Downstream-biased"),
      if (!is.null(go_up))   go_up   |> dplyr::arrange(p.adjust) |>
        dplyr::slice_head(n = 20) |> dplyr::mutate(direction_label = "Upstream-biased")
    )
    if (nrow(go_combined) > 0) {
      go_combined <- go_combined |>
        dplyr::group_by(direction_label) |>
        dplyr::arrange(gene_ratio_num, .by_group = TRUE) |>
        dplyr::mutate(term_ordered = factor(paste0(direction_label, "__", Description),
                                            levels = unique(paste0(direction_label,
                                                                   "__", Description)))) |>
        dplyr::ungroup()
      max_sig <- max(go_combined$neg_log10_padj, na.rm = TRUE)
      p_comb <- ggplot(go_combined,
                       aes(x = gene_ratio_num, y = term_ordered,
                           size = Count, color = neg_log10_padj)) +
        geom_point() +
        scale_y_discrete(labels = function(x) sub("^.*?__", "", x)) +
        facet_wrap(~ direction_label, scales = "free_y", ncol = 2) +
        scale_color_gradient2(low = "#7B9EC9", mid = "#AA5588", high = "#C0392B",
                              midpoint = max_sig / 2, limits = c(0, max_sig),
                              name = expression(-log[10](p[adj]))) +
        scale_size_continuous(name = "Gene count", range = c(2, 8)) +
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
        labs(title = "GO Biological Process enrichment",
             subtitle = paste0("Top 20 terms per direction | Background: ",
                               scales::comma(length(all_genes)), " genes"),
             x = "Gene ratio", y = NULL) +
        theme_thesis() +
        theme(axis.text.y = element_text(size = 8),
              strip.text  = element_text(face = "bold", size = 11),
              panel.spacing = unit(1.5, "lines"))
      save_fig(p_comb, "10_GO_enrichment_combined", w = 16, h = 10)
    }
  }
}

# --------------------------------------------------------------------------- #
# Write significant pairs TSV
# --------------------------------------------------------------------------- #

sig_pairs <- df |>
  filter(sig) |>
  arrange(desc(abs_effect), padj)
write.table(sig_pairs,
            file.path(opt$outdir, "significant_pairs.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("  significant_pairs.tsv written: ", nrow(sig_pairs), " rows")

# --------------------------------------------------------------------------- #
# Summary
# --------------------------------------------------------------------------- #

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY\n")
cat(strrep("=", 60), "\n")
cat(sprintf("  Total intron pairs:           %s\n", scales::comma(nrow(df))))
cat(sprintf("  Significant (FDR<%.2f, |e|>=%.2f): %s (%.1f%%)\n",
            opt$fdr, opt$`effect-threshold`,
            scales::comma(n_sig), 100 * n_sig / nrow(df)))
cat(sprintf("    Downstream-biased:          %s\n", scales::comma(n_down)))
cat(sprintf("    Upstream-biased:            %s\n", scales::comma(n_up)))
cat(sprintf("  Unique genes (all):           %s\n",
            scales::comma(length(unique(df$gene_symbol)))))
cat(sprintf("  Unique genes (significant):   %s\n",
            scales::comma(length(unique(sig_pairs$gene_symbol)))))
cat(sprintf("  Figures saved to:             %s/\n", opt$outdir))
cat(strrep("=", 60), "\n")