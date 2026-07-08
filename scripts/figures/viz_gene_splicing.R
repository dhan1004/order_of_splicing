#!/usr/bin/env Rscript
# =============================================================================
# plot_gene_splicing_order.R
#
# Draws a gene's exon-intron structure and annotates it with pairwise
# splicing order arcs from your actual data files.
#
# Input: significant_pairs.tsv  (pre-filtered, has gene_symbol,
#   direction, binom_padj columns; produced by viz_order_enrichment*.R)
#
# Gene structure is reconstructed directly from the intron coordinates in
# either file — no separate annotation file is needed.
#
# Usage (command line):
#   Rscript plot_gene_splicing_order.R \
#       --gene     NBPF10 \
#       --input    /path/to/significant_pairs.tsv \
#       --outdir   ./figures
#
# Columns used from significant_pairs.tsv:
#   chr, gene_symbol, intron1_start, intron1_end, intron2_start, intron2_end,
#   direction   ("downstream_biased" | "upstream_biased" | "not_significant")
#   binom_padj  (used to scale arc line weight by significance)
#   fraction_downstream, upstream, downstream, total
#
# Columns used from splicing_order_pooled.tsv:
#   chr, gene_id, intron1_start, intron1_end, intron2_start, intron2_end,
#   upstream, downstream, total, fraction_downstream
#   (gene_symbol added via GENE_BED or org.Hs.eg.db, same as your other scripts)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggforce)   # geom_bezier — install: install.packages("ggforce")
  library(scales)
})

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

GENE_NAME  <- "FUS"   # ← change to any gene_symbol in your data

# Path to either significant_pairs.tsv, change
INPUT_FILE <- "/users/dhan30/scratch/splicing_order/figures/0407/significant_pairs.tsv"

# Output directory (created if needed)
OUT_DIR <- "figures/gene_structures"

# For pooled TSV only: significance thresholds (ignored if sig pairs file supplied)
FDR_THRESH    <- 0.05
EFFECT_THRESH <- 0.25   # |fraction_downstream - 0.5|
MIN_READS     <- 10

# Gene BED for symbol mapping (needed only with pooled TSV if gene_id is ENST/ENSG).
GENE_BED <- "/users/dhan30/reference/gencode_genes.bed"

# Command-line override
args <- commandArgs(trailingOnly = TRUE)
for (i in seq_along(args)) {
  if (args[i] == "--gene"   && i < length(args)) GENE_NAME  <- args[i + 1]
  if (args[i] == "--input"  && i < length(args)) INPUT_FILE <- args[i + 1]
  if (args[i] == "--outdir" && i < length(args)) OUT_DIR    <- args[i + 1]
}

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: LOAD + NORMALISE DATA
# ─────────────────────────────────────────────────────────────────────────────

message("Reading: ", INPUT_FILE)
raw <- read.delim(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)
message("  Rows: ", nrow(raw), "  Cols: ", paste(colnames(raw), collapse = ", "))

# ── Normalize column names ────────────────────────────────────────────────────
rename_col <- function(df, old, new) {
  if (old %in% names(df) && !new %in% names(df)) names(df)[names(df) == old] <- new
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

# ── Add gene_symbol if missing (pooled TSV path) ─────────────────────────────
if (!"gene_symbol" %in% names(raw)) {
  raw$gene_id <- sub("\\.\\d+$", "", raw$gene_id)

  gene_symbols <- NULL

  # coordinate-based BED lookup
  if (!is.null(GENE_BED) && file.exists(GENE_BED)) {
    message("  Mapping gene symbols via BED: ", GENE_BED)
    genes <- tryCatch({
      g <- read.table(GENE_BED, header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "", comment.char = "#")
      colnames(g)[1:4] <- c("chr", "start", "end", "gene_symbol")
      g$start <- as.integer(g$start); g$end <- as.integer(g$end)
      g[!is.na(g$start), ]
    }, error = function(e) { message("  BED load failed: ", e$message); NULL })

    if (!is.null(genes)) {
      genes_by_chr <- split(genes, genes$chr)
      gene_symbols <- mapply(function(ch, ps, pe) {
        g <- genes_by_chr[[ch]]
        if (is.null(g)) return(NA_character_)
        hits <- g[g$start < pe & g$end > ps, "gene_symbol", drop = TRUE]
        if (length(hits) > 0) hits[1] else NA_character_
      }, raw$chr, raw$intron1_start, raw$intron2_end,
      SIMPLIFY = TRUE, USE.NAMES = FALSE)
      message("  Mapped: ", sum(!is.na(gene_symbols)), "/", nrow(raw))
    }
  }

  raw$gene_symbol <- if (!is.null(gene_symbols))
    ifelse(is.na(gene_symbols) | gene_symbols == "", raw$gene_id, gene_symbols)
  else
    raw$gene_id
}

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: SUBSET TO TARGET GENE
# ─────────────────────────────────────────────────────────────────────────────

# significant_pairs.tsv already contains only sig pairs — use all of them
df_gene_sig <- raw %>%
  filter(gene_symbol == GENE_NAME,
         direction %in% c("downstream_biased", "upstream_biased"))

# All pairs for this gene (to reconstruct full intron set)
df_gene_all <- raw %>% filter(gene_symbol == GENE_NAME)

if (nrow(df_gene_all) == 0)
  stop("Gene '", GENE_NAME, "' not found in data.\n",
       "Sample gene symbols: ",
       paste(head(unique(raw$gene_symbol), 12), collapse = ", "))

message("\nGene: ", GENE_NAME)
message("  Total pairs in gene:           ", nrow(df_gene_all))
message("  Significant pairs to annotate: ", nrow(df_gene_sig))
if (nrow(df_gene_sig) == 0)
  message("  (No significant pairs found — plotting gene structure only)")

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: RECONSTRUCT GENE MODEL FROM COORDINATES
# ─────────────────────────────────────────────────────────────────────────────
# Each TSV row is an intron pair. Collect all unique individual introns,
# then infer exon positions as the gaps between consecutive introns.

introns_tbl <- df_gene_all %>%
  select(intron1_start, intron1_end, intron2_start, intron2_end) %>%
  pivot_longer(everything(),
               names_to  = c("which", ".value"),
               names_pattern = "(intron[12])_(start|end)") %>%
  distinct(start, end) %>%
  arrange(start) %>%
  mutate(
    number = row_number(),
    mid    = (start + end) / 2,
    length = end - start
  )

message("  Unique introns reconstructed: ", nrow(introns_tbl))

# Terminal exon stubs (display-only flanking regions)
EXON_STUB <- 200

exons_internal <- if (nrow(introns_tbl) >= 2) {
  tibble(
    start  = introns_tbl$end[-nrow(introns_tbl)],
    end    = introns_tbl$start[-1],
    number = seq_len(nrow(introns_tbl) - 1) + 1L
  ) %>% filter(end > start)
} else {
  tibble(start = numeric(), end = numeric(), number = integer())
}

exons_all <- bind_rows(
  tibble(start = min(introns_tbl$start) - EXON_STUB,
         end   = min(introns_tbl$start),
         number = 1L),
  exons_internal,
  tibble(start = max(introns_tbl$end),
         end   = max(introns_tbl$end) + EXON_STUB,
         number = nrow(introns_tbl) + 1L)
) %>%
  arrange(start) %>%
  mutate(
    mid = (start + end) / 2,
    # Constitutive = internal exon flanked on both sides by observed introns.
    # Terminal stubs (first and last exon) are treated as non-constitutive.
    constitutive = number > 1 & number < max(number)
  )

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: BUILD ARC DATA
# ─────────────────────────────────────────────────────────────────────────────
# Arrow direction encodes which intron is spliced SECOND (i.e., later):
#   downstream_biased = intron1 spliced first  = arrow points right (to intron2)
#   upstream_biased   = intron2 spliced first  = arrow points left  (to intron1)
# Arcs alternate above/below to reduce clutter; span > 1 gets extra height.

COL_DOWNSTREAM <- "#E05C5C"
COL_UPSTREAM   <- "#4B8BBE"
COL_EXON       <- "#7a5195"

build_arcs <- function(pairs_df, introns) {
  if (nrow(pairs_df) == 0) return(NULL)

  # Map each pair to intron indices
  idx_map <- introns %>% select(start, end, number)

  pairs_indexed <- pairs_df %>%
    left_join(idx_map %>% rename(i1_num = number),
              by = c("intron1_start" = "start", "intron1_end" = "end")) %>%
    left_join(idx_map %>% rename(i2_num = number),
              by = c("intron2_start" = "start", "intron2_end" = "end")) %>%
    filter(!is.na(i1_num), !is.na(i2_num)) %>%
    mutate(span = abs(i2_num - i1_num)) %>%
    arrange(span, i1_num)

  if (nrow(pairs_indexed) == 0) {
    warning("No pairs could be matched to intron indices. ",
            "Check that intron coordinates in the pairs table match those in ",
            "the full data used to reconstruct intron positions.")
    return(NULL)
  }

  arc_list <- vector("list", nrow(pairs_indexed))

  for (k in seq_len(nrow(pairs_indexed))) {
    row   <- pairs_indexed[k, ]
    i_mid <- introns$mid[introns$number == row$i1_num]
    j_mid <- introns$mid[introns$number == row$i2_num]
    if (length(i_mid) == 0 || length(j_mid) == 0) next

    above  <- (k %% 2 == 1)
    dir_y  <- if (above) 1 else -1
    height <- dir_y * (0.40 + 0.18 * (row$span - 1))

    # Bezier from "spliced first" intron → "spliced second" intron
    x_from <- if (row$direction == "downstream_biased") i_mid else j_mid
    x_to   <- if (row$direction == "downstream_biased") j_mid else i_mid

    padj_val <- if ("binom_padj" %in% names(row) && !is.na(row$binom_padj))
      row$binom_padj else NA_real_

    arc_list[[k]] <- tibble(
      x           = c(x_from, (x_from + x_to) / 2, x_to),
      y           = c(0,       height,               0),
      group       = k,
      direction   = row$direction,
      span        = row$span,
      above       = above,
      x_tip       = x_to,
      padj        = padj_val,
      fd          = row$fraction_downstream,
      total_reads = row$total_reads,
      # Adjacent pairs (span == 1) get high alpha; long-range get low alpha
      constitutive = row$span == 1
    )
  }
  bind_rows(arc_list)
}

arcs <- build_arcs(df_gene_sig, introns_tbl)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: TABLE OUTPUT
# ─────────────────────────────────────────────────────────────────────────────
# One row per significant intron pair for the target gene, sorted by intron
# index. Includes intron indices, coordinates, lengths, read counts,
# fraction_downstream, effect size, FDR, and a plain-English splicing_order
# column ("I1 first" / "I2 first") for quick reading.

idx_map_for_table <- introns_tbl %>% select(start, end, number)

pair_table <- df_gene_sig %>%
  left_join(idx_map_for_table %>% rename(intron_i = number),
            by = c("intron1_start" = "start", "intron1_end" = "end")) %>%
  left_join(idx_map_for_table %>% rename(intron_j = number),
            by = c("intron2_start" = "start", "intron2_end" = "end")) %>%
  mutate(
    coords        = paste0(chr, ":", format(intron1_start, big.mark = ","),
                           "-", format(intron2_end,   big.mark = ",")),
    effect        = fraction_downstream - 0.5,
    # "spliced first" = the intron with the larger fraction
    splicing_order = case_when(
      direction == "downstream_biased" ~ paste0("I", intron_i, " first"),
      direction == "upstream_biased"   ~ paste0("I", intron_j, " first"),
      TRUE ~ "no bias"
    ),
    direction_label = recode(direction,
      downstream_biased = "Downstream-biased",
      upstream_biased   = "Upstream-biased"
    )
  ) %>%
  arrange(intron_i, intron_j) %>%
  select(
    intron_i, intron_j,
    coords,
    intron1_length, intron2_length,
    upstream_count, downstream_count, total_reads,
    fraction_downstream, effect,
    any_of("binom_padj"),
    splicing_order, direction_label
  )

tsv_out <- file.path(OUT_DIR, paste0(GENE_NAME, "_pair_table.tsv"))
write.table(pair_table, tsv_out, sep = "\t", row.names = FALSE, quote = FALSE)
message("Pair table written: ", tsv_out, "  (", nrow(pair_table), " rows)")

# Also print to console for interactive use
cat("\n── Intron pair splicing order table: ", GENE_NAME, " ──\n", sep = "")
pair_table_print <- pair_table %>%
  mutate(across(where(is.numeric), ~ signif(., 4)),
         across(everything(), as.character))
print(as.data.frame(pair_table_print), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: PLOT
# ─────────────────────────────────────────────────────────────────────────────

x_range    <- range(c(exons_all$start, exons_all$end))
x_pad      <- diff(x_range) * 0.03
x_limits   <- c(x_range[1] - x_pad, x_range[2] + x_pad)
arc_extent <- if (!is.null(arcs) && nrow(arcs) > 0) max(abs(arcs$y)) + 0.12 else 0.55
y_limits   <- c(-arc_extent - 0.10, arc_extent + 0.10)

# Layer order (bottom → top):
#   1. backbone + strand arrow
#   2. arcs (Bézier curves)
#   3. arrowhead points on arcs
#   4. intron midpoint ticks, gene structure drawn ON TOP of arcs
#   5. exon boxes
#   6. exon labels (below axis)
#   7. intron labels (above backbone)

p <- ggplot() +

  # ── 1. Backbone ────────────────────────────────────────────────────────────
  annotate("segment",
           x = min(exons_all$start), xend = max(exons_all$end),
           y = 0, yend = 0,
           colour = "black", linewidth = 0.4) +

  # Strand arrow (3′ end)
  annotate("segment",
           x    = max(exons_all$end) - diff(x_range) * 0.015,
           xend = max(exons_all$end),
           y = 0, yend = 0,
           colour = "grey40", linewidth = 0.9,
           arrow = arrow(length = unit(5, "pt"), type = "closed"))

# ── 2 + 3. Arc layers (drawn before gene structure) ──────────────────────────
if (!is.null(arcs) && nrow(arcs) > 0) {
  for (g in unique(arcs$group)) {
    ag    <- filter(arcs, group == g)
    col_g <- if (ag$direction[1] == "downstream_biased") COL_DOWNSTREAM else COL_UPSTREAM

    # Line weight scaled by sqrt(total_reads), clamped to [0.3, 1.5]
    max_reads <- max(arcs$total_reads, na.rm = TRUE)
    lwd_g <- pmax(0.3, pmin(1.5, 0.3 + 1.2 * sqrt(ag$total_reads[1]) / sqrt(max_reads)))

    # Alpha: higher for constitutive pairs, lower for pairs involving terminal introns
    alpha_g <- if (isTRUE(ag$constitutive[1])) 0.4 else 0.4

    p <- p +
      geom_bezier(data = ag,
                  aes(x = x, y = y, group = group),
                  colour = col_g, linewidth = lwd_g, alpha = alpha_g,
                  lineend = "round", show.legend = FALSE) +
      # Arrowhead point at tip — also drawn before exon boxes
      annotate("point",
               x = ag$x_tip[1], y = 0,
               colour = col_g, size = 2.0, shape = 17, alpha = alpha_g)
  }
}

p <- p +
  # ── 4. Intron midpoint ticks ────────────────────────────────────────────────
  # geom_segment(data = introns_tbl,
  #              aes(x = mid, xend = mid, y = -0.07, yend = 0.07),
  #              colour = "grey65", linewidth = 0.5) +

  # ── 5. Exon boxes ───────────────────────────────────────────────────────────
  geom_rect(data = exons_all,
            aes(xmin = start, xmax = end, ymin = -0.28, ymax = 0.28),
            fill = COL_EXON, colour = COL_EXON, linewidth = 0.3) +

  # ── 6. Exon labels (below axis) ─────────────────────────────────────────────
  geom_text(data = exons_all,
            aes(x = mid, y = -0.5, label = number),
            size = 2.6, colour = COL_EXON, fontface = "bold", vjust = 1) +

  # ── 7. Intron labels (just above backbone, on top of everything) ─────────────
  geom_text(data = introns_tbl,
            aes(x = mid, y = 0.5, label = number),
            size = 2.4, colour = "black", fontface = "bold", vjust = 0)

# ── 8. Dummy layer for legend ────────────────────────────────────────────────
leg_df <- tibble(
  x    = c(NA_real_, NA_real_),
  y    = c(NA_real_, NA_real_),
  bias = factor(c("Downstream-biased", "Upstream-biased"))
)

p <- p +
  geom_line(data = leg_df,
            aes(x = x, y = y, colour = bias),
            linewidth = 1.3, na.rm = TRUE) +
  scale_colour_manual(
    values = c("Downstream-biased" = COL_DOWNSTREAM,
               "Upstream-biased"   = COL_UPSTREAM),
    name = NULL, drop = FALSE
  ) +
  coord_cartesian(xlim = x_limits, ylim = y_limits, clip = "off") +
  scale_x_continuous(
    labels = function(x) paste0(round(x / 1e3, 1), " kb"),
    expand = c(0, 0)
  ) +
  labs(
    title    = GENE_NAME,
    x = "Genomic position",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line.y      = element_blank(),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.line.x      = element_line(colour = "grey70", linewidth = 0.4),
    axis.text.x      = element_text(colour = "grey40", size = 8),
    axis.ticks.x     = element_line(colour = "grey70"),
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(colour = "grey40", size = 8.5, hjust = 0,
                                    margin = margin(b = 8)),
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    legend.key.width = unit(26, "pt"),
    plot.margin      = margin(10, 16, 6, 16)
  )

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: SAVE PLOT
# ─────────────────────────────────────────────────────────────────────────────

n_introns  <- nrow(introns_tbl)
fig_width  <- 10
fig_height <- 4.0

pdf_out <- file.path(OUT_DIR, paste0(GENE_NAME, "_splicing_order.pdf"))
png_out <- file.path(OUT_DIR, paste0(GENE_NAME, "_splicing_order.png"))

ggsave(pdf_out, plot = p, width = fig_width, height = fig_height, device = cairo_pdf)
ggsave(png_out, plot = p, width = fig_width, height = fig_height, dpi = 300)

message("\nSaved:\n  ", pdf_out, "\n  ", png_out)
print(p)