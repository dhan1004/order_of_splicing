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
  # Rscript scripts/figures/viz_gene_splicing.R \
  #     --gene     Pol1ra \
  #     --input    results/ortholog/mouse_subset_pairs.tsv \
  #     --outdir   ./figures
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

# Point at the species-appropriate introns BED. Switch this (and GENE_BED,
# if used for symbol mapping) when plotting mouse.
#   human: /users/dhan30/reference/hg38.gencode.basic.v43.introns.bed.gz
#   mouse: /users/dhan30/reference/mm39.gencode.basic.vM36.introns.bed.gz
INTRON_BED <- "/users/dhan30/reference/mm39.gencode.basic.vM36.introns.bed.gz"
 
# The introns BED has no gene symbol. We resolve gene -> transcripts by
# matching the gene's observed pair coordinates against BED introns, then
# pick the longest transcript overlapping them. If your BED name field
# encodes the gene symbol directly, set BED_NAME_HAS_SYMBOL <- TRUE and
# adjust parse_bed_name() below.
BED_NAME_HAS_SYMBOL <- TRUE

GENE_NAME  <- "Polr1a"   # ← change to any gene_symbol in your data

# Path to either significant_pairs.tsv, change
INPUT_FILE <- "results/ortholog/mouse_subset_pairs.tsv"

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

# ── Resolve gene_symbol ──────────────────────────────────────────────────────
# 1. ortholog subset files carry mouse_symbol (or human_symbol); alias it.
if (!"gene_symbol" %in% names(raw) && "mouse_symbol" %in% names(raw))
  raw$gene_symbol <- raw$mouse_symbol
if (!"gene_symbol" %in% names(raw) && "human_symbol" %in% names(raw))
  raw$gene_symbol <- raw$human_symbol

# 2. pooled TSV path: only if we STILL have no symbol AND a gene_id exists.
if (!"gene_symbol" %in% names(raw) && "gene_id" %in% names(raw)) {
  raw$gene_id <- sub("\\.\\d+$", "", raw$gene_id)

  gene_symbols <- NULL
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
  else raw$gene_id
}

# 3. still nothing? fail loudly instead of subsetting a closure later.
if (!"gene_symbol" %in% names(raw))
  stop("No gene_symbol / mouse_symbol / human_symbol / gene_id column in ",
       INPUT_FILE, ".\n  Columns present: ",
       paste(colnames(raw), collapse = ", "))
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

# The backbone now comes from GENCODE. Each BED row is one intron:
#   col1 chr, col2 start(0-based), col3 end, col4 name=<tx>_intron_<n>,
#   col5 score(ignored), col6 strand.
# intron_num in the name is 1-based, transcription-order (your documented
# convention): on minus strand it DEcreASES as genomic coordinate increases.
 
suppressPackageStartupMessages({
  library(readr)
})
 
message("  Loading GENCODE introns BED: ", INTRON_BED)
stopifnot(file.exists(INTRON_BED))
 
bed <- read_tsv(
  INTRON_BED,
  col_names = c("chr", "start", "end", "name", "score", "strand"),
  col_types = cols(chr = "c", start = "i", end = "i",
                   name = "c", score = "c", strand = "c"),
  comment = "#", progress = FALSE
)
 
# name format (UCSC/bedparse style):
#   <transcript>_intron_<num>_<count>_<chr>_<start1based>_<f|r>
# e.g. ENSMUST00000070533.5_intron_1_0_chr1_3287192_r
# The intron number is the field IMMEDIATELY after "_intron_", NOT the tail.
# The transcript base may itself contain no "_intron_", so anchor on the first.
parse_bed_name <- function(nm) {
  m <- regmatches(nm, regexec("^(.*?)_intron_(\\d+)_", nm))
  tx  <- vapply(m, function(x) if (length(x) == 3) x[2] else NA_character_, "")
  num <- vapply(m, function(x) if (length(x) == 3) as.integer(x[3]) else NA_integer_, 0L)
  list(tx = tx, intron_num = num)
}
pn <- parse_bed_name(bed$name)
bed$transcript <- pn$tx
bed$intron_num <- pn$intron_num
bed <- bed[!is.na(bed$transcript), ]

# ── Resolve this gene to a set of candidate transcripts ─────────────────────
# We don't have a symbol column in the BED, so intersect the gene's observed
# pair introns (from df_gene_all) with BED introns to find which transcripts
# this gene lives on, then keep the LONGEST (max genomic span ≈ UCSC default).
gene_chr <- unique(df_gene_all$chr)
if (length(gene_chr) != 1) {
  # pick the modal chr defensively
  gene_chr <- names(sort(table(df_gene_all$chr), decreasing = TRUE))[1]
}
 
obs_introns <- df_gene_all %>%
  transmute(a1 = intron1_start, b1 = intron1_end,
            a2 = intron2_start, b2 = intron2_end) %>%
  { bind_rows(select(., start = a1, end = b1),
              select(., start = a2, end = b2)) } %>%
  distinct()
 
bed_chr <- bed %>% filter(chr == gene_chr)
 
# A transcript "belongs" to this gene if any of its BED introns coincides
# (exact coordinate match, allowing the BED 0-based vs data convention offset)
# with an observed intron. Try exact match first; if nothing matches, retry
# with a ±1 tolerance on start to absorb 0-/1-based differences.
match_tx <- function(tol) {
  bed_chr %>%
    rowwise() %>%
    mutate(hit = any(abs(start - obs_introns$start) <= tol &
                     abs(end   - obs_introns$end)   <= tol)) %>%
    ungroup() %>%
    filter(hit) %>%
    pull(transcript) %>%
    unique()
}
cand_tx <- match_tx(0L)
if (length(cand_tx) == 0) cand_tx <- match_tx(1L)
if (length(cand_tx) == 0)
  stop("No GENCODE transcript on ", gene_chr,
       " matched any observed intron for ", GENE_NAME,
       ". Check INTRON_BED build/coordinate convention.")
 
# Longest transcript = widest (min start, max end) across its introns.
tx_span <- bed_chr %>%
  filter(transcript %in% cand_tx) %>%
  group_by(transcript) %>%
  summarise(lo = min(start), hi = max(end),
            n_introns = n(), .groups = "drop") %>%
  mutate(span = hi - lo) %>%
  arrange(desc(span), desc(n_introns))
 
chosen_tx <- tx_span$transcript[1]
message("  Gene ", GENE_NAME, " -> ", length(cand_tx),
        " candidate transcript(s); using longest: ", chosen_tx,
        " (", tx_span$n_introns[1], " introns, span ",
        tx_span$span[1], " bp)")
 
# ── Authoritative intron table from GENCODE, ordered by transcription ───────
# Keep BED intron_num as the TRUE index. Sort by genomic start for plotting,
# but preserve intron_num for labels so numbering matches the browser and is
# correct on minus-strand transcripts.
introns_tbl <- bed_chr %>%
  filter(transcript == chosen_tx) %>%
  transmute(start, end, strand,
            intron_num,                 # true GENCODE index (1-based, tx-order)
            mid = (start + end) / 2,
            length = end - start) %>%
  arrange(start) %>%
  mutate(number = intron_num)           # 'number' drives arc join + labels
 
gene_strand <- introns_tbl$strand[1]
 
message("  Backbone introns from GENCODE: ", nrow(introns_tbl),
        " (strand ", gene_strand, ")")
 
# ── Exons = gaps between consecutive GENCODE introns (+ terminal stubs) ──────
EXON_STUB <- 200
 
introns_sorted <- introns_tbl %>% arrange(start)
 
exons_internal <- if (nrow(introns_sorted) >= 2) {
  tibble(
    start = introns_sorted$end[-nrow(introns_sorted)],
    end   = introns_sorted$start[-1]
  ) %>% filter(end > start)
} else tibble(start = numeric(), end = numeric())
 
# Exon numbering in transcription order: on + strand exon 1 is leftmost;
# on - strand exon 1 is rightmost. We number left->right then flip labels
# for minus strand so they read like the browser.
exons_all <- bind_rows(
  tibble(start = min(introns_sorted$start) - EXON_STUB,
         end   = min(introns_sorted$start)),
  exons_internal,
  tibble(start = max(introns_sorted$end),
         end   = max(introns_sorted$end) + EXON_STUB)
) %>%
  arrange(start) %>%
  mutate(
    mid          = (start + end) / 2,
    genomic_rank = row_number(),
    number       = if (identical(gene_strand, "-"))
                     (n() - genomic_rank + 1L) else genomic_rank,
    constitutive = genomic_rank > 1 & genomic_rank < n()
  )
 
# ── Warn if observed pairs reference introns absent from chosen transcript ──
# (These arcs will be dropped by the coordinate join in build_arcs — surfacing
# it here prevents silent loss.)
obs_keys <- unique(paste(obs_introns$start, obs_introns$end))
bed_keys <- unique(paste(introns_tbl$start, introns_tbl$end))
n_unmatched <- sum(!obs_keys %in% bed_keys)
if (n_unmatched > 0)
  message("  NOTE: ", n_unmatched, " observed intron(s) not in chosen ",
          "transcript ", chosen_tx, "; their arcs will be skipped. ",
          "If this is high, the gene's reads may favour a different isoform ",
          "than the longest transcript.")

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