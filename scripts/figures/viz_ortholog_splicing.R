#!/usr/bin/env Rscript
# =============================================================================
# viz_ortholog_splicing.R
#
# Draws mouse + human splicing-order gene structures for a RANGE of ortholog
# pairs, stacked (mouse on top, human below), one pair per PDF page.
# Designed to be driven by a SLURM array: each task handles a row range of the
# common-genes TSV and emits ONE multi-page PDF.
#
# Common-genes TSV columns (tab-sep, no header assumed -- see COMMON_HAS_HEADER):
#   db_key   human_symbol   mouse_symbol
#   e.g.  51806712   ALDH1L1   Aldh1l1
#
# Per-species intron TSVs already carry their own gene_symbol column
# (human -> gene_symbol matches human_symbol; mouse -> gene_symbol matches
#  mouse_symbol). No BED lookup needed.
#
# Usage:
#   Rscript viz_ortholog_splicing.R \
#     --common   common_genes.tsv \
#     --human    human_introns.tsv \
#     --mouse    mouse_introns.tsv \
#     --outdir   figures/ortho \
#     --start    1 \
#     --end      50
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# ── Defaults ──────────────────────────────────────────────────────────────────
COMMON_FILE <- NULL
HUMAN_FILE  <- NULL
MOUSE_FILE  <- NULL
OUT_DIR     <- "figures/ortho"
START_ROW   <- 1L
END_ROW     <- NA_integer_
COMMON_HAS_HEADER <- FALSE
BLURB_FILE <- NULL
# Column names of the intron TSVs that identify the gene.
# Leave NULL to auto-detect (gene_symbol, human_symbol, mouse_symbol, gene_id, ...).
HUMAN_SYMBOL_COL <- NULL
MOUSE_SYMBOL_COL <- NULL

# Output format: "png" (grid of pairs per file) or "pdf" (multi-page).
OUT_FORMAT <- "png"
PER_FILE   <- 4L     # ortholog pairs per PNG (ignored for pdf, which is 1/page)
DPI        <- 150L

# ── CLI ───────────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag) { i <- which(args == flag); if (length(i) && i < length(args)) args[i + 1] else NULL }
if (!is.null(v <- getarg("--common")))  COMMON_FILE <- v
if (!is.null(v <- getarg("--human")))   HUMAN_FILE  <- v
if (!is.null(v <- getarg("--mouse")))   MOUSE_FILE  <- v
if (!is.null(v <- getarg("--outdir")))  OUT_DIR     <- v
if (!is.null(v <- getarg("--start")))   START_ROW   <- as.integer(v)
if (!is.null(v <- getarg("--end")))     END_ROW     <- as.integer(v)
if ("--common-header" %in% args)        COMMON_HAS_HEADER <- TRUE
if (!is.null(v <- getarg("--human-symbol-col"))) HUMAN_SYMBOL_COL <- v
if (!is.null(v <- getarg("--mouse-symbol-col"))) MOUSE_SYMBOL_COL <- v
if (!is.null(v <- getarg("--format")))   OUT_FORMAT <- tolower(v)
if (!is.null(v <- getarg("--per-png")))  PER_FILE   <- as.integer(v)
if (!is.null(v <- getarg("--dpi")))      DPI        <- as.integer(v)
if (!is.null(v <- getarg("--blurbs")))  BLURB_FILE <- v   # add in CLI block
if (!OUT_FORMAT %in% c("png", "pdf"))
  stop("--format must be 'png' or 'pdf' (got '", OUT_FORMAT, "')")

stopifnot(!is.null(COMMON_FILE), !is.null(HUMAN_FILE), !is.null(MOUSE_FILE))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Palette ───────────────────────────────────────────────────────────────────
COL_DOWNSTREAM <- "#E05C5C"
COL_UPSTREAM   <- "#4B8BBE"
COL_EXON       <- "#6a7b8c"
COL_BACKBONE   <- "grey75"

# ── Load intron TSVs ONCE (shared across all pairs in this task) ──────────────
message("Reading human introns: ", HUMAN_FILE)
human_raw <- read.delim(HUMAN_FILE, stringsAsFactors = FALSE, check.names = FALSE)
message("Reading mouse introns: ", MOUSE_FILE)
mouse_raw <- read.delim(MOUSE_FILE, stringsAsFactors = FALSE, check.names = FALSE)

# Normalise the columns each species table needs.
prep <- function(df, symcol) {
  ren <- function(d, old, new) { if (old %in% names(d) && !new %in% names(d)) names(d)[names(d) == old] <- new; d }
  df <- ren(df, "upstream",   "upstream_count")
  df <- ren(df, "downstream", "downstream_count")
  df <- ren(df, "total",      "total_reads")
  if (!"fraction_downstream" %in% names(df) &&
      all(c("downstream_count", "total_reads") %in% names(df)))
    df$fraction_downstream <- df$downstream_count / df$total_reads
  if (!"intron1_length" %in% names(df))
    df$intron1_length <- df$intron1_end - df$intron1_start
  if (!"intron2_length" %in% names(df))
    df$intron2_length <- df$intron2_end - df$intron2_start
  # Resolve the gene-symbol column and normalise it to "gene_symbol".
  #  1. If the caller passed an explicit column that exists, use it.
  #  2. Else auto-detect among common spellings.
  #  3. Else fail loudly, showing the columns actually present.
  if (!is.null(symcol) && symcol %in% names(df)) {
    if (symcol != "gene_symbol") df$gene_symbol <- df[[symcol]]
  } else if ("gene_symbol" %in% names(df)) {
    # already correct
  } else {
    candidates <- c("human_symbol", "mouse_symbol", "symbol",
                    "gene_name", "geneName", "gene", "gene_id")
    hit <- candidates[candidates %in% names(df)][1]
    if (is.na(hit))
      stop("No gene-symbol column found. Looked for: gene_symbol, ",
           paste(candidates, collapse = ", "),
           ".\n  Columns present: ", paste(names(df), collapse = ", "),
           "\n  Set HUMAN_SYMBOL_COL / MOUSE_SYMBOL_COL (or --human-symbol-col / ",
           "--mouse-symbol-col) to the right name.")
    message("  Auto-detected gene-symbol column: '", hit, "'")
    df$gene_symbol <- df[[hit]]
  }
  df
}
human_raw <- prep(human_raw, HUMAN_SYMBOL_COL)
mouse_raw <- prep(mouse_raw, MOUSE_SYMBOL_COL)

# Pre-split by gene for fast per-gene subsetting.
human_by_gene <- split(human_raw, human_raw$gene_symbol)
mouse_by_gene <- split(mouse_raw, mouse_raw$gene_symbol)

# ── Gene-function blurbs (optional) ──────────────────────────────────────────
blurb_lookup <- function(symbol, species) ""   # default no-op
if (!is.null(BLURB_FILE) && file.exists(BLURB_FILE)) {
  bl <- read.delim(BLURB_FILE, stringsAsFactors = FALSE)
  bl$key <- paste(tolower(bl$species), toupper(bl$symbol), sep = "|")
  bmap <- setNames(bl$blurb, bl$key)
  blurb_lookup <- function(symbol, species) {
    v <- bmap[[paste(tolower(species), toupper(symbol), sep = "|")]]
    if (is.null(v) || is.na(v)) "" else v
  }
}

# wrap long text so it fits the panel; trim to ~1-2 lines
wrap_blurb <- function(txt, width = 90, max_lines = 2) {
  if (!nzchar(txt)) return("")
  w <- strwrap(txt, width = width)
  if (length(w) > max_lines) { w <- w[seq_len(max_lines)]; w[max_lines] <- paste0(w[max_lines], " \u2026") }
  paste(w, collapse = "\n")
}

# ── Load common-genes list, slice to this task's range ────────────────────────
common <- read.delim(COMMON_FILE, header = COMMON_HAS_HEADER,
                     stringsAsFactors = FALSE, check.names = FALSE)
if (!COMMON_HAS_HEADER) {
  # positional: db_key, human_symbol, mouse_symbol
  names(common)[1:3] <- c("db_key", "human_symbol", "mouse_symbol")
} else {
  # accept common header spellings
  lc <- tolower(names(common))
  names(common)[grepl("human", lc)] <- "human_symbol"
  names(common)[grepl("mouse", lc)] <- "mouse_symbol"
}

n_total <- nrow(common)
if (is.na(END_ROW)) END_ROW <- n_total
START_ROW <- max(1L, START_ROW)
END_ROW   <- min(n_total, END_ROW)
if (START_ROW > END_ROW) { message("Empty range; nothing to do."); quit(status = 0) }
slice <- common[START_ROW:END_ROW, , drop = FALSE]
message(sprintf("Plotting pairs %d..%d of %d", START_ROW, END_ROW, n_total))

# =============================================================================
# Build one gene panel from that gene's intron-pair rows.
# Reconstructs the exon/intron backbone from intron coordinates and draws a
# significance-coloured arc per informative pair.
# =============================================================================
build_gene_panel <- function(gdf, gene_name, species_label, blurb="") {

  if (is.null(gdf) || nrow(gdf) == 0) {
    return(
      ggplot() + theme_void() +
        annotate("text", x = 0, y = 0,
                 label = paste0(species_label, ": ", gene_name,
                                " \u2014 no data"),
                 size = 3.2, colour = "grey50") +
        coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))
    )
  }

  # Assemble the unique intron set from both members of each pair.
  introns <- bind_rows(
    gdf %>% transmute(start = intron1_start, end = intron1_end),
    gdf %>% transmute(start = intron2_start, end = intron2_end)
  ) %>%
    distinct() %>%
    arrange(start) %>%
    mutate(number = row_number(), mid = (start + end) / 2)

  gene_lo <- min(introns$start)
  gene_hi <- max(introns$end)

  # Exons = the gaps flanking/between introns (approximate model).
  bnds <- sort(unique(c(introns$start, introns$end)))
  exons <- tibble(start = c(gene_lo, bnds), end = c(bnds, gene_hi)) %>%
    filter(end > start) %>%
    anti_join(introns, by = c("start", "end")) %>%
    mutate(mid = (start + end) / 2, number = row_number())

  # Direction per pair (fallback if no 'direction' column present).
  if (!"direction" %in% names(gdf)) {
    gdf <- gdf %>% mutate(
      direction = case_when(
        fraction_downstream > 0.5 ~ "downstream_biased",
        fraction_downstream < 0.5 ~ "upstream_biased",
        TRUE ~ "not_significant"
      )
    )
  }

  arcs <- gdf %>%
    mutate(
      x1  = (intron1_start + intron1_end) / 2,
      x2  = (intron2_start + intron2_end) / 2,
      col = case_when(
        direction == "downstream_biased" ~ COL_DOWNSTREAM,
        direction == "upstream_biased"   ~ COL_UPSTREAM,
        TRUE                             ~ "grey70"
      )
    )

  x_pad <- (gene_hi - gene_lo) * 0.02
  p <- ggplot() +
    # backbone
    annotate("segment", x = gene_lo, xend = gene_hi, y = 0, yend = 0,
             colour = COL_BACKBONE, linewidth = 1) +
    # exon blocks
    geom_rect(data = exons,
              aes(xmin = start, xmax = end, ymin = -0.25, ymax = 0.25),
              fill = COL_EXON, colour = NA)

  # arcs: simple quadratic curves via geom_curve, coloured by direction
  if (nrow(arcs) > 0) {
    p <- p + geom_curve(
      data = arcs,
      aes(x = x1, xend = x2, y = 0.28, yend = 0.28),
      curvature = -0.3, linewidth = 0.5, colour = arcs$col, alpha = 0.8
    )
  }

  p +
    geom_text(data = introns, aes(x = mid, y = 0.42, label = number),
              size = 2.1, colour = "black") +
    coord_cartesian(xlim = c(gene_lo - x_pad, gene_hi + x_pad),
                    ylim = c(-0.6, 1.0), clip = "off") +
    scale_x_continuous(labels = function(x) paste0(round(x / 1e3, 1), " kb"),
                       expand = c(0, 0)) +
    labs(title = paste0(species_label, ": ", gene_name), 
      subtitle = wrap_blurb(blurb, width = 180),
      x = NULL, y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      axis.line.y  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x  = element_line(colour = "grey70", linewidth = 0.4),
      axis.text.x  = element_text(colour = "grey40", size = 7),
      plot.title   = element_text(face = "bold", size = 10, hjust = 0),
      plot.margin  = margin(6, 12, 4, 12),
      plot.subtitle = element_text(size = 6.5, colour = "grey35",
                                   lineheight = 0.95, hjust = 0)
    )
}

# ── Legend strip (built once, appended under each pair) ───────────────────────
legend_df <- tibble(
  x = c(NA_real_, NA_real_), y = c(NA_real_, NA_real_),
  bias = factor(c("Downstream-biased", "Upstream-biased"))
)
legend_plot <- ggplot(legend_df, aes(x, y, colour = bias)) +
  # geom_point needs only one obs per group -> no "single observation" warning
  geom_point(size = 2.4, na.rm = TRUE) +
  scale_colour_manual(values = c("Downstream-biased" = COL_DOWNSTREAM,
                                 "Upstream-biased"   = COL_UPSTREAM),
                      name = NULL, drop = FALSE) +
  guides(colour = guide_legend(override.aes = list(shape = 15, size = 4))) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# =============================================================================
# Build one stacked panel (mouse over human) per ortholog pair.
# =============================================================================
build_pair_panel <- function(hs, ms) {
  hdf <- human_by_gene[[hs]]
  mdf <- mouse_by_gene[[ms]]
  if (is.null(hdf) && is.null(mdf)) return(NULL)   # signal: skip

  p_mouse <- build_gene_panel(mdf, ms, "Mouse", blurb_lookup(ms, "mouse"))
  p_human <- build_gene_panel(hdf, hs, "Human", blurb_lookup(hs, "human"))

  (p_mouse / p_human) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = paste0(hs, " / ", ms),
      theme = theme(plot.title = element_text(face = "bold", size = 12))
    )
}

# Collect panels + labels for the pairs that have data.
panels <- list(); labels <- character(0)
n_skip <- 0L
for (i in seq_len(nrow(slice))) {
  pnl <- build_pair_panel(slice$human_symbol[i], slice$mouse_symbol[i])
  if (is.null(pnl)) { n_skip <- n_skip + 1L; next }
  panels[[length(panels) + 1L]] <- pnl
  labels <- c(labels, paste0(slice$human_symbol[i], "/", slice$mouse_symbol[i]))
}
n_ok <- length(panels)

if (n_ok == 0L) {
  message(sprintf("No pairs with data in range %d..%d; nothing written.",
                  START_ROW, END_ROW))
  quit(status = 0)
}

# ── PDF: one pair per page (legend strip under each) ──────────────────────────
if (OUT_FORMAT == "pdf") {
  out_pdf <- file.path(OUT_DIR,
    sprintf("ortholog_splicing_%06d_%06d.pdf", START_ROW, END_ROW))
  message("Writing PDF: ", out_pdf)
  cairo_pdf(out_pdf, width = 10, height = 6, onefile = TRUE)
  for (pnl in panels) print(pnl / legend_plot + plot_layout(heights = c(1, 0.08)))
  dev.off()
  message(sprintf("Done. %d pages written, %d pairs skipped (no data).",
                  n_ok, n_skip))

# ── PNG: grid of PER_FILE pairs per image ─────────────────────────────────────
} else {
  per <- max(1L, PER_FILE)
  n_files <- ceiling(n_ok / per)
  # Per-pair panel height in inches; total canvas grows with rows.
  pair_h <- 3.4
  fig_w  <- 10
  message(sprintf("Writing %d PNG(s), up to %d pairs each (%d pairs, %d skipped).",
                  n_files, per, n_ok, n_skip))
  for (f in seq_len(n_files)) {
    lo <- (f - 1L) * per + 1L
    hi <- min(f * per, n_ok)
    grp <- panels[lo:hi]
    n_in <- length(grp)

    # Stack this file's pairs vertically, add one shared legend at the bottom.
    combined <- wrap_plots(grp, ncol = 1) / legend_plot +
      plot_layout(heights = c(rep(1, n_in), 0.12))

    # Absolute row-range this PNG covers, for a stable filename.
    row_lo <- START_ROW + (lo - 1L)
    row_hi <- START_ROW + (hi - 1L)
    out_png <- file.path(OUT_DIR,
      sprintf("ortholog_splicing_%06d_%06d.png", row_lo, row_hi))

    ggsave(out_png, combined, width = fig_w,
           height = pair_h * n_in + 0.6, dpi = DPI, bg = "white",
           limitsize = FALSE)
    message("  saved: ", basename(out_png), "  (", n_in, " pairs)")
  }
  message(sprintf("Done. %d pairs across %d PNG(s), %d skipped (no data).",
                  n_ok, n_files, n_skip))
}