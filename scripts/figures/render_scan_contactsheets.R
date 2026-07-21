#!/usr/bin/env Rscript
# =============================================================================
# render_scan_contactsheets.R
#
# Exploratory cross-species scan. For every distinct ortholog-matched gene,
# draws the human panel above the mouse panel (via plot_gene_splicing from
# plot_gene_splicing_lib.R) and packs GENES_PER_PAGE genes per page into
# multi-page PDF contact sheets. Genes are sorted by cross-species concordance
# so the most-conserved-order genes appear on the first pages.
#
# Inputs (produced by build_ortholog_subset.py):
#   ortholog_matched_genes.tsv   db_key, human_symbol, mouse_symbol
#   human_subset_pairs.tsv       human significant_pairs rows (has gene_symbol)
#   mouse_subset_pairs.tsv       mouse rows + mouse_symbol column (GTF-bridged)
#
# The mouse subset file has no `gene_symbol` column (the original mouse
# significant_pairs.tsv had none); the builder added `mouse_symbol`, so we
# alias it to gene_symbol here before plotting.
#
# Usage:
#   Rscript scripts/figures/render_scan_contactsheets.R \
#     results/ortholog/ortholog_matched_genes.tsv \
#     results/ortholog/human_subset_pairs.tsv \
#     results/ortholog/mouse_subset_pairs.tsv \
#     figures/ortholog/scan_out

# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
})
source("scripts/figures/plot_gene_splicing_lib.R")

args      <- commandArgs(trailingOnly = TRUE)
matched_f <- args[1]
human_f   <- args[2]
mouse_f   <- args[3]
outdir    <- if (length(args) >= 4) args[4] else "scan_out"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

GENES_PER_PAGE <- 6   # 6 genes x 2 species panels = 12 panels/page
PNG_PAGES      <- TRUE # also write one PNG per contact-sheet page (page_001.png ...)
PNG_DPI        <- 150  # 150 is plenty for on-screen scanning; 300 for print

matched <- read.delim(matched_f, stringsAsFactors = FALSE)
human   <- normalise_pairs(read.delim(human_f, stringsAsFactors = FALSE,
                                      check.names = FALSE))
mouse   <- normalise_pairs(read.delim(mouse_f, stringsAsFactors = FALSE,
                                      check.names = FALSE))

# mouse subset uses `mouse_symbol`; alias to gene_symbol for the plot fn
if (!"gene_symbol" %in% names(mouse) && "mouse_symbol" %in% names(mouse))
  mouse$gene_symbol <- mouse$mouse_symbol

# ── dedup paralogs: one row per distinct human gene ─────────────────────────
gene_pairs <- matched %>%
  distinct(human_symbol, mouse_symbol) %>%
  distinct(human_symbol, .keep_all = TRUE)   # first mouse partner per human gene

# ── concordance score for sort order ────────────────────────────────────────
# mean signed effect (fraction_downstream - 0.5) per gene, each species;
# concordance = product (same sign & large magnitude in both -> high).
gene_effect <- function(df, sym) {
  sub <- df[df$gene_symbol == sym, "fraction_downstream"]
  if (length(sub) == 0) return(NA_real_)
  mean(sub, na.rm = TRUE) - 0.5
}
gene_pairs$h_eff <- vapply(gene_pairs$human_symbol,
                           function(s) gene_effect(human, s), numeric(1))
gene_pairs$m_eff <- vapply(gene_pairs$mouse_symbol,
                           function(s) gene_effect(mouse, s), numeric(1))
gene_pairs$concordance <- gene_pairs$h_eff * gene_pairs$m_eff
gene_pairs <- gene_pairs %>%
  arrange(desc(concordance))   # concordant (both same direction) first

n_genes <- nrow(gene_pairs)
message("distinct genes to render: ", n_genes)

# ── render loop: one PDF, GENES_PER_PAGE genes per page ─────────────────────
pdf_out <- file.path(outdir, "ortholog_scan_contactsheet.pdf")
# page height scales with genes per page (2 panels each)
pdf(pdf_out, width = 11, height = 2.2 * GENES_PER_PAGE, onefile = TRUE)

page_h <- 2.2 * GENES_PER_PAGE
page_panels <- list()
page_idx    <- 0
flush_page <- function(panels) {
  if (length(panels) == 0) return(invisible())
  page_idx <<- page_idx + 1
  sheet <- wrap_plots(panels, ncol = 1)
  print(sheet)                                   # -> open PDF device
  if (PNG_PAGES) {
    png_page <- file.path(outdir, sprintf("page_%03d.png", page_idx))
    ggplot2::ggsave(png_page, sheet, width = 11, height = page_h,
                    dpi = PNG_DPI, limitsize = FALSE)
  }
}

rendered <- 0
for (i in seq_len(n_genes)) {
  hs <- gene_pairs$human_symbol[i]
  ms <- gene_pairs$mouse_symbol[i]

  ph <- tryCatch(plot_gene_splicing(human, hs, "human"),
                 error = function(e) NULL)
  pm <- tryCatch(plot_gene_splicing(mouse, ms, "mouse"),
                 error = function(e) NULL)
  if (is.null(ph) && is.null(pm)) next

  # stack this gene's two species panels
  pair_plot <- if (!is.null(ph) && !is.null(pm)) ph / pm
               else if (!is.null(ph)) ph else pm
  page_panels[[length(page_panels) + 1]] <- pair_plot
  rendered <- rendered + 1

  if (length(page_panels) >= GENES_PER_PAGE) {
    flush_page(page_panels)
    page_panels <- list()
  }
}
flush_page(page_panels)   # trailing partial page
dev.off()

message("rendered ", rendered, " genes -> ", pdf_out)
# also write the sorted gene list so you can map page position -> gene
write.table(gene_pairs, file.path(outdir, "scan_gene_order.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
message("gene order/concordance -> ",
        file.path(outdir, "scan_gene_order.tsv"))