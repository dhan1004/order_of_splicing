#!/bin/bash
#SBATCH --job-name=ortho_plots
#SBATCH --array=0-119          # <- set to ceil(N_pairs / CHUNK) - 1  (see below)
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/ortho_%A_%a.out
#SBATCH --error=logs/ortho_%A_%a.err

# =============================================================================
# run_ortholog_plots.sh
#
# Parallelizes viz_ortholog_splicing.R over a SLURM array. Each array task
# handles a contiguous CHUNK of rows from the common-genes TSV and writes ONE
# multi-page PDF (one page per ortholog pair).
#
# --- Sizing the array -------------------------------------------------------
#   N_pairs = number of lines in COMMON_TSV
#   CHUNK   = pairs per task (below)
#   ntasks  = ceil(N_pairs / CHUNK)
#   Set   #SBATCH --array=0-(ntasks-1)
#
#   Example: 6200 pairs, CHUNK=50  ->  ceil(6200/50)=124 tasks -> --array=0-123
#   Keep ntasks < 1001 (Oscar MaxArraySize). If it would exceed, raise CHUNK.
#
#   Quick helper to compute it before submitting:
#     N=$(wc -l < common_genes.tsv); CHUNK=50
#     echo "--array=0-$(( (N + CHUNK - 1)/CHUNK - 1 ))"
#
# --- Submit -----------------------------------------------------------------
#   mkdir -p logs figures/ortho
#   sbatch run_ortholog_plots.sh
# =============================================================================

set -euo pipefail

# ── Paths (edit these) ────────────────────────────────────────────────────────
COMMON_TSV=/users/dhan30/splicing_order/results/ortholog/ortholog_matched_genes.tsv
HUMAN_TSV=/users/dhan30/splicing_order/results/ortholog/human_subset_pairs.tsv
MOUSE_TSV=/users/dhan30/splicing_order/results/ortholog/mouse_subset_pairs.tsv
SCRIPT=/users/dhan30/splicing_order/scripts/figures/viz_ortholog_splicing.R
OUT_DIR=/users/dhan30/splicing_order/figures/ortholog
BLURB_TSV=/users/dhan30/splicing_order/results/ortholog/gene_blurbs.tsv

CHUNK=50       # pairs handled per task
FORMAT=png     # png (grid) or pdf (multi-page)
PER_PNG=4      # ortholog pairs per PNG (ignored for pdf)
DPI=150

# ── Environment ───────────────────────────────────────────────────────────────
module load miniforge3/25.3.0-3
eval "$(conda shell.bash hook)"
conda activate order_env
export PATH="/users/dhan30/.conda/envs/order_env/bin:$PATH"
module load r/4.5.1

mkdir -p "$OUT_DIR" logs

# ── Compute this task's row range (1-based, inclusive) ────────────────────────
TASK=${SLURM_ARRAY_TASK_ID:-0}
START=$(( TASK * CHUNK + 1 ))
END=$(( START + CHUNK - 1 ))

echo "[task $TASK] rows ${START}..${END}  (chunk=${CHUNK})"

Rscript "$SCRIPT" \
  --common "$COMMON_TSV" \
  --human  "$HUMAN_TSV" \
  --mouse  "$MOUSE_TSV" \
  --blurbs "$BLURB_TSV" \
  --outdir "$OUT_DIR" \
  --start  "$START" \
  --end    "$END" \
  --format "$FORMAT" \
  --per-png "$PER_PNG" \
  --dpi    "$DPI"

echo "[task $TASK] done"