#!/usr/bin/env python3
"""
coverage_bias_check.py

Tests whether the global downstream-first splicing bias could be an artifact
of 3' read coverage bias in the RNA-seq data.

For each intron pair, computes its normalized position within the transcript
(0 = 5' end, 1 = 3' end). Then writes a TSV with the position added, for
plotting in R.

Logic:
  - intron1 is always the upstream (5') intron, intron2 is always downstream
  - Strand is inferred from genomic coordinates: if intron1_start < intron2_start,
    gene is on the plus strand; otherwise minus strand
  - Gene extent is computed from the data itself (min/max coords per gene_id)
  - Minus-strand positions are flipped so 0 = 5' end for all genes

Usage:
    python coverage_bias_check.py \
        --input /users/dhan30/scratch/splicing_order/merged/0407/splicing_order_pooled.tsv \
        --output /users/dhan30/scratch/splicing_order/coverage_bias/pooled_with_position.tsv \
        --min-reads 10

Author: Squash
"""

import argparse
import pandas as pd
import numpy as np
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Add normalized transcript position to splicing order TSV")
    parser.add_argument("--input", required=True,
                        help="Path to pooled splicing order TSV")
    parser.add_argument("--output", required=True,
                        help="Output TSV path (input + normalized_position column)")
    parser.add_argument("--min-reads", type=int, default=10,
                        help="Minimum total reads per pair [default: 10]")
    return parser.parse_args()

def main():
    args = parse_args()

    print(f"Reading: {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    print(f"  Loaded {len(df):,} rows")
    print(f"  Columns: {list(df.columns)}")

    # Apply minimum read filter
    before = len(df)
    df = df[df["total"] >= args.min_reads].copy()
    print(f"  After min_reads >= {args.min_reads}: {len(df):,} rows "
          f"(dropped {before - len(df):,})")

    # -------------------------------------------------------------------------
    # Infer strand from relative positions of intron1 vs intron2
    # Plus strand:  intron1_start < intron2_start (transcription left→right)
    # Minus strand: intron1_start > intron2_start (transcription right→left)
    # -------------------------------------------------------------------------
    df["strand"] = np.where(df["intron1_start"] < df["intron2_start"], "+", "-")

    plus_pct = (df["strand"] == "+").mean() * 100
    minus_pct = (df["strand"] == "-").mean() * 100
    print(f"\n  Strand inference:")
    print(f"    Plus-strand pairs:  {plus_pct:.1f}%")
    print(f"    Minus-strand pairs: {minus_pct:.1f}%")

    # -------------------------------------------------------------------------
    # Compute gene extent per gene_id from the data itself
    # Use all intron coordinates to define the gene span
    # -------------------------------------------------------------------------
    print("\n  Computing gene extents from intron coordinates...")

    gene_extent = df.groupby("gene_id").agg(
        gene_start=("intron1_start", "min"),
        gene_end=("intron2_end", "max")
    ).reset_index()

    # For minus-strand genes the min start is actually the 3' end
    # but we just need the span — the flip handles direction
    df = df.merge(gene_extent, on="gene_id", how="left")
    df["transcript_length"] = df["gene_end"] - df["gene_start"]

    # Guard against zero-length (shouldn't happen but just in case)
    df = df[df["transcript_length"] > 0].copy()

    # -------------------------------------------------------------------------
    # Compute midpoint of the intron pair, then normalized genomic position
    # -------------------------------------------------------------------------
    df["pair_midpoint"] = (df["intron1_start"] + df["intron2_end"]) / 2

    df["raw_position"] = (
        (df["pair_midpoint"] - df["gene_start"]) / df["transcript_length"]
    )

    # Flip minus-strand so that 0 = 5' end for all genes
    df["normalized_position"] = np.where(
        df["strand"] == "-",
        1.0 - df["raw_position"],
        df["raw_position"]
    )

    # Sanity check — should be in [0, 1]
    out_of_range = ((df["normalized_position"] < 0) |
                    (df["normalized_position"] > 1)).sum()
    if out_of_range > 0:
        print(f"  WARNING: {out_of_range} pairs have normalized_position outside [0,1]")

    print(f"\n  normalized_position summary:")
    print(f"    min:    {df['normalized_position'].min():.4f}")
    print(f"    median: {df['normalized_position'].median():.4f}")
    print(f"    max:    {df['normalized_position'].max():.4f}")
    print(f"    mean:   {df['normalized_position'].mean():.4f}")

    # -------------------------------------------------------------------------
    # Quick Pearson correlation check (printed here, also done in R)
    # -------------------------------------------------------------------------
    from scipy.stats import pearsonr
    r, p = pearsonr(df["normalized_position"], df["fraction_downstream"])
    print(f"\n  Pearson r (position vs fraction_downstream): {r:.4f}  p = {p:.2e}")
    if abs(r) < 0.05:
        print("  → No meaningful positional trend. Coverage bias unlikely.")
    else:
        print("  → Non-trivial correlation detected. Inspect the plot carefully.")

    # -------------------------------------------------------------------------
    # Write output
    # -------------------------------------------------------------------------
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    # Drop helper columns not needed downstream
    out_cols = [c for c in df.columns if c not in
                ["gene_start", "gene_end", "transcript_length",
                 "pair_midpoint", "raw_position"]]
    df[out_cols].to_csv(args.output, sep="\t", index=False)
    print(f"\n  Output written to: {args.output}")
    print(f"  Rows: {len(df):,}")

if __name__ == "__main__":
    main()