#!/usr/bin/env python3
"""
merge_and_split.py

Merges all pairwise splicing order TSVs from the 50-sample run,
pools read counts for intron pairs seen in multiple samples,
optionally filters by intron length, and splits into N chunks
for parallel structure analysis.

Usage:
    python3 scripts/merge_and_split.py \
        --results-dir /users/dhan30/scratch/splicing_order/results_subset50 \
        --output-dir  /users/dhan30/scratch/splicing_order/merged/0414 \
        --n-chunks    50 \
        --max-intron-length 10000
"""

import argparse
import glob
import os
import pandas as pd
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results-dir',       required=True)
    parser.add_argument('--output-dir',        required=True)
    parser.add_argument('--n-chunks',          type=int, default=20)
    parser.add_argument('--min-reads',         type=int, default=10,
                        help='Min pooled read support after merging')
    parser.add_argument('--max-intron-length', type=int, default=None,
                        help='Max length (bp) for BOTH introns in a pair. '
                             'Default: no filter')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load all splicing order TSVs
    pattern = os.path.join(args.results_dir, '*', '*_pairwise_splicing_order.tsv')
    files = sorted(glob.glob(pattern))
    print(f"Found {len(files)} TSV files")

    dfs = []
    for f in files:
        sample_id = os.path.basename(os.path.dirname(f))
        df = pd.read_csv(f, sep='\t')
        df['sample_id'] = sample_id
        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)
    print(f"Total rows before processing: {len(merged):,}")

    # Compute intron lengths
    merged['intron1_length'] = merged['intron1_end'] - merged['intron1_start']
    merged['intron2_length'] = merged['intron2_end'] - merged['intron2_start']

    print(f"\nIntron length distribution (upstream intron):")
    for p in [25, 50, 75, 90, 95, 99]:
        print(f"  p{p}: {merged['intron1_length'].quantile(p/100):.0f} bp")

    # Optional intron length filter
    if args.max_intron_length:
        before = len(merged)
        merged = merged[
            (merged['intron1_length'] <= args.max_intron_length) &
            (merged['intron2_length'] <= args.max_intron_length)
        ]
        after = len(merged)
        print(f"\nAfter max_intron_length <= {args.max_intron_length}bp: "
              f"{after:,} rows ({100*after/before:.1f}% retained)")
    else:
        print(f"\nNo intron length filter applied.")

    # Pool read counts across samples for the same intron pair
    group_cols = ['chr', 'gene_id', 'intron1_start', 'intron1_end',
                  'intron2_start', 'intron2_end',
                  'intron1_length', 'intron2_length']

    pooled = merged.groupby(group_cols).apply(
        lambda g: pd.Series({
            'upstream':   g['upstream'].sum(),
            'downstream': g['downstream'].sum(),
            'total':      g['total'].sum(),
            'n_samples':  len(g),
        })
    ).reset_index()

    pooled['fraction_downstream'] = pooled['downstream'] / pooled['total']

    before = len(pooled)
    pooled = pooled[pooled['total'] >= args.min_reads].reset_index(drop=True)
    print(f"\nAfter pooling and min_reads>={args.min_reads}: "
          f"{len(pooled):,} unique intron pairs "
          f"(dropped {before - len(pooled):,} low-coverage pairs)")
    print(f"  Pairs seen in >1 sample: {(pooled['n_samples'] > 1).sum():,}")
    print(f"  fraction_downstream mean: {pooled['fraction_downstream'].mean():.3f}")
    print(f"  fraction_downstream std:  {pooled['fraction_downstream'].std():.3f}")

    merged_path = os.path.join(args.output_dir, 'splicing_order_pooled.tsv')
    pooled.to_csv(merged_path, sep='\t', index=False)
    print(f"\nPooled TSV written to {merged_path}")

    # Shuffle then split into chunks
    # Shuffle so each chunk gets a mix of long/short introns — avoids
    # one chunk being all huge introns and running much slower than others
    pooled = pooled.sample(frac=1, random_state=42).reset_index(drop=True)

    chunk_size = int(np.ceil(len(pooled) / args.n_chunks))
    chunk_dir = os.path.join(args.output_dir, 'chunks')
    os.makedirs(chunk_dir, exist_ok=True)

    actual_chunks = 0
    for i in range(args.n_chunks):
        chunk = pooled.iloc[i * chunk_size : (i + 1) * chunk_size]
        if len(chunk) == 0:
            continue
        chunk_path = os.path.join(chunk_dir, f'chunk_{i:03d}.tsv')
        chunk.to_csv(chunk_path, sep='\t', index=False)
        actual_chunks += 1

    print(f"Split into {actual_chunks} chunks (~{chunk_size:,} pairs each)")
    print(f"Chunks written to {chunk_dir}")
    print(f"\nNext steps:")
    print(f"  1. Update run_structure_array.sh: #SBATCH --array=0-{actual_chunks-1}")
    print(f"  2. sbatch run_structure_array.sh")


if __name__ == '__main__':
    main()