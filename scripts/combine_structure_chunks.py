#!/usr/bin/env python3
"""
combine_structure_chunks.py

Combines all structure chunk TSVs into one final output file.
Run this after the structure array job completes.

Usage:
    python3 combine_structure_chunks.py \
        --chunk-dir  /users/dhan30/scratch/splicing_order/merged/sig_structure_chunks \
        --output     /users/dhan30/scratch/splicing_order/merged/sig_structure_features_final.tsv
"""

import argparse
import glob
import os
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chunk-dir', required=True)
    parser.add_argument('--output',    required=True)
    args = parser.parse_args()

    files = sorted(glob.glob(os.path.join(args.chunk_dir, 'structure_chunk_*.tsv')))
    print(f"Found {len(files)} chunk files")

    dfs = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        dfs.append(df)
        print(f"  {os.path.basename(f)}: {len(df):,} rows")

    combined = pd.concat(dfs, ignore_index=True)
    print(f"\nTotal rows: {len(combined):,}")
    print(f"Columns: {list(combined.columns)}")

    # Sanity check: confirm diff_ columns exist and have the right sign convention
    diff_cols = [c for c in combined.columns if c.startswith('diff_')]
    print(f"\ndiff_* columns (downstream - upstream): {diff_cols}")

    combined.to_csv(args.output, sep='\t', index=False)
    print(f"\nWritten to {args.output}")

    try:
        from scipy.stats import pearsonr
        print("\nQuick correlation preview (diff_* vs fraction_downstream):")
        for col in diff_cols:
            valid = combined[['fraction_downstream', col]].dropna()
            if len(valid) > 100:
                r, p = pearsonr(valid['fraction_downstream'], valid[col])
                print(f"  {col}: r={r:.3f}, p={p:.2e}, n={len(valid):,}")
    except ImportError:
        pass

if __name__ == '__main__':
    main()