#!/usr/bin/env python3
"""
pool_pairwise_spanning.py

Pools the per-sample {gsm}_pairwise_spanning.tsv files (from
pairwise_splice_site_spanning.py) into pairwise_spanning_pooled.tsv.

Counts are SUMMED per (retained intron, spliced intron, direction) across
samples, then fractions recomputed from the pooled sums (coverage-weighted,
never averaged — same rule as merge_and_split.py and pool_spanning.py).

The grouping key includes `direction` so an intron pair that shows both
orders in different transcripts is kept as two rows (retained-when-downstream
vs retained-when-upstream). For the frozen-intermediate question you usually
want the rows where the intron of interest is the RETAINED (later-spliced) one,
which is already what each row encodes.

Usage:
  python3 pool_pairwise_spanning.py \
      --results-dir /users/dhan30/scratch/data/results \
      --output-dir  results/merged \
      --min-over    10
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--min-over", type=int, default=10,
                    help="min pooled overlapping reads at BOTH sites for the "
                         "filtered summary (default: 10)")
    args = ap.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    pattern = os.path.join(args.results_dir, "*", "*_pairwise_spanning.tsv")
    files = sorted(glob.glob(pattern))
    print(f"Found {len(files)} pairwise spanning TSV files")
    if not files:
        raise SystemExit("No pairwise spanning TSVs found — did the hook run?")

    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f, sep="\t")
        except Exception as e:
            print(f"  WARNING skip {f}: {e}")
            continue
        if not df.empty:
            dfs.append(df)
    merged = pd.concat(dfs, ignore_index=True)
    print(f"Total per-sample pair rows: {len(merged):,}")

    count_cols = ["fivep_span", "fivep_over", "threep_span", "threep_over"]
    for c in count_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0).astype(int)

    key = ["chr", "gene_id",
           "retained_intron_start", "retained_intron_end", "retained_strand",
           "spliced_intron_start", "spliced_intron_end", "direction"]

    pooled = (merged.groupby(key)
              .apply(lambda g: pd.Series({
                  "fivep_span":  g["fivep_span"].sum(),
                  "fivep_over":  g["fivep_over"].sum(),
                  "threep_span": g["threep_span"].sum(),
                  "threep_over": g["threep_over"].sum(),
                  "n_samples":   len(g),
              }))
              .reset_index())

    pooled["fivep_frac"] = np.where(
        pooled["fivep_over"] > 0,
        pooled["fivep_span"] / pooled["fivep_over"], np.nan)
    pooled["threep_frac"] = np.where(
        pooled["threep_over"] > 0,
        pooled["threep_span"] / pooled["threep_over"], np.nan)
    pooled["span_ratio_3p_over_5p"] = np.where(
        (pooled["fivep_frac"] > 0) & pooled["fivep_frac"].notna()
        & pooled["threep_frac"].notna(),
        pooled["threep_frac"] / pooled["fivep_frac"], np.nan)

    pooled = pooled.sort_values(["chr", "retained_intron_start",
                                 "retained_intron_end"])
    out_all = os.path.join(args.output_dir, "pairwise_spanning_pooled.tsv")
    pooled.to_csv(out_all, sep="\t", index=False)
    print(f"\nWrote all pairs: {out_all}  ({len(pooled):,} pair-direction rows)")

    keep = ((pooled["fivep_over"] >= args.min_over)
            & (pooled["threep_over"] >= args.min_over))
    filt = pooled[keep].copy()
    out_filt = os.path.join(args.output_dir, "pairwise_spanning_pooled_filtered.tsv")
    filt.to_csv(out_filt, sep="\t", index=False)
    print(f"Wrote powered pairs (both over >= {args.min_over}): "
          f"{out_filt}  ({len(filt):,} rows)")

    if len(filt):
        med = filt["span_ratio_3p_over_5p"].median()
        frac = (filt["threep_frac"] > filt["fivep_frac"]).mean()
        print("\n--- Retained-intron spanning summary (powered pairs) ---")
        print(f"  median 3'ss/5'ss span ratio: {med:.3f}")
        print(f"  fraction with 3'ss > 5'ss coverage: {frac:.3f}")
        print("  (>1 / >0.5 consistent with step-1 intermediates frozen "
              "before the second step on the later-spliced intron)")


if __name__ == "__main__":
    main()