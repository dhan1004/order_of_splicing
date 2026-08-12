#!/usr/bin/env python3
"""
pool_spanning.py

Pool all GSM*_pairwise_spanning.tsv files, compute the frozen-intermediate
metric per retained intron by POOLING RAW COUNTS (never averaging per-row
fractions), and emit summary stats + a plotting-ready table.

Frozen signature (PI hypothesis / cola-seq):
  A retained intron frozen between step 1 and step 2 has its 5'ss already
  cleaved (fivep_span absent) but 3'ss intact (threep_span present).
  frozen_score = threep_total / (threep_total + fivep_total), bounded [0,1].
  High frozen_score  -> consistent with step-1-frozen intermediate.

Input : /users/dhan30/scratch/data/results/GSM*/GSM*_pairwise_spanning.tsv
Output: pooled_spanning_by_intron.tsv   (one row per retained intron, pooled)
        spanning_summary_stats.txt
        spanning_pooling_manifest.tsv    (per-sample row counts, for QC)
"""

import argparse
import glob
import os
import sys
import gzip
import pandas as pd
import numpy as np

RESULTS_DIR = "/users/dhan30/scratch/data/results"

# columns that define a unique retained intron (the join key later)
INTRON_KEY = ["chr", "gene_id", "retained_intron_start",
              "retained_intron_end", "retained_strand"]

# raw count columns we pool
COUNT_COLS = ["fivep_span", "fivep_over", "threep_span", "threep_over"]


def find_files(results_dir):
    pat = os.path.join(results_dir, "GSM*", "GSM*_pairwise_spanning.tsv")
    files = sorted(glob.glob(pat))
    return files


def load_one(path):
    """Load a single sample TSV, coercing count cols to numeric.
    Per-row fraction/ratio columns are IGNORED on purpose — we recompute
    from pooled raw counts."""
    opener = gzip.open if path.endswith(".gz") else open
    df = pd.read_csv(path, sep="\t", dtype={"chr": str, "gene_id": str,
                                            "retained_strand": str,
                                            "direction": str})
    # keep only key + raw counts (+ direction/is_first if present, carried later)
    for c in COUNT_COLS:
        if c not in df.columns:
            df[c] = 0
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    # sample id from filename
    gsm = os.path.basename(path).split("_pairwise_spanning")[0]
    df["sample"] = gsm
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default=RESULTS_DIR)
    ap.add_argument("--out-prefix", default="spanning")
    ap.add_argument("--min-total-span", type=int, default=10,
                    help="min (fivep_span+threep_span) pooled per intron to "
                         "keep it in the analysis set")
    args = ap.parse_args()

    files = find_files(args.results_dir)
    if not files:
        sys.exit(f"No files matched under {args.results_dir}")
    print(f"Found {len(files)} sample files", file=sys.stderr)

    frames = []
    manifest = []
    for i, f in enumerate(files, 1):
        try:
            df = load_one(f)
        except Exception as e:
            print(f"  WARN skipping {f}: {e}", file=sys.stderr)
            manifest.append({"sample": os.path.basename(f), "rows": 0,
                             "status": f"error:{e}"})
            continue
        manifest.append({"sample": df["sample"].iat[0] if len(df) else
                         os.path.basename(f), "rows": len(df), "status": "ok"})
        frames.append(df)
        if i % 100 == 0:
            print(f"  loaded {i}/{len(files)}", file=sys.stderr)

    allrows = pd.concat(frames, ignore_index=True)
    print(f"Concatenated {len(allrows):,} rows from {len(frames)} samples",
          file=sys.stderr)

    # --- POOL raw counts across all samples, per retained intron ---
    # also carry the modal direction (upstream/downstream) for context
    agg = (allrows
           .groupby(INTRON_KEY, as_index=False)[COUNT_COLS]
           .sum())

    # attach n_samples the intron was seen in, and modal direction
    nsamp = (allrows.groupby(INTRON_KEY)["sample"].nunique()
             .rename("n_samples").reset_index())
    agg = agg.merge(nsamp, on=INTRON_KEY, how="left")

    if "direction" in allrows.columns:
        # modal direction per intron
        modedir = (allrows.groupby(INTRON_KEY)["direction"]
                   .agg(lambda s: s.mode().iat[0] if len(s.mode()) else "NA")
                   .rename("direction").reset_index())
        agg = agg.merge(modedir, on=INTRON_KEY, how="left")

    # --- frozen metric from POOLED counts ---
    total_span = agg["fivep_span"] + agg["threep_span"]
    agg["total_span"] = total_span
    agg["frozen_score"] = np.where(total_span > 0,
                                   agg["threep_span"] / total_span, np.nan)
    # log ratio (guarded), useful for modelling / diagnostics
    agg["log2_3p_over_5p"] = np.log2(
        (agg["threep_span"] + 0.5) / (agg["fivep_span"] + 0.5))
    agg["retained_intron_length"] = (agg["retained_intron_end"]
                                     - agg["retained_intron_start"]).abs()

    # analysis set: coverage floor
    keep = agg["total_span"] >= args.min_total_span
    analysis = agg[keep].copy()

    # --- write tables ---
    out_pool = f"{args.out_prefix}_pooled_by_intron.tsv"
    agg.to_csv(out_pool, sep="\t", index=False, float_format="%.4f")
    out_analysis = f"{args.out_prefix}_analysis_set.tsv"
    analysis.to_csv(out_analysis, sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(manifest).to_csv(f"{args.out_prefix}_pooling_manifest.tsv",
                                  sep="\t", index=False)

    # --- summary stats ---
    lines = []
    def add(s): lines.append(s); print(s)
    add("=" * 60)
    add("PAIRWISE SPANNING — POOLED SUMMARY")
    add("=" * 60)
    add(f"samples pooled           : {len(frames)}")
    add(f"total input rows         : {len(allrows):,}")
    add(f"unique retained introns  : {len(agg):,}")
    add(f"  with >=1 span          : {(agg['total_span']>0).sum():,}")
    add(f"analysis set (>= {args.min_total_span} span): {len(analysis):,}")
    add("")
    add("Pooled span totals (analysis set):")
    add(f"  5'ss span total        : {int(analysis['fivep_span'].sum()):,}")
    add(f"  3'ss span total        : {int(analysis['threep_span'].sum()):,}")
    gt = analysis['threep_span'].sum()
    g5 = analysis['fivep_span'].sum()
    add(f"  global 3p/(3p+5p)      : {gt/(gt+g5):.4f}")
    add("")
    add("frozen_score distribution (analysis set):")
    add(analysis["frozen_score"].describe().to_string())
    add("")
    # fraction of introns that are 3'ss-dominant (frozen-like)
    fz = (analysis["frozen_score"] >= 0.9).mean()
    add(f"introns with frozen_score >= 0.90 : {fz:.3f}")
    fzc = (analysis["frozen_score"] <= 0.1).mean()
    add(f"introns with frozen_score <= 0.10 : {fzc:.3f}")
    add("")
    if "direction" in analysis.columns:
        add("frozen_score by direction (mean):")
        add(analysis.groupby("direction")["frozen_score"].agg(
            ["mean", "median", "count"]).to_string())

    with open(f"{args.out_prefix}_summary_stats.txt", "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(f"\nWrote:\n  {out_pool}\n  {out_analysis}\n"
          f"  {args.out_prefix}_summary_stats.txt\n"
          f"  {args.out_prefix}_pooling_manifest.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()