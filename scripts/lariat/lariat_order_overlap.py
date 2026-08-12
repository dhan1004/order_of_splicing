#!/usr/bin/env python3
"""
lariat_order_overlap.py

Intersect DBR1-KO lariat data with pairwise splicing-order calls.

Two framings, written to two outputs:

  1. PAIRWISE  -> per-intron-per-pair rows (introns NOT deduplicated; an intron
     appears once per pair it participates in). Each row is tagged with the
     intron's role in that pair (first_spliced / second_spliced) and whether it
     carries lariat support. Because the same intron repeats across its pairs,
     these rows are NOT independent -- use a model with a random effect / clustered
     SEs on `intron_key` when you test (see the R snippet printed at the end).

  2. GENE-ORDINAL -> one row per intron. An intron's position is scored by the
     fraction of its tested pairs in which it splices FIRST (a Copeland-style
     "win fraction"), then ranked within gene. This tolerates intransitive
     (cyclic) pairwise calls and missing pairs, unlike a strict topological sort.

Join key: (chrom, strand, 5'ss). 5'ss = intronX_start on '+' strand,
intronX_end on '-' strand. Lariat `fivep_pos` is already the true 5'ss.
A +/- TOL nt tolerance absorbs BED half-open vs read-derived off-by-one.

Usage:
python scripts/lariat/lariat_order_overlap.py \
    --lariat  data/lariat/all_datasets_HEK293T_lariat_reads_merged.txt \
    --order   data/significant_pairs.tsv \
    --outdir  results/lariat_overlap \
    --tol 1
"""

import argparse
import os
import sys
from collections import defaultdict

import pandas as pd
import numpy as np


# ----------------------------------------------------------------------
# lariat parsing
# ----------------------------------------------------------------------
def load_lariat(path, tol):
    """Return dict: (chrom, strand, rounded_5pss_bucket) -> lariat read count,
    plus a set-based lookup that tolerates +/- tol via bucket expansion."""
    df = pd.read_csv(path, sep="\t", dtype={"chrom": str})
    needed = {"chrom", "strand", "fivep_pos", "read_id"}
    missing = needed - set(df.columns)
    if missing:
        sys.exit(f"[lariat] missing columns: {missing}")

    # distinct reads per (chrom, strand, fivep_pos)
    grp = (
        df.groupby(["chrom", "strand", "fivep_pos"])["read_id"]
        .nunique()
        .reset_index(name="lariat_reads")
    )

    # Build a tolerant lookup: for each observed 5'ss, register it at every
    # offset in [-tol, tol]. Store max reads if collisions.
    lut = defaultdict(int)
    for chrom, strand, pos, n in grp.itertuples(index=False):
        for d in range(-tol, tol + 1):
            key = (chrom, strand, int(pos) + d)
            lut[key] = max(lut[key], int(n))
    return lut, grp


def fivep_of_intron(row, which):
    """5'ss genomic coord for intron1 or intron2 given strand."""
    strand = row["_strand"]
    if which == 1:
        return row["intron1_start"] if strand == "+" else row["intron1_end"]
    else:
        return row["intron2_start"] if strand == "+" else row["intron2_end"]


# ----------------------------------------------------------------------
# order table parsing
# ----------------------------------------------------------------------
def load_order(path):
    df = pd.read_csv(path, sep="\t", dtype={"chr": str})
    # strand isn't a column in the sample you showed; derive it.
    # For a genuine adjacent pair, intron1 is upstream in transcript order.
    # On '+' strand transcript coords increase; intron1_start < intron2_start.
    # On '-' strand transcript order runs opposite to genomic, so intron1 sits
    # at HIGHER genomic coords than intron2. Use that to infer strand.
    if "strand" in df.columns:
        df["_strand"] = df["strand"]
    else:
        df["_strand"] = np.where(
            df["intron1_start"] < df["intron2_start"], "+", "-"
        )
    return df


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    print("main")
    ap = argparse.ArgumentParser()
    ap.add_argument("--lariat", required=True)
    ap.add_argument("--order", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tol", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lut, lariat_grp = load_lariat(args.lariat, args.tol)
    order = load_order(args.order)

    print(f"[info] {len(order)} pairs, {lariat_grp['lariat_reads'].sum()} "
          f"lariat reads across {len(lariat_grp)} distinct 5'ss")

    # -----------------------------------------------------------------
    # Determine, per pair, which intron splices FIRST.
    # direction == 'downstream_biased' -> downstream (intron2) splices first.
    # direction == 'upstream_biased'   -> upstream  (intron1) splices first.
    # fraction_downstream > 0.5 == downstream first (consistent w/ direction).
    # -----------------------------------------------------------------
    def first_is_intron2(r):
        if "direction" in r and isinstance(r["direction"], str):
            return r["direction"] == "downstream_biased"
        return r["fraction_downstream"] > 0.5

    # -----------------------------------------------------------------
    # FRAMING 1: pairwise, per-intron-per-pair (no dedup)
    # -----------------------------------------------------------------
    rows = []
    for _, r in order.iterrows():
        d2_first = first_is_intron2(r)
        for which in (1, 2):
            fp = int(fivep_of_intron(r, which))
            key = (str(r["chr"]), r["_strand"], fp)
            lreads = lut.get(key, 0)
            # role: is THIS intron the first-spliced member of the pair?
            if which == 2:
                is_first = d2_first
                ilen = r["intron2_length"]
            else:
                is_first = not d2_first
                ilen = r["intron1_length"]
            rows.append({
                "gene_id": r["gene_id"],
                "gene_symbol": r.get("gene_symbol", ""),
                "chr": r["chr"],
                "strand": r["_strand"],
                "intron_key": f'{r["chr"]}:{r["_strand"]}:{fp}',
                "fivep_pos": fp,
                "which_in_pair": "downstream" if which == 2 else "upstream",
                "first_spliced": bool(is_first),
                "intron_length": ilen,
                "pair_total_reads": r["total"],
                "is_first_pair": (which == 1) and (r["intron1_start"] == r.get("intron1_start")),
                "lariat_reads": lreads,
                "lariat_present": lreads > 0,
                "source_pair": f'{r["gene_id"]}:{r["intron1_start"]}-{r["intron2_start"]}',
            })
    pairwise = pd.DataFrame(rows)
    pairwise_path = os.path.join(args.outdir, "pairwise_intron_lariat.tsv")
    pairwise.to_csv(pairwise_path, sep="\t", index=False)
    print(f"[out] {pairwise_path}  ({len(pairwise)} intron-in-pair rows)")

    # quick 2x2: first vs second spliced x lariat present
    ct = pd.crosstab(pairwise["first_spliced"], pairwise["lariat_present"])
    print("\n[pairwise] first_spliced x lariat_present (rows not independent!):")
    print(ct)

    # -----------------------------------------------------------------
    # FRAMING 2: gene-ordinal via win-fraction (Copeland score)
    # For each intron (unique 5'ss within a gene), tally across all pairs it is
    # in: how many it splices first (win) vs second (loss). score = wins/(wins+losses).
    # Higher score == earlier in order. Rank within gene by score.
    # -----------------------------------------------------------------
    wins = defaultdict(int)
    losses = defaultdict(int)
    meta = {}  # intron_key -> (gene_id, symbol, chr, strand, fivep, best_len, lariat)
    for _, r in order.iterrows():
        d2_first = first_is_intron2(r)
        for which in (1, 2):
            fp = int(fivep_of_intron(r, which))
            key = (r["gene_id"], str(r["chr"]), r["_strand"], fp)
            this_first = d2_first if which == 2 else (not d2_first)
            if this_first:
                wins[key] += 1
            else:
                losses[key] += 1
            lreads = lut.get((str(r["chr"]), r["_strand"], fp), 0)
            ilen = r["intron2_length"] if which == 2 else r["intron1_length"]
            if key not in meta:
                meta[key] = {
                    "gene_id": r["gene_id"],
                    "gene_symbol": r.get("gene_symbol", ""),
                    "chr": r["chr"],
                    "strand": r["_strand"],
                    "fivep_pos": fp,
                    "intron_length": ilen,
                    "lariat_reads": lreads,
                }

    grows = []
    for key, m in meta.items():
        w, l = wins[key], losses[key]
        n = w + l
        m = dict(m)
        m["n_pairs"] = n
        m["wins_first"] = w
        m["losses_second"] = l
        m["win_fraction"] = w / n if n else np.nan  # 1.0 == always first
        m["lariat_present"] = m["lariat_reads"] > 0
        grows.append(m)
    gene_ord = pd.DataFrame(grows)

    # rank within gene: rank 1 == earliest (highest win_fraction).
    gene_ord["order_rank"] = (
        gene_ord.groupby("gene_id")["win_fraction"]
        .rank(ascending=False, method="dense")
    )
    gene_ord["n_introns_in_gene"] = (
        gene_ord.groupby("gene_id")["fivep_pos"].transform("count")
    )
    # normalized position in [0,1]: 0 == first, 1 == last. Comparable across genes.
    gene_ord["norm_position"] = np.where(
        gene_ord["n_introns_in_gene"] > 1,
        (gene_ord["order_rank"] - 1) / (gene_ord["n_introns_in_gene"] - 1),
        0.0,
    )
    gene_ord = gene_ord.sort_values(["gene_id", "order_rank"])
    gene_path = os.path.join(args.outdir, "gene_ordinal_intron_lariat.tsv")
    gene_ord.to_csv(gene_path, sep="\t", index=False)
    print(f"\n[out] {gene_path}  ({len(gene_ord)} unique introns)")

    # lariat rate by rank (only multi-intron genes)
    multi = gene_ord[gene_ord["n_introns_in_gene"] > 1]
    by_rank = (
        multi.groupby("order_rank")["lariat_present"]
        .agg(["mean", "count"])
        .head(10)
    )
    print("\n[gene-ordinal] lariat rate by order_rank (1 = earliest-spliced):")
    print(by_rank)

    # -----------------------------------------------------------------
    # R modeling snippet (mixed model handles the pairwise non-independence)
    # -----------------------------------------------------------------
    print("""
------------------------------------------------------------------
Test the position effect in R.

PAIRWISE (accounts for repeated introns via random effect):
  library(lme4)
  d <- read.delim("pairwise_intron_lariat.tsv")
  m <- glmer(lariat_present ~ first_spliced + log(intron_length) +
             log(pair_total_reads) + (1 | intron_key),
             data = d, family = binomial)
  summary(m)

GENE-ORDINAL (one row per intron, correlate position with lariat):
  g <- read.delim("gene_ordinal_intron_lariat.tsv")
  g <- subset(g, n_introns_in_gene > 1)
  m2 <- glm(lariat_present ~ norm_position + log(intron_length),
            data = g, family = binomial)
  summary(m2)
  # or continuous: lariat_reads ~ norm_position (negbin, MASS::glm.nb)
------------------------------------------------------------------""")


if __name__ == "__main__":
    main()