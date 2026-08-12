#!/usr/bin/env python3
"""
join_spanning_order.py

Join the pooled frozen metric (from pool_spanning.py) to splicing order
(significant_pairs.tsv) and prepare model-ready tables for BOTH framings:

  (A) Pairwise  : each retained intron is the first- or second-spliced member
                  of an informative pair. Frozen prediction: the LATER-spliced
                  (retained/pre-mRNA) intron shows higher frozen_score.
  (B) Gene-ordinal: each retained intron gets a Copeland-derived norm_position
                  [0,1]. Frozen prediction: frozen_score rises with norm_position.

Both framings split on is_first_pair (first introns splice fast, shouldn't
freeze -> internal negative control).

Join key tolerance: BED half-open vs read-derived coords can differ by +-1 nt,
so we match retained-intron coordinates within JOIN_TOL. Minus-strand
separations use abs().

Outputs:
  frozen_pairwise.tsv       (for glmer: cbind(threep,fivep) ~ role + covars + (1|intron))
  frozen_gene_ordinal.tsv   (for glm:   cbind(threep,fivep) ~ norm_position + covars)
"""

import argparse
import sys
import numpy as np
import pandas as pd

JOIN_TOL = 1  # +-1 nt tolerance for BED half-open vs read-derived coords


def load_pooled(path):
    d = pd.read_csv(path, sep="\t", dtype={"chr": str, "gene_id": str,
                                           "retained_strand": str})
    need = {"chr", "gene_id", "retained_intron_start", "retained_intron_end",
            "fivep_span", "threep_span"}
    missing = need - set(d.columns)
    if missing:
        sys.exit(f"pooled file missing columns: {missing}")
    if "frozen_score" not in d.columns:
        tot = d["fivep_span"] + d["threep_span"]
        d["frozen_score"] = np.where(tot > 0, d["threep_span"] / tot, np.nan)
    if "retained_intron_length" not in d.columns:
        d["retained_intron_length"] = (d["retained_intron_end"]
                                       - d["retained_intron_start"]).abs()
    if "total_span" not in d.columns:
        d["total_span"] = d["fivep_span"] + d["threep_span"]
    return d


def load_order(path):
    """significant_pairs.tsv. Expected to identify each informative pair with
    two introns and which spliced first. We stay tolerant to column naming and
    resolve the fields we actually need."""
    o = pd.read_csv(path, sep="\t", dtype=str)
    return o


def build_tol_key(df, chrom, start, end, strand=None):
    """Return a merge-friendly integer key list at +-JOIN_TOL by expanding
    each interval into candidate (chr, start_bucket) — simpler: we round
    coords into JOIN_TOL-wide buckets so half-open/read offsets collide."""
    # bucket = floor(coord / (JOIN_TOL... )) is fragile; instead do an explicit
    # tolerance merge below. This helper kept for clarity of intent.
    raise NotImplementedError


def tol_merge(pooled, order, o_chr, o_start, o_end, o_strand):
    """Merge pooled frozen introns onto an order-table intron column set,
    allowing +-JOIN_TOL on start AND end. Implemented by generating the small
    set of (start+ds, end+de) offset variants on the ORDER side and merging
    exact — cheap because tolerance is +-1."""
    keyp = ["chr", "retained_intron_start", "retained_intron_end"]
    frames = []
    for ds in range(-JOIN_TOL, JOIN_TOL + 1):
        for de in range(-JOIN_TOL, JOIN_TOL + 1):
            oo = order.copy()
            oo["_start"] = pd.to_numeric(oo[o_start], errors="coerce") + ds
            oo["_end"]   = pd.to_numeric(oo[o_end],   errors="coerce") + de
            oo["chr"] = oo[o_chr].astype(str)
            m = pooled.merge(
                oo, left_on=keyp, right_on=["chr", "_start", "_end"],
                how="inner", suffixes=("", "_ord"))
            frames.append(m)
    if not frames:
        return pooled.iloc[0:0]
    out = pd.concat(frames, ignore_index=True)
    # a pooled intron may match several offset variants -> dedup on the
    # pooled-intron identity + the order-row identity
    out = out.drop_duplicates()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pooled", required=True,
                    help="spanning_analysis_set.tsv from pool_spanning.py")
    ap.add_argument("--order", required=True,
                    help="significant_pairs.tsv")
    ap.add_argument("--out-prefix", default="frozen")
    # column-name knobs for significant_pairs.tsv (override if yours differ)
    ap.add_argument("--o-chr", default="chr")
    ap.add_argument("--o-strand", default="strand")
    # the intron-1 / intron-2 coordinate columns in significant_pairs
    ap.add_argument("--i1-start", default="intron1_start")
    ap.add_argument("--i1-end",   default="intron1_end")
    ap.add_argument("--i2-start", default="intron2_start")
    ap.add_argument("--i2-end",   default="intron2_end")
    ap.add_argument("--fraction-downstream", default="fraction_downstream")
    ap.add_argument("--is-first-pair", default="is_first_pair")
    ap.add_argument("--norm-position", default="norm_position",
                    help="gene-ordinal Copeland position; if absent it is "
                         "derived per gene from fraction_downstream ranks")
    args = ap.parse_args()

    pooled = load_pooled(args.pooled)
    order  = load_order(args.order)

    print(f"pooled introns: {len(pooled):,}", file=sys.stderr)
    print(f"order rows    : {len(order):,}", file=sys.stderr)
    print(f"order columns : {list(order.columns)}", file=sys.stderr)

    # ---------------- (A) PAIRWISE framing --------------------------------
    # For each informative pair, the retained (pre-mRNA-side) intron is one
    # member. We match the pooled intron to BOTH intron slots and record
    # whether it is the first- or second-spliced member of that pair.
    #
    # fraction_downstream = P(downstream spliced first). If retained intron is
    # the downstream one: it is "first" when fraction_downstream high.
    # We label role = 'second_spliced' (i.e. later -> expected frozen) vs
    # 'first_spliced' based on the pair call. Because significant_pairs schemas
    # vary, we derive role from fraction_downstream + which slot matched.
    def match_slot(sstart, send, slot):
        m = tol_merge(pooled, order, args.o_chr, sstart, send, args.o_strand)
        m["matched_slot"] = slot
        return m

    m1 = match_slot(args.i1_start, args.i1_end, "intron1")
    m2 = match_slot(args.i2_start, args.i2_end, "intron2")
    pw = pd.concat([m1, m2], ignore_index=True)

    if len(pw) == 0:
        print("WARN: zero pairwise matches. Check --i1/--i2 column names "
              "against significant_pairs.tsv header printed above.",
              file=sys.stderr)
    else:
        fd_col = args.fraction_downstream
        if fd_col in pw.columns:
            pw[fd_col] = pd.to_numeric(pw[fd_col], errors="coerce")
            # Convention: intron1 = upstream, intron2 = downstream (adjust if
            # your table differs). Retained intron is 'later-spliced' when the
            # OTHER intron tends to splice first.
            #   matched intron2 (downstream) is later-spliced when fd is LOW
            #   matched intron1 (upstream)   is later-spliced when fd is HIGH
            later = np.where(
                pw["matched_slot"].eq("intron2"), pw[fd_col] < 0.5,
                pw[fd_col] > 0.5)
            pw["retained_is_later_spliced"] = later
            pw["role"] = np.where(later, "later_spliced", "earlier_spliced")
        for c in (args.is_first_pair,):
            if c in pw.columns:
                pw[c] = pw[c].astype(str).str.lower().isin(["true", "1", "yes"])
        keep_cols = (["chr", "gene_id", "retained_intron_start",
                      "retained_intron_end", "retained_strand",
                      "fivep_span", "threep_span", "total_span",
                      "frozen_score", "retained_intron_length",
                      "matched_slot", "role", "retained_is_later_spliced"]
                     + [c for c in (fd_col, args.is_first_pair)
                        if c in pw.columns])
        keep_cols = [c for c in keep_cols if c in pw.columns]
        pw_out = pw[keep_cols].drop_duplicates()
        pw_out.to_csv(f"{args.out_prefix}_pairwise.tsv", sep="\t",
                      index=False, float_format="%.4f")
        print(f"pairwise matched rows: {len(pw_out):,} -> "
              f"{args.out_prefix}_pairwise.tsv", file=sys.stderr)

    # ---------------- (B) GENE-ORDINAL framing ----------------------------
    # One row per retained intron with a norm_position in [0,1]. Prefer an
    # existing norm_position column on the order table; else derive per gene.
    go = None
    if args.norm_position in order.columns:
        npos = order[[args.o_chr, args.i1_start, args.i1_end,
                      args.norm_position]].copy()
        npos.columns = ["chr", "_s", "_e", "norm_position"]
        go = tol_merge(pooled, npos.rename(columns={"_s": args.i1_start,
                                                    "_e": args.i1_end}),
                       "chr", args.i1_start, args.i1_end, args.o_strand)
    else:
        print("norm_position not in order table; falling back to per-gene "
              "rank of frozen introns by genomic position within gene.",
              file=sys.stderr)
        go = pooled.copy()
        go = go.sort_values(["gene_id", "retained_intron_start"])
        go["rank_in_gene"] = go.groupby("gene_id").cumcount()
        n = go.groupby("gene_id")["rank_in_gene"].transform("max").clip(lower=1)
        go["norm_position"] = go["rank_in_gene"] / n

    go_cols = ["chr", "gene_id", "retained_intron_start", "retained_intron_end",
               "retained_strand", "fivep_span", "threep_span", "total_span",
               "frozen_score", "retained_intron_length", "norm_position"]
    go_cols = [c for c in go_cols if c in go.columns]
    go_out = go[go_cols].drop_duplicates()
    go_out.to_csv(f"{args.out_prefix}_gene_ordinal.tsv", sep="\t",
                  index=False, float_format="%.4f")
    print(f"gene-ordinal rows: {len(go_out):,} -> "
          f"{args.out_prefix}_gene_ordinal.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()