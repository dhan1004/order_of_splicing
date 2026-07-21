#!/usr/bin/env python3
"""
Ortholog-matched gene subset from human + mouse significant_pairs.tsv.

Keyspace problem this solves:
  - human significant_pairs.tsv : has gene_id (ENST) AND gene_symbol -> use symbol
  - mouse  significant_pairs.tsv: gene_id is a versioned ENSMUST transcript id,
                                  NO symbol -> must bridge ENSMUST -> gene_name
  - MGI .rpt                    : keyed on SYMBOL, grouped by 'DB Class Key'
                                  (the homology group). NO Ensembl ids.

Chain:
  mouse ENSMUST --(GENCODE GTF: transcript_id -> gene_name)--> mouse symbol
  mouse symbol  --(MGI: symbol -> DB Class Key, mouse rows)--> homology key
  homology key  --(MGI: DB Class Key -> symbol, human rows)--> human symbol
  human symbol  == human gene_symbol from the human pairs file

Matching orthologs by DB Class Key (not string equality) is the correct MGI
approach: human & mouse rows sharing a key are curated orthologs, which handles
case differences and symbol drift.

Fetch inputs once on Oscar:
  curl -sL https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt \
    -o HOM_MouseHumanSequence.rpt
  # GTF already unzipped for STAR, e.g. gencode.vM36.annotation.gtf

Usage:
python scripts/supp/build_ortholog_subset.py \
  --human-pairs data/significant_pairs.tsv \
  --mouse-pairs results/mouse/significant_pairs.tsv \
  --gtf /users/dhan30/reference/gencode.vM36.basic.annotation.gtf \
  --ortholog-rpt /users/dhan30/reference/HOM_MouseHumanSequence.rpt \
  --out-prefix results/ortholog/
"""
import argparse
import re
import pandas as pd

TX_RE = re.compile(r'transcript_id "([^"]+)"')
GN_RE = re.compile(r'gene_name "([^"]+)"')


def strip_ver(s):
    return s.str.replace(r"\.\d+$", "", regex=True)


def build_tx2sym(gtf_path):
    """ENSMUST (unversioned) -> mouse gene symbol, from a GENCODE GTF."""
    seen = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if 'transcript_id "' not in line:
                continue
            tx = TX_RE.search(line)
            gn = GN_RE.search(line)
            if not (tx and gn):
                continue
            txid = tx.group(1).split(".")[0]
            if txid not in seen:
                seen[txid] = gn.group(1)
    return pd.DataFrame(
        {"mouse_tx": list(seen.keys()), "mouse_symbol": list(seen.values())}
    )


def build_ortholog_map(rpt_path):
    """DB Class Key -> (human_symbol, mouse_symbol), one row per H/M pairing."""
    df = pd.read_csv(rpt_path, sep="\t", dtype=str)
    key, org, sym = "DB Class Key", "Common Organism Name", "Symbol"
    df = df[[key, org, sym]].dropna(subset=[key, sym])
    h = (df[df[org].str.contains("human", case=False, na=False)][[key, sym]]
         .rename(columns={sym: "human_symbol"}))
    m = (df[df[org].str.contains("mouse", case=False, na=False)][[key, sym]]
         .rename(columns={sym: "mouse_symbol"}))
    return (h.merge(m, on=key)
              .rename(columns={key: "db_key"})
              .drop_duplicates())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--human-pairs", required=True)
    ap.add_argument("--mouse-pairs", required=True)
    ap.add_argument("--gtf", required=True, help="GENCODE vM36 GTF (unzipped)")
    ap.add_argument("--ortholog-rpt", required=True)
    ap.add_argument("--out-prefix", default="")
    args = ap.parse_args()

    human_df = pd.read_csv(args.human_pairs, sep="\t", dtype=str)
    mouse_df = pd.read_csv(args.mouse_pairs, sep="\t", dtype=str)

    # ---- mouse: attach a symbol via the GTF bridge ----
    mouse_df["mouse_tx"] = strip_ver(mouse_df["gene_id"])
    tx2sym = build_tx2sym(args.gtf)
    mouse_df = mouse_df.merge(tx2sym, on="mouse_tx", how="left")
    unmapped = mouse_df["mouse_symbol"].isna().sum()

    mouse_syms = set(mouse_df["mouse_symbol"].dropna().str.upper())
    human_syms = set(human_df["gene_symbol"].dropna().str.upper())

    # ---- ortholog groups present on BOTH sides ----
    omap = build_ortholog_map(args.ortholog_rpt)
    omap["h_ok"] = omap["human_symbol"].str.upper().isin(human_syms)
    omap["m_ok"] = omap["mouse_symbol"].str.upper().isin(mouse_syms)
    matched = (omap[omap["h_ok"] & omap["m_ok"]]
               [["db_key", "human_symbol", "mouse_symbol"]]
               .drop_duplicates())
    matched.to_csv(f"{args.out_prefix}ortholog_matched_genes.tsv",
                   sep="\t", index=False)

    # ---- subset each pair file by qualifying symbols (case-insensitive) ----
    h_keep = set(matched["human_symbol"].str.upper())
    m_keep = set(matched["mouse_symbol"].str.upper())
    (human_df[human_df["gene_symbol"].str.upper().isin(h_keep)]
        .to_csv(f"{args.out_prefix}human_subset_pairs.tsv", sep="\t", index=False))
    (mouse_df[mouse_df["mouse_symbol"].str.upper().isin(m_keep)]
        .to_csv(f"{args.out_prefix}mouse_subset_pairs.tsv", sep="\t", index=False))

    print(f"mouse transcripts unmapped to symbol: {unmapped}")
    print(f"human symbols with pairs: {len(human_syms)}")
    print(f"mouse symbols with pairs: {len(mouse_syms)}")
    print(f"ortholog-matched qualifying groups: {len(matched)}")


if __name__ == "__main__":
    main()