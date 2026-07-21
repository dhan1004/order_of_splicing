#!/usr/bin/env python3
# scripts/supp/gene_blurb.py
# Usage:
#   python scripts/supp/gene_blurb.py \
#     --common results/ortholog/ortholog_matched_genes.tsv \
#     --out    results/ortholog/gene_blurbs.tsv
import argparse, time, requests, pandas as pd

def query(symbols, species):
    # species: "human" (9606) or "mouse" (10090)
    out = {}
    url = "https://mygene.info/v3/query"
    for i in range(0, len(symbols), 100):          # batch
        batch = symbols[i:i+100]
        r = requests.post(url, data={
            "q": ",".join(batch),
            "scopes": "symbol",
            "fields": "symbol,name,summary",
            "species": species,
        })
        r.raise_for_status()
        for hit in r.json():
            sym = hit.get("query")
            if sym and sym not in out:
                summ = hit.get("summary") or hit.get("name") or ""
                out[sym] = summ
        time.sleep(0.3)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--common", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    c = pd.read_csv(a.common, sep="\t", dtype=str)
    # header-detect: your file may be positional
    if "human_symbol" not in c.columns:
        c.columns = ["db_key", "human_symbol", "mouse_symbol"] + list(c.columns[3:])
    h = query(sorted(c["human_symbol"].dropna().unique()), "human")
    m = query(sorted(c["mouse_symbol"].dropna().unique()), "mouse")
    rows = []
    for sym, txt in h.items(): rows.append(("human", sym, txt))
    for sym, txt in m.items(): rows.append(("mouse", sym, txt))
    pd.DataFrame(rows, columns=["species","symbol","blurb"]).to_csv(
        a.out, sep="\t", index=False)
    print(f"Wrote {len(rows)} blurbs to {a.out}")

if __name__ == "__main__":
    main()