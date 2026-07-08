#!/usr/bin/env python3
import gzip
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

INTRON_BED = "/users/dhan30/reference/hg38.gencode.basic.v43.introns.bed.gz"
POOLED_TSV = "/users/dhan30/splicing_order/data/splicing_order_pooled.tsv"
OUT_TSV    = "/users/dhan30/splicing_order/data/splicing_order_pooled_adjacent.tsv"
OUT_FIG    = "/users/dhan30/splicing_order/adjacent_fraction_downstream.png"

# --- Load annotation ---
print("Loading intron BED...", flush=True)
transcript_intron_starts = defaultdict(set)

with gzip.open(INTRON_BED, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 4:
            continue
        start = int(parts[1])
        name  = parts[3]
        transcript_id = name.split('_intron_')[0] if '_intron_' in name else name
        transcript_intron_starts[transcript_id].add(start)

print(f"Loaded {len(transcript_intron_starts):,} transcripts", flush=True)

# --- Load pooled TSV ---
print("Loading pooled TSV...", flush=True)
df = pd.read_csv(POOLED_TSV, sep='\t')
print(f"Total pairs: {len(df):,}", flush=True)

# --- Filter ---
def is_adjacent(row):
    tid = row['gene_id']
    starts = transcript_intron_starts.get(tid)
    if starts is None:
        starts = transcript_intron_starts.get(tid.rsplit('.', 1)[0])
    if starts is None:
        return False
    return not any(row['intron1_end'] < s < row['intron2_start'] for s in starts)

print("Filtering...", flush=True)
df_adj = df[df.apply(is_adjacent, axis=1)].copy()
print(f"Adjacent pairs: {len(df_adj):,} ({100*len(df_adj)/len(df):.1f}%)", flush=True)

df_adj.to_csv(OUT_TSV, sep='\t', index=False)
print(f"Written to {OUT_TSV}", flush=True)

# --- Plot ---
mean_val = df_adj['fraction_downstream'].mean()

fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(df_adj['fraction_downstream'], bins=50, color='steelblue', edgecolor='white', linewidth=0.3)
ax.axvline(0.5,      color='black', linestyle='dashed', linewidth=1.2, label='Null (0.5)')
ax.axvline(mean_val, color='#E64B35', linestyle='solid', linewidth=1.2,
           label=f'Mean = {mean_val:.3f}')
ax.set_xlabel('Fraction Downstream Spliced First')
ax.set_ylabel('Number of Intron Pairs')
ax.set_title(f'Splicing order: adjacent pairs only (n = {len(df_adj):,})')
ax.legend()
plt.tight_layout()
fig.savefig(OUT_FIG, dpi=150)
print(f"Figure saved to {OUT_FIG}", flush=True)