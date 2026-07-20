# Co-transcriptional Splicing Order — Analysis Pipeline

## Setup

Create conda environment using yml file:
```bash
conda env create -f environment/order_env.yml
```
Miniforge needs to be loaded in order to access conda commands and environment, but this
should all be taken care of in the scripts themselves.

Reference files expected at `/users/dhan30/reference/`:
- `hg38.fa` — genome FASTA
- `hg38.gencode.basic.v43.introns.bed.gz` — intron annotations
---

## Step-by-Step

### Step 1 - GEO Metadata Filtering
**Script:** `scripts/preprocess/filter_gse_minimal_metadata.py`
 
Reads a local HDF5 file of GEO minimal XML metadata for all human high-throughput sequencing studies. Applies regular-expression filters to find HEK293 RNA-seq studies that use total RNA (not polyA-selected), are paired-end, and have reads ≥80 bp. Excludes single-cell studies.
 
**Input:** `human_expression_profiling_by_HT_seq_gds_metadata.hdf5`  
**Output:** `human_GSE_candidates_hek293.txt.gz`, `human_expression_profiling_by_HT_seq_gds_metadata_indices_hek293.txt`

---

### Step 2 - Parse Sample Metadata
**Script:** `scripts/preprocess/parse_gse_candidate_metadata.py`
 
For each candidate GSE, extracts GSM sample IDs and their associated SRR/SRX run accessions using the SRA–GSM accession table. Produces two output files: a TSV of sample titles and characteristics, and an HDF5 file with structured metadata strings used downstream.
 
**Input:** candidate GSE list (Step 1), `SRA_GSM_accessions.txt.gz`  
**Outputs:**
- `human_gene_perturbation_GSE_candidates_sample_title_characteristics_hek293.txt.gz` - TSV with columns `gse_id`, `gsm_id`, `sample_title`, `sample_characteristics` 
- `human_gene_perturbation_GSE_candidates_sample_metadata_hek293.hdf5` - structured metadata strings per GSE, including GSM IDs and their SRR/SRX accessions; used to build file for part 3
- `gsm_sra_list_for_pipeline_hek293.tsv` - file fed to the downstream pipeline (columns: `gsm_id`, `sra_ids`)

---

### Step 3 - Per-Sample Processing
**Wrapper:** `scripts/run_filtering_all.sh`  
**Per-sample script:** `scripts/order_gse_filtering.sh`
 
Iterates over the GSM/SRR manifest and runs the per-sample pipeline to determine splicing order for each sample. The per-sample script has checkpoints so it can resume from any completed stage (raw fastq -> trimmed -> aligned BAM -> informative pairs -> splicing order).
 
Per-sample steps:
1. **Checkpoints** for check if final outputs, alignments, or trimmed files exist
2. **Download** reads with `fastq-dl` (SRR accessions, concatenating multiple runs) and **trim** with `fastp` (adapter detection, paired-end)
3. **Align** with STAR (two-pass mode, hg38/GENCODE v43)
4. **Extract informative pairs**: junction-spanning reads paired with intron-overlapping reads (Kim et al. 2017 framework), using `samtools` + `bedtools`
5. **Classify splicing order**: for each adjacent intron pair, count upstream-first vs. downstream-first informative pairs; compute `fraction_downstream` using `analyze_splicing_order.py`
**Output per sample:** `{GSM_ID}_pairwise_splicing_order.tsv`  
Columns include: `chr`, `gene_id`, `intron1_start`, `intron1_end`, `intron2_start`, `intron2_end`, `upstream`, `downstream`, `total`, `fraction_downstream`
 
**Note:** Raw and trimmed FASTQs and the STAR BAM are deleted after each stage
 
---

### Step 4 - Pool samples and split into chunks

```
python scripts/merge_and_split.py \
    --results-dir /users/dhan30/scratch/data/results \
    --output-dir  results/merged \
    --n-chunks    100 \
    --min-reads   10 \
    --max-intron-length 10000
```

Pools read counts across samples for each unique intron pair, then splits
into N chunks for parallel structure analysis.

Output: `splicing_order_pooled.tsv` + `chunks/chunk_000.tsv` ...

Use `splicing_order_pooled.tsv` as input to `filter_significant_pairs.R` to
get tsv of significant, strongly biased intron pairs:

```
Rscript scripts/filter_significant_pairs.R \
    --input  results/merged/splicing_order_pooled.tsv \
    --outdir results/merged/
```

Output:
significant_pairs.tsv   filtered to sig pairs only (input for structure array)
ene_summary.tsv        per-gene counts of sig/total pairs
all_pairs_tested.tsv    full table with binom_p, binom_padj, effect, direction

---

### Step 5 - RNA secondary structure prediction

Update the array bound in `run_structure_array.sh` then run

```
sbatch scripts/structure/run_structure_array.sh
```

Each task runs `analyze_transcript.py` on one chunk using ViennaRNA/RNAfold
(sliding window, per-nt accessibility and H-bond counts).

Output: `results/merged/sig_structure_chunks_*/structure_chunk_*.tsv`

---

### Step 6 - Combine structure chunks

```
python scripts/structure/combine_structure_chunks.py \
    --chunk-dir results/merged/sig_structure_chunks_0414/ \
    --output    results/merged/sig_structure_features_final.tsv
```

---

### Step 7 - Figures

All figure scripts live in `scripts/figures/`. Each is self-contained and
reads from the pooled/structure TSVs. Edit the `INPUT_TSV` / `OUT_DIR`
variables at the top of each script, or pass CLI flags.

| Script | Figures produced |
|---|---|
| `viz_order_enrichment_new.R` | Significance testing, volcano, GO enrichment, length/coverage distributions |
| `viz_bonds.R` | H-bond and accessibility correlation figures |
| `viz_splicing_order.R` | Global `fraction_downstream` distributions |
| `viz_gene_splicing.R` | Per-gene intron-level splicing order diagrams |

```bash
Rscript scripts/figures/viz_order_enrichment_new.R \
    --input  results/merged/splicing_order_pooled.tsv \
    --outdir results/figures/
```