#!/bin/bash
#SBATCH --job-name=structure_array
#SBATCH --array=0-99          # matches --n-chunks in merge_and_split.py (0-indexed)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00       # 24h per chunk — adjust if chunks are large
#SBATCH --output=out_files/structure_%A_%a.out
#SBATCH --error=out_files/structure_%A_%a.err
#SBATCH --partition=batch

# -----------------------------------------------------------------------
# Runs analyze_transcript.py on one chunk of the merged splicing order TSV.
# Each array task processes one chunk file independently.
# Combine outputs with combine_structure_chunks.py after all tasks finish.
# -----------------------------------------------------------------------

CHUNK_DIR="/users/dhan30/scratch/splicing_order/merged/sig_chunks"   # new chunk dir
OUTPUT_DIR="/users/dhan30/scratch/splicing_order/merged/sig_structure_chunks_0414"  # new output dir
SCRIPT_DIR="/users/dhan30/splicing_order/scripts"
GENOME_FASTA="/users/dhan30/reference/hg38.fa"

mkdir -p "$OUTPUT_DIR"
mkdir -p out_files

CHUNK_IDX=$(printf "%03d" "$SLURM_ARRAY_TASK_ID")
CHUNK_FILE="${CHUNK_DIR}/chunk_${CHUNK_IDX}.tsv"
OUTPUT_FILE="${OUTPUT_DIR}/structure_chunk_${CHUNK_IDX}.tsv"

if [ ! -f "$CHUNK_FILE" ]; then
    echo "ERROR: chunk file not found: $CHUNK_FILE"
    exit 1
fi

# skip if already done
if [ -f "$OUTPUT_FILE" ] && [ -s "$OUTPUT_FILE" ]; then
    line_count=$(wc -l < "$OUTPUT_FILE")
    if [ "$line_count" -gt 1 ]; then
        echo "Output already exists with $((line_count - 1)) rows, skipping."
        exit 0
    fi
fi

echo "$(date) | Array task $SLURM_ARRAY_TASK_ID | Processing $CHUNK_FILE"
n_pairs=$(tail -n +2 "$CHUNK_FILE" | wc -l)
echo "$(date) | Pairs in this chunk: $n_pairs"

# load environment
module load miniforge3/25.3.0-3
eval "$(conda shell.bash hook)"
conda activate order_env
export PATH="/users/dhan30/.conda/envs/order_env/bin:$PATH"

# Run structure analysis
python3 "${SCRIPT_DIR}/analyze_transcript.py" \
    "$CHUNK_FILE" \
    "$GENOME_FASTA" \
    -o "$OUTPUT_FILE" \
    || { echo "ERROR: analyze_transcript.py failed for chunk $CHUNK_IDX"; exit 1; }

echo "$(date) | Chunk $CHUNK_IDX complete"
n_out=$(tail -n +2 "$OUTPUT_FILE" | wc -l)
echo "$(date) | Output rows: $n_out / $n_pairs"

conda deactivate

