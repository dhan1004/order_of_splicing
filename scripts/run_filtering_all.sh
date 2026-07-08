#!/bin/bash
#SBATCH --job-name=order_splice
#SBATCH --array=1-1000%100
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G 
#SBATCH --time=36:00:00
#SBATCH --output=out_files/order_genome_match_%A_%a.out
#SBATCH --error=out_files/order_genome_match_%A_%a.err
#SBATCH --partition=batch

# Configuration
INPUT_LIST_TSV="/users/dhan30/scratch/splicing_order/data/gsm_sra_list_paired_only.tsv"
THREADS=16 
BASE_OUT_DIR="/users/dhan30/scratch/splicing_order/results"
PIPELINE_SCRIPT="/users/dhan30/scratch/splicing_order/scripts/order_gse_filtering.sh"

# Job Info
printf "$(date +'%d/%b/%Y %H:%M:%S') | Running on node: $(hostname)\n"
printf "$(date +'%d/%b/%Y %H:%M:%S') | SLURM assigned node: $SLURM_NODELIST\n"
printf "$(date +'%d/%b/%Y %H:%M:%S') | Job ID: $SLURM_JOB_ID\n"
printf "$(date +'%d/%b/%Y %H:%M:%S') | Array Task ID: $SLURM_ARRAY_TASK_ID\n"

# Get sample info for this array task
LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${LINE_NUM}p" "$INPUT_LIST_TSV")

gsm_id=$(echo "$LINE" | cut -f1 | xargs)
sra_ids=$(echo "$LINE" | cut -f2 | tr -d '\r' | xargs)

if [ -z "$gsm_id" ] || [ -z "$sra_ids" ]; then
    echo "ERROR: Empty GSM ID or SRA IDs for array task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

printf "$(date +'%d/%b/%Y %H:%M:%S') | Processing GSM: $gsm_id\n"
printf "$(date +'%d/%b/%Y %H:%M:%S') | SRA IDs: $sra_ids\n"

# output directory
sample_out_dir="${BASE_OUT_DIR}/${gsm_id}"
mkdir -p "$sample_out_dir"

# run pipeline
log_file="${sample_out_dir}/${gsm_id}_pipeline.log"

printf "$(date +'%d/%b/%Y %H:%M:%S') | Starting pipeline for $gsm_id\n"
bash "$PIPELINE_SCRIPT" "$THREADS" "$sample_out_dir" "$gsm_id" "$sra_ids" > "$log_file" 2>&1

if [ $? -eq 0 ]; then
    printf "$(date +'%d/%b/%Y %H:%M:%S') | SUCCESS: Pipeline completed for $gsm_id\n"
else
    printf "$(date +'%d/%b/%Y %H:%M:%S') | FAILED: Pipeline failed for $gsm_id (Check: $log_file)\n"
    exit 1
fi