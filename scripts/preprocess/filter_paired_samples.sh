#!/bin/bash
set -uo pipefail

INPUT_TSV="/users/dhan30/splicing_order/data/mouse_gsm_sra_list_for_pipeline_3t3.tsv"
OUTPUT_TSV="/users/dhan30/splicing_order/data/mouse_gsm_sra_list_paired_only.tsv"
FAILED_TSV="/users/dhan30/splicing_order/data/mouse_gsm_check_failed.tsv"
TEMP_DIR="/users/dhan30/scratch/splicing_order/data/temp_check"

mkdir -p "$TEMP_DIR"
module load miniconda3/23.11.0s
source /oscar/runtime/software/external/miniconda3/23.11.0/etc/profile.d/conda.sh
conda activate order_env

head -n 1 "$INPUT_TSV" > "$OUTPUT_TSV"
: > "$FAILED_TSV"

total=0; paired=0; single=0; failed=0

while IFS=$'\t' read -r gsm_id sra_ids; do
    gsm_id=$(echo "$gsm_id" | xargs)
    sra_ids=$(echo "$sra_ids" | tr -d '\r' | xargs)
    [ -z "$gsm_id" ] || [ -z "$sra_ids" ] && continue

    total=$((total + 1))
    first_srr=$(echo "$sra_ids" | cut -d',' -f1)
    info="$TEMP_DIR/fastq-run-info.tsv"

    # CRITICAL: clear before every call, not after
    rm -f "$info" "$TEMP_DIR/fastq-run-mergers.tsv"

    echo -n "Checking $gsm_id ($first_srr)... "
    fastq-dl --only-download-metadata --outdir "$TEMP_DIR" --accession "$first_srr" >/dev/null 2>&1

    if [ ! -s "$info" ]; then
        echo "FAILED (no metadata)"
        echo -e "${gsm_id}\t${sra_ids}" >> "$FAILED_TSV"
        failed=$((failed + 1))
        continue
    fi

    # Parse the library_layout column by name, check data rows only
    layout=$(awk -F'\t' '
        NR==1 { for (i=1; i<=NF; i++) if ($i == "library_layout") col=i; next }
        NR==2 && col { print toupper($col) }
    ' "$info")

    if [ "$layout" == "PAIRED" ]; then
        echo "PAIRED-END"
        echo -e "${gsm_id}\t${sra_ids}" >> "$OUTPUT_TSV"
        paired=$((paired + 1))
    elif [ -n "$layout" ]; then
        echo "single-end (skipping)"
        single=$((single + 1))
    else
        echo "FAILED (no library_layout column)"
        echo -e "${gsm_id}\t${sra_ids}" >> "$FAILED_TSV"
        failed=$((failed + 1))
    fi
done < <(tail -n +2 "$INPUT_TSV")

rm -rf "$TEMP_DIR"

echo ""
echo "Total:  $total"
echo "Paired: $paired"
echo "Single: $single"
echo "Failed: $failed  -> $FAILED_TSV"
echo "Output: $OUTPUT_TSV"