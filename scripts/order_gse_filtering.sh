#!/bin/bash

# SLURM directives (if running as standalone job)
#SBATCH --job-name=splicing_pipeline
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --time=36:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8

################################################################################
# CONFIGURATION
################################################################################

# Inputs
threads=$1
out_dir=$2
gsm_id=$3
srr_id_input=$4

printf "$(date +'%d/%b/%Y %H:%M:%S') | Processing sample %s...\n" "${gsm_id}"
printf "$(date +'%d/%b/%Y %H:%M:%S') | Skip structure analysis: %s\n" "${SKIP_STRUCTURE}"

# Raw & trimmed files (PAIRED-END ONLY)
paired_read_one_file="${out_dir}/${gsm_id}_all_reads_1.fq.gz"
paired_read_two_file="${out_dir}/${gsm_id}_all_reads_2.fq.gz"
paired_trimmed_one_file="${out_dir}/${gsm_id}_trimmed_1.fq.gz"
paired_trimmed_two_file="${out_dir}/${gsm_id}_trimmed_2.fq.gz"

# STAR alignment
star_out_prefix="${out_dir}/${gsm_id}_STAR_"
aligned_bam_file="${star_out_prefix}Aligned.sortedByCoord.out.bam"

# Intermediate files
TEMP_PREFIX="temp_${gsm_id}"
GSM_PREFIX="${gsm_id}"
ALL_PAIRED_SAM="${out_dir}/${GSM_PREFIX}_all_paired_reads.sam"
PAIRED_JUNCS_SAM="${out_dir}/${GSM_PREFIX}_junction_paired_reads.sam"
INTRON_OVERLAP_BED="${out_dir}/${GSM_PREFIX}_intron_reads.bed"
SORTED_BAM_PREFIX="${out_dir}/${GSM_PREFIX}_intermediate_reads_sorted"

# Final outputs
informative_pairs_bam="${out_dir}/${gsm_id}_informative_pairs.bam"
splicing_order_output="${out_dir}/${gsm_id}_pairwise_splicing_order.tsv"

# Reference files
reference_genome_dir="/users/dhan30/reference/mm39"
reference_genome_fasta="/users/dhan30/reference/mm39.fa"
intron_bed_file="/users/dhan30/reference/mm39.gencode.basic.vM36.introns.bed.gz"

# Script directory
SCRIPT_DIR="/users/dhan30/splicing_order/scripts"

# Ensure output directory exists
mkdir -p "${out_dir}"

################################################################################
# ENVIRONMENT SETUP
################################################################################

printf "$(date +'%d/%b/%Y %H:%M:%S') | Loading conda environment...\n"
module load miniforge3/25.3.0-3
eval "$(conda shell.bash hook)"
conda activate order_env
export PATH="/users/dhan30/.conda/envs/order_env/bin:$PATH"

################################################################################
# CHECKPOINT 1: Check if final outputs already exist
################################################################################

if [ -f "$splicing_order_output" ] && [ -s "$splicing_order_output" ]; then
    line_count=$(wc -l < "$splicing_order_output")
    
    if [ "$line_count" -gt 1 ]; then
        printf "$(date +'%d/%b/%Y %H:%M:%S') | Splicing order output exists with $((line_count - 1)) intron pairs.\n"
        # Continue to structure analysis if not done
        skip_to_structure=true
    else
        printf "$(date +'%d/%b/%Y %H:%M:%S') | Output file is empty (header only). Reprocessing...\n"
        rm -f "$splicing_order_output"
        skip_to_structure=false
    fi
else
    skip_to_structure=false
fi

################################################################################
# CHECKPOINT 2: Check if alignment exists  
################################################################################

if [ "$skip_to_structure" = false ]; then
    if [ -f "$aligned_bam_file" ] && [ -s "$aligned_bam_file" ]; then
        printf "$(date +'%d/%b/%Y %H:%M:%S') | STAR alignment file found. Skipping to informative pairs extraction...\n"
        skip_to_extraction=true
    else
        skip_to_extraction=false
    fi
fi

################################################################################
# CHECKPOINT 3: Check if trimmed files exist
################################################################################

if [ "$skip_to_structure" = false ] && [ "$skip_to_extraction" = false ]; then
    if [ -f "$paired_trimmed_one_file" ] && [ -f "$paired_trimmed_two_file" ]; then
        printf "$(date +'%d/%b/%Y %H:%M:%S') | Trimmed files found. Skipping download and trimming...\n"
        skip_to_alignment=true
    else
        skip_to_alignment=false
    fi
fi

################################################################################
# DOWNLOAD & TRIM (only if we need to align)
################################################################################

if [ "$skip_to_structure" = false ] && [ "$skip_to_extraction" = false ] && [ "$skip_to_alignment" = false ]; then

    # CHECKPOINT: Check if raw reads exist
    if [ -f "$paired_read_one_file" ] && [ -f "$paired_read_two_file" ] && [ -s "$paired_read_one_file" ]; then
        printf "$(date +'%d/%b/%Y %H:%M:%S') | Raw paired reads found. Skipping download...\n"
    else
        printf "$(date +'%d/%b/%Y %H:%M:%S') | Downloading reads for %s...\n" "${gsm_id}"
        IFS=',' read -r -a srr_id_list <<< "$srr_id_input"

        for srr in "${srr_id_list[@]}"
        do
            printf "$(date +'%d/%b/%Y %H:%M:%S') |   Downloading %s...\n" "${srr}"
            
            fastq-dl -F --outdir "${out_dir}" --cpus "${threads}" --accession "${srr}" ||
                { printf "ERROR: fastq-dl failed for %s\n" "${srr}"; exit 1; }

            srr_output_prefix="${out_dir}/${srr}"

            # Create paired-end read files by concatenating if they exist
            if [ -f "${srr_output_prefix}_1.fastq.gz" ] && [ -f "${srr_output_prefix}_2.fastq.gz" ]; then
                cat "${srr_output_prefix}_1.fastq.gz" >> "$paired_read_one_file"
                cat "${srr_output_prefix}_2.fastq.gz" >> "$paired_read_two_file"
                rm "${srr_output_prefix}_1.fastq.gz" "${srr_output_prefix}_2.fastq.gz"
            else
                printf "ERROR: Expected paired-end data for %s but didn't find _1 and _2 files.\n" "${srr}"
                printf "This pipeline only supports paired-end RNA-seq data.\n"
                
                # Log failure reason
                echo "FAILURE: Single-end data detected for ${srr}" > "${out_dir}/FAILED_REASON.txt"
                echo "Date: $(date)" >> "${out_dir}/FAILED_REASON.txt"
                echo "SRA ID: ${srr}" >> "${out_dir}/FAILED_REASON.txt"
                
                # Delete everything except logs
                find "${out_dir}" -type f ! -name "*.log" ! -name "FAILED_REASON.txt" -delete
                
                exit 1
            fi

            rm -f "${out_dir}/fastq-run-info.tsv" "${out_dir}/fastq-run-mergers.tsv"
        done

        printf "$(date +'%d/%b/%Y %H:%M:%S') | Download complete for %s.\n" "${gsm_id}"
    fi

    # Validate we have paired-end data
    if [ ! -f "$paired_read_one_file" ] || [ ! -f "$paired_read_two_file" ]; then
        printf "ERROR: Paired-end read files not found for %s\n" "${gsm_id}"
        printf "This pipeline requires paired-end RNA-seq data.\n"
        printf "Cleaning up incomplete data for sample %s...\n" "${gsm_id}"
        
        # Clean up the output directory
        rm -rf "${out_dir}"
        
        exit 1
    fi

    # Trimming with fastp
    
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Trimming reads for %s...\n" "${gsm_id}"
    fastp_report_json="${out_dir}/${gsm_id}_fastp.json"
    fastp_report_html="${out_dir}/${gsm_id}_fastp.html"

    fastp \
        --in1 "$paired_read_one_file" \
        --in2 "$paired_read_two_file" \
        --out1 "$paired_trimmed_one_file" \
        --out2 "$paired_trimmed_two_file" \
        --detect_adapter_for_pe \
        --thread "$threads" \
        --json "$fastp_report_json" \
        --html "$fastp_report_html" ||
        { printf "ERROR: fastp failed for %s\n" "${gsm_id}"; 
          rm -f "$paired_trimmed_one_file" "$paired_trimmed_two_file"; 
          exit 1; }

    # Verify trimmed files are not corrupted
    if ! gzip -t "$paired_trimmed_one_file" 2>/dev/null; then
        printf "ERROR: Corrupted trimmed file 1 for %s\n" "${gsm_id}"
        rm -f "$paired_trimmed_one_file" "$paired_trimmed_two_file"
        exit 1
    fi
    if ! gzip -t "$paired_trimmed_two_file" 2>/dev/null; then
        printf "ERROR: Corrupted trimmed file 2 for %s\n" "${gsm_id}"
        rm -f "$paired_trimmed_one_file" "$paired_trimmed_two_file"
        exit 1
    fi

    rm "$paired_read_one_file" "$paired_read_two_file"
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Trimming complete for %s.\n" "${gsm_id}"
    
fi

################################################################################
# ALIGNMENT (only if we don't have STAR BAM)
################################################################################

if [ "$skip_to_structure" = false ] && [ "$skip_to_extraction" = false ]; then
    # Check if trimmed files exist before attempting alignment
    if [ ! -f "$paired_trimmed_one_file" ] || [ ! -f "$paired_trimmed_two_file" ]; then
        printf "ERROR: Cannot run alignment - trimmed files not found\n"
        exit 1
    fi
    
    # STAR Alignment
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Aligning reads for %s using STAR...\n" "${gsm_id}"

    STAR --genomeDir "$reference_genome_dir" \
        --readFilesIn "$paired_trimmed_one_file" "$paired_trimmed_two_file" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$star_out_prefix" \
        --outSAMtype BAM SortedByCoordinate \
        --runThreadN "$threads" \
        --twopassMode Basic ||
        { printf "ERROR: STAR alignment failed for %s\n" "${gsm_id}"; exit 1; }

    if [ ! -f "$aligned_bam_file" ]; then
        printf "ERROR: STAR alignment failed for %s. BAM file not found.\n" "${gsm_id}"
        exit 1
    fi

    printf "$(date +'%d/%b/%Y %H:%M:%S') | STAR alignment complete for %s.\n" "${gsm_id}"
    
    # Clean up trimmed fastq files after successful alignment
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Cleaning up trimmed fastq files...\n"
    rm -f "$paired_trimmed_one_file" "$paired_trimmed_two_file"
    
    # Clean up STAR intermediate files
    rm -f "${star_out_prefix}"Log.out "${star_out_prefix}"Log.progress.out
    rm -f "${star_out_prefix}"SJ.out.tab
fi  # End of alignment block

################################################################################
# EXTRACT INFORMATIVE PAIRS (pairs spanning exon-intron junction and intron overlaps)
################################################################################

if [ "$skip_to_structure" = false ]; then
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Extracting informative pairs...\n"

    # convert to name-sorted SAM
    printf "$(date +'%d/%b/%Y %H:%M:%S') |   Converting to name-sorted SAM...\n"
    samtools sort -n -O SAM -@ "$threads" -o "${ALL_PAIRED_SAM}" "${aligned_bam_file}"

    # extract junction-containing pairs
    printf "$(date +'%d/%b/%Y %H:%M:%S') |   Extracting junction pairs...\n"
    awk 'BEGIN {OFS="\t"}
        /^@/ {print; next}
        $6 ~ /N/ {junction_reads[$1]=1}
        END {
            while ((getline < "'"${ALL_PAIRED_SAM}"'") > 0) {
                if ($0 ~ /^@/) continue
                if ($1 in junction_reads) print
            }
        }' "${ALL_PAIRED_SAM}" > "${PAIRED_JUNCS_SAM}"

    rm "${ALL_PAIRED_SAM}"

    # find intron-overlapping reads
    printf "$(date +'%d/%b/%Y %H:%M:%S') |   Finding intron overlaps (f=0.05)...\n"
    samtools view -bh "${PAIRED_JUNCS_SAM}" | \
        bedtools bamtobed -split -i stdin | \
        bedtools intersect -f 0.05 -a stdin -b "${intron_bed_file}" > "${INTRON_OVERLAP_BED}"

    # extract pairs where one has junction and one overlaps intron
    printf "$(date +'%d/%b/%Y %H:%M:%S') |   Extracting informative pairs...\n"
    awk '{print $4}' "${INTRON_OVERLAP_BED}" | \
        sed 's/\/[12]$//' | \
        sort -u > "${out_dir}/${TEMP_PREFIX}_intron_read_names.txt"

    (grep "^@" "${PAIRED_JUNCS_SAM}"; \
        grep -wFf "${out_dir}/${TEMP_PREFIX}_intron_read_names.txt" "${PAIRED_JUNCS_SAM}") | \
        samtools view -bh - | \
        samtools fixmate -m - - | \
        samtools sort -@ "$threads" -o "${SORTED_BAM_PREFIX}.bam" - ||
        { printf "ERROR: Failed to create sorted BAM\n"; exit 1; }

    # remove duplicates
    printf "$(date +'%d/%b/%Y %H:%M:%S') |   Removing duplicates...\n"
    samtools markdup -r -@ "$threads" "${SORTED_BAM_PREFIX}.bam" "${informative_pairs_bam}"
    samtools index "${informative_pairs_bam}"

    # cleanup intermediate files to save space
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Cleaning up intermediate files...\n"
    rm -f "${PAIRED_JUNCS_SAM}" "${INTRON_OVERLAP_BED}" 
    rm -f "${SORTED_BAM_PREFIX}.bam" "${out_dir}/${TEMP_PREFIX}_intron_read_names.txt"

    # delete STAR alignment BAM (big files)
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Deleting STAR alignment to save disk space...\n"
    rm -f "${aligned_bam_file}" "${aligned_bam_file}.bai"

    printf "$(date +'%d/%b/%Y %H:%M:%S') | Informative pairs extraction complete.\n"
    printf "  Total reads: $(samtools view -c ${informative_pairs_bam})\n"
    printf "  Read pairs: $(($(samtools view -c ${informative_pairs_bam}) / 2))\n"
fi

################################################################################
# ANALYZE SPLICING ORDER
################################################################################

if [ "$skip_to_structure" = false ]; then
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Analyzing splicing order...\n"

    printf "$(date +'%d/%b/%Y %H:%M:%S') | Ensuring conda environment is active...\n"
    module load miniforge3/25.3.0-3
    eval "$(conda shell.bash hook)"
    conda activate order_env
    export PATH="/users/dhan30/.conda/envs/order_env/bin:$PATH"

    printf "$(date +'%d/%b/%Y %H:%M:%S') | Active environment: $CONDA_DEFAULT_ENV\n"
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Python location: $(which python3)\n"

    python3 ${SCRIPT_DIR}/analyze_splicing_order.py \
        --bam "${informative_pairs_bam}" \
        --intron-bed "${intron_bed_file}" \
        --output "${splicing_order_output}" \
        --min-reads 10 \
        --min-mapq 10 \
        --tolerance 10 \
        --max-intron-length 15000 ||
        { printf "ERROR: Splicing order analysis failed\n"; exit 1; }

     # delete informative pairs BAM to save disk space, just keep tsv
    printf "$(date +'%d/%b/%Y %H:%M:%S') | Deleting informative pairs BAM to save disk...\n"
    rm -f "${informative_pairs_bam}" "${informative_pairs_bam}.bai"

    printf "$(date +'%d/%b/%Y %H:%M:%S') | Splicing order analysis complete.\n"

    n_pairs=$(tail -n +2 "${splicing_order_output}" | wc -l)
    # cleanup
    rm -f "$paired_read_one_file" "$paired_read_two_file" 2>/dev/null
    rm -f "$paired_trimmed_one_file" "$paired_trimmed_two_file" 2>/dev/null
    rm -f "${star_out_prefix}"*.out "${star_out_prefix}"*.tab 2>/dev/null
    rm -f "${out_dir}"/temp_* 2>/dev/null

    conda deactivate
    exit 0

    printf "$(date +'%d/%b/%Y %H:%M:%S') | Splicing order analysis complete.\n"
fi 

################################################################################
# FINAL SUMMARY
################################################################################

printf "$(date +'%d/%b/%Y %H:%M:%S') | Sample %s processing complete!\n" "${gsm_id}"
printf "$(date +'%d/%b/%Y %H:%M:%S') | Splicing order results: %s\n" "${splicing_order_output}"
printf "$(date +'%d/%b/%Y %H:%M:%S') | Structure features: %s\n" "${structure_features_output}"

# Count results
if [ -f "${splicing_order_output}" ]; then
    n_pairs=$(tail -n +2 "${splicing_order_output}" | wc -l)
    printf "  Intron pairs analyzed: %s\n" "$n_pairs"
fi

# Cleanup intermediate files to save space
printf "$(date +'%d/%b/%Y %H:%M:%S') | Final cleanup...\n"
rm -f "$paired_read_one_file" "$paired_read_two_file" 2>/dev/null
rm -f "$paired_trimmed_one_file" "$paired_trimmed_two_file" 2>/dev/null
rm -f "${star_out_prefix}"*.out "${star_out_prefix}"*.tab 2>/dev/null
rm -f "${out_dir}"/temp_* 2>/dev/null

printf "$(date +'%d/%b/%Y %H:%M:%S') | All processing complete for %s!\n" "${gsm_id}"

conda deactivate