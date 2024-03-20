#!/bin/bash

SECONDS=0

# Define paths
FASTQ_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNApipeline_CPT/CPT_trimmed_fastq_files_v0"
TRIMMED_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNApipeline_CPT/CPT_trimmed_fastq_files_v1"
FASTQC_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNApipeline_CPT/CPT_fastqc_reports_v1"
ADAPTERS_PATH="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNApipeline_CPT/adapters/TruSeq3-PE.fa"

# Create directories if they don't exist
mkdir -p "$TRIMMED_DIR" "$FASTQC_DIR"

# Define a processing function
process_sample() {
    sample=$(basename "$1" _paired_R1.fastq.gz)
    echo "Processing sample: $sample"

    # Run Trimmomatic
    trimmomatic PE -phred33 \
    "${FASTQ_DIR}/${sample}_paired_R1.fastq.gz" \
    "${FASTQ_DIR}/${sample}_paired_R2.fastq.gz" \
    "${TRIMMED_DIR}/${sample}_paired2_R1.fastq.gz" \
    "${TRIMMED_DIR}/${sample}_unpaired2_R1.fastq.gz" \
    "${TRIMMED_DIR}/${sample}_paired2_R2.fastq.gz" \
    "${TRIMMED_DIR}/${sample}_unpaired2_R2.fastq.gz" \
    ILLUMINACLIP:"${ADAPTERS_PATH}":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    # Run FastQC
    fastqc -o "${FASTQC_DIR}" \
    "${TRIMMED_DIR}/${sample}_paired2_R1.fastq.gz" \
    "${TRIMMED_DIR}/${sample}_paired2_R2.fastq.gz"
}

# Export the function and necessary variables to be used by GNU Parallel
export -f process_sample
export FASTQ_DIR TRIMMED_DIR FASTQC_DIR ADAPTERS_PATH

# Find FASTQ files and run the processing in parallel
find "$FASTQ_DIR" -name "*_paired_R1.fastq.gz" | parallel --jobs 12 process_sample

echo "Processing complete."

duration=$SECONDS
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."
