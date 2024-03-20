#!/bin/bash

SECONDS=0

# Input, output, and FastQC directories
INPUT_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNApipeline_CPT/CPT_trimmed_fastq_files_test"
OUTPUT_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNApipeline_CPT/CPT_trimmed_fastq_files_v2"
FASTQC_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNApipeline_CPT/CPT_fastqc_reports_test"
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FASTQC_DIR"

# Trimmomatic and FastQC commands wrapped in a function
process_files() {
  file=$1
  INPUT_DIR=$2
  OUTPUT_DIR=$3
  FASTQC_DIR=$4
  PHRED="-phred33"
  OPTIONS="HEADCROP:3"
  
  base=$(basename "$file" _paired2_R1.fastq.gz)
  forward_in="$file"
  reverse_in="${INPUT_DIR}/${base}_paired2_R2.fastq.gz"
  forward_paired_out="${OUTPUT_DIR}/${base}_paired3_R1.fastq.gz"
  forward_unpaired_out="${OUTPUT_DIR}/${base}_unpaired3_R1.fastq.gz"
  reverse_paired_out="${OUTPUT_DIR}/${base}_paired3_R2.fastq.gz"
  reverse_unpaired_out="${OUTPUT_DIR}/${base}_unpaired3_R2.fastq.gz"

  # Run Trimmomatic
  trimmomatic PE \
    $PHRED \
    $forward_in $reverse_in \
    $forward_paired_out $forward_unpaired_out \
    $reverse_paired_out $reverse_unpaired_out \
    $OPTIONS

  # FastQC quality check for the trimmed files
  fastqc -o "$FASTQC_DIR" "$forward_paired_out" "$reverse_paired_out" "$forward_unpaired_out" "$reverse_unpaired_out"
}

export -f process_files

# Use GNU Parallel to run the process_files function in parallel
find "$INPUT_DIR" -name "*_paired2_R1.fastq.gz" | parallel -j 12 process_files {} "$INPUT_DIR" "$OUTPUT_DIR" "$FASTQC_DIR"

echo "Trimming and FastQC analysis complete."

duration=$SECONDS
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."

