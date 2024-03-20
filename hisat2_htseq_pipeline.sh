#!/bin/bash

# Define paths
GENOME_FA="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNAseq_pipeline/ferret_genome/GCF_011764305.1_ASM1176430v1.1_genomic.fna"
GENOME_INDEX="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNAseq_pipeline/hisat2_output/ferret_genome_index"
ANNOTATION_FILE="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNAseq_pipeline/ferret_genome/genomic.gtf"
FASTQ_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNAseq_pipeline/arun_trimmed_fastq_files"
OUTPUT_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNAseq_pipeline/hisat2_output_unpaired"
ALIGN_DIR="${OUTPUT_DIR}/alignments"
COUNTS_DIR="${OUTPUT_DIR}/counts"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$ALIGN_DIR" "$COUNTS_DIR"

# Step 1: Index the genome with HISAT2
hisat2-build "$GENOME_FA" "$GENOME_INDEX"

# Step 2: Align reads with HISAT2 for each pair of FASTQ files
for file in "$FASTQ_DIR"/*_R1_unpaired.fastq; do
    base=$(basename "$file" _R1_unpaired.fastq)
    echo "Processing $base"
    hisat2 -x "$GENOME_INDEX" -1 "$file" -2 "${FASTQ_DIR}/${base}_R2_unpaired.fastq" \
           -S "${ALIGN_DIR}/${base}.sam"
    
    # Step 3: Convert SAM to BAM, sort, and index with Samtools
    samtools view -bS "${ALIGN_DIR}/${base}.sam" | samtools sort -o "${ALIGN_DIR}/${base}_sorted.bam"
    samtools index "${ALIGN_DIR}/${base}_sorted.bam"
    
    # Step 4: Count reads with HTSeq
    htseq-count -f bam -r pos -s no -t exon -i gene_id "${ALIGN_DIR}/${base}_sorted.bam" "$ANNOTATION_FILE" > "${COUNTS_DIR}/${base}_counts.txt"
done

echo "Pipeline completed."

