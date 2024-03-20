#!/bin/bash

# Set variables
GENOME_URL="http://example.com/ferret_genome.fa"
GENOME_FA="GCA_011764305.2_ASM1176430v1.1_genomic.fna"
GENOME_INDEX="ferret_genome_index"
ANNOTATION_URL="http://example.com/ferret_annotation.gff"
ANNOTATION_FILE= "ferret_genomic_annotation.gff" #"ferret_annotation.gff"
FASTQ_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNAseq_pipeline/arun_trimmed_fastq_files"
OUTPUT_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNAseq_pipeline/hisat2_output"
ALIGN_DIR="${OUTPUT_DIR}/alignments"
COUNTS_DIR="${OUTPUT_DIR}/counts"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$ALIGN_DIR" "$COUNTS_DIR"

# Step 1: Download the Ferret Genome
wget -O "${OUTPUT_DIR}/${GENOME_FA}" "$GENOME_URL"

# Step 2: Index the genome with HISAT2
hisat2-build "${OUTPUT_DIR}/${GENOME_FA}" "${OUTPUT_DIR}/${GENOME_INDEX}"

# Download the annotation file
wget -O "${OUTPUT_DIR}/${ANNOTATION_FILE}" "$ANNOTATION_URL"

# Step 3: Align reads with HISAT2 for each pair of FASTQ files
for file in "${FASTQ_DIR}"/*_1.fastq.gz; do
    base=$(basename "$file" _1.fastq.gz)
    echo "Processing $base"
    hisat2 -x "${OUTPUT_DIR}/${GENOME_INDEX}" -1 "$file" -2 "${FASTQ_DIR}/${base}_2.fastq.gz" \
           -S "${ALIGN_DIR}/${base}.sam"
    
    # Step 4: Convert SAM to BAM, sort, and index with Samtools
    samtools view -bS "${ALIGN_DIR}/${base}.sam" | samtools sort - -o "${ALIGN_DIR}/${base}_sorted.bam"
    samtools index "${ALIGN_DIR}/${base}_sorted.bam"
    
    # Step 5: Count reads with HTSeq
    htseq-count -f bam -r pos -s no -t exon -i gene_id "${ALIGN_DIR}/${base}_sorted.bam" "${OUTPUT_DIR}/${ANNOTATION_FILE}" > "${COUNTS_DIR}/${base}_counts.txt"
done

echo "Pipeline completed."
