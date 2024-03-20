#!/bin/bash

SECONDS=0

# Directories for downloads and processing
WORK_DIR="/Users/seidmuhie/Documents/02052024_RNAseqDataFolder/RNApipeline_CPT"
GENOME_DIR="${WORK_DIR}/grch38"
ANNOTATION_DIR="${WORK_DIR}/grch38_annotation"
FASTQ_DIR="${WORK_DIR}/CPT_trimmed_fastq_files_v2"
OUTPUT_DIR="${WORK_DIR}/CPT_hisat2_output"
ALIGN_DIR="${OUTPUT_DIR}/CPT_alignments"
COUNTS_DIR="${OUTPUT_DIR}/CPT_counts"

# Create necessary directories
mkdir -p "$GENOME_DIR" "$ANNOTATION_DIR" "$FASTQ_DIR" "$OUTPUT_DIR" "$ALIGN_DIR" "$COUNTS_DIR"

# URLs for GRCh38 primary assembly and GENCODE annotation
GENOME_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz"
ANNOTATION_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz"

# Download GRCh38 primary assembly and GENCODE annotation
wget -c $GENOME_URL -O "$GENOME_DIR/GRCh38.primary_assembly.genome.fa.gz"
wget -c $ANNOTATION_URL -O "$ANNOTATION_DIR/gencode.v38.primary_assembly.annotation.gtf.gz"

# Decompress downloaded files
gunzip -c "$GENOME_DIR/GRCh38.primary_assembly.genome.fa.gz" > "$GENOME_DIR/GRCh38.primary_assembly.genome.fa"
gunzip -c "$ANNOTATION_DIR/gencode.v38.primary_assembly.annotation.gtf.gz" > "$ANNOTATION_DIR/gencode.v38.primary_assembly.annotation.gtf"

# Index the genome with HISAT2 using 12 cores
GENOME_FA="${GENOME_DIR}/GRCh38.primary_assembly.genome.fa"
GENOME_INDEX="${GENOME_DIR}/grch38_genome_index"
hisat2-build -p 12 "${GENOME_FA}" "${GENOME_INDEX}"

# RNA-seq alignment and analysis
for file in "${FASTQ_DIR}"/*_paired3_R1.fastq.gz; do
    base=$(basename "$file" "_paired3_R1.fastq.gz") # Corrected basename extraction
    echo "Processing $base"
    hisat2 -p 12 -x "${GENOME_INDEX}" -1 "$file" -2 "${FASTQ_DIR}/${base}_paired3_R2.fastq.gz" \
           -S "${ALIGN_DIR}/${base}.sam"
    
    # Convert SAM to BAM, sort, and index with Samtools
    samtools view -@ 12 -bS "${ALIGN_DIR}/${base}.sam" | samtools sort -@ 12 - -o "${ALIGN_DIR}/${base}_sorted.bam"
    samtools index "${ALIGN_DIR}/${base}_sorted.bam"
    
    # Count reads with HTSeq
    htseq-count -f bam -r pos -s no -t exon -i gene_id "${ALIGN_DIR}/${base}_sorted.bam" "${ANNOTATION_DIR}/gencode.v38.primary_assembly.annotation.gtf" > "${COUNTS_DIR}/${base}_counts.txt"
done

echo "Pipeline completed."

duration=$SECONDS
echo "$((duration / 60)) minutes and $((duration % 60)) seconds elapsed."

