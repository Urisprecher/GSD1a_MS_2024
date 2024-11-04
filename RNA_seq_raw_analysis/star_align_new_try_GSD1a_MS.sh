#!/bin/bash
#variables and dirs
STAR_INDEX="custom_star_index" #locate index
INPUT_DIR="star_fastq_unzipped"
OUTPUT_DIR="star_aligned" 
THREADS=8 #threads to use.
#create out dir
mkdir -p "$OUTPUT_DIR"

#run over paired-end fastq files
for file in "$INPUT_DIR"/*_trimmed_forward.fastq; do
    #samples name
    sample=$(basename "$file" _trimmed_forward.fastq)

    # Run STAR alignment
    STAR \
    	--runThreadN "$THREADS" \
    	--genomeDir "$STAR_INDEX"   \
    	--readFilesIn  "$INPUT_DIR/${sample}_trimmed_forward.fastq" "$INPUT_DIR/${sample}_trimmed_reverse.fastq"  \
    	--outFileNamePrefix "$OUTPUT_DIR/${sample}_" \
    	--outSAMtype BAM Unsorted \
    	--quantMode TranscriptomeSAM GeneCounts

 
done
    	




