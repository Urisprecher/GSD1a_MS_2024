#!/bin/bash
#input output dirs
input_folder="trim_data"
output_folder="bam_files_ready"

mkdir -p "$output_folder"

#find R1 trimmed fastqs
for forward_read in "$input_folder"/*_forward_trimmed.fastq.gz
do

    base_name=$(basename "$forward_read" _forward_trimmed.fastq.gz)

    # find corresponding file
    reverse_read="${input_folder}/${base_name}_reverse_trimmed.fastq.gz"

    #run Bowtie2 to align paired end files 
    bowtie2 -x GRCh38_index -1 "$forward_read" -2 "$reverse_read" -S "$output_folder/${base_name}.sam"

    # sam to bam, sorting and indexing
    samtools view -bS "$output_folder/${base_name}.sam" | samtools sort -o "$output_folder/${base_name}.sorted.bam"
    samtools index "$output_folder/${base_name}.sorted.bam"

    # qc
    samtools flagstat "$output_folder/${base_name}.sorted.bam" > "$output_folder/${base_name}_mapping_stats.txt"

    #removw intermediate SAM files
    rm "$output_folder/${base_name}.sam"

    echo "Alignment and post-processing completed for $base_name"
done



