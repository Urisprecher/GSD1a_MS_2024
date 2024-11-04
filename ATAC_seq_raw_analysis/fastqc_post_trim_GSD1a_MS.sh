#!/bin/bash
#input output dirs
input_dir="trimmed_fastq_results"
output_dir="/media/urisp/Elements/uri_rna_feb24/fastqc_trim_results"
mkdir -p "$output_dir"

#run fastqc-each fastq in input
for fastq_file in "$input_dir"/*.fastq.gz; do
    filename=$(basename -- "$fastq_file")
    filename_no_ext="${filename%.*}"


    fastqc "$fastq_file" -o "$output_dir"
done
##############
