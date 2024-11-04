#!/bin/bash
#input and output directories
input_dir="fastq_all"
output_dir="/media/urisp/Elements/uri_rna_feb24/fastqc_results"
mkdir -p "$output_dir"

#qc for each fastq file in input directory
for fastq_file in "$input_dir"/*.fq.gz; do
    filename=$(basename -- "$fastq_file")
    filename_no_ext="${filename%.*}"

    fastqc "$fastq_file" -o "$output_dir"
done
