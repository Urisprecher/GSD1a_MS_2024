##!/bin/bash
#input and output dirs 
input_folder="fastq_all_true"
output_folder="trim_data"


mkdir -p "$output_folder"

#find all R1 fastqs
for forward_read in "$input_folder"/*_R1.fastq.gz
do
# names R1
    base_name=$(basename "$forward_read" _R1.fastq.gz)

    #find corresponding R2
    reverse_read="${input_folder}/${base_name}_R2.fastq.gz"

    #run
    java -jar trimmomatic-0.39.jar PE "$forward_read" "$reverse_read" \
        "$output_folder/${base_name}_forward_trimmed.fastq.gz" "$output_folder/${base_name}_forward_unpaired.fastq.gz" \
        "$output_folder/${base_name}_reverse_trimmed.fastq.gz" "$output_folder/${base_name}_reverse_unpaired.fastq.gz" \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo "Trimmomatic completed for $base_name"
done


