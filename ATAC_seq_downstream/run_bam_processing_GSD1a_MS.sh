#!/bin/bash
#dirs
input_folder="bam_files_ready"
output_folder="peak_calling_data"
#mito fa file
mito_reference="mitochondrial_reference.fa"
mkdir -p "$output_folder"

#remove mtDNA
for bam_file in "$input_folder"/*.bam
do

    base_name=$(basename "$bam_file" .bam)

    #temp bed
    bed_file="$output_folder/${base_name}_mito_regions.bed"
    samtools faidx "$mito_reference"
    grep -v '^>' "${mito_reference}.fai" | awk -v OFS='\t' '{print $1, 0, $2}' > "$bed_file"

    #remove mtDNA reads
    samtools view -b -h -L "$bed_file" "$bam_file" > "$output_folder/${base_name}_no_mt.bam"

    # remove bed file
    rm "$bed_file"

    echo "mitochondrial reads removed for $base_name"

    # PCR-duplicated reads removal using Picard 
    edited_bam_file="$output_folder/${base_name}_no_mt.bam"
    base_name_no_mt=$(basename "$edited_bam_file" _no_mt.bam)
    
    java -jar picard.jar MarkDuplicates INPUT="$edited_bam_file" OUTPUT="$output_folder/${base_name_no_mt}_no_dups.bam" METRICS_FILE="$output_folder/${base_name_no_mt}_dups_metrics.txt" REMOVE_DUPLICATES=true

    echo "PCR-duplicates removed for $base_name_no_mt"


    #qc
    samtools flagstat "$output_folder/${base_name_no_mt}_no_dups.bam" > "$output_folder/${base_name_no_mt}_alignment_stats.txt"

    #mapping qc
    samtools view -H "$output_folder/${base_name_no_mt}_no_dups.bam" 
    samtools view "$output_folder/${base_name_no_mt}_no_dups.bam" | head -n 10 | cut -f 5 > "$output_folder/${base_name_no_mt}_mapping_quality.txt"


    samtools depth "$output_folder/${base_name_no_mt}_no_dups.bam" | awk '{sum+=$3} END {print "average Coverage:", sum/NR}' > "$output_folder/${base_name_no_mt}_coverage_stats.txt"

    echo "alignment statistics and assessment completed for $base_name_no_mt"
done
