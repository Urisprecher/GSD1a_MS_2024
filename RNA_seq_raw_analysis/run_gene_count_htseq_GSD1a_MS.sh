#!/bin/bash
#input output dirs
input_dir="dedup_bams_ready_for_count"
output_dir="htseq_counts_ready"
mkdir -p "$output_dir"
#annotation file-verify location
annotation_file="Homo_sapiens.GRCh38.112.gtf"
#track running jobs
running_jobs=0

#run each bam file
for bam_file in "$input_dir"/*.bam; do
    if [ -f "$bam_file" ]; then
        
        output_count="$output_dir/$(basename "$bam_file" .bam).counts.txt"
       
        #run
        htseq-count -f bam -r name -s no -i gene_id -t exon "$bam_file" "$annotation_file" > "$output_count"
       
    else
        echo "bam file not found: $bam_file"
    fi
done
