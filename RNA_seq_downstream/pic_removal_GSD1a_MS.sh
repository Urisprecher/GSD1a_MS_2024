#!/bin/bash
#input dir with aligned and sorted bams
input_dir="merged_bam_all"
#output dir for deduplicated bams+qc
output_dir="Res_dedup_ready_for_count"
qc_dir="Res_ready-qc"
#create dirs
mkdir -p "$output_dir"
mkdir -p "$qc_dir"

#picard tools dir - pre install picard and locate 
picard_dir="picard_analysis"
#run each bam
for bam_file in "$input_dir"/*.bam; do
    #base name of file
    filename=$(basename -- "$bam_file")
    sample_name="${filename%.bam}"

    output_bam="$output_dir/${sample_name}_dedup.bam"

    #run
    java -jar "$picard_dir/picard.jar" MarkDuplicates \
        INPUT="$bam_file" \
        OUTPUT="$output_bam" \
        METRICS_FILE="$output_dir/${sample_name}_duplication_metrics.txt" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true

    #qc for deduplicated bam
    qc_metrics="$qc_dir/${sample_name}_deduplication_qc.txt"
    samtools flagstat "$output_bam" > "$qc_metrics"
done
