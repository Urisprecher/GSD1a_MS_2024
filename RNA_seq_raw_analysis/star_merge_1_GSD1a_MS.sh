#!/bin/bash
#input output qc dirs
input_dir="star_Res_aligned"
output_dir="star_Res_merged_sorted"
qc_dir="star_Res_sorted_merged-qc"
#create output and qc dirs
mkdir -p "$output_dir"
mkdir -p "$qc_dir"
##run over each unique sample identifier and merge aligned bams
for sample_id in $(ls "$input_dir"/*_Aligned.out_sorted.bam | sed -n 's/.*\(s[0-9][0-9]*_mRNA_[T0-9]\+\)_Aligned.out_sorted.bam/\1/p' | sort -u); do
    echo "merging bam files for sample $sample_id...!"


    merged_bam="$output_dir/${sample_id}_merged.bam"


    samtools merge "$merged_bam" $(ls "$input_dir"/*"${sample_id}"*_Aligned.out_sorted.bam)

    #index
    samtools index "$merged_bam"

    # qc
    output_qc="$qc_dir/${sample_id}_alignment_qc.txt"
    samtools flagstat "$merged_bam" > "$output_qc"

    echo "merging and qc completed for sample $sample_id.:)"
done

echo "All samples processed successfully!"

