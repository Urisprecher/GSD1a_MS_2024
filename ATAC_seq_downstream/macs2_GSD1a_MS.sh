#!/bin/bash

#dirs
MACS2_PATH="macs2_dir"  
INPUT_BAM_FOLDER="peak_calling_data"   
OUTPUT_FOLDER="ready_peaks_all   


mkdir -p $OUTPUT_FOLDER

# peak call for each bam in input folder
for BAM_FILE in $INPUT_BAM_FOLDER/*.bam; do
    FILENAME=$(basename "$BAM_FILE" .bam)
   
    #RUN MACS2
    $MACS2_PATH callpeak -t "$BAM_FILE" -f BAMPE --outdir "$OUTPUT_FOLDER/$FILENAME" -n "$FILENAME"_peaks
   
   
    echo "finished processing $FILENAME"
done

echo "peak calling for all files completed!"

