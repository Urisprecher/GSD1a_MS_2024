#input output dirs
input_dir="fastq_all"
output_dir="/media/urisp/Elements/uri_rna_feb24/trimmed_fastq_results"
mkdir -p "$output_dir"

#parameters
trimmomatic_path="trimmomatic-0.39.jar"
adapter_file="adapters/TruSeq3-PE-2.fa"
min_length="36"

#iterating over forward fastq files
for forward_file in "$input_dir"/*_1.fq.gz; do
    #coresponding reverse files
    reverse_file="${forward_file/_1/_2}"

    #file name
    filename=$(basename -- "$forward_file")
    filename_no_ext="${filename%_1.fq.gz}"

    #output  files path
    output_forward="$output_dir/${filename_no_ext}_trimmed_forward.fastq.gz"
    output_reverse="$output_dir/${filename_no_ext}_trimmed_reverse.fastq.gz"
    output_forward_unpaired="$output_dir/${filename_no_ext}_unpaired_forward.fastq.gz"
    output_reverse_unpaired="$output_dir/${filename_no_ext}_unpaired_reverse.fastq.gz"

    #run
    java -jar "$trimmomatic_path" PE -phred33 "$forward_file" "$reverse_file" \
        "$output_forward" "$output_forward_unpaired" \
        "$output_reverse" "$output_reverse_unpaired" \
        ILLUMINACLIP:"$adapter_file":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:"$min_length"
done


