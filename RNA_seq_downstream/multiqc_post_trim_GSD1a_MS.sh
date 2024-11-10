#input  output dirs
input_dir="/media/urisp/Elements/uri_rna_feb24/fastqc_trim_results"
output_dir="multiqc_res_trim"
mkdir -p "$output_dir"

#run 
multiqc "$input_dir" -o "$output_dir"
