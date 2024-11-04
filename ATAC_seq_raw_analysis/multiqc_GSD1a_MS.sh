#input output dirs
input_dir="/media/urisp/Elements/uri_rna_feb24/fastqc_results_multi"
output_dir="multiqc_res"
mkdir -p "$output_dir"
#run
multiqc "$input_dir" -o "$output_dir"


