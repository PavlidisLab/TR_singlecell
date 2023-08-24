#!/bin/bash

# This script takes an input file of IDS, and creates a copy of the respective
# preprocessing script save for replacing the count normalization step to CPM

input_file=$1
input_dir="/home/amorin/Projects/TR_singlecell/R/preprocessing_scripts/"
output_dir="/home/amorin/Projects/TR_singlecell/R/preprocessing_scripts/CPM"

mkdir -p "$output_dir"


while IFS=$'\t' read -r id species; do
    
    input_filename="${input_dir}${id}.R"
    output_filename="${output_dir}/${id}_cpm.R"
    
    # Load the content of the input file
    content=$(<"$input_filename")

    modified1="${content//_clean_mat_and_meta.RDS/_clean_mat_and_meta_CPM.RDS}"
    modified2="${modified1//_RSR_allrank.tsv/_RSR_allrank_CPM.tsv}"
    modified3="${modified2//_NA_mat.tsv/_NA_mat_CPM.tsv}"
    modified4="${modified3//Seurat::LogNormalize(., verbose = FALSE)/Seurat::NormalizeData(., normalization.method = "RC", scale.factor = 1e6, verbose = FALSE)}"

    # Save the modified content to the output file
    echo "$modified4" > "$output_filename"
    
done < "$input_file"