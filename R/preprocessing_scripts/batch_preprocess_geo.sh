#!/bin/bash

input_file=$1
preprocess_dir="/home/amorin/Projects/TR_singlecell/R/preprocessing_scripts/"

while IFS=$'\t' read -r id species; do

    r_script="${preprocess_dir}/${id}.R"
    
    echo "Beginning $id"
    
    Rscript --vanilla "$r_script"
    
done < "$input_file"