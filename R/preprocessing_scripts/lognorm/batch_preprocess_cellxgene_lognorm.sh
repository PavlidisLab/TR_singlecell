#!/bin/bash

# Input file assumes that the first column is the ID and the second column is
# the species. Iterates over the fields of the input file, pre-processing the
# associated data and generating an aggregate coexpression network.


sc_dir="/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets/"
input_file=$1
preprocess="/home/amorin/Projects/TR_singlecell/R/preprocessing_scripts/lognorm/preprocess_cellxgene_lognorm.R"


while IFS=$'\t' read -r id species url; do

  echo "Beginning $id"

  Rscript --vanilla "$preprocess" "$id" "$species"

done < "$input_file"
