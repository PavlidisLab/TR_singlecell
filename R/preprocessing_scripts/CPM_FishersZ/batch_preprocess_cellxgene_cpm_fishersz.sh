#!/bin/bash

# Preprocess all cell x gene datasets provided in the input file.
# Input file is assumed to be a .tsv of ID in first column and species in second.

set -euo pipefail

input_file="$1"
ncore=4  # Number of parallel xargs processes

# Preprocessing script
preprocess="/home/amorin/Projects/TR_singlecell/R/preprocessing_scripts/CPM_FishersZ/preprocess_cellxgene_cpm_fishersz.R"


if [ ! -f "$input_file" ]; then
    echo "Error: input file $input_file not found." >&2
    exit 1
fi


cut -f 1,2 "$input_file" | xargs -n 2 -P "$ncore" -- Rscript "$preprocess"
