library(WGCNA)
library(tidyverse)
library(Matrix)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")

id <- "GSE216019"
sc_dir <- paste0("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datases/Human/", id)
dat_path <- file.path(sc_dir, paste0(id, ".RDS"))
out_path <- paste0("/space/scratch/amorin/R_objects/", id, "_mat_and_meta.RDS")
pc <- read.delim("/home/amorin/Data/Metadata/ensembl_human_protein_coding_105.tsv", stringsAsFactors = FALSE)


dat <- readRDS(dat_path)

meta <- dat@meta.data %>% 
  dplyr::rename(Cell_type = cell_type) %>% 
  rownames_to_column(var = "ID")

mat_norm <- as.matrix(GetAssayData(dat, slot = "data"))
mat_counts <- as.matrix(GetAssayData(dat, slot = "counts"))

avg_norm <- rowMeans(mat_norm)
avg_counts <- rowMeans(mat_counts)

avg_df <- data.frame(ID = rownames(mat_norm),
                     Count = avg_counts, 
                     Norm = avg_norm)
