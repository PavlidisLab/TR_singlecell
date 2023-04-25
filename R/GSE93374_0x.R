library(WGCNA)
library(tidyverse)
library(Matrix)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")

sc_dir <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/mice/has_celltype_metadata/GSE93374"
dat_path <- file.path(sc_dir, "GSE93374_Merged_all_020816_BatchCorrected_LNtransformed_doubletsremoved_Data.txt")
meta_path <- file.path(sc_dir, "GSE93374_cell_metadata.txt")
out_path <- "/space/scratch/amorin/R_objects/GSE93374_mat_and_meta.RDS"

pc <- read.delim("/home/amorin/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)


if (!file.exists(out_path)) {
  
  dat <- read.delim(dat_path)
  meta <- read.delim(meta_path)
  
  dat <- as.matrix(dat)
  
  # Get common IDs (different between data and metadata...)
  meta_ids <- str_replace(meta$X1.ID, ".*_", "")
  dat_ids <- str_replace(colnames(dat), ".*_", "")
  
  # Remove handful of duplicate IDs
  meta_ids_dup <- duplicated(meta_ids)
  dat_ids_dup <- duplicated(dat_ids)
  
  dat <- dat[, !dat_ids_dup]
  meta <- meta[!meta_ids_dup, ]
  
  colnames(dat) <- str_replace(colnames(dat), ".*_", "")
  
  meta <- meta %>% 
    mutate(ID = str_replace(X1.ID, ".*_", "")) %>% 
    filter(ID %in% intersect(colnames(dat), ID)) %>% 
    dplyr::rename(Cell_type = X7.clust_all)
  
  
  common_genes <- intersect(pc$Symbol, rownames(dat))

  mat <- dat[common_genes, meta$ID]
  
  saveRDS(list(mat, meta), file = out_path)
  
} else {
  
  dat <- readRDS(out_path)
  meta <- dat[[2]]
  mat <- dat[[1]]
  
}


stopifnot(identical(colnames(mat), meta$ID))

rsr1 <- all_RSR_aggregate1(mat, meta)
saveRDS(rsr1, file = "/space/scratch/amorin/R_objects/GSE93374_RSR1.RDS")

z1 <- all_zscore_aggregate(mat, meta)
saveRDS(z1, file = "/space/scratch/amorin/R_objects/GSE93374_Z1.RDS")

rsr2 <- all_RSR_aggregate2(mat, meta)
saveRDS(rsr2, file = "/space/scratch/amorin/R_objects/GSE93374_RSR2.RDS")
