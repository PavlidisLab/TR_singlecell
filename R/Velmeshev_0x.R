library(WGCNA)
library(tidyverse)
library(Matrix)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")

sc_dir <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/jules_garreau_sc_datasets/Velmeshev/"
dat_path <- file.path(sc_dir, "matrix.mtx")
meta_path <- file.path(sc_dir, "meta.txt")
genes_path <- file.path(sc_dir, "genes.tsv")
out_path <- "/space/scratch/amorin/R_objects/Velmeshev_mat_and_meta.RDS"
pc <- read.delim("/home/amorin/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)


if (!file.exists(out_path)) {
  
  dat <- as.matrix(Matrix::readMM(dat_path))
  meta <- read.delim(meta_path)
  genes <- read.delim(genes_path, header = FALSE)
  
  meta <- meta %>% 
    dplyr::rename(ID = cell, Cell_type = cluster)
  
  mat <- dat
  colnames(mat) <- meta$ID
  rownames(mat) <- genes$V2
  
  stopifnot(identical(colnames(mat), meta$ID))
  
  common_genes <- intersect(pc$Symbol, genes$V2)
  
  mat <- mat[common_genes, meta$ID]

  zero_genes <- names(which(apply(mat, 1, function(x) sum(x != 0)) == 0))
  keep_genes <- unique(setdiff(rownames(mat), zero_genes))
  mat <- mat[keep_genes, ]
  
  saveRDS(list(mat, meta), file = out_path)
  
  rm(dat)
  gc()
  
} else {
  
  dat <- readRDS(out_path)
  meta <- dat[[2]]
  mat <- dat[[1]]
  
}


rsr1 <- all_RSR_aggregate1(mat, meta)
saveRDS(rsr1, file = "/space/scratch/amorin/R_objects/Velmeshev_RSR1.RDS")

z1 <- all_zscore_aggregate(mat, meta)
saveRDS(z1, file = "/space/scratch/amorin/R_objects/Velmeshev_Z1.RDS")

rsr2 <- all_RSR_aggregate2(mat, meta)
saveRDS(rsr2, file = "/space/scratch/amorin/R_objects/Velmeshev_RSR2.RDS")