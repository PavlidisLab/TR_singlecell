library(WGCNA)
library(tidyverse)
library(Matrix)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")

sc_dir <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/jules_garreau_sc_datasets/ROSMAP"
dat_path <- file.path(sc_dir, "sc_exprmats_rosmap.rds")
out_path <- "/space/scratch/amorin/R_objects/ROSMAP_mat_and_meta.RDS"

pc <- read.delim("/home/amorin/Data/Metadata/ensembl_human_protein_coding_105.tsv", stringsAsFactors = FALSE)


if (!file.exists(out_path)) {
  
  dat <- readRDS(dat_path)
  
  # Filling out mat with separated cell types
  mat <- as.matrix(dat$raw$ctmat)
  
  for (i in names(dat$exprmats)) {
    ct_mat <- as.matrix(dat$exprmats[[i]])
    mat[rownames(ct_mat), colnames(ct_mat)] <- ct_mat
  }

  stopifnot(identical(
    as.matrix(dat$exprmats$astrocyte), 
    mat[rownames(dat$exprmats$astrocyte), colnames(dat$exprmats$astrocyte)]))
 
  # Extract metadata
  meta <- dat$raw$samples %>% 
    dplyr::rename(Cell_type = cell_type, ID = sample)

  stopifnot(all(colnames(mat) %in% meta$ID))

  # Ensembl IDs to gene symbols
  common_ens_ids <- intersect(pc$Gene_ID, rownames(mat))
  
  common_genes <- data.frame(
    Ensembl_ID = common_ens_ids,
    Symbol = pc$Symbol[match(common_ens_ids, pc$Gene_ID)])

  mat <- mat[common_genes$Ensembl_ID, meta$ID]
  rownames(mat) <- common_genes$Symbol
  
  zero_genes <- names(which(apply(mat, 1, function(x) sum(x != 0)) == 0))
  keep_genes <- unique(setdiff(rownames(mat), zero_genes))
  keep_genes <- keep_genes[keep_genes != ""]
  
  mat <- mat[keep_genes, ]

  saveRDS(list(mat, meta), file = out_path)
  
  rm(dat)
  gc()
  
} else {
  
  dat <- readRDS(out_path)
  meta <- dat[[2]]
  mat <- dat[[1]]
  
}


stopifnot(identical(colnames(mat), meta$ID))

rsr1 <- all_RSR_aggregate1(mat, meta)
saveRDS(rsr1, file = "/space/scratch/amorin/R_objects/ROSMAP_RSR1.RDS")

z1 <- all_zscore_aggregate(mat, meta)
saveRDS(z1, file = "/space/scratch/amorin/R_objects/ROSMAP_Z1.RDS")

rsr2 <- all_RSR_aggregate2(mat, meta)
saveRDS(rsr2, file = "/space/scratch/amorin/R_objects/ROSMAP_RSR2.RDS")
