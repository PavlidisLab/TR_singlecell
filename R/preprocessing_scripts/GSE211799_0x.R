library(WGCNA)
library(tidyverse)
library(Matrix)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")

id <- "GSE211799"
sc_dir <- paste0("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datases/Mouse/", id)
dat_path <- file.path(sc_dir, paste0(id, ".RDS"))
out_path <- paste0("/space/scratch/amorin/R_objects/", id, "_mat_and_meta.RDS")
pc <- read.delim("/home/amorin/Data/Metadata/ensembl_mouse_protein_coding_105.tsv", stringsAsFactors = FALSE)


if (!file.exists(out_path)) {
  
  dat <- readRDS(dat_path)
  
  meta <- dat@meta.data %>% 
    dplyr::rename(Cell_type = cell_type) %>% 
    rownames_to_column(var = "ID")
  
  mat <- as.matrix(GetAssayData(dat, slot = "data"))
  
  stopifnot(identical(colnames(mat), meta$ID))
  
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


rsr1 <- all_RSR_aggregate1(mat, meta)
saveRDS(rsr1, file = paste0("/space/scratch/amorin/R_objects/", id, "_RSR1.RDS"))

z1 <- all_zscore_aggregate(mat, meta)
saveRDS(z1, file = paste0("/space/scratch/amorin/R_objects/", id, "_Z1.RDS"))

rsr2 <- all_RSR_aggregate2(mat, meta)
saveRDS(rsr2, file = paste0("/space/scratch/amorin/R_objects/", id, "_RSR2.RDS"))