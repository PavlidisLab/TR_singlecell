## Generate gene-gene correlation matrix by cell type
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(future)
library(parallel)
library(WGCNA)
source("R/utils/functions.R")
source("R/00_config.R")

# Output RDS
cor_ct_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_cormat_celltype.RDS"

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

cts <- unique(sdat$Cell_type)

# Subset of genes to look at
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

# Min count of non-zero expressing cells to keep gene for correlating
min_count <- 20


# Loop over TFs and cor with every other gene

cor_ct <- lapply(cts, function(x) {
  
  message(paste(x, Sys.time()))
  
  # Isolate cells from cell type. If not enough cells meet the minimum non-zero
  # threshold, set all counts to NA (producing an NA cor)
  
  cells <- sdat$Cell_type == x
  mat <- t(as.matrix(sdat@assays$RNA@data[, cells]))
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) <= min_count
  mat[, na_genes] <- NA
  
  cmat <- WGCNA::cor(x = mat,
                     y = mat[, tfs],
                     nThreads = ncore,
                     use = "pairwise.complete.obs")
})

names(cor_ct) <- cts


saveRDS(cor_ct, cor_ct_path)
