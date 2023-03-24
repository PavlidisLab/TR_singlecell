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

# Prepare count mat
mat <- t(as.matrix(sdat@assays$RNA@data))


# Loop over TFs and cor with every other gene

cor_ct <- lapply(cts, function(x) {
  
  message(paste(x, Sys.time()))
  
  cells <- sdat$Cell_type == x
  
  cmat <- WGCNA::cor(x = mat[cells, ],
                     nThreads = ncore,
                     use = "pairwise.complete.obs")
})

names(cor_ct) <- cts


saveRDS(cor_ct, cor_ct_path)
