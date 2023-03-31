## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(parallel)
source("R/utils/functions.R")
source("R/00_config.R")

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Load correlation matrices generated per cell type and across all cells
cor_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_celltype.RDS")

tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

min_count <- 20


## Trying a post-hoc sweep of min count filter


# ct <- unique(sdat$Cell_type)[1]
# cor_mat <- cor_ct[[ct]]
# expr_mat <- sdat@assays$RNA@counts[, sdat$Cell_type == ct]
# cor_mat[1:5, 1:5]
# expr_mat[1:5, 1:5]
# i <- 1
# j <- 500
# sum(expr_mat[i, ] > 0)
# sum(expr_mat[j, ] > 0)
# sum(expr_mat[i, ] > 0 & expr_mat[j, ] > 0) >= min_count



genes <- rownames(sdat)
thresh_mat <- matrix(0, nrow = length(genes), ncol = length(tfs))
rownames(thresh_mat) <- genes
colnames(thresh_mat) <- tfs


thresh_l <- mclapply(unique(sdat$Cell_type), function(ct) {
  
  expr_mat <- sdat@assays$RNA@counts[, sdat$Cell_type == ct]
  
  message(ct, Sys.time())
  
  for (i in genes) {
    for (j in tfs) {
      if (sum(expr_mat[i, ] > 0 & expr_mat[j, ] > 0) >= min_count) {
        thresh_mat[i, j] <- 1
      }
    }
  }
  
  return(thresh_mat)
  
}, mc.cores = 8)

names(thresh_l) <- unique(sdat$Cell_type)

saveRDS(thresh_l, "/space/scratch/amorin/R_objects/Hochgerner2022_cormat_celltype_mincounts.RDS")
