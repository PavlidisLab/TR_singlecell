## Generate the aggregated cell type correlation as in Harris et al., 2021.
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(future)
library(parallel)
library(WGCNA)
source("R/utils/functions.R")
source("R/00_config.R")

# Output RDS
out_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_binthresh_celltype_cor.RDS"

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Cell types to loop over
cts <- unique(sdat$Cell_type)

# Min count of non-zero expressing cells to keep gene for correlating
min_count <- 20

# Init aggregated matrix
amat <- matrix(0, nrow = nrow(sdat), ncol = nrow(sdat))
rownames(amat) <- colnames(amat) <- rownames(sdat)


# Loop over TFs and cor with every other gene

thresh <- floor(nrow(sdat) * 0.005)
thresh <- thresh + 1  # +1 for the moment to account for keeping self cor


for (i in cts) {
  
  message(paste(i, Sys.time()))
  
  # Isolate cells from cell type. If not enough cells meet the minimum non-zero
  # threshold, set all counts to NA (producing an NA cor)
  
  cells <- sdat$Cell_type == i
  mat <- t(as.matrix(sdat@assays$RNA@data[, cells]))
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) <= min_count
  mat[, na_genes] <- NA
  
  # Pearson cor mat
  cmat <- WGCNA::corAndPvalue(x = mat, use = "pairwise.complete.obs", alternative = "greater")
  
  # Threshold by pval
  
  binmat <- mclapply(1:ncol(cmat$p), function(x) {
    
    vec <- cmat$p[, x]
    passes <- names(sort(vec)[1:thresh])
    vec <- ifelse(names(vec) %in% passes, 1, 0)
    return(vec)
  
  }, mc.cores = ncore)
  
  binmat <- do.call(cbind, binmat)

  # Running addition to aggregate matrix
  amat <- amat + binmat
  
  # Clean up
  gc()
  
}


# Rank of sum of ranks
# amat <- colrank_mat(-amat)

# Save out
saveRDS(amat, out_path)
