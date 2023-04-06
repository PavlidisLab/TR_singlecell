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
out_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_aggregate_celltype_cor.RDS"

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

for (i in cts[1]) {
  
  message(paste(i, Sys.time()))
  
  # Isolate cells from cell type. If not enough cells meet the minimum non-zero
  # threshold, set all counts to NA (producing an NA cor)
  
  cells <- sdat$Cell_type == i
  mat <- t(as.matrix(sdat@assays$RNA@data[, cells]))
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) <= min_count
  mat[, na_genes] <- NA
  
  # Pearson cor mat
  cmat <- WGCNA::cor(x = mat, nThreads = ncore, use = "pairwise.complete.obs")
  
  # Convert to ranks and set NA ranks to mean of rank matrix
  rmat <- colrank_mat(cmat)
  rmat <- na_to_mean(rmat)
  
  # Running addition to aggregate matrix
  amat <- amat + rmat
  
  # Clean up
  gc()
  
}


# Rank of sum of ranks
amat <- colrank_mat(-amat)

# Save out
saveRDS(amat, out_path)
