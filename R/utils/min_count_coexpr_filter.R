## Implementing a check to see if both genes meet a minimum threshold of cells
## with non-zero expression in order to be considered for correlation
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(parallel)
library(microbenchmark)
source("R/utils/functions.R")
source("R/00_config.R")

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# TFs of interest
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

# How many cells must have non-zero expression for pairs of genes
min_count <- 20

# Demonstration on a single cell type
ct <- unique(sdat$Cell_type)[1]
expr_mat <- GetAssayData(object = sdat, slot = "counts")
expr_mat <- as.matrix(expr_mat[, sdat$Cell_type == ct, drop = FALSE])
dim(expr_mat)  #  21074 309


# Init with 0s the binary/mask matrix: 1 if gene pair meets threshold, 0 otherwise.
# Assumes expr_mat rows are genes and columns are cells.

build_mask <- function(expr_mat, subset_genes) {
  
  stopifnot(is.matrix(expr_mat))
  mat <- matrix(0, nrow = nrow(expr_mat), ncol = length(subset_genes))
  rownames(mat) <- rownames(expr_mat)
  colnames(mat) <- subset_genes
  
  return(mat)
}



# Original naive implementation

function_1 <- function(expr_mat, subset_genes, min_count) {
  
  thresh_mat <- build_mask(expr_mat, subset_genes)
  
  # Get indices to prevent costly name look-up
  subset_ix <- which(rownames(expr_mat) %in% subset_genes)
  
  for (i in 1:nrow(thresh_mat)) {
    for (j in 1:length(subset_ix)) {
      if (sum(expr_mat[i, ] > 0 & expr_mat[subset_ix[j], ] > 0) >= min_count) {
        thresh_mat[i, j] <- 1
      }
    }
  }
  
  return(thresh_mat)
}



# Eric's suggestion

function_2 <- function(expr_mat, subset_genes, min_count) {
  
  thresh_mat <- build_mask(expr_mat, subset_genes)
  bin_mat <- expr_mat > 0
  
  # Get indices to prevent costly name look-up
  subset_ix <- which(rownames(expr_mat) %in% subset_genes)
  
  for (i in 1:nrow(thresh_mat)) {
    for (j in 1:length(subset_ix)) {
      if (sum(bin_mat[i, ] & bin_mat[subset_ix[j], ]) >= min_count) {
        thresh_mat[i, j] <- 1
      }
    }
  }
  
  return(thresh_mat)
}


# Matrix multiplication: 0/1 binarization means that product will yield sum of
# non-zero expressing cells between genes, which can then threshold.

function_3 <- function(expr_mat, subset_genes, min_count) {
  
  # Binary matrix as 1/0 rather than T/F 
  # https://stackoverflow.com/questions/30943167/replace-logical-values-true-false-with-numeric-1-0
  bin_mat <- (expr_mat > 0) * 1

  # [Genes, Cells] %*% [Cells, Sub_genes] = [Genes, Sub_genes]
  thresh_mat <- crossprod(t(bin_mat), t(bin_mat[subset_genes, ]))
  thresh_mat <- (thresh_mat >= min_count) * 1
  
  return(thresh_mat)
}



tt1 <- function_1(expr_mat, tfs, min_count)
tt2 <- function_2(expr_mat, tfs, min_count)
tt3 <- function_3(expr_mat, tfs, min_count)


identical(tt1, tt2, tt3)


res <- microbenchmark::microbenchmark(
  F1 = function_1(expr_mat, tfs, min_count),
  F2 = function_2(expr_mat, tfs, min_count),
  F3 = function_3(expr_mat, tfs, min_count),
  times = 10
)



###
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
