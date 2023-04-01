## Implementing a check to see if both genes meet a minimum threshold of cells
## with non-zero expression in order to be considered for correlation
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(parallel)
library(microbench)
source("R/utils/functions.R")
source("R/00_config.R")

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Load correlation matrices generated per cell type and across all cells
cor_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_celltype.RDS")

# Genes and TFs of interest
genes <- rownames(sdat)
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

# How many cells must have non-zero expression for pairs of genes
min_count <- 20

# Demonstration on a single cell type
ct <- unique(sdat$Cell_type)[1]
expr_mat <- GetAssayData(object = sdat, slot = "counts")
expr_mat <- as.matrix(expr_mat[, sdat$Cell_type == ct, drop = FALSE])
dim(expr_mat)  #  21074 309


# Construct the binary/mask matrix: 1 if gene pair meets threshold, 0 otherwise.
# Assumes expr_mat rows are genes and columns are cells.

build_mask <- function(expr_mat, subset_genes) {
  
  stopifnot(is.matrix(expr_mat))
  mat <- matrix(0, nrow = nrow(expr_mat), ncol = length(subset_genes))
  rownames(mat) <- rownames(expr_mat)
  colnames(mat) <- subset_genes
  
  return(mat)
}


thresh_mat <- build_mask(expr_mat, tfs)
# thresh_mat <- matrix(0, nrow = length(genes), ncol = length(tfs))
# rownames(thresh_mat) <- genes
# colnames(thresh_mat) <- tfs


# Each function takes in a transcript count matrix (expr_mat) and returns XXX


# Original naive implementation


function_1 <- function(expr_mat, subset_genes, min_count) {
  
  thresh_mat <- build_mask(expr_mat, subset_genes)
  
  for (i in 1:nrow(expr_mat)) {
    for (j in 1:length(tfs)) {
      if (sum(expr_mat[i, ] > 0 & expr_mat[j, ] > 0) >= min_count) {
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
  
  for (i in 1:nrow(thresh_mat)) {
    for (j in 1:nrow(thresh_mat)) {
      if (sum((bin_mat[i, ] + bin_mat[j, ]) == 2) >= min_count) {
        thresh_mat[i, j] <- 1
      }
    }
  }
  
  
  return(thresh_mat)
}

# eric's suggestion for NA matching
expr_mat <- as.matrix(sdat@assays$RNA@counts[, sdat$Cell_type == sdat$Cell_type[1]])
min_count <- 2
bin_mat <- expr_mat > 0
thresh_mat <- matrix(0, nrow = nrow(expr_mat), ncol = nrow(expr_mat))
rownames(thresh_mat) <- colnames(thresh_mat) <- rownames(expr_mat)


# a <- rep(0, 10)
# a[c(1, 3)] <- 1
# b <- rep(0, 10)
# b[c(1, 4)] <- 1
# a + b
# sum((a + b) == 2) >= min_count

tt <- rowSums(bin_mat)

i <- 6
j <- 9

# sum(bin_mat[i, ])
# sum(bin_mat[j, ])
# sum(bin_mat[i, ] + bin_mat[j, ] == 2)

for (i in 1:nrow(thresh_mat)) {
  for (j in 1:nrow(thresh_mat)) {
    if (sum((bin_mat[i, ] + bin_mat[j, ]) == 2) >= min_count) {
      thresh_mat[i, j] <- 1
    }
  }
}


sum(thresh_mat)







microbenchmark::microbenchmark(
  function_1(expr_mat, tfs, min_count),
  times = 2
)





# Angus's suggestion


#for each gene get number of cells expressing the gene 
#cells_expressing <- rowSums(expresission_matrix>0)
# for testing I will just random gen the vector this would produce

# gene length vector of cell counts with reads > 0 
cells_expressing <- floor(runif(n = 10, min=0, max = 20))

ngenes <- length(cells_expressing)

# make df to use a mask, first replicate each gene's cell count
# to make a ngene x ngene matrix
both_genes_pass_thresh_mask <- matrix(replicate(ngenes,cells_expressing),nrow= ngenes)

# function to use in Map
and_pass <- function(x,y, threshold){
  out <- (y>threshold)&(x>threshold)
  return(out)
}

# set threshold both genes must pass
thresh <- 5
both_genes_pass_thresh_mask[] <- Map(and_pass, both_genes_pass_thresh_mask,
                                     y = cells_expressing, threshold = thresh)

# the above produces a vector we turn it back into a matrix
both_genes_pass_thresh_mask <- matrix(both_genes_pass_thresh_mask,nrow= ngenes)



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


function1 <- function() {
  print("hello")
}

function2 <- function() {
  print("hello")
}

res <- microbenchmark::microbenchmark(function1, function2, times = 1e3)
boxplot(res)
