## Sample genes from experiments and calculate the topk intersect across 
## experiments to generate a null. Save out the result as a list.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 1000 
n_samps <- 1000

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)



check_k <- function(vec_sort, k) {
  
  if (vec_sort[k] == vec_sort[k - 1]) {
    tie_start <- head(sort(table(vec_sort), decreasing = TRUE))[1]
    k <- sum(vec_sort > as.numeric(names(tie_start)), na.rm = TRUE)
  }
  
  return(k)
}



topk_sort <- function(vec, k, check_k_arg = TRUE) {
  vec_sort <- sort(vec, decreasing = TRUE)
  if (check_k_arg) k <- check_k(vec_sort, k)
  names(vec_sort[1:k])
}



topk_intersect <- function(vec1, vec2) length(intersect(vec1, vec2))



colwise_topk_intersect <- function(mat, k = 1000, check_k_arg = TRUE) {
  
  col_list <- asplit(mat, 2)
  topk_list <- lapply(col_list, topk_sort, k = k, check_k_arg = check_k_arg)
  topk_mat <- outer(topk_list, topk_list, Vectorize(topk_intersect))
  
  return(topk_mat)
}




# NOTE: The following implementation is very slow when calling multiple times,
# as each dataset is loaded each call. This could be made faster by generating
# the sampled genes first and loading the sampled genes for each dataset once. 
# This can quickly into memory issues however as each sampled matrix needs to be
# held in memory.


# For the given ids, sample a single gene that is measured from each experiment.
# Return a dataframe of the unique pairs and their topk intersect.

sample_topk_intersect <- function(ids, genes, msr_mat, k = 1000) {
  
  # Sample one gene that is measured in each data set
  sample_genes <- lapply(ids, function(x) {
    msr_genes <- msr_mat[msr_mat[, x] == 1, x]
    names(sample(msr_genes, 1))
  })
  
  # Load the sampled gene for each experiment and bind into a matrix
  sample_mat <- lapply(1:length(ids), function(x) {
    load_agg_mat_list(ids = ids[x], genes = genes, sub_genes = sample_genes[[x]])[[1]]
  })
  sample_mat <- do.call(cbind, sample_mat)
  
  # Get topk overlap between sampled genes
  sample_topk <- colwise_topk_intersect(sample_mat)
  sample_df <- mat_to_df(sample_topk, symmetric = TRUE, value_name = "Topk")
  
  return(sample_df)
}



set.seed(5)


if (!file.exists(null_topk_hg_path)) {
  
  topk_hg <- lapply(1:n_samps, function(x) {
    
    message(paste("Human sample #", x, Sys.time()))
    
    sample_topk_intersect(ids = ids_hg,
                          genes = pc_hg$Symbol,
                          msr_mat = msr_hg,
                          k = k)
  })
  
  saveRDS(topk_hg, null_topk_hg_path)
  
} 



if (!file.exists(null_topk_mm_path)) {
  
  topk_mm <- lapply(1:n_samps, function(x) {
    
    message(paste("Mouse sample #", x, Sys.time()))
    
    sample_topk_intersect(ids = ids_mm,
                          genes = pc_mm$Symbol,
                          msr_mat = msr_mm,
                          k = k)
  })
  
  saveRDS(topk_mm, null_topk_mm_path)
  
} 
