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
source("R/utils/functions.R")
source("R/00_config.R")

k <- 1000 
n_samps <- 1000
null_topk_hg_path <- "/space/scratch/amorin/R_objects/sampled_TF_null_topk_intersect_human.RDS"
null_topk_hg_path <- "/space/scratch/amorin/R_objects/sampled_TF_null_topk_intersect_mouse.RDS"


# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# For each dataset, load the subset gene x TF aggregation matrix 

if (!file.exists(agg_tf_hg_path)) {
  agg_tf_hg <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
  saveRDS(agg_tf_hg, agg_tf_hg_path)
} else {
  agg_tf_hg <- readRDS(agg_tf_hg_path)
}


if (!file.exists(agg_tf_mm_path)) {
  agg_tf_mm <- load_agg_mat_list(ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)
  saveRDS(agg_tf_mm, agg_tf_mm_path)
} else {
  agg_tf_mm <- readRDS(agg_tf_mm_path)
}



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



# 

sample_topk_intersect <- function(agg_l, genes, msr_mat, k = 1000) {
  
  ids <- names(agg_l)
  
  # Sample one gene that is measured in each data set
  sample_genes <- unlist(lapply(ids, function(x) {
    msr_gene <- names(which(msr_mat[genes, x] == 1))
    sample(msr_gene, 1)
  }))
  
  # Bind sampled genes into a matrix
  sample_mat <- lapply(1:length(sample_genes), function(x) {
    agg_l[[x]][, sample_genes[x]]
  })
  sample_mat <- do.call(cbind, sample_mat)
  colnames(sample_mat) <- paste0(ids, "_", sample_genes)
  
  # Get topk overlap between sampled genes
  sample_topk <- colwise_topk_intersect(sample_mat)
  sample_df <- mat_to_df(sample_topk, symmetric = TRUE, value_name = "Topk")
  
  return(sample_df)
}



set.seed(5)



if (!file.exists(null_topk_hg_path)) {
  
  topk_hg <- lapply(1:n_samps, function(x) {
    
    message(paste("Human sample #", x, Sys.time()))
    
    sample_topk_intersect(agg_l = agg_tf_hg,
                          genes = tfs_hg$Symbol,
                          msr_mat = msr_hg,
                          k = k)
  })
  
  saveRDS(topk_hg, null_topk_hg_path)
  
} 



if (!file.exists(null_topk_mm_path)) {
  
  topk_mm <- lapply(1:n_samps, function(x) {
    
    message(paste("Mouse sample #", x, Sys.time()))
    
    sample_topk_intersect(agg_l = agg_tf_mm,
                          genes = tfs_mm$Symbol,
                          msr_mat = msr_mm,
                          k = k)
  })
  
  saveRDS(topk_mm, null_topk_mm_path)
  
} 
