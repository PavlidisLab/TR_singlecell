## TODO
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



# TODO: order of sort + check


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



# Binarize matrix such that top k and bottom k is 1, everything else 0

binarize_topk_btmk <- function(mat, k = 1000) {
  
  bin_mat <- apply(mat, 2, function(x) {
    sort_values <- sort(x, decreasing = TRUE)
    topk <- sort_values[k]
    btmk <- sort_values[length(x) - k + 1]
    ifelse(x >= topk | x <= btmk, 1, 0)
  })
  
  return(bin_mat)
}



colwise_jaccard <- function(mat) {
  
  jaccard <- function(vec1, vec2) {
    sum(vec1 & vec2, na.rm = TRUE) / sum(vec1 | vec2, na.rm = TRUE)
  }
  
  col_list <- asplit(mat, 2)
  jacc_mat <- outer(col_list, col_list, Vectorize(jaccard))
  return(jacc_mat)
}


# Convert similarity matrices into a dataframe of unique pairs and values

get_similarity_pair_df <- function(cor_mat, topk_mat, jacc_mat) {
  
  df <-
    mat_to_df(cor_mat, symmetric = TRUE, value_name = "Scor") %>%
    cbind(
      Topk = mat_to_df(topk_mat, symmetric = TRUE)$Value,
      Jaccard = mat_to_df(jacc_mat, symmetric = TRUE)$Value
    )
  
  return(df)
}



get_all_similarity <- function(agg_l, 
                               msr_mat, 
                               genes, 
                               k = 1000, 
                               check_k_arg = TRUE) {

  sim_l <- lapply(genes, function(x) {

    message(paste(x, Sys.time()))

    # Gene matrix with only measured experiments and prevent self gene inflating overlap
    gene_mat <- gene_vec_to_mat(agg_l, x)
    gene_mat <- gene_mat[setdiff(rownames(gene_mat), x), ]
    gene_mat <- gene_mat[, which(msr_mat[x, ] == 1), drop = FALSE]

    if (ncol(gene_mat) < 2) {  # need more than one experiment for pairing
      return(NA)
    }

    # Cor, Topk, and Jaccard of top and bottom k
    gene_cor <- colwise_cor(gene_mat, cor_method = "spearman")
    gene_topk <- colwise_topk_intersect(gene_mat, k = k, check_k_arg = check_k_arg)
    gene_jacc <- colwise_jaccard(binarize_topk_btmk(gene_mat, k = k))

    get_similarity_pair_df(gene_cor, gene_topk, gene_jacc)
  
  })

  names(sim_l) <- genes
  return(sim_l)
}



#
# ------------------------------------------------------------------------------


if (!file.exists(tf_sim_hg_path)) {
  
  sim_hg <- get_all_similarity(agg_l = agg_tf_hg, 
                               msr_mat = msr_hg, 
                               genes = tfs_hg$Symbol, 
                               k = k)
  
  sim_hg <- sim_hg[!is.na(sim_hg)]
  
  saveRDS(sim_hg, tf_sim_hg_path)
  
}



if (!file.exists(tf_sim_mm_path)) {
  
  sim_mm <- get_all_similarity(agg_l = agg_tf_mm, 
                               msr_mat = msr_mm, 
                               genes = tfs_mm$Symbol, 
                               k = k)
  
  sim_mm <- sim_mm[!is.na(sim_mm)]
  
  saveRDS(sim_mm, tf_sim_mm_path)
}
