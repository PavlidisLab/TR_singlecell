## Generate various metrics of similarity for TF vectors across aggregate
## coexpression networks
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
force_resave <- TRUE

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



# For a given set of input genes and list of aggregate coexpression matrices,
# calculate multiple metrics of similarity within the provided genes across networks. 
# agg_l: A list of aggregate coexpression matrices (gene x gene numeric matrix)
# msr_mat: binary gene x experiment matrix tracking if a gene was measured
# k: an integer of the top/bottom elements to select
# check_k_arg: logical controls if k should be reduced to the largest non-tied element
# returns: a list for each gene containing a dataframe of all unique experiment
# pairs in which the given gene was measured, and three measures of similarity:
# 1) The spearman correlation between the gene vectors of the experiment pair;
# 2) The size of the top k intersect between the experiments;
# 3) The Jaccard of the top and bottom k elements between the experiments
# ------------------------------------------------------------------------------



get_all_similarity <- function(agg_l, 
                               msr_mat, 
                               genes, 
                               k = 1000, 
                               check_k_arg = TRUE) {
  
  stopifnot(genes %in% colnames(agg_l[[1]]))

  sim_l <- lapply(genes, function(x) {

    message(paste(x, Sys.time()))

    # For the given gene, create a matrix from the experiments that it is measured
    # in. Remove the gene itself to prevent inflated overlap
    gene_mat <- gene_vec_to_mat(agg_l, x)
    gene_mat <- gene_mat[setdiff(rownames(gene_mat), x), ]
    gene_mat <- gene_mat[, which(msr_mat[x, ] == 1), drop = FALSE]

    if (ncol(gene_mat) < 2) {  # need more than one experiment for pairing
      return(NA)
    }

    # Experiment x experiment similarity matrices
    gene_cor <- colwise_cor(gene_mat, cor_method = "spearman")
    gene_topk <- colwise_topk_intersect(gene_mat, k = k, check_k_arg = check_k_arg)
    gene_jacc <- colwise_jaccard(binarize_topk_btmk(gene_mat, k = k))

    # Data frame of unique dataset pairs and their similarities
    get_similarity_pair_df(gene_cor, gene_topk, gene_jacc)
    
  })

  names(sim_l) <- genes
  return(sim_l)
}



# Run and save
# ------------------------------------------------------------------------------


if (!file.exists(tf_sim_hg_path) || force_resave) {
  
  sim_hg <- get_all_similarity(agg_l = agg_tf_hg, 
                               msr_mat = msr_hg, 
                               genes = tfs_hg$Symbol, 
                               k = k)
  
  sim_hg <- sim_hg[!is.na(sim_hg)]
  
  saveRDS(sim_hg, tf_sim_hg_path)
  
}



if (!file.exists(tf_sim_mm_path) || force_resave) {
  
  sim_mm <- get_all_similarity(agg_l = agg_tf_mm, 
                               msr_mat = msr_mm, 
                               genes = tfs_mm$Symbol, 
                               k = k)
  
  sim_mm <- sim_mm[!is.na(sim_mm)]
  
  saveRDS(sim_mm, tf_sim_mm_path)
}
