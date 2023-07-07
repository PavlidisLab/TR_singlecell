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

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Load aggregate matrix into list
agg_hg <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_mm <- load_agg_mat_list(ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

stopifnot(all(unlist(lapply(agg_hg, function(x) identical(rownames(x), pc_hg$Symbol)))))
stopifnot(all(unlist(lapply(agg_mm, function(x) identical(rownames(x), pc_mm$Symbol)))))

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)



# Functions
# ------------------------------------------------------------------------------



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




get_all_similarity <- function(agg_l, msr_mat, genes, k = 1000) {
  
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
    gene_topk <- colwise_topk_intersect(gene_mat, k = k)
    gene_jacc <- colwise_jaccard(binarize_topk_btmk(gene_mat, k = k))
    
    list(
      Sim_df = get_similarity_pair_df(gene_cor, gene_topk, gene_jacc),
      N_exp = ncol(gene_mat)
    )

  })
  
  names(sim_l) <- genes
  
  return(sim_l)
}



#
# ------------------------------------------------------------------------------


if (!file.exists(tf_sim_hg_path)) {
  sim_hg <- get_all_similarity(agg_l = agg_hg, msr_mat = msr_hg, genes = tfs_hg$Symbol, k = 1000)
  sim_hg <- sim_hg[!is.na(sim_hg)]
  saveRDS(sim_hg, tf_sim_hg_path)
} else {
  sim_hg <- readRDS(tf_sim_hg_path)
}



if (!file.exists(tf_sim_mm_path)) {
  sim_mm <- get_all_similarity(agg_l = agg_mm, msr_mat = msr_mm, genes = tfs_mm$Symbol, k = 1000)
  sim_mm <- sim_mm[!is.na(sim_mm)]
  saveRDS(sim_mm, tf_sim_mm_path)
} else {
  sim_mm <- readRDS(tf_sim_mm_path)
}
