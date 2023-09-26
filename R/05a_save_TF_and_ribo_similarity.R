## Generate various metrics of similarity for TF and L/S ribosomal vectors 
## across aggregate coexpression networks
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

k <- 200

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

# Ribosomal genes
sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Loading the subset TF and ribosomal aggregate matrices

agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

tfs_hg <- filter(tfs_hg, Symbol %in% colnames(agg_tf_hg[[1]]))
tfs_mm <- filter(tfs_mm, Symbol %in% colnames(agg_tf_mm[[1]]))

agg_ribo_hg <- load_or_generate_agg(path = agg_ribo_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = ribo_genes$Symbol_hg)
agg_ribo_mm <- load_or_generate_agg(path = agg_ribo_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = ribo_genes$Symbol_mm)



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
                               k, 
                               check_k_arg = TRUE) {
  
  stopifnot(genes %in% colnames(agg_l[[1]]))

  sim_l <- lapply(genes, function(x) {

    message(paste(x, Sys.time()))

    # For the given gene, create a matrix from the experiments that it is measured
    # in. Remove the gene itself to prevent inflated overlap
    gene_mat <- gene_vec_to_mat(agg_l, x)
    gene_mat <- gene_mat[setdiff(rownames(gene_mat), x), ]
    gene_mat <- subset_to_measured(gene_mat, msr_mat = msr_mat, gene = x)

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



save_all_similarity <- function(path, 
                                agg_l, 
                                msr_mat, 
                                genes, 
                                k, 
                                check_k_arg = TRUE,
                                force_resave = FALSE) {
  
  if (!file.exists(path) || force_resave) {
    
    sim <- get_all_similarity(agg_l = agg_l, 
                              msr_mat = msr_mat, 
                              genes = genes, 
                              k = k,
                              check_k_arg = check_k_arg)
    
    sim <- sim[!is.na(sim)]
    saveRDS(sim, path)
    return(invisible(NULL))
    
  }
}



# Run and save
# ------------------------------------------------------------------------------


# Human TFs
save_all_similarity(path = sim_tf_hg_path, 
                    agg_l = agg_tf_hg, 
                    msr_mat = msr_hg, 
                    genes = tfs_hg$Symbol, 
                    k = k)


# Mouse TFs
save_all_similarity(path = sim_tf_mm_path, 
                    agg_l = agg_tf_mm, 
                    msr_mat = msr_mm, 
                    genes = tfs_mm$Symbol, 
                    k = k)


# Human ribo
save_all_similarity(path = sim_ribo_hg_path, 
                    agg_l = agg_ribo_hg, 
                    msr_mat = msr_hg, 
                    genes = ribo_genes$Symbol_hg, 
                    k = k)


# Mouse ribo
save_all_similarity(path = sim_ribo_mm_path, 
                    agg_l = agg_ribo_mm, 
                    msr_mat = msr_mm, 
                    genes = ribo_genes$Symbol_mm, 
                    k = k)
