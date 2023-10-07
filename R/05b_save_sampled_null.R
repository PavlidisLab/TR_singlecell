## Sample TFs across experiments and calculate their similarity to generate
## a null. Save out the results as a list.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 1000
n_samps <- 1000
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
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)



# For a given set of input genes and list of aggregate coexpression matrices, 
# sample one gene that is measured in each experiment from the list of input
# genes, and calculate the similarity between these sampled genes
# agg_l: A list of aggregate coexpression matrices (gene x gene numeric matrix)
# msr_mat: binary gene x experiment matrix tracking if a gene was measured
# k: an integer of the top elements to select
# check_k_arg: logical controls if k should be reduced to the largest non-tied element
# returns: a list for each gene containing a dataframe of all unique experiment
# pairs in which the given gene was measured, and their top k intersect
# ------------------------------------------------------------------------------


calc_sample_similarity <- function(agg_l,
                                   genes,
                                   msr_mat,
                                   k,
                                   check_k_arg = TRUE) {
  
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
  
  # Experiment x experiment similarity matrices
  cor_mat <- colwise_cor(sample_mat, cor_method = "spearman")
  topk_mat <- colwise_topk_intersect(sample_mat, k = k)
  bottomk_mat <- colwise_topk_intersect(sample_mat, k = k, decreasing = FALSE)
  jacc_mat <- colwise_jaccard(sample_mat, k = k)
  
  # Data frame of unique dataset pairs and their similarities
  sample_df <- get_similarity_pair_df(cor_mat, topk_mat, bottomk_mat, jacc_mat)
  
  return(sample_df)
}



set.seed(5)



sim_null_hg_path <- paste0("/space/scratch/amorin/R_objects/similarity_null_hg_k=", k, ".RDS")
sim_null_mm_path <- paste0("/space/scratch/amorin/R_objects/similarity_null_mm_k=", k, ".RDS")




if (!file.exists(sim_null_hg_path) || force_resave) {
  
  sim_hg <- lapply(1:n_samps, function(x) {
    
    message(paste("Human sample #", x, Sys.time()))
    
    calc_sample_similarity(agg_l = agg_tf_hg,
                           genes = tfs_hg$Symbol,
                           msr_mat = msr_hg,
                           k = k)
  })
  
  saveRDS(sim_hg, sim_null_hg_path)
  
} 



if (!file.exists(sim_null_mm_path) || force_resave) {
  
  sim_mm <- lapply(1:n_samps, function(x) {
    
    message(paste("Mouse sample #", x, Sys.time()))
    
    calc_sample_similarity(agg_l = agg_tf_mm,
                           genes = tfs_mm$Symbol,
                           msr_mat = msr_mm,
                           k = k)
  })
  
  saveRDS(sim_mm, sim_null_mm_path)
  
} 
