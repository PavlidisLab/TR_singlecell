## Sample TFs across experiments and calculate their similarity to generate
## a null. Save out the results as a list.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 200
n_samps <- 1000
force_resave <- FALSE
set.seed(5)

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes, TFs, and L/S ribo genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
tfs_hg <- read.delim(tfs_hg_path, stringsAsFactors = FALSE)
tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)
ribo_genes <- read.delim(ribo_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Loading the subset TF and ribosomal aggregate matrices
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)
agg_ribo_hg <- load_or_generate_agg(path = agg_ribo_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = ribo_genes$Symbol_hg)
agg_ribo_mm <- load_or_generate_agg(path = agg_ribo_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = ribo_genes$Symbol_mm)

# Output paths of similarity lists
# TODO: pathing in config
sim_null_hg_path <- paste0("/space/scratch/amorin/R_objects/TRsc/similarity_null_hg_k=", k, ".RDS")
sim_null_mm_path <- paste0("/space/scratch/amorin/R_objects/TRsc/similarity_null_mm_k=", k, ".RDS")
sim_tf_hg_path <- paste0("/space/scratch/amorin/R_objects/TRsc/similarity_TF_hg_k=", k, ".RDS")
sim_tf_mm_path <- paste0("/space/scratch/amorin/R_objects/TRsc/similarity_TF_mm_k=", k, ".RDS")
sim_ribo_hg_path <- paste0("/space/scratch/amorin/R_objects/TRsc/similarity_ribo_hg_k=", k, ".RDS")
sim_ribo_mm_path <- paste0("/space/scratch/amorin/R_objects/TRsc/similarity_ribo_mm_k=", k, ".RDS")



# For a given set of input genes and list of aggregate coexpression matrices,
# calculate multiple metrics of similarity within the provided genes across networks. 
# agg_l: A list of aggregate coexpression matrices (gene x gene numeric matrix)
# msr_mat: binary gene x experiment matrix tracking if a gene was measured
# k: an integer of the top/bottom elements to select
# returns: a list for each gene containing a dataframe of all unique experiment
# pairs in which the given gene was measured, and three measures of similarity:
# 1) The spearman correlation between the gene vectors of the experiment pair;
# 2) The size of the top k intersect between the experiments;
# 3) The size of the bottom k intersect between the experiments;
# 4) The Jaccard of the top and bottom k elements between the experiments
# ------------------------------------------------------------------------------




# Generates the similarity matrices for the input gene x exp matrix and format
# the resulting symmetric exp x exp similarity matrices into a summary df of
# the unique pairs.

gen_similarity_df <- function(gene_mat, k) {
  
  # Experiment x experiment similarity matrices
  cor_mat <- colwise_cor(gene_mat, cor_method = "spearman")
  topk_mat <- colwise_topk_intersect(gene_mat, k = k)
  bottomk_mat <- colwise_topk_intersect(gene_mat, k = k, decreasing = FALSE)
  jacc_mat <- colwise_jaccard(gene_mat, k = k)
    
  # Extract unique elements of similarity matrices into a single long df
  df <-
    mat_to_df(cor_mat, symmetric = TRUE, value_name = "Scor") %>%
    mutate(
      Topk = mat_to_df(topk_mat, symmetric = TRUE)$Value,
      Bottomk = mat_to_df(bottomk_mat, symmetric = TRUE)$Value,
      Jaccard = mat_to_df(jacc_mat, symmetric = TRUE)$Value
    )
  
  return(df)
}



# Main function to retrieve the gene list of summary dfs of exp similarity pairs

calc_gene_similarity <- function(agg_l,
                                 msr_mat,
                                 genes,
                                 k,
                                 ncores = 1) {
  
  stopifnot(genes %in% colnames(agg_l[[1]]))
  
  sim_l <- mclapply(genes, function(x) {
    
    message(paste(x, Sys.time()))
    
    # Bind the gene profiles into a single matrix. Remove self to prevent inflated overlap
    gene_mat <- gene_vec_to_mat(agg_l, gene = x, msr_mat = msr_mat)
    gene_mat <- gene_mat[setdiff(rownames(gene_mat), x), , drop = FALSE]

    # need more than one experiment for pairing
    if (is.null(gene_mat) || ncol(gene_mat) < 2) {  
      return(NA)
    }
    
    # Data frame of unique dataset pairs and their similarities
    gen_similarity_df(gene_mat, k)
    
  }, mc.cores = ncores)
  
  names(sim_l) <- genes
  sim_l <- sim_l[!is.na(sim_l)]
  
  return(sim_l)
}




# This samples from the input genes to generate null similarities 

calc_sample_similarity <- function(agg_l,
                                   msr_mat,
                                   genes,
                                   k, 
                                   n_samps,
                                   ncores = 1) {
  
  ids <- names(agg_l)
  
  sim_l <- mclapply(1:n_samps, function(i) {
    
    message(paste("Sample #", i, Sys.time()))
    
    # Sample one gene that is measured in each data set
    sample_genes <- vapply(ids, function(x) {
      msr_gene <- names(which(msr_mat[genes, x] == 1))
      sample(msr_gene, 1)
    }, FUN.VALUE = character(1))
    
    # Bind sampled genes into a matrix
    sample_mat <- lapply(1:length(sample_genes), function(x) {
      agg_l[[x]][, sample_genes[x]]
    })
    sample_mat <- do.call(cbind, sample_mat)
    colnames(sample_mat) <- paste0(ids, "_", sample_genes)
    
    # Data frame of unique dataset pairs and their similarities
    gen_similarity_df(sample_mat, k)
    
  }, mc.cores = ncores)
  
  return(sim_l)
 
}




# Handle saving RDS output of either similarity function if it does not exist 

save_similarity_results <- function(path, 
                                    fun, 
                                    args, 
                                    force_resave = FALSE) {
  
  if (!file.exists(path) || force_resave) {
    
    result <- do.call(fun, args)
    
    if (!is.null(result)) {
      saveRDS(result, path)
    }
    
    return(invisible(NULL))
  }
}



# Run and save
# ------------------------------------------------------------------------------



# Human TFs
save_similarity_results(
  path = sim_tf_hg_path,
  fun = calc_gene_similarity,
  args = list(
    agg_l = agg_tf_hg,
    msr_mat = msr_hg,
    genes = tfs_hg$Symbol,
    k = k,
    ncores = ncore
  ),
  force_resave = force_resave
)



# Mouse TFs
save_similarity_results(
  path = sim_tf_mm_path,
  fun = calc_gene_similarity,
  args = list(
    agg_l = agg_tf_mm,
    msr_mat = msr_mm,
    genes = tfs_mm$Symbol,
    k = k,
    ncores = ncore
  ),
  force_resave = force_resave
)



# Human ribo
save_similarity_results(
  path = sim_ribo_hg_path,
  fun = calc_gene_similarity,
  args = list(
    agg_l = agg_ribo_hg,
    msr_mat = msr_hg,
    genes = ribo_genes$Symbol_hg,
    k = k,
    ncores = ncore
  ),
  force_resave = force_resave
)



# Mouse ribo
save_similarity_results(
  path = sim_ribo_mm_path,
  fun = calc_gene_similarity,
  args = list(
    agg_l = agg_ribo_mm,
    msr_mat = msr_mm,
    genes = ribo_genes$Symbol_mm,
    k = k,
    ncores = ncore
  ),
  force_resave = force_resave
)


# Human null
save_similarity_results(
  path = sim_null_hg_path,
  fun = calc_sample_similarity,
  args = list(
    agg_l = agg_tf_hg,
    msr_mat = msr_hg,
    genes = tfs_hg$Symbol,
    k = k,
    n_samps = n_samps,
    ncores = ncore
  ),
  force_resave = force_resave
)



# Mouse null
save_similarity_results(
  path = sim_null_mm_path,
  fun = calc_sample_similarity,
  args = list(
    agg_l = agg_tf_mm,
    msr_mat = msr_mm,
    genes = tfs_mm$Symbol,
    k = k,
    n_samps = n_samps,
    ncores = ncore
  ),
  force_resave = force_resave
)
