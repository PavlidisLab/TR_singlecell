## For each TF generate a dataframe ranking genes by their aggregate coexpression
## with the given TF, and save out as a list
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
msr_mat_hg_path <- "/space/scratch/amorin/R_objects/binary_measurement_matrix_hg_lognorm.RDS"
msr_mat_mm_path <- "/space/scratch/amorin/R_objects/binary_measurement_matrix_mm_lognorm.RDS"
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# For each dataset, load the subset gene x TF aggregation matrix 
agg_tf_hg_path <- "/space/scratch/amorin/R_objects/TF_agg_mat_list_hg_lognorm.RDS"
agg_tf_mm_path <- "/space/scratch/amorin/R_objects/TF_agg_mat_list_mm_lognorm.RDS"
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol, pattern = "_RSR_allrank.tsv")
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol, pattern = "_RSR_allrank.tsv")



# TODO:
# ------------------------------------------------------------------------------




# TODO:

get_topk_count <- function(gene_mat, k, check_k_arg = TRUE) {
  
  bin_l <- lapply(1:ncol(gene_mat), function(x) {
    
    vec_sort <- sort(gene_mat[, x], decreasing = TRUE)
    if (check_k_arg) k <- check_k(vec_sort, k)
    gene_mat[, x] >= vec_sort[k]
    
  })
  
  bin_mat <- do.call(cbind, bin_l)
  k_count <- rowSums(bin_mat, na.rm = TRUE)
  
  return(k_count)
}



# TODO: make explicit that proportion is count of datasets in which both
# genes are measured, not just TF

# For each gene in genes, generate a summary dataframe that includes:
# TODO: update
# 1) The gene ranking of coexpressed partners (average aggregate rank)
# 2) The count of times the gene pair was mutually measured
# 3) The count of NAs (no mutual measurement)
# 4) For the best correlation partner, what was its best rank across datasets?
# 5) And in what dataset did it achieve this best rank?

all_rank_summary <- function(agg_l, 
                             msr_mat, 
                             genes, 
                             k = 1000) {
  
  summ_l <- lapply(genes, function(x) {
    
    message(paste(x))
    
    gene_mat <- gene_vec_to_mat(agg_l, x)
    gene_mat[x, ] <- NA
    gene_mat <- subset_to_measured(gene_mat, msr_mat = msr_mat, gene = x)
    
    if (length(gene_mat) == 0) {
      return(NA)
    }
    
    # Co-measurement between TF and genes
    comsr <- vapply(rownames(gene_mat), function(y) {
      sum(msr_mat[x,] & msr_mat[y,])
    }, integer(1))
    
    
    colrank_gene_mat <- colrank_mat(gene_mat, ties_arg = "min")
    
    avg_rsr <- rowMeans(gene_mat, na.rm = TRUE)
    # avg_colrank <- rowMeans(colrank_gene_mat, na.rm = TRUE)
    
    rank_rsr <- rank(-avg_rsr, ties.method = "min")
    # rank_colrank <- rank(avg_colrank, ties.method = "min")
    
    best_rank <- apply(colrank_gene_mat, 1, min)
    
    topk_count <- get_topk_count(gene_mat, k)
    topk_prop <- round(topk_count / comsr, 3)
    
    data.frame(
      Symbol = rownames(gene_mat),
      N_comeasured = comsr,
      Avg_RSR = avg_rsr,
      # Avg_colrank = avg_colrank,
      Rank_RSR = rank_rsr,
      # Rank_colrank = rank_colrank,
      Best_rank = best_rank,
      Topk_count = topk_count,
      Topk_proportion = topk_prop)
    
  })
  
  names(summ_l) <- genes
  
  return(summ_l)
}




rank_tf_hg_path <- "/space/scratch/amorin/R_objects/TF_agg_ranking_hg_lognorm.RDS"
rank_tf_mm_path <- "/space/scratch/amorin/R_objects/TF_agg_ranking_mm_lognorm.RDS"



if (!file.exists(rank_tf_hg_path) || force_resave) {
  
  summ_hg <- all_rank_summary(agg_l = agg_tf_hg, 
                              msr_mat = msr_hg, 
                              genes = tfs_hg$Symbol,
                              k = k)
  
  summ_hg <- summ_hg[!is.na(summ_hg)]
  
  saveRDS(summ_hg, rank_tf_hg_path)
  
}



if (!file.exists(rank_tf_mm_path) || force_resave) {
  
  summ_mm <- all_rank_summary(agg_l = agg_tf_mm, 
                              msr_mat = msr_mm, 
                              genes = tfs_mm$Symbol,
                              k = k)
  
  summ_mm <- summ_mm[!is.na(summ_mm)]
  
  saveRDS(summ_mm, rank_tf_mm_path)
  
}
