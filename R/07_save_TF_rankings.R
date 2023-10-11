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
agg_ribo_hg <- load_or_generate_agg(path = agg_ribo_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = ribo_genes$Symbol_hg)
agg_ribo_mm <- load_or_generate_agg(path = agg_ribo_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = ribo_genes$Symbol_mm)

tfs_hg <- filter(tfs_hg, Symbol %in% colnames(agg_tf_hg[[1]]))
tfs_mm <- filter(tfs_mm, Symbol %in% colnames(agg_tf_mm[[1]]))


# Functions
# ------------------------------------------------------------------------------



# Given a gene x experiment matrix, return a count vector of the times that each
# gene was in the top k of values for each column/experiment of gene_mat

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




# agg_l: list of gene x TF RSR matrices for each experiment
# msr_mat: gene x experiment binary matrix indicating gene measurement
# genes: the list of genes (assumed TFs) to generate a summary df for
# k: integer used as the index for the cut-off of the "top" of each list
# verbose: declare time for each calculated

# returns a list of dataframes for each gene in genes that contains:
# 1) Symbol: all genes in the rows of the RSR matrices of agg_l
# 2) N_comeasured: How many datasets the TF and gene were both measured in
# 3) Avg_RSR: Average RSR across experiments after imputation
# 4) Rank_RSR: Rank of the average RSR (higher RSR = lower/better rank)
# 5) Best_rank: The single best rank a gene obtained across experiments
# 6) Topk_count: Count of times a gene was in the top K across experiments
# 7) Topk_proportion: Topk_count divided by the count of comeasured experiments

# NOTE: When a TF-gene is not co-measured, its RSR is imputed to the median
# of the entire gene x experiment RSR matrix for the given TF (in practice 
# ~0.51). This is done to set these non-measured genes to the "middle of the 
# pack" as some will have an artificially high RSR (0.8+) for a given dataset 
# because ties (NA values) were encountered early in the ranking due to shallow
# measurement. This imputation down weights this artifact and still preserves 
# the bottom of the list for the most consistently negative TF-gene pairs.


all_rank_summary <- function(agg_l, 
                             msr_mat, 
                             genes, 
                             k = 1000,
                             verbose = TRUE) {
  
  summ_l <- lapply(genes, function(x) {
    
    if (verbose) message(paste(x, Sys.time()))
    
    # Gene x experiment matrix for given TF, and TF itself set to NA
    gene_mat <- gene_vec_to_mat(agg_l, x)
    gene_mat[x, ] <- NA
    gene_mat <- subset_to_measured(gene_mat, msr_mat = msr_mat, gene = x)
    
    if (length(gene_mat) == 0) {
      return(NA)
    }

    # Logical matrix of whether TF-genes were co-measured in an experiment
    comsr <- lapply(rownames(gene_mat), function(y) {
      msr_mat[x, colnames(gene_mat)] & msr_mat[y, colnames(gene_mat)]
    })
    comsr <- do.call(rbind, comsr)
    rownames(comsr) <- rownames(gene_mat)
    
    # Count of times the TF-gene were co-measured across experiments
    comsr_sum <- rowSums(comsr)
    
    # When a TF-gene was not co-measured, impute to the median NA
    med <- median(gene_mat, na.rm = TRUE)
    gene_mat[comsr == FALSE] <- med

    # Get a gene's ranking within each experiment and select the best
    colrank_gene_mat <- colrank_mat(gene_mat, ties_arg = "min")
    best_rank <- apply(colrank_gene_mat, 1, min)
    
    # Average the RSR and rank (higher RSR = more important = lower/better rank)
    avg_rsr <- rowMeans(gene_mat, na.rm = TRUE)
    rank_rsr <- rank(-avg_rsr, ties.method = "min")

    # Count/proportion of times a gene was in the top K across experiments  
    topk_count <- get_topk_count(gene_mat, k)
    topk_prop <- round(topk_count / comsr_sum, 3)
    
    # Organize summary df
    data.frame(
      Symbol = rownames(gene_mat),
      N_comeasured = comsr_sum,
      Avg_RSR = avg_rsr,
      Rank_RSR = rank_rsr,
      Best_rank = best_rank,
      Topk_count = topk_count,
      Topk_proportion = topk_prop)
    
  })
  
  names(summ_l) <- genes
  
  return(summ_l)
}



# Run and save
# ------------------------------------------------------------------------------



if (!file.exists(rank_tf_hg_path) || force_resave) {
  
  tf_hg <- all_rank_summary(agg_l = agg_tf_hg,
                            msr_mat = msr_hg,
                            genes = tfs_hg$Symbol,
                            k = k)
  
  tf_hg <- tf_hg[!is.na(tf_hg)]
  
  saveRDS(tf_hg, rank_tf_hg_path)

}



if (!file.exists(rank_tf_mm_path) || force_resave) {
  
  tf_mm <- all_rank_summary(agg_l = agg_tf_mm,
                            msr_mat = msr_mm,
                            genes = tfs_mm$Symbol,
                            k = k)
  
  tf_mm <- tf_mm[!is.na(tf_mm)]
  
  saveRDS(tf_mm, rank_tf_mm_path)
  
}



if (!file.exists(rank_ribo_hg_path) || force_resave) {
  
  ribo_hg <- all_rank_summary(agg_l = agg_ribo_hg, 
                              msr_mat = msr_hg, 
                              genes = ribo_genes$Symbol_hg,
                              k = k)
  
  ribo_hg <- ribo_hg[!is.na(ribo_hg)]
  
  saveRDS(ribo_hg, rank_ribo_hg_path)
  
}



if (!file.exists(rank_ribo_mm_path) || force_resave) {
  
  ribo_mm <- all_rank_summary(agg_l = agg_ribo_mm, 
                              msr_mat = msr_mm, 
                              genes = ribo_genes$Symbol_mm,
                              k = k)
  
  ribo_mm <- ribo_mm[!is.na(ribo_mm)]
  
  saveRDS(ribo_mm, rank_ribo_mm_path)
  
}
