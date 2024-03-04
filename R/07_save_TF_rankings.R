## For each TF generate a dataframe ranking genes by their aggregate coexpression
## with the given TF, and save out as a list
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 500
force_resave <- TRUE

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes, TFs, and L/S ribo genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)
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


# Functions
# ------------------------------------------------------------------------------


# Assumes gene mat is a gene x experiment matrix of a given TF's coexpression
# profiles. Returns a logical gene x experiment matrix tracking if the given
# TF was co-measured with each gene in each experiment

gen_comsr_mat <- function(tf, gene_mat, msr_mat) {
  
  exps <- colnames(gene_mat)
  genes <- rownames(gene_mat)
  
  tf_msr <- matrix(rep(msr_mat[tf, exps], each = length(genes)), 
                   nrow = length(genes))

  gene_msr <- msr_mat[genes, exps, drop = FALSE]
  
  
  comsr <- tf_msr & gene_msr
  
  return(comsr)
}



# Given a gene x experiment matrix, return a count vector of the times that each
# gene was in the top k of values for each column/experiment of gene_mat

calc_topk_count <- function(gene_mat, k, check_k_arg = TRUE) {
  
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

gen_rank_list <- function(agg_l,
                          msr_mat,
                          genes,
                          k,
                          ncores = 1,
                          verbose = TRUE) {
  
  rank_l <- mclapply(genes, function(x) {
    
    if (verbose) message(paste(x, Sys.time()))
    
    # Gene x experiment aggregate coexpression matrix for current TF
    gene_mat <- gene_vec_to_mat(agg_l, gene = x, msr_mat = msr_mat)
    gene_mat[x, ] <- NA  # prevent self cor from being #1 rank
    
    if (length(gene_mat) == 0) {
      return(NA)
    }

    # Get count of times each gene was co-measured with the TF
    comsr <- gen_comsr_mat(x, gene_mat, msr_mat)
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

    # Count of times a gene was in the top K across experiments  
    topk_count <- calc_topk_count(gene_mat, k)
    topk_var <- paste0("Top", k, "_count")
    
    # Organize summary df
    rank_df <- data.frame(
      Symbol = rownames(gene_mat),
      Rank_aggr_coexpr = rank_rsr,
      N_comeasured = comsr_sum,
      Avg_aggr_coexpr = avg_rsr,
      Rank_single_best = best_rank,
      Topk_count = topk_count
    )
    
    colnames(rank_df)[colnames(rank_df) == "Topk_count"] <- paste0("Top", k, "_count")
    rank_df <- arrange(rank_df, Rank_aggr_coexpr)
    
    return(rank_df)
  })
  
  names(rank_l) <- genes
  rank_l <- rank_l[!is.na(rank_l)]
  
  return(rank_l)
}



# Gene is assumed to be an ortho gene found in pc_df. Retrieve and join the 
# corresponding gene rank dfs from the mouse and human rank lists. Re-rank
# average RSR just using ortho genes

join_ortho_ranks <- function(pc_ortho,
                             rank_hg,
                             rank_mm,
                             ncores = 1) {
  
  # Ortho TFs with rankings in both species
  tf_ortho <- filter(pc_ortho, 
                     Symbol_hg %in% names(rank_tf_hg) & 
                     Symbol_mm %in% names(rank_tf_mm))
  
  # Create a list of ranked dfs for each ortho TF
  rank_ortho <- mclapply(tf_ortho$Symbol_hg, function(x) {
    
    gene_ortho <- filter(pc_ortho, Symbol_hg == x | Symbol_mm == x)
    
    # Prepare species specific rankings, re-ranking after filtering for ortho
    df_hg <-
      left_join(rank_hg[[gene_ortho$Symbol_hg]],
                pc_ortho[, c("Symbol_hg", "ID")],
                by = c("Symbol" = "Symbol_hg")) %>%
      filter(!is.na(ID)) %>%
      mutate(Rank_RSR = rank(-Avg_RSR, ties.method = "min"))
  
    df_mm <-
      left_join(rank_mm[[gene_ortho$Symbol_mm]],
                pc_ortho[, c("Symbol_mm", "ID")],
                by = c("Symbol" = "Symbol_mm")) %>%
      filter(!is.na(ID)) %>%
      mutate(Rank_RSR = rank(-Avg_RSR, ties.method = "min"))
  
    # Join and create a rank that combines each species ranking
    df_ortho <-
      left_join(df_hg, df_mm,
                by = "ID",
                suffix = c("_hg", "_mm")) %>%
      filter((!is.na(Avg_RSR_hg) & !is.na(Avg_RSR_mm))) %>% 
      dplyr::select(-c(ID)) %>% 
      mutate(Rank_aggr_coexpr_ortho = rank(Rank_RSR_hg * Rank_RSR_mm)) %>% 
      arrange(Rank_aggr_coexpr_ortho) %>% 
      relocate(Symbol_hg, Symbol_mm, Rank_aggr_coexpr_ortho, Rank_RSR_hg, Rank_RSR_mm)
      
    })
    
    names(rank_ortho) <- tf_ortho$Symbol_hg
    return(rank_ortho)
}



# Run and save
# ------------------------------------------------------------------------------


# Human TFs
save_function_results(
  path = rank_tf_hg_path,
  fun = gen_rank_list,
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
save_function_results(
  path = rank_tf_mm_path,
  fun = gen_rank_list,
  args = list(
    agg_l = agg_tf_mm,
    msr_mat = msr_mm,
    genes = tfs_mm$Symbol,
    k = k,
    ncores = ncore
  ),
  force_resave = force_resave
)



# Join ortho
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

save_function_results(
  path = rank_tf_ortho_path,
  fun = join_ortho_ranks,
  args = list(
    pc_ortho = pc_ortho,
    rank_hg = rank_tf_hg,
    rank_mm = rank_tf_mm,
    ncores = ncore
  ),
  force_resave = force_resave
)



# Human ribo
save_function_results(
  path = rank_ribo_hg_path,
  fun = gen_rank_list,
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
save_function_results(
  path = rank_ribo_mm_path,
  fun = gen_rank_list,
  args = list(
    agg_l = agg_ribo_mm,
    msr_mat = msr_mm,
    genes = ribo_genes$Symbol_mm,
    k = k,
    ncores = ncore
  ),
  force_resave = force_resave
)
