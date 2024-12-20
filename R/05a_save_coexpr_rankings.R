## For each TF generate a dataframe ranking genes by their aggregate coexpression
## with the given TF, and save out as a list
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 200
force_resave <- FALSE

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes, TFs, and L/S ribo genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path, stringsAsFactors = FALSE)
tfs_hg <- read.delim(tfs_hg_path, stringsAsFactors = FALSE)
tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)
ribo_genes <- read.delim(ribo_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Comeasure matrices tracking count of experiments gene pairs were comeasured
comsr_hg <- readRDS(comsr_mat_hg_path)
comsr_mm <- readRDS(comsr_mat_mm_path)

# Loading the subset TF and ribosomal aggregate matrices
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)
agg_ribo_hg <- load_or_generate_agg(path = agg_ribo_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = ribo_genes$Symbol_hg)
agg_ribo_mm <- load_or_generate_agg(path = agg_ribo_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = ribo_genes$Symbol_mm)


# Functions
# ------------------------------------------------------------------------------


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



# agg_l: list of gene x TF aggr coexpr matrices for each experiment
# msr_mat: gene x experiment binary matrix indicating gene measurement
# comsr_mat: gene x gene matrix count gene pair comeasurement across experiments
# genes: the list of genes to generate a rank summary df for each
# k: integer used as the index for the cut-off of the "top" of each list
# verbose: declare time for each calculated

# returns a list of dataframes for each gene in genes that contains:
# 1) Symbol: all genes in the rows of the aggr coexpr matrices of agg_l
# 2) N_comeasured: How many datasets the TF and gene were both measured in
# 3) Avg_aggr_coexpr: Average aggr coexpr across experiments after imputation
# 4) Rank_aggr_coexpr: Rank of the aggr coexpr (higher coexpr = lower rank)
# 5) Best_rank: The single best rank a gene obtained across experiments
# 6) Topk_count: Count of times a gene was in the top K across experiments
# 7) Topk_proportion: Topk_count divided by the count of comeasured experiments

# --- NOTE about median imputation: 
# When a TF-gene is not co-measured, its aggr coexpr is imputed to the median
# of the entire gene x experiment aggr coexpr matrix for the given TF (in 
# practice ~0.51). This is done to set these non-measured genes to the "middle 
# of the pack" as some will have an artificially high aggr coexpr (0.8+) for a 
# given dataset because ties (NA values) were encountered early in the ranking.
# This imputation down weights this artifact and still preserves the bottom of
# the list for the most consistently negative TF-gene pairs.

gen_rank_list <- function(agg_l,
                          msr_mat,
                          comsr_mat,
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

    # When a TF-gene was not co-measured, impute to the median NA
    msr_mat <- msr_mat[, colnames(gene_mat)]
    med <- median(gene_mat, na.rm = TRUE)
    gene_mat[msr_mat == 0] <- med
    
    # Get count of times each gene was co-measured with the TF
    comsr_sum <- comsr_mat[x, rownames(gene_mat)]
    
    # Get a gene's ranking within each experiment and select the best
    colrank_gene_mat <- colrank_mat(gene_mat, ties_arg = "min")
    best_rank <- apply(colrank_gene_mat, 1, min)
    
    # Average the aggr coexpr and rank (lower rank = positive correlation)
    avg_aggr_coexpr <- rowMeans(gene_mat, na.rm = TRUE)
    rank_aggr_coexpr <- rank(-avg_aggr_coexpr, ties.method = "min")

    # Count of times a gene was in the top K across experiments  
    topk_count <- calc_topk_count(gene_mat, k)
    topk_var <- paste0("Top", k, "_count")
    
    # Organize summary df
    rank_df <- data.frame(
      Symbol = rownames(gene_mat),
      Rank_aggr_coexpr = rank_aggr_coexpr,
      N_comeasured = comsr_sum,
      Avg_aggr_coexpr = avg_aggr_coexpr,
      Rank_single_best = best_rank,
      Topk_count = topk_count
    )
    
    # Actual k value in col name and set TF's self stats (self-cor) to NAs
    rank_df <- rank_df %>% 
      dplyr::rename(!!(paste0("Top", k, "_count")) := Topk_count) %>% 
      arrange(Rank_aggr_coexpr) 
    
    rank_df[rank_df$Symbol == x, 
            setdiff(colnames(rank_df), c("Symbol", "N_comeasured"))] <- NA
      
    return(rank_df)
  })
  
  names(rank_l) <- genes
  rank_l <- rank_l[!is.na(rank_l)]
  
  return(rank_l)
}



# Gene is assumed to be an ortho gene found in pc_df. Retrieve and join the 
# corresponding gene rank dfs from the mouse and human rank lists. Re-rank
# average aggr coexpr just using ortho genes

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
      mutate(
        Rank_aggr_coexpr = rank(-Avg_aggr_coexpr, 
                                ties.method = "min", 
                                na.last = "keep"))
  
    df_mm <-
      left_join(rank_mm[[gene_ortho$Symbol_mm]],
                pc_ortho[, c("Symbol_mm", "ID")],
                by = c("Symbol" = "Symbol_mm")) %>%
      filter(!is.na(ID)) %>%
      mutate(Rank_aggr_coexpr = rank(-Avg_aggr_coexpr, 
                                     ties.method = "min", 
                                     na.last = "keep"))
  
    # Join and create a rank that combines each species ranking
    df_ortho <-
      left_join(df_hg, df_mm,
                by = "ID",
                suffix = c("_hg", "_mm")) %>%
      dplyr::select(-c(ID)) %>% 
      mutate(
        Rank_aggr_coexpr_ortho = rank(
          Rank_aggr_coexpr_hg * Rank_aggr_coexpr_mm,
          ties.method = "min",
          na.last = "keep"
        )
      ) %>% 
      arrange(Rank_aggr_coexpr_ortho) %>% 
      relocate(
        Symbol_hg,
        Symbol_mm,
        Rank_aggr_coexpr_ortho,
        Rank_aggr_coexpr_hg,
        Rank_aggr_coexpr_mm
      )
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
    comsr_mat = comsr_hg,
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
    comsr_mat = comsr_mm,
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
    comsr_mat = comsr_hg,
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
    comsr_mat = comsr_mm,
    genes = ribo_genes$Symbol_mm,
    k = k,
    ncores = ncore
  ),
  force_resave = force_resave
)
