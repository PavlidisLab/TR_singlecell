## 
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
library(aggtools)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 200

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

# # Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
# msr_mm <- readRDS(msr_mat_mm_path)

# Comeasure matrices tracking count of experiments gene pairs were comeasured
comsr_hg <- readRDS(comsr_mat_hg_path)
# comsr_mm <- readRDS(comsr_mat_mm_path)

# Loading the subset TF and ribosomal aggregate matrices
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
# agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)


# Final rankings
rank_tf_hg <- readRDS(rank_tf_hg_path)



# Count of times a TF was measured
# n_tf_msr <- msr_hg[tfs_hg$Symbol]


tf <- "E2F8"
tf_mat <- gene_vec_to_mat(agg_tf_hg, tf, msr_hg)
tf_mat[tf, ] <- NA  # prevent self cor from being #1 rank
# rank_tf_hg[[tf]] %>% head
tf_ids <- colnames(tf_mat)


# When a TF-gene was not co-measured, impute to the median NA
msr_mat <- msr_hg[, tf_ids]
med <- median(tf_mat, na.rm = TRUE)
tf_mat_imp <- tf_mat
tf_mat_imp[msr_mat == 0] <- med

# Average the aggr coexpr and rank (lower rank = positive correlation)
avg_aggr_coexpr <- rowMeans(tf_mat_imp, na.rm = TRUE)



# Individual experiments that are most similar


tt1 <- cbind(tf_mat_imp, Avg = avg_aggr_coexpr)
tt2 <- colwise_topk_intersect(tt1, k = k)
# tt3 <- mat_to_df(tt2, value_name = "Topk")
tt3 <- tt2[setdiff(rownames(tt2), "Avg"), "Avg"]
head(sort(tt3, decreasing = TRUE))

# Sample procedure
# NOTE: Sampling from the full imputed matrix. I tested imputing NAs from only
# the sampled experiments, which did not make a difference, so using full for
# simplicity!


steps <- 2:(length(tf_ids) - 1)
n_iter <- 100


step <- 14


# Use full matrix (impute from global)
topk_l1 <- lapply(1:n_iter, function(iter) {
  
  sample_ids <- sample(tf_ids, size = step, replace = FALSE)
  sample_avg <- rowMeans(tf_mat_imp[, sample_ids], na.rm = TRUE)
  
  topk_intersect(
    topk_sort(sample_avg, k = k),
    topk_sort(avg_aggr_coexpr, k = k)
  )
  
})


hist(unlist(topk_l1))
summary(unlist(topk_l1))




# Impute from sampled vectors only
# topk_l2 <- lapply(1:n_iter, function(iter) {
# 
#   sample_ids <- sample(tf_ids, size = step, replace = FALSE)
#   sample_msr_mat <- msr_hg[, sample_ids]
#   sample_tf_mat <- tf_mat[, sample_ids]
#   sample_med <- median(sample_tf_mat, na.rm = TRUE)
#   sample_tf_mat[sample_msr_mat == 0] <- sample_med
# 
#   # Average the aggr coexpr and rank (lower rank = positive correlation)
#   sample_avg <- rowMeans(sample_tf_mat, na.rm = TRUE)
# 
#   topk_intersect(
#     topk_sort(sample_avg, k = k),
#     topk_sort(avg_aggr_coexpr, k = k)
#   )
# 
# })
# 
# 
# hist(unlist(topk_l2))
# summary(unlist(topk_l2))



# check_df <- data.frame(coexpr = avg_aggr_coexpr) %>% 
#   rownames_to_column(var = "Symbol") %>%
#   left_join(., rank_tf_hg[[tf]], by = "Symbol")
  
