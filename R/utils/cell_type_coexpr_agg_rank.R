## TODO: fix inconsistent dims of cor->rank
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(parallel)
source("R/utils/functions.R")
source("R/00_config.R")

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Load correlation matrices generated per cell type and across all cells
cor_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_celltype.RDS")
cor_all <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_all.RDS")

# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)


#

get_cor_all_df <- function(cmat, tf, rm_tf = TRUE) {
  
  if (rm_tf) {
    cmat <- cmat[setdiff(rownames(cmat), tf), ]
  }
  
  df <- 
    data.frame(Cor_all = cmat[, tf]) %>%
    rownames_to_column(var = "Symbol") %>%
    mutate(
      Cor_all_abs = abs(Cor_all),
      Rank_cor_all = rank(-Cor_all, ties.method = "min"),
      Rank_cor_all_abs = rank(-Cor_all_abs, ties.method = "min")
    )
  
  return(df)
}


# Given a list of cor matrices per cell type and a specified TF, get that TFs
# gene cors in a gene x cell type matrix

tf_by_ct_cmat <- function(cor_list, tf, rm_tf = TRUE) {
  
  stopifnot(identical(rownames(cor_list[[1]]), rownames(cor_list[[2]])))
  
  cor_tf <- do.call(cbind, lapply(cor_ct, function(x) x[tf, ]))
  
  if (rm_tf) {
    cor_tf <- cor_tf[rownames(cor_tf) != tf, ]
  }
  
  return(cor_tf)
}


# NA counts

# na_count_mat <- function(cor_ct_list, ncores = ncore) {
#   
#   na_l <- mclapply(cor_ct_list, function(mat) {
#     apply(mat, 1, function(x) sum(is.na(x)))
#   }, mc.cores = ncore)
#   
#   mat <- do.call(rbind, na_l)
#   
#   return(mat)
# }


# For the given TF, return a df of the count of times the cor for each gene-TF 
# pair was NA across cell types

tf_na_count <- function(cor_ct_list, tf) {
  
  na_l <- lapply(cor_ct_list, function(x) {
    is.na(x[tf, ])
  })
  
  count_na <- data.frame(Count_NA = rowSums(do.call(cbind, na_l))) %>% 
    rownames_to_column(var = "Symbol")
  
  return(count_na)
}


# NA over everything
tt <- apply(simplify2array(cor_ct), 1:2, function(x) sum(is.na(x)))



# Return a df of the max cor for each cell type for the given TR. NAs removed

max_cor_df <- function(cor_ct_list, tf) {
  
  cor_max <- lapply(names(cor_ct_list), function(x) {
    
    cor_mat <- cor_ct_list[[x]]
    
    vec <- cor_mat[tf, setdiff(colnames(cor_mat), tf)]
    
    if (all(is.na(vec))) {
      return(NA)
    }
    
    data.frame(
      Cell_type = x,
      Value = max(vec, na.rm = TRUE),
      Symbol = names(vec)[which.max(vec)])
  })
  
  cor_max <- cor_max[!is.na(cor_max)]
  
  cor_max <- data.frame(do.call(rbind, cor_max)) %>% 
    arrange(desc(Value))
  
  return(cor_max)
}






###


cor_ct_nato0 <- lapply(cor_ct, function(x) {
  x[is.na(x)] <- 0
  return(x)
})


cor_ct_rank <- mclapply(cor_ct, function(x) {
  rank_cormat(t(x))
}, mc.cores = ncore)


cor_ct_rank_nato0 <- mclapply(cor_ct_nato0, function(x) {
  rank_cormat(t(x))
}, mc.cores = ncore)


# Ensure max cor equals min (best) rank

assertthat::are_equal(
  names(which.max(cor_ct[[1]][1, ])),
  names(which.min(cor_ct_rank[[1]][, 1])))

assertthat::are_equal(
  names(which.max(cor_ct_nato0[[1]][1, ])),
  names(which.min(cor_ct_rank_nato0[[1]][, 1])))

# NAs ranked last, want min non-NA ranking

assertthat::are_equal(
  names(which.min(cor_ct[[1]][1, ])),
  names(which.max(cor_ct_rank[[1]][, 1])))

assertthat::are_equal(
  names(which.min(cor_ct_nato0[[1]][1, ])),
  names(which.max(cor_ct_rank_nato0[[1]][, 1])))



# Aggregate cell type rank matrices
# https://stackoverflow.com/questions/42628385/sum-list-of-matrices-with-nas
# Sum ignores NAs (treated as 0), but converts sum==0 (all NAs) to NA
# Mean_nona divides by the number of non-NA observations
# Mean_all divides by the total number of cell types


cor_ct_rank_sum <- apply(simplify2array(cor_ct_rank), 1:2, sum, na.rm = TRUE)
cor_ct_rank_nato0_sum <- apply(simplify2array(cor_ct_rank_nato0), 1:2, sum, na.rm = TRUE)
cor_ct_rank_0tona_sum <- cor_ct_rank_sum
cor_ct_rank_0tona_sum[cor_ct_rank_0tona_sum == 0] <- NA

cor_ct_rank_mean <- apply(simplify2array(cor_ct_rank), 1:2, mean, na.rm = TRUE)
cor_ct_rank_nato0_mean <- apply(simplify2array(cor_ct_rank_nato0), 1:2, mean, na.rm = TRUE)
cor_ct_rank_all_mean <- cor_ct_rank_sum / length(cor_ct_rank)


assertthat::are_equal(
  sum(sapply(cor_ct_rank, function(x) x[1, 1]), na.rm = TRUE),
  cor_ct_rank_sum[1, 1])


assertthat::are_equal(
  mean(sapply(cor_ct_rank, function(x) x[1, 1]), na.rm = TRUE),
  cor_ct_rank_mean[1, 1])


# Then convert these aggregations to ranks. Mean_all and sum have same rank.

cor_agg_sum <- rank_cormat(-cor_ct_rank_sum)
cor_agg_nato0_sum <- rank_cormat(-cor_ct_rank_nato0_sum)
cor_agg_0tona_sum <- rank_cormat(-cor_ct_rank_0tona_sum)

cor_agg_mean <- rank_cormat(-cor_ct_rank_mean)
cor_agg_nato0_mean <- rank_cormat(-cor_ct_rank_nato0_mean)
cor_agg_all_mean <- rank_cormat(-cor_ct_rank_all_mean)

cor_agg_nato0_final <- agg_cor_ct(cor_ct, zero_nas = TRUE, ncores = ncore)
cor_agg_keepna_final <- agg_cor_ct(cor_ct, zero_nas = FALSE, ncores = ncore)



# cor_ct[[1]][1:5, 1:5]
# cor_ct_nato0[[1]][1:5, 1:5]
# 
# cor_ct_rank[[1]][1:5, 1:5]
# cor_ct_rank_nato0[[1]][1:5, 1:5]

# cor_ct_rank_sum[1:5, 1:5]
# cor_ct_rank_nato0_sum[1:5, 1:5]
# cor_ct_rank_0tona_sum[1:5, 1:5]
# cor_ct_rank_mean[1:5, 1:5]
# cor_ct_rank_nato0_mean[1:5, 1:5]
# cor_ct_rank_all_mean[1:5, 1:5]

# cor_agg_sum[1:5, 1:5]
# cor_agg_nato0_sum[1:5, 1:5]
# cor_agg_0tona_sum[1:5, 1:5]
# cor_agg_mean[1:5, 1:5]
# cor_agg_nato0_mean[1:5, 1:5]
# cor_agg_all_mean[1:5, 1:5]

# cor_agg_nato0_final[1:5, 1:5]
# cor_agg_keepna_final[1:5, 1:5]



# Focus on single TF and organize gene rankings
# ------------------------------------------------------------------------------


tf <- "Ascl1"


# Which TF-gene pair had max cor in each cell type
max_cor <- max_cor_df(cor_ct, tf)


# Gene rankings by different aggregations in single data frame
rank_df <- data.frame(
  Symbol = rownames(cor_agg_sum),
  Sum = cor_ct_rank_sum[, tf],
  Sum_nato0 = cor_ct_rank_nato0_sum[, tf],
  Sum_0tona = cor_ct_rank_0tona_sum[, tf],
  Mean = cor_ct_rank_mean[, tf],
  Mean_nato0 = cor_ct_rank_nato0_mean[, tf],
  Mean_all = cor_ct_rank_all_mean[, tf],
  Agg_sum = cor_agg_sum[, tf],
  Agg_sum_nato0 = cor_agg_nato0_sum[, tf],
  Agg_sum_0tona = cor_agg_0tona_sum[, tf],
  Agg_mean = cor_agg_mean[, tf],
  Agg_mean_nato0 = cor_agg_nato0_mean[, tf],
  Agg_mean_all = cor_agg_all_mean[, tf],
  Agg_final_nato0 = cor_agg_nato0_final[, tf],
  Agg_final_keepna = cor_agg_keepna_final[, tf]
) %>%
  left_join(tf_na_count(cor_ct, tf), by = "Symbol") %>% 
  arrange(Agg_final_nato0)


# Correlation of the ranks and the count of NAs
rank_df_cor <- cor(select_if(rank_df, is.numeric), use = "pairwise.complete.obs")


# Gene NA count. When equal to number of cell types, it means that those genes
# and the TF were never expressed in the same cell type. The minimum NA count
# reflects the count of cell types in which the TF was not detected. Genes
# equal to this count have a non-NA cor in every cell type with TF detection.

hist(rank_df$Count_NA, breaks = 100)
sum(rank_df$Count_NA == n_distinct(sdat$Cell_type))
sum(rank_df$Count_NA == min(rank_df$Count_NA))

filter(rank_df, Count_NA < 120) %>% head()
filter(rank_df, Count_NA == min(rank_df$Count_NA)) %>% head()


#
# ------------------------------------------------------------------------------


gene2 <- "Ank2"
ct <- "VGLUT2-40-Otp_Pnoc"

sort(sapply(cor_ct, function(x) x[tf, gene2]), decreasing = TRUE)
sort(sapply(cor_ct_rank, function(x) x[gene2, tf]))
head(sort(cor_ct[[ct]][tf, ], decreasing = TRUE))


# By cell type
sdat_sub <- subset(sdat, idents = ct)
plot_scatter(sdat_sub, tf, gene2, slot = "data", jitter = TRUE)
plot_scatter(sdat_sub, tf, gene2, slot = "data", jitter = FALSE)
plot_scatter(sdat_sub, tf, gene2, slot = "counts", jitter = TRUE)
plot_scatter(sdat_sub, tf, gene2, slot = "counts", jitter = FALSE)
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "data")[c(tf, gene2), ])))
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "counts")[c(tf, gene2), ])))

sum(sdat@assays$RNA@counts[gene2, sdat$Cell_type == ct] > 0)

sum(sdat@assays$RNA@counts[tf, sdat$Cell_type == ct] > 0)

sum(sdat@assays$RNA@counts[tf, sdat$Cell_type == ct] > 0 & 
      sdat@assays$RNA@counts[gene2, sdat$Cell_type == ct] > 0)

# Across all cells
plot_scatter(sdat, tf, gene2, slot = "data", jitter = TRUE)
plot_scatter(sdat, tf, gene2, slot = "data", jitter = FALSE)
plot_scatter(sdat, tf, gene2, slot = "counts", jitter = TRUE)
plot_scatter(sdat, tf, gene2, slot = "counts", jitter = FALSE)
cor_all[tf, gene2]
cor(t(as.matrix(GetAssayData(object = sdat, slot = "counts")[c(tf, gene2), ])))


###

rank_df <- rank_l$Mouse[[tf]] %>% 
  dplyr::select(Symbol, Rank_integrated, Curated_target) %>% 
  left_join(rank_df, by = "Symbol") %>% 
  filter(Symbol != tf)


keep_cols <- c("Rank_integrated", "Agg_final_nato0", "Agg_final_keepna")


pr_df <- all_perf_df(rank_df, keep_cols, label_col = "Curated_target", measure = "PR")
auprc <- all_au_perf(rank_df, keep_cols, label_col = "Curated_target", measure = "AUPRC")

roc_df <- all_perf_df(rank_df, keep_cols, label_col = "Curated_target", measure = "ROC")
auroc <- all_au_perf(rank_df, keep_cols, label_col = "Curated_target", measure = "AUROC")

cols <- c(rep("lightgrey", length(keep_cols)))
names(cols) <- keep_cols
cols["Rank_integrated"] <- "black"


plot_perf(df = roc_df, auc_l = auroc, measure = "ROC", cols = cols, title = tf, ncol_legend = 1)
plot_perf(df = pr_df, auc_l = auprc, measure = "PR", cols = cols, title = tf, ncol_legend = 1)




###


# Inspect top negative cor
gene2 <- names(which.min(a[gene1,]))
sdat_sub <- subset(sdat, idents = names(cor_ct)[ct_ix])
plot_scatter(sdat_sub, gene1, gene2, slot = "data")
plot_scatter(sdat_sub, gene1, gene2, slot = "counts")
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "counts")[c(gene1, gene2), ])))
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "data")[c(gene1, gene2), ])))

