## 
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


# Convert cor matrix to column-wise (cell-type) ranks such that 1 is best

rank_cormat <- function(cmat) {
  
  rank_mat <- apply(-cmat, 2, rank, ties.method = "min", na.last = "keep")
  
  return(rank_mat)
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


##


cor_rank_ct <- mclapply(cor_ct, function(x) {
  rank_cormat(t(x))
}, mc.cores = ncore)



# https://stackoverflow.com/questions/42628385/sum-list-of-matrices-with-nas
# Sum ignores NAs (treated as 0), but converts sum==0 (all NAs) to NA
# Mean_nona divides by the number of non-NA observations
# Mean_all divides by the total number of cell types


cor_rank_sum <- apply(simplify2array(cor_rank_ct), 1:2, sum, na.rm = TRUE)
cor_rank_sum[cor_rank_sum == 0] <- NA

cor_rank_mean_nona <- apply(simplify2array(cor_rank_ct), 1:2, mean, na.rm = TRUE)
cor_rank_mean_all <- cor_rank_sum / length(cor_rank_ct)


assertthat::are_equal(
  sum(sapply(cor_rank_ct, function(x) x[1, 1]), na.rm = TRUE),
  cor_rank_sum[1, 1])


assertthat::are_equal(
  mean(sapply(cor_rank_ct, function(x) x["Xkr4", "Ascl1"]), na.rm = TRUE),
  cor_rank_mean_nona["Xkr4", "Ascl1"])


# Then convert these aggregations to ranks. Mean_all and sum same rank.

cor_agg_sum <- rank_cormat(-cor_rank_sum)
cor_agg_mean_nona <- rank_cormat(-cor_rank_mean_nona)
cor_agg_mean_all <- rank_cormat(-cor_rank_mean_all)
# identical(cor_agg_sum[, tf], cor_agg_mean_all[, tf])



# Focus on single TF and organize gene rankings


tf <- "Ascl1"


rank_df <- data.frame(
  Symbol = rownames(cor_agg_sum),
  Sum_raw = cor_rank_sum[, tf],
  Sum_rank = cor_agg_sum[, tf],
  Mean_raw_nona = cor_rank_mean_nona[, tf],
  Mean_rank_nona = cor_agg_mean_nona[, tf],
  Mean_raw_all = cor_rank_mean_all[, tf],
  Mean_rank_all = cor_agg_mean_all[, tf]
  ) %>%
  left_join(tf_na_count(cor_ct, tf), by = "Symbol")


cor(select_if(rank_df, is.numeric), use = "pairwise.complete.obs")



filter(rank_df, Count_NA < 120) %>% arrange(Mean_rank_nona) %>% head()


gene2 <- "Stk17b"
ct <- "VGLUT1-31-Igfn1_Ndst4"

sort(sapply(cor_ct, function(x) x[tf, gene2]), decreasing = TRUE)
head(sort(cor_ct[[ct]][tf, ], decreasing = TRUE))


# By cell type
sdat_sub <- subset(sdat, idents = ct)
plot_scatter(sdat_sub, tf, gene2, slot = "data", jitter = TRUE)
plot_scatter(sdat_sub, tf, gene2, slot = "data", jitter = FALSE)
plot_scatter(sdat_sub, tf, gene2, slot = "counts", jitter = TRUE)
plot_scatter(sdat_sub, tf, gene2, slot = "counts", jitter = FALSE)
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "data")[c(tf, gene2), ])))
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "counts")[c(tf, gene2), ])))


# Across all cells
plot_scatter(sdat, tf, gene2, slot = "data", jitter = TRUE)
plot_scatter(sdat, tf, gene2, slot = "data", jitter = FALSE)
plot_scatter(sdat, tf, gene2, slot = "counts", jitter = TRUE)
plot_scatter(sdat, tf, gene2, slot = "counts", jitter = FALSE)
cor_all[tf, gene2]
cor(t(as.matrix(GetAssayData(object = sdat, slot = "counts")[c(tf, gene2), ])))


# Which TF-gene pair had max cor in each cell type
max_cor <- max_cor_df(cor_ct, tf)



###

ct_ix <- 7
a <- cor_ct[[ct_ix]]
b <- cor_rank_ct[[ct_ix]]


gene1 <- "Ascl1"
not_na <- which(!is.na(a[gene1, ]))

# Ensure max cor equals min (best) rank
assertthat::assert_that(
  names(which.max(a[gene1,])) == names(which.min(b[, gene1])))

# NAs ranked last, want min non-NA ranking
assertthat::assert_that(
  names(which.min(a[gene1,])) == names(which.max(b[not_na, gene1])))


summary(a[gene1, ])
a[1:5, 1:5]
b[1:5, 1:5]
head(sort(a[gene1, ], decreasing = TRUE))
head(sort(a[gene1, ]))
head(sort(b[not_na, gene1], decreasing = TRUE))
head(sort(b[, gene1]))



# Inspect top negative cor
gene2 <- names(which.min(a[gene1,]))
sdat_sub <- subset(sdat, idents = names(cor_ct)[ct_ix])
plot_scatter(sdat_sub, gene1, gene2, slot = "data")
plot_scatter(sdat_sub, gene1, gene2, slot = "counts")
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "counts")[c(gene1, gene2), ])))
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "data")[c(gene1, gene2), ])))

