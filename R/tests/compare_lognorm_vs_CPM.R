## Examine differences between results generated from log norm vs CPM norm
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(Seurat)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 200

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# Measurement matrices used for finding when a gene was measured in an experiment
msr_cpm_hg <- readRDS(msr_mat_hg_path)
msr_cpm_mm <- readRDS(msr_mat_mm_path)
msr_ln_hg <- readRDS("/space/scratch/amorin/R_objects/binary_measurement_matrix_hg_lognorm.RDS")
msr_ln_mm <- readRDS("/space/scratch/amorin/R_objects/binary_measurement_matrix_mm_lognorm.RDS")

stopifnot(identical(rownames(msr_cpm_hg), rownames(msr_ln_hg)),
          identical(colnames(msr_cpm_hg), colnames(msr_ln_hg)))

stopifnot(identical(rownames(msr_cpm_mm), rownames(msr_ln_mm)),
          identical(colnames(msr_cpm_mm), colnames(msr_ln_mm)))

# TF Similarity objects
sim_cpm_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_TF_hg_k=", k, ".RDS"))
sim_cpm_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_TF_mm_k=", k, ".RDS"))
sim_ln_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_TF_hg_k=", k, "_lognorm.RDS"))
sim_ln_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_TF_mm_k=", k, "_lognorm.RDS"))

# TF rankings
rank_cpm_hg <- readRDS(rank_tf_hg_path)
rank_cpm_mm <- readRDS(rank_tf_mm_path)
rank_ln_hg <- readRDS("/space/scratch/amorin/R_objects/TF_agg_ranking_hg_lognorm.RDS")
rank_ln_mm <- readRDS("/space/scratch/amorin/R_objects/TF_agg_ranking_mm_lognorm.RDS")

stopifnot(identical(names(rank_cpm_hg), names(rank_ln_hg)))
stopifnot(identical(names(rank_cpm_mm), names(rank_ln_mm)))
stopifnot(identical(rank_cpm_hg[[1]]$Symbol, rank_ln_hg[[1]]$Symbol))
stopifnot(identical(rank_cpm_mm[[1]]$Symbol, rank_ln_mm[[1]]$Symbol))

# Focus on human aggregate lists due to size
agg_cpm_hg <- readRDS("/space/scratch/amorin/R_objects/TF_agg_mat_list_hg.RDS")
agg_ln_hg <- readRDS("/space/scratch/amorin/R_objects/TF_agg_mat_list_hg_lognorm.RDS")


# 1. Differences in binary measurement
# ------------------------------------------------------------------------------


get_msr_diff <- function(msr_cpm, msr_ln) {
  
  msr_diff <- msr_cpm - msr_ln 
  msr_diff_df <- data.frame(which(msr_diff != 0, arr.ind = TRUE))
  msr_diff_df$row <- rownames(msr_cpm)[msr_diff_df$row]
  msr_diff_df$col <- colnames(msr_cpm)[msr_diff_df$col]
  rownames(msr_diff_df) <- NULL
  colnames(msr_diff_df) <- c("Gene", "Experiment")
  msr_diff_df$Diff <- msr_diff[msr_diff != 0]
  
  return(msr_diff_df)
}



msr_diff_hg <- get_msr_diff(msr_cpm_hg, msr_ln_hg)
msr_diff_mm <- get_msr_diff(msr_cpm_mm, msr_ln_mm)



# 2. Examining TF similarity
# ------------------------------------------------------------------------------


# Returns a list of dataframes of summary stats for each TF's similarity

get_summary_df <- function(sim_l, msr_mat) {
  
  stats <- c("Scor", "Topk", "Bottomk", "Jaccard")
  
  df_l <- lapply(stats, function(stat) {
    
    lapply(sim_l, function(x) summary(x[[stat]])) %>% 
      do.call(rbind, .) %>%
      as.data.frame() %>% 
      rownames_to_column(var = "Symbol") %>% 
      mutate(N_exp = rowSums(msr_mat[names(sim_l), ])) %>%
      arrange(desc(Median)) %>% 
      mutate(Symbol = factor(Symbol, levels = unique(Symbol)))
  })
  
  names(df_l) <- stats
  return(df_l)
}



summ_cpm_hg <- get_summary_df(sim_cpm_hg, msr_cpm_hg)
summ_cpm_mm <- get_summary_df(sim_cpm_mm, msr_cpm_mm)
summ_ln_hg <- get_summary_df(sim_ln_hg, msr_ln_hg)
summ_ln_mm <- get_summary_df(sim_ln_mm, msr_ln_mm)


sim_df_hg <- left_join(summ_cpm_hg$Topk,
                       summ_ln_hg$Topk,
                       by = "Symbol",
                       suffix = c("_CPM", "_LN")) %>%
  mutate(Diff = Mean_CPM - Mean_LN)


sim_df_mm <- left_join(summ_cpm_mm$Topk,
                       summ_ln_mm$Topk,
                       by = "Symbol",
                       suffix = c("_CPM", "_LN")) %>%
  mutate(Diff = Mean_CPM - Mean_LN)


qplot(df = sim_df_hg, yvar = "Mean_CPM", xvar = "Mean_LN") + 
  xlab(paste0("Mean LN topk=", k)) + 
  ylab(paste0("Mean CPM topk=", k)) +
  ggtitle("Human")

qplot(df = sim_df_mm, yvar = "Mean_CPM", xvar = "Mean_LN") + 
  xlab(paste0("Mean LN topk=", k)) + 
  ylab(paste0("Mean CPM topk=", k)) +
  ggtitle("Mouse")

cor(sim_df_hg$Mean_CPM, sim_df_hg$Mean_LN, method = "spearman")
cor(sim_df_mm$Mean_CPM, sim_df_mm$Mean_LN, method = "spearman")


plot_hist(sim_df_hg, stat_col = "Diff")
plot_hist(sim_df_mm, stat_col = "Diff")



# 3) Examining aggregate rankings
# ------------------------------------------------------------------------------


# Cor of aggregate ranking

get_rank_cor_df <- function(cpm_l, lognorm_l, ncores = 1) {
  
  symbols <- intersect(names(cpm_l), names(lognorm_l))
  
  cor_l <- mclapply(symbols, function(x) {
    cor(cpm_l[[x]]$Avg_RSR, lognorm_l[[x]]$Avg_RSR, 
        use = "pairwise.complete.obs", method = "spearman")
  }, mc.cores = ncores)
  
  cor_df <- data.frame(Symbol = symbols, Cor = unlist(cor_l))
  
  return(cor_df)
}



rank_cor_hg <- get_rank_cor_df(rank_cpm_hg, rank_ln_hg, ncore)
rank_cor_mm <- get_rank_cor_df(rank_cpm_mm, rank_ln_mm, ncore)

plot_hist(rank_cor_hg, stat_col = "Cor")
plot_hist(rank_cor_mm, stat_col = "Cor")


# Top/bottom k between rankings

get_topk_df <- function(cpm_l, lognorm_l, ncores = 1) {
  
  symbols <- intersect(names(cpm_l), names(lognorm_l))
  
  topk_l <- mclapply(symbols, function(x) {
    
    vec1 <- setNames(cpm_l[[x]]$Avg_RSR, cpm_l[[x]]$Symbol)
    vec2 <- setNames(lognorm_l[[x]]$Avg_RSR, lognorm_l[[x]]$Symbol)
    
    top_k200 <- topk_intersect(topk_sort(vec1, k = 200), 
                               topk_sort(vec2, k = 200))
    
    top_k1000 <- topk_intersect(topk_sort(vec1, k = 1000), 
                                topk_sort(vec2, k = 1000))
    
    bottom_k200 <- topk_intersect(topk_sort(vec1, k = 200, decreasing = FALSE), 
                                  topk_sort(vec2, k = 200, decreasing = FALSE))
    
    bottom_k1000 <- topk_intersect(topk_sort(vec1, k = 1000, decreasing = FALSE),
                                   topk_sort(vec2, k = 1000, decreasing = FALSE))
    
    data.frame(
      Symbol = x,
      Top_K200 = top_k200,
      Top_K1000 = top_k1000,
      Bottom_K200 = bottom_k200,
      Bottom_K1000 = bottom_k1000)

  }, mc.cores = ncores)
  
  topk_df <- do.call(rbind, topk_l)
  return(topk_df)
}


topk_df_hg <- get_topk_df(rank_cpm_hg, rank_ln_hg, ncores = ncore)
topk_df_mm <- get_topk_df(rank_cpm_mm, rank_ln_mm, ncores = ncore)



qplot(topk_df_hg, xvar = "Top_K200", yvar = "Top_K1000")
qplot(topk_df_mm, xvar = "Top_K200", yvar = "Top_K1000")
qplot(topk_df_hg, xvar = "Top_K200", yvar = "Top_K1000")
qplot(topk_df_mm, xvar = "Top_K200", yvar = "Top_K1000")

qplot(topk_df_hg, xvar = "Top_K200", yvar = "Bottom_K200")
qplot(topk_df_mm, xvar = "Top_K200", yvar = "Bottom_K200")
qplot(topk_df_hg, xvar = "Top_K1000", yvar = "Bottom_K1000")
qplot(topk_df_mm, xvar = "Top_K1000", yvar = "Bottom_K1000")

plot_hist(topk_df_hg, stat_col = "Top_K200")
plot_hist(topk_df_mm, stat_col = "Top_K200")
plot_hist(topk_df_hg, stat_col = "Top_K1000")
plot_hist(topk_df_mm, stat_col = "Top_K1000")

plot_hist(topk_df_hg, stat_col = "Bottom_K200")
plot_hist(topk_df_mm, stat_col = "Bottom_K200")
plot_hist(topk_df_hg, stat_col = "Bottom_K1000")
plot_hist(topk_df_mm, stat_col = "Bottom_K1000")



# Examine difference in ranks for a given TF

gene <- "PAX6"

diff_rank_df <- left_join(rank_cpm_hg[[gene]], 
                          rank_ln_hg[[gene]],
                          by = "Symbol",
                          suffix = c("_CPM", "_LN")) %>% 
  mutate(Diff = Rank_RSR_CPM - Rank_RSR_LN)



qplot(diff_rank_df, xvar = "Avg_RSR_CPM", yvar = "Avg_RSR_LN")


# Explore at level of aggregate matrices

mat_cpm <- subset_to_measured(gene_vec_to_mat(agg_cpm_hg, gene), msr_cpm_hg, gene)
mat_ln <- subset_to_measured(gene_vec_to_mat(agg_ln_hg, gene), msr_ln_hg, gene)

comsr <- lapply(rownames(mat_cpm), function(y) {
  mat_cpm[gene, colnames(mat_cpm)] & msr_cpm_hg[y, colnames(mat_cpm)]
})
comsr <- do.call(rbind, comsr)
rownames(comsr) <- rownames(mat_cpm)

# When a TF-gene was not co-measured, impute to the median NA
med_cpm <- median(mat_cpm, na.rm = TRUE)
med_ln <- median(mat_ln, na.rm = TRUE)

mat_cpm[comsr == FALSE] <- med_cpm
mat_ln[comsr == FALSE] <- med_ln

arrange(rank_cpm_hg[[gene]], Rank_RSR) %>% head
arrange(rank_ln_hg[[gene]], Rank_RSR) %>% head

head(sort(rowMeans(mat_cpm), decreasing = TRUE), 10)
head(sort(rowMeans(mat_ln), decreasing = TRUE), 10)
