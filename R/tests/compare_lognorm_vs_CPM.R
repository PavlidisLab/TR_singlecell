## Examine differences between results generated from log norm vs CPM norm
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(Seurat)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 1000

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



# 2. 
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



# Most diff

diff_rank_df <- left_join(rank_cpm_hg$YBX1, 
                          rank_ln_hg$YBX1,
                          by = "Symbol",
                          suffix = c("_CPM", "_LN")) %>% 
  mutate(Diff = Rank_RSR_CPM - Rank_RSR_LN)


# Cor of aggregate ranking


get_rank_cor_df <- function(cpm_l, lognorm_l, ncore = 1) {
  
  symbols <- intersect(names(cpm_l), names(lognorm_l))
  
  cor_l <- mclapply(symbols, function(x) {
    cor(cpm_l[[x]]$Avg_RSR, lognorm_l[[x]]$Avg_RSR, 
        use = "pairwise.complete.obs", method = "spearman")
  }, mc.cores = ncore)
  
  cor_df <- data.frame(Symbol = symbols, Cor = unlist(cor_l))
  
  return(cor_df)
}



rank_cor_hg <- get_rank_cor_df(rank_cpm_hg, rank_ln_hg, ncore)
rank_cor_mm <- get_rank_cor_df(rank_cpm_mm, rank_ln_mm, ncore)

plot_hist(rank_cor_hg, stat_col = "Cor")
plot_hist(rank_cor_mm, stat_col = "Cor")

