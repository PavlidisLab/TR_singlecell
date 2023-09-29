## Examine differences between results generated from log norm vs CPM norm
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(Seurat)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

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
sim_cpm_hg <- readRDS(sim_tf_hg_path)
sim_cpm_mm <- readRDS(sim_tf_mm_path)
sim_ln_hg <- readRDS("/space/scratch/amorin/R_objects/similarity_TF_hg_lognorm.RDS")
sim_ln_mm <- readRDS("/space/scratch/amorin/R_objects/similarity_TF_mm_lognorm.RDS")


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


sim_df_hg <- left_join(summ_cpm_hg$Topk, summ_ln_hg$Topk, by = "Symbol", suffix = c("_CPM", "_LN"))
