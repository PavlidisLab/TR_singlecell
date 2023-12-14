## Save out a list of the AUC performances of aggregated coexpression's ability
## to recover curated targets, reversing the lists to prioritize negative coexpr
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

n_samps <- 1000

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

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Saved list RDS of the summaries
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)


# Flip sign of Average_RSR (used as importance score) to make the bottom of
# the rankings (most consistent negative coexpr) more important
# ------------------------------------------------------------------------------


flip_avgrsr <- function(rank_df) {
  rank_df$Avg_RSR <- -rank_df$Avg_RSR
  return(rank_df)
}



rank_tf_hg <- lapply(rank_tf_hg, flip_avgrsr)
rank_tf_mm <- lapply(rank_tf_mm, flip_avgrsr)



# Get ortho-matched symbols of TFs with available data, as well as all targets
# which are used for null
# ------------------------------------------------------------------------------


# Human

ortho_tf_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_hg %in% names(rank_tf_hg)) %>% 
  pull(Symbol_hg)

ortho_target_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_hg %in% rownames(rank_tf_hg[[1]])) %>% 
  pull(Symbol_hg)

tf_hg <- union(
  intersect(names(rank_tf_hg), str_to_upper(curated$TF_Symbol)),
  ortho_tf_hg)

target_hg <- union(
  intersect(rownames(rank_tf_hg[[1]]), str_to_upper(curated$Target_Symbol)),
  ortho_target_hg)


# Mouse

ortho_tf_mm <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_mm %in% names(rank_tf_mm)) %>% 
  pull(Symbol_mm)

ortho_target_mm <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_mm %in% rownames(rank_tf_mm[[1]])) %>% 
  pull(Symbol_mm)

tf_mm <- union(
  intersect(names(rank_tf_mm), str_to_title(curated$TF_Symbol)),
  ortho_tf_mm)

target_mm <- union(
  intersect(rownames(rank_tf_mm[[1]]), str_to_title(curated$Target_Symbol)),
  ortho_target_mm)



# Run and save
# ------------------------------------------------------------------------------


set.seed(5)


# Human 

save_curated_auc_list(path = rev_coexpr_auc_hg_path,
                      tfs = tf_hg,
                      rank_l = rank_tf_hg,
                      score_col = "Avg_RSR",
                      curated_df = curated,
                      label_all = target_hg,
                      ortho_df = pc_ortho,
                      pc_df = pc_hg,
                      species = "Human",
                      n_samps = n_samps,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)


# Mouse 

save_curated_auc_list(path = rev_coexpr_auc_mm_path,
                      tfs = tf_mm,
                      rank_l = rank_tf_mm,
                      score_col = "Avg_RSR",
                      curated_df = curated,
                      label_all = target_mm,
                      ortho_df = pc_ortho,
                      pc_df = pc_mm,
                      species = "Mouse",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)
