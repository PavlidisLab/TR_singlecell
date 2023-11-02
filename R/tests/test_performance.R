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

# Curated TFs with ChIP-seq and all targets for null
tfs_curated_hg <- intersect(tfs_hg$Symbol, str_to_upper(curated$TF_Symbol))
tfs_curated_mm <- intersect(tfs_mm$Symbol, str_to_title(curated$TF_Symbol))
targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))


# coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)



#
# ------------------------------------------------------------------------------



tf <- "Pax6"
rank_l <- rank_tf_mm
rank_df <- rank_l[[tf]]
score_col <- "Avg_RSR"
curated_df <- curated
label_all <- targets_curated_mm
pc_df <- pc_mm
species <- "Mouse"
n_samps <- n_samps
ncores <- 8


# Reversing ranks
# rank_df$Avg_RSR <- -rank_df$Avg_RSR


set.seed(5)

tx <- curated_obs_and_null_auc(
  tf = tf,
  rank_df = rank_df,
  score_col = score_col,
  curated_df = curated,
  label_all = label_all,
  pc_df = pc_df,
  species = species,
  n_samps = n_samps,
  ncores = ncores
)



#
# ------------------------------------------------------------------------------


#

tfs <- tfs_curated_hg[1:3]
rank_l <- rank_tf_hg
score_col <- "Avg_RSR"
curated_df <- curated
label_all <- targets_curated_hg
pc_df <- pc_hg
species <- "Human"
n_samps <- n_samps
ncores <- 8



tx <- curated_obs_and_null_auc_list(
  tfs = tfs,
  rank_l = rank_l,
  score_col = score_col,
  curated_df = curated,
  label_all = label_all,
  pc_df = pc_df,
  species = species,
  n_samps = n_samps,
  ncores = ncores,
  verbose = TRUE
)
