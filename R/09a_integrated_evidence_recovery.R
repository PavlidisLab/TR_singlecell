## TODO: point to missing genes
## TODO: missing genes excluded from rank df makes a difference?
## TODO: common process for curated and top k evidence for get_rank_df
## TODO: remove tf?
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

# For each dataset, load the subset gene x TF aggregation matrix 
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

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



# 
# ------------------------------------------------------------------------------
# TODO: Better way to assign species-specific argument since already have to 
# specifiy species?



avg_vs_ind_recover_curated_hg <- get_colwise_curated_auc_list(
  tfs = tfs_curated_hg,
  agg_l = agg_tf_hg,
  msr_mat = msr_hg,
  curated_df = curated,
  pc_df = pc_hg,
  species = "Human",
  ncores = 8,
  verbose = TRUE
)



avg_vs_ind_recover_curated_mm <- get_colwise_curated_auc_list(
  tfs = tfs_curated_mm,
  agg_l = agg_tf_mm,
  msr_mat = msr_mm,
  curated_df = curated,
  pc_df = pc_mm,
  species = "Mouse",
  ncores = 8,
  verbose = TRUE
)



saveRDS(avg_vs_ind_recover_curated_hg, avg_vs_ind_recover_curated_hg_path)
saveRDS(avg_vs_ind_recover_curated_mm, avg_vs_ind_recover_curated_mm_path)
