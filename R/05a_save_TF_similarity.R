## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Load aggregate matrix into list
agg_hg <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol[1:3])
agg_mm <- load_agg_mat_list(ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol[1:3])

stopifnot(all(unlist(lapply(agg_hg, function(x) identical(rownames(x), pc_hg$Symbol)))))
stopifnot(all(unlist(lapply(agg_mm, function(x) identical(rownames(x), pc_mm$Symbol)))))

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)



#
# ------------------------------------------------------------------------------


if (!file.exists(tf_sim_hg_path)) {
  sim_hg <- get_all_similarity(agg_l = agg_hg, msr_mat = msr_hg, genes = tfs_hg$Symbol, k = 1000)
  sim_hg <- sim_hg[!is.na(sim_hg)]
  saveRDS(sim_hg, tf_sim_hg_path)
}



if (!file.exists(tf_sim_mm_path)) {
  sim_mm <- get_all_similarity(agg_l = agg_mm, msr_mat = msr_mm, genes = tfs_mm$Symbol, k = 1000)
  sim_mm <- sim_mm[!is.na(sim_mm)]
  saveRDS(sim_mm, tf_sim_mm_path)
}
