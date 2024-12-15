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
n_tf_msr <- msr_hg[tfs_hg$Symbol]


tf_mat <- agg_tf_hg$E2F8


n_iter <- 100



