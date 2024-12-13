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
# msr_hg <- readRDS(msr_mat_hg_path)
# msr_mm <- readRDS(msr_mat_mm_path)
# 
# # Comeasure matrices tracking count of experiments gene pairs were comeasured
# comsr_hg <- readRDS(comsr_mat_hg_path)
# comsr_mm <- readRDS(comsr_mat_mm_path)

# Loading the subset TF and ribosomal aggregate matrices
# agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
# agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)



# Generate steps of the number of cells to subsample from the total number of cells.
# Steps are exponentially decided by: (prev_step / ratio) where ratio is assumed
# to be less than 1.0.

dynamic_steps <- function(max_cells, min_cells, max_steps, ratio = 0.5) {
  
  steps <- c(min_cells)
  
  while (tail(steps, 1) < max_cells) {
    next_step <- ceiling(tail(steps, 1) / ratio)
    
    if (next_step > max_cells || length(steps) >= max_steps) {
      break
    }
    
    steps <- c(steps, next_step)
  }
  
  steps <- unique(sort(c(steps, max_cells)))
  
  return(steps)
}



iterations_per_step <- function(max_steps, min_iter = 100) {
  
  iter <- rep(min_iter, max_steps)
  iter[max_steps] <- 1  # last step is full count of cells, no need to iter.
  
  return(iter)
  
}




id <- "GSE180928"
dat <- load_dat_list(id)[[1]]
n_cells <- count(dat$Meta, Cell_type, sort = TRUE)

ct <- "Ependymal"
ct_ids <- filter(dat$Meta, Cell_type == ct)$ID
max_cells <- filter(n_cells, Cell_type == ct)$n
min_cells <- 100


# Isolate cell type count matrix and prep for correlation
ct_mat <- prepare_celltype_mat(mat = dat$Mat, 
                              meta = dat$Meta, 
                              cell_type = ct,
                              pc_df = pc_hg)


# Full correlation matrix
cor_mat <- sparse_pcor(ct_mat)


n
steps <- dynamic_steps(max_cells = max_cells, min_cells = min_cells)
iters <- iter


