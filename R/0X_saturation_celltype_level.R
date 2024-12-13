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


# TODO: into function that does not collide with aggtools::prepare_celltype_mat
# (removing genes in advance speeds up)

ct_mat <- t(dat$Mat[, ct_ids])
keep_genes <- names(which(colSums(ct_mat > 0) >= 20))
keep_tfs <- intersect(tfs_hg$Symbol, keep_genes)
ct_mat <- ct_mat[, keep_genes]



# Full correlation matrix
cor_mat <- sparse_pcor(ct_mat)


n_steps <- 10
min_iter <- 10
steps <- dynamic_steps(max_cells = max_cells, min_cells = min_cells, max_steps = n_steps)
iters <- iterations_per_step(max_steps = n_steps, min_iter = min_iter)


# sample_cells <- function(cell_ids, steps, iters) {
#   
#   
# }



topk_at_iter <- lapply(1:iters[1], function(x) {
  
  
  sample_ids <- sample(ct_ids, steps[4], replace = FALSE)
  ct_mat_sub <- ct_mat[sample_ids, ]
  cor_mat_sub <- sparse_pcor(ct_mat_sub)
  
  # NA cors to 0 to allow overlap and remove self cor
  cor_mat_sub[is.na(cor_mat_sub)] <- 0  
  diag(cor_mat) <- diag(cor_mat_sub) <- 0 
  
  # Only consider TFs for overlap
  topk <- pair_colwise_topk(mat1 = cor_mat[, keep_tfs], 
                            mat2 = cor_mat_sub[, keep_tfs], 
                            k = k, 
                            ncores = ncore)
  
  gc(verbose = FALSE)
  
  return(topk)
  
})


topk_at_iter <- bind_rows(topk_at_iter)


generate_summ_df <- function(topk_l) {
  
  summ <- sapply(topk_l, summary)
  summ <- t(summ)
  summ <- as.data.frame(summ) %>% rownames_to_column(var = "Symbol")
  summ <- arrange(summ, Mean)
  summ$Symbol <- factor(summ$Symbol, levels = unique(summ$Symbol))
  
  return(summ)
}




summ_df <- generate_summ_df(topk_at_iter) 


ggplot(summ_df, aes(x = Symbol, y = Mean)) +
  geom_crossbar(aes(x = Symbol, ymin = `1st Qu.`, ymax = `3rd Qu.`)) +
  geom_point(aes(x = Symbol, y = Mean), shape = 21, colour = "firebrick") +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())





##

# sample_ids <- sample(ct_ids, steps[1], replace = FALSE)
# ct_mat_sub <- ct_mat[sample_ids, ]
# cor_mat_sub <- sparse_pcor(ct_mat_sub)
# 
# # NA cors to 0 to allow overlap and remove self cor
# cor_mat_sub[is.na(cor_mat_sub)] <- 0  
# diag(cor_mat) <- diag(cor_mat_sub) <- 0 
# 
# # Only consider TFs for overlap
# topk <- pair_colwise_topk(mat1 = cor_mat[, keep_tfs], 
#                           mat2 = cor_mat_sub[, keep_tfs], 
#                           k = k, 
#                           ncores = ncore)