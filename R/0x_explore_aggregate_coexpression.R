##
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(googlesheets4)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# 
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

#
pc_ortho <- read.delim(pc_ortho_path)

# Ranked targets from paper
evidence_l <- readRDS(evidence_path)
tfs_mm <- names(evidence_l$Mouse)
tfs_hg <- names(evidence_l$Human)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID
ids <- union(ids_hg, ids_mm)

# Load aggregate matrix into list
# agg_files <- list.files(amat_dir, pattern = "_RSR_allrank.RDS", recursive = TRUE, full.names = TRUE)
# agg_l <- lapply(agg_files, readRDS)
agg_hg <- lapply(ids_hg, function(x) readRDS(file.path(amat_dir, x, paste0(x, "_RSR_allrank.RDS"))))
names(agg_hg) <- ids_hg
agg_hg <- lapply(agg_hg, lowertri_to_symm)
agg_mm <- lapply(ids_mm, function(x) readRDS(file.path(amat_dir, x, paste0(x, "_RSR_allrank.RDS"))))
names(agg_mm) <- ids_mm
agg_mm <- lapply(agg_mm, lowertri_to_symm)

gc()

# Load mat and meta
# dat_files <- list.files(amat_dir, pattern = "_clean_mat_and_meta.RDS", recursive = TRUE, full.names = TRUE)
# dat_l <- lapply(dat_files, readRDS)
# dat_mm <- lapply(ids_mm, function(x) readRDS(file.path(amat_dir, x, paste0(x, "_clean_mat_and_meta.RDS")))
# names(dat_mm) <- ids_mm
# dat_hg <- lapply(ids_hg, function(x) readRDS(file.path(amat_dir, x, paste0(x, "_clean_mat_and_meta.RDS")))
# names(dat_hg) <- ids_hg


genes_hg <- rownames(agg_hg[[1]])
genes_mm <- rownames(agg_mm[[1]])



gene_vec_to_mat <- function(agg_l, gene) {
  do.call(cbind, lapply(agg_l, function(x) x[gene, ]))
}


gene <- "RPL3"




gene_mat <- gene_vec_to_mat(agg_hg, gene)


# Correlation of gene's rankings across studies
gene_rank_cor <- cor(gene_mat, method = "spearman", use = "pairwise.complete.obs")


# For the given gene, retrieve the rank position of its similarity for the 
# other studies


rank_similarity_matrix <- function(agg_l, gene, ncores = 1) {
  
  ids <- names(agg_l)
  genes <- rownames(agg_l[[1]])
  
  stopifnot(length(ids) > 0, gene %in% genes)
  
  rank_cor_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  rownames(rank_cor_mat) <- colnames(rank_cor_mat) <- ids
  
  for (i in ids) {
    for (j in ids) {
      if (i == j) {
        next
      }
      rank_cor <-  WGCNA::cor(x = agg_l[[i]][, gene], y = agg_l[[j]], nThreads = ncores) 
      rank_cor <- sort(rank_cor[1, ], decreasing = TRUE)
      rank_cor_mat[i, j] <- which(names(rank_cor) == gene)
    }
  }
  
  return(rank_cor_mat)
}



rmat_hg <- rank_similarity_matrix(agg_hg, gene, ncores = ncore)
# rmat_mm <- rank_similarity_matrix(agg_mm, gene, ncores = ncore)




# Using ribosomal genes as a sanity check
lapply(agg_hg, function(x) head(sort(x[, "RPL3"], decreasing = TRUE)))


# List of max pairs for each data set

max_pair_l <- function(agg_l, ncores = 1) {
  mclapply(agg_l, function(x) {
    diag(x) <- NA
    which(x == max(x, na.rm = TRUE), arr.ind = TRUE)
  }, mc.cores = ncores)
}


max_pair_hg <- max_pair_l(agg_hg, ncore)
max_pair_mm <- max_pair_l(agg_mm, ncore)


dat_files <- list.files(amat_dir, pattern = "_clean_mat_and_meta.RDS", recursive = TRUE, full.names = TRUE)
# dat_l <- lapply(dat_files, readRDS)
ct_cors <- all_celltype_cor(mat, meta, gene1 = gene1, gene2 = gene2)


