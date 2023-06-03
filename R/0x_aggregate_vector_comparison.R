## This script generates metrics of similarity between a query gene vector from 
## one aggregated correlation matrix and all other gene vectors from another
## network. The main idea idea is to find the rank position of the gene of
## interest with itself in another network.
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID
ids <- union(ids_hg, ids_mm)

# Load aggregate matrix into list
agg_hg <- lapply(ids_hg, function(x) lowertri_to_symm(readRDS(file.path(amat_dir, x, paste0(x, "_RSR_allrank.RDS")))))
names(agg_hg) <- ids_hg
agg_mm <- lapply(ids_mm, function(x) lowertri_to_symm(readRDS(file.path(amat_dir, x, paste0(x, "_RSR_allrank.RDS")))))
names(agg_mm) <- ids_mm
gc()

#
sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
pc_ortho <- read.delim(pc_ortho_path)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))

#
genes_hg <- rownames(agg_hg[[1]])
genes_mm <- rownames(agg_mm[[1]])




gene <- "RPL13"  # "RPL13"


gene_mat <- gene_vec_to_mat(agg_hg, gene)


# Correlation of gene's rankings across studies

gene_rank_cor <- cor(gene_mat, method = "spearman", use = "pairwise.complete.obs")





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