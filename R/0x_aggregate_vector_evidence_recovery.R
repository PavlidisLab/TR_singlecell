## TODO
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

# Ranked targets from paper
evidence_l <- readRDS(evidence_path)
tfs_mm <- names(evidence_l$Mouse)
tfs_hg <- names(evidence_l$Human)

#
genes_hg <- rownames(agg_hg[[1]])
genes_mm <- rownames(agg_mm[[1]])



# Functions
# ------------------------------------------------------------------------------



which_ranked_auprc <- function(rank_vec, 
                               rank_mat, 
                               gene, 
                               k = 1000,
                               ncores = 1) {
  
  genes <- rownames(rank_mat)
  x <- sort(rank_vec, decreasing = TRUE)
  x <- x[names(x) != gene]
  
  auprc_l <- mclapply(genes, function(g) {
    
    y <- sort(rank_mat[, g], decreasing = TRUE)[1:(k+1)]
    y <- y[names(y) != g]
    
    rank_df <- data.frame(
      Symbol = names(x),
      Rank = x,
      Label = names(x) %in% names(y)
    )
    
    if (all(rank_df$Label)) return(1)
    if (all(!rank_df$Label)) return(0)
    
    get_au_perf(rank_df, label_col = "Label", score_col = "Rank", measure = "AUPRC")
    
  }, mc.cores = ncores)
  names(auprc_l) <- genes
  
  rank_auprc <- sort(unlist(auprc_l), decreasing = TRUE)
  rank_ix <- which(names(rank_auprc) == gene)
  
  return(list(Rank = rank_ix, List = rank_auprc))
}



rank_similarity_matrix <- function(agg_l, 
                                   gene, 
                                   msr,  # AUPRC|Intersect|Cor
                                   k = 1000,
                                   ncores = 1) {
  
  ids <- names(agg_l)
  genes <- rownames(agg_l[[1]])
  
  stopifnot(length(ids) > 0, gene %in% genes)
  
  rank_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  rownames(rank_mat) <- colnames(rank_mat) <- ids
  
  for (i in ids) {
    for (j in ids) {
      if (i == j) {
        next
      }
      
      rank_ix <- if (msr == "AUPRC") {
        
        which_ranked_auprc(
          rank_vec = agg_l[[i]][, gene],
          rank_mat = agg_l[[j]],
          gene = gene,
          k = k,
          ncores = ncores)$Rank
        
      } else if (msr == "Intersect") {
        
        which_ranked_intersect(
          rank_vec = agg_l[[i]][, gene],
          rank_mat = agg_l[[j]],
          gene = gene,
          k = k,
          ncores = ncores)$Rank
        
      } else {
        
        which_ranked_cor(
          rank_vec = agg_l[[i]][, gene],
          rank_mat = agg_l[[j]],
          gene = gene,
          ncores = ncores)$Rank
      }
      
      rank_mat[i, j] <- rank_ix
    }
  }
  
  return(rank_mat)
}



# 
# ------------------------------------------------------------------------------



