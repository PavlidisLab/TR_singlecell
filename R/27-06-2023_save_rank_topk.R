library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 1000 

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Loading genes of interest
# sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
# lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
# pc_ortho <- read.delim(pc_ortho_path)
# ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

# Genes to focus/subset on when loading aggregate coexpression matrices
subset_hg <- NULL
subset_mm <- NULL

# Load aggregate matrix into list
# agg_hg <- load_agg_mat_list(ids = ids_hg, sub_genes = subset_hg)
# agg_mm <- load_agg_mat_list(ids = ids_mm, sub_genes = subset_mm)

# genes_hg <- rownames(agg_hg[[1]])
# genes_mm <- rownames(agg_mm[[1]])

# TODO: this needs to be replaced with upstream ordering
# agg_hg <- lapply(agg_hg, function(x) x[genes_hg, ])
# agg_mm <- lapply(agg_mm, function(x) x[genes_mm, ])


# stopifnot(all(unlist(lapply(agg_hg, function(x) identical(rownames(x), genes_hg)))))
# stopifnot(all(unlist(lapply(agg_mm, function(x) identical(rownames(x), genes_mm)))))


out_dir <- "~/scratch/R_objects/27-06-2023/"

test_genes_hg <- c(str_to_upper(tfs), "RPL3", "RPL19", "EGR1", "FOS", "FOSB")
test_genes_mm <- c(str_to_title(tfs), "Rpl3", "Rpl19", "Egr1", "Fos", "Fosb")



query_gene_rank_topk <- function(query_vec,
                                 subject_mat,
                                 gene,
                                 k = 1000,
                                 ncores = 1) {
  
  # genes <- rownames(subject_mat)
  genes <- colnames(subject_mat)
  
  # stopifnot(gene %in% genes, identical(names(query_vec), genes))
  
  query_topk <- topk_sort(query_vec, k)
  
  topk_l <- mclapply(genes, function(x) {
    subject_topk <- topk_sort(subject_mat[, x], k)
    topk_intersect(query_topk, subject_topk)
  }, mc.cores = ncores)
  
  names(topk_l) <- genes
  topk_sort <- sort(unlist(topk_l), decreasing = TRUE)
  rank_ix <- which(names(topk_sort) == gene)
  
  return(list(Rank = rank_ix, Topk = topk_sort))
  # return(rank_ix)
}



query_gene_rank_topk_all <- function(ids,
                                     gene,
                                     k = 1000,
                                     ncores = 1) {
  
  rank_mat <- matrix(1, nrow = length(ids), ncol = length(ids))
  rownames(rank_mat) <- colnames(rank_mat) <- ids
  
  for (i in ids) {
    
    query_vec  <-  load_agg_mat_list(ids = i, sub_genes = gene)[[1]][, 1]
    message(paste(i, Sys.time()))
    
    for (j in ids) {
      if (i == j) next
      
      subject_mat <- load_agg_mat_list(ids = j)[[1]]
      
      rank_mat[i, j] <- query_gene_rank_topk(query_vec,
                                             subject_mat,
                                             gene,
                                             k,
                                             ncores)$Rank
      
      rm(subject_mat)
      gc(verbose = FALSE)
    }
  }
  
  return(rank_mat)
}




for (gene in test_genes_hg) {
  
  message(paste("Beginning", gene, Sys.time()))
  
  rank_mat <- query_gene_rank_topk_all(ids_hg, gene, k, ncore)
  
  saveRDS(rank_mat, file = paste0(out_dir, gene, ".RDS"))
  
}



for (gene in test_genes_mm) {
  
  message(paste("Beginning", gene, Sys.time()))
  
  rank_mat <- query_gene_rank_topk_all(ids_hg, gene, k, ncore)
  
  saveRDS(rank_mat, file = paste0(out_dir, gene, ".RDS"))
  
}