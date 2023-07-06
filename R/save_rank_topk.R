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

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

# Loading genes of interest
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)


out_dir <- "~/scratch/R_objects/04-07-2023/"


test_genes_hg <- c(str_to_upper(tfs), "RPL3", "RPL19", "EGR1", "FOS", "FOSB")
test_genes_hg <- setdiff(test_genes_hg, "ASCL1")
test_genes_mm <- c(str_to_title(tfs), "Rpl3", "Rpl19", "Egr1", "Fos", "Fosb")



query_gene_rank_topk <- function(query_vec,
                                 subject_mat,
                                 query_gene,
                                 k = 1000,
                                 ncores = 1) {
  
  # genes <- rownames(subject_mat)
  genes <- colnames(subject_mat)
  
  # stopifnot(query_gene %in% genes, identical(names(query_vec), genes))
  
  query_topk <- topk_sort(query_vec, k)
  
  topk_l <- mclapply(genes, function(x) {
    subject_topk <- topk_sort(subject_mat[, x], k)
    topk_intersect(query_topk, subject_topk)
  }, mc.cores = ncores)
  
  names(topk_l) <- genes
  topk_sort <- sort(unlist(topk_l), decreasing = TRUE)
  rank_ix <- which(names(topk_sort) == query_gene)
  
  return(list(Rank = rank_ix, Topk = topk_sort))
}



# ids <- ids_hg
# i <- ids[1]
# j <- ids[2]
# query_gene <- "ASCL1"
# genes <- pc_hg$Symbol
# k <- 1000
# msr_mat <- msr_hg
# ncores <- 8


query_gene_rank_topk_all <- function(ids,
                                     query_gene,
                                     genes,
                                     k = 1000,
                                     msr_mat,
                                     ncores = 1) {
  
  rank_mat <- matrix(NA, nrow = length(ids), ncol = length(ids))
  rownames(rank_mat) <- colnames(rank_mat) <- ids
  diag(rank_mat) <- 1
  
  for (i in ids) {
    
    query_vec  <-  load_agg_mat_list(ids = i, sub_genes = query_gene, genes = genes)[[1]][, 1]
    message(paste(i, Sys.time()))
    
    for (j in ids) {
      
      # Skip pairing the same experiments or if query is not measured
      if (i == j || any(msr_mat[query_gene, c(i, j)] == 0)) {
        next
      } 
      
      subject_mat <- load_agg_mat_list(ids = j, genes = genes)[[1]]
      
      rank_mat[i, j] <- query_gene_rank_topk(query_vec = query_vec,
                                             subject_mat = subject_mat,
                                             query_gene = query_gene,
                                             k = k,
                                             ncores = ncores)$Rank
      
      rm(subject_mat)
      gc(verbose = FALSE)
    }
  }
  
  return(rank_mat)
}




for (gene in test_genes_hg) {
  
  message(paste("======== Beginning", gene, Sys.time()))
  
  rank_mat <- query_gene_rank_topk_all(ids = ids_hg, 
                                       query_gene = gene,
                                       genes = pc_hg$Symbol,
                                       k = k, 
                                       msr_mat = msr_hg,
                                       ncores = ncore)
  
  saveRDS(rank_mat, file = paste0(out_dir, gene, ".RDS"))
  
}



for (gene in test_genes_mm) {
  
  message(paste("======== Beginning", gene, Sys.time()))
  
  rank_mat <- query_gene_rank_topk_all(ids = ids_mm, 
                                       query_gene = gene,
                                       genes = pc_mm$Symbol,
                                       k = k, 
                                       msr_mat = msr_mm,
                                       ncores = ncore)
  
  saveRDS(rank_mat, file = paste0(out_dir, gene, ".RDS"))
  
}
