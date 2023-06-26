## Testing speed of implementations for functions for comparing vector similarity

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
source("R/utils/vector_comparison_functions.R")
source("R/00_config.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID[1:3]
ids_mm <- filter(sc_meta, Species == "Mouse")$ID[1:3]

# Genes to focus/subset on when loading aggregate coexpression matrices
subset_hg <- NULL    # c("RPL3", "EGR1")
subset_mm <- NULL    # c("Rpl3", "Egr1")

# Load aggregate matrix into list
agg_hg <- load_agg_mat_list(ids = ids_hg, sub_genes = subset_hg)
agg_mm <- load_agg_mat_list(ids = ids_mm, sub_genes = subset_mm)

genes_hg <- rownames(agg_hg[[1]])
genes_mm <- rownames(agg_mm[[1]])
stopifnot(all(unlist(lapply(agg_hg, function(x) identical(rownames(x), genes_hg)))))
stopifnot(all(unlist(lapply(agg_mm, function(x) identical(rownames(x), genes_mm)))))


gene <- "EGR1"



# The outer() version is slower than the nested loop, as can't skip diagonal

colwise_topk_auprc <- function(mat, k = 1000) {
  
  
  col_list <- lapply(1:ncol(mat), function(x) {
    
    list(Scores = sort(mat[, x], decreasing = TRUE),
         Labels = topk_sort(mat[, x], k))
    
  })
  
  auprc_mat <- outer(lapply(col_list, `[[`, "Scores"), 
                     lapply(col_list, `[[`, "Labels"), 
                     Vectorize(vec_auprc))
  
  colnames(auprc_mat) <- rownames(auprc_mat) <- colnames(mat)
  
  return(auprc_mat)
}



# Nested loop

colwise_topk_auprc2 <- function(mat, k = 1000) {
  
  auprc_mat <- matrix(1, nrow = ncol(mat), ncol = ncol(mat))
  colnames(auprc_mat) <- rownames(auprc_mat) <- colnames(mat)
  
  for (i in 1:nrow(auprc_mat)) {
    for (j in 1:ncol(auprc_mat)) {
      if (i == j) next
      scores <- sort(mat[, i], decreasing = TRUE)
      labels <- topk_sort(mat[, j], k)
      auprc_mat[i, j] <- vec_auprc(scores, labels)
    }
  }
  
  return(auprc_mat)
}




# Directly access list instead of binding vectors

colwise_topk_auprc3 <- function(agg_l, gene, k = 1000) {
  
  auprc_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  colnames(auprc_mat) <- rownames(auprc_mat) <- names(agg_l)
  
  for (i in 1:nrow(auprc_mat)) {
    for (j in 1:ncol(auprc_mat)) {
      if (i == j) next
      scores <- sort(agg_l[[i]][, gene], decreasing = TRUE)
      labels <- topk_sort(agg_l[[j]][, gene], k)
      auprc_mat[i, j] <- vec_auprc(scores, labels)
    }
  }
  
  return(auprc_mat)
}



res1 <- microbenchmark::microbenchmark(
  F1 = gene_vec_to_mat(agg_hg, gene) %>% colwise_topk_auprc(),
  F2 = gene_vec_to_mat(agg_hg, gene) %>% colwise_topk_auprc2(),
  F3 = colwise_topk_auprc3(agg_hg, gene),
  times = 5)





# Nested loop slower than one pass over columns into list for input to outer

colwise_topk_intersect2 <- function(mat, k = 1000) {

  ids <- colnames(mat)
  topk_mat <- matrix(NA, nrow = ncol(mat), ncol = ncol(mat))
  diag(topk_mat) <- k
  rownames(topk_mat) <- colnames(topk_mat) <- ids

  for (i in ids) {
    vec1 <- topk_sort(mat[, i], k)
    for (j in ids) {
      if (i == j) next
      vec2 <- topk_sort(mat[, j], k)
      topk_mat[i, j] <- topk_intersect(vec1, vec2)
    }
  }

  return(topk_mat)
}



res2 <- microbenchmark::microbenchmark(
  F1 = gene_vec_to_mat(agg_hg, gene) %>% colwise_topk_intersect(),
  F2 = gene_vec_to_mat(agg_hg, gene) %>% colwise_topk_intersect2(),
  times = 5)





query_gene_rank_topk <- function(query_vec,
                                 subject_mat,
                                 gene,
                                 k = 1000,
                                 ncores = 1) {
  
  genes <- rownames(subject_mat)
  
  stopifnot(gene %in% genes, identical(names(query_vec), genes))
  
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



query_gene_rank_topk2 <- function(query_vec,
                                  subject_mat,
                                  gene,
                                  k = 1000,
                                  ncores = 1) {
  
  genes <- rownames(subject_mat)
  
  stopifnot(gene %in% genes, identical(names(query_vec), genes))
  
  query_topk <- topk_sort(query_vec, k)
  
  topk_l <- mclapply(genes, function(x) {
    subject_topk <- topk_sort(subject_mat[, x], k)
    topk_intersect(query_topk, subject_topk)
  }, mc.cores = ncores)
  
  names(topk_l) <- genes
  topk_sort <- sort(unlist(topk_l), decreasing = TRUE)
  rank_ix <- which(names(topk_sort) == gene)
  
  # return(list(Rank = rank_ix, Topk = topk_sort))
  return(rank_ix)
}





query_gene_rank_topk_all <- function(agg_l,
                                     gene,
                                     k = 1000,
                                     ncores = 1) {
  ids <- names(agg_l)
  genes <- rownames(agg_l[[1]])
  
  stopifnot(length(ids) > 0, gene %in% genes)
  
  rank_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  rownames(rank_mat) <- colnames(rank_mat) <- ids
  
  for (i in ids) {
    for (j in ids) {
      if (i == j) next
      query_vec <- agg_l[[i]][, gene]
      subject_mat <- agg_l[[j]]
      rank_mat[i, j] <- query_gene_rank_topk(query_vec, subject_mat, gene, k, ncores)$Rank
    }
  }

  return(rank_mat)
}


query_gene_rank_topk_all2 <- function(agg_l,
                                      gene,
                                      k = 1000,
                                      ncores = 1) {
  ids <- names(agg_l)
  genes <- rownames(agg_l[[1]])
  
  stopifnot(length(ids) > 0, gene %in% genes)
  
  rank_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  rownames(rank_mat) <- colnames(rank_mat) <- ids
  
  for (i in ids) {
    for (j in ids) {
      if (i == j) next
      query_vec <- agg_l[[i]][, gene]
      subject_mat <- agg_l[[j]]
      rank_mat[i, j] <- query_gene_rank_topk2(query_vec, subject_mat, gene, k, ncores)
    }
  }
  
  return(rank_mat)
}




res3 <- microbenchmark::microbenchmark(
  F1 = query_gene_rank_topk_all(agg_hg, gene, ncores = ncore),
  F2 = query_gene_rank_topk_all2(agg_hg, gene, ncores = ncore),
  times = 3
)



F1 = query_gene_rank_topk_all(agg_hg, gene, ncores = ncore) 
F2 = query_gene_rank_cor_all(agg_hg, gene, ncores = ncore)
