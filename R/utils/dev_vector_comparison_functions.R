# rank_vec is a gene ranking vector (1=best)
# rank_mat is a gene x gene (1=best)
# gene is a pcoding gene of interest corresponding to rank_vec
# k is cutoff for the top genes to use as labels for every column in rank_mat
# ncores is number of cores
# 
# Calculate the size of the intersect between the topk genes in rank_vec and
# the top k genes of every column of rank_mat (topk). Return a list of these 
# sorted topks, as well as the rank (1=best) of gene's topk with its name-matched 
# column in rank_mat

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




# rank_vec is a gene ranking vector (1=best)
# rank_mat is a gene x gene (1=best)
# gene is a pcoding gene of interest corresponding to rank_vec
# cor_method is the type of correlation to be fed to WGCNA::cor
# ncores is number of cores
# 
# Calculate the correlation between rank_vec and every column of rank_mat, and
# return a list of the sorted correlations, as well as the rank (1=best) of 
# gene's cor with its name-matched column in rank_mat

query_gene_rank_cor <- function(query_vec,
                                subject_mat,
                                gene,
                                cor_method = "spearman",
                                ncores = 1) {
  
  genes <- rownames(subject_mat)
  
  stopifnot(gene %in% genes, identical(names(query_vec), genes))
  
  rank_cor <-  WGCNA::cor(x = query_vec, 
                          y = subject_mat,
                          method = cor_method,
                          nThreads = ncores)
  
  rank_cor <- sort(rank_cor[1, ], decreasing = TRUE)
  rank_ix <- which(names(rank_cor) == gene)
  
  return(list(Rank = rank_ix, List = rank_cor))
}




# rank_vec is a gene ranking vector (1=best)
# rank_mat is a gene x gene (1=best)
# gene is a pcoding gene of interest corresponding to rank_vec
# k is cutoff for the top genes to use as labels for every column in rank_mat
# ncores is number of cores
# 
# Calculate the area under the precision recall curve (AUPRC) using rank_vec and
# the top k genes of every column of rank_mat. Return a list of the sorted 
# AUPRCs, as well as the rank (1=best) of gene's AUPRC with its name-matched 
# column in rank_mat

query_gene_rank_auprc <- function(query_vec,
                                  subject_mat,
                                  gene,
                                  k = 1000,
                                  ncores = 1) {
  
  genes <- rownames(subject_mat)
  
  stopifnot(gene %in% genes, identical(names(query_vec), genes))
  
  scores <- sort(query_vec, decreasing = TRUE)

  auprc_l <- mclapply(genes, function(x) {
    labels <- topk_sort(subject_mat[, x], k)
    vec_auprc(scores, labels)
  }, mc.cores = ncores)
  
  names(auprc_l) <- genes
  
  rank_auprc <- sort(unlist(auprc_l), decreasing = TRUE)
  rank_ix <- which(names(rank_auprc) == gene)
  
  return(list(Rank = rank_ix, List = rank_auprc))
}




# agg_l is a list of aggregate gene-gene rank matrices
# gene is a pcoding gene of interest
# msr is the type of similarity to be calculated
# cor_method is the type of correlation to be fed to WGCNA::cor
# k is cutoff for the top genes to use as labels for every column in rank_mat
# ncores is number of cores
# 
# Return an n x n matrix where n is equal to the count of matrices in agg_l. Each
# element corresponds to the rank (1=best) of the given gene's similarity rank
# across all pairs of experiments in agg_l.

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




query_gene_rank_cor_all <- function(agg_l,
                                    gene,
                                    cor_method = "spearman",
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
      rank_mat[i, j] <- query_gene_rank_cor(query_vec, subject_mat, gene, cor_method, ncores)$Rank
    }
  }
  
  return(rank_mat)
}




# TODO:

get_colwise_performance_df <- function(score_mat, labels, ncores = 1) {
  
  perf_df_l <- mclapply(colnames(score_mat), function(x) {
    
    score_vec <- sort(score_mat[, x], decreasing = TRUE)
    label_vec <- names(score_vec) %in% labels
    perf_df <- get_performance_df(score_vec, label_vec, measure = "both")
    perf_df$ID <- x
    return(perf_df)
    
  }, mc.cores = ncores)
  
  perf_df_all <- do.call(rbind, perf_df_l)
  
  return(perf_df_all)
}
