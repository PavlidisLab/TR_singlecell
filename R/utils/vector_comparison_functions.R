## Current house of functions dedicated to calculating similarity between
## ranked vectors

# TODO: consider removal of self gene (and k+1) # rank_vec <- rank_vec[names(rank_vec) != gene]
# TODO: row vs col gene in rank function when using a non-symmetrix matrix


source("R/utils/functions.R")


# TODO:

topk_sort <- function(vec, k) names(sort(vec, decreasing = TRUE)[1:k])


# TODO:

topk_intersect <- function(vec1, vec2) length(intersect(vec1, vec2))



# TODO:

vec_auprc <- function(vec1, vec2) {
  
  rank_df <- data.frame(Score = vec1, Label = names(vec1) %in% vec2)
  if (all(rank_df$Label)) return(1)
  if (all(!rank_df$Label)) return(0)
  get_au_perf(rank_df, label_col = "Label", score_col = "Score", measure = "AUPRC")
}


  

# Binarize matrix such that top k and bottom k is 1, everything else 0

binarize_topk_btmk <- function(mat, k = 1000) {
  
  bin_mat <- apply(mat, 2, function(x) {
    sort_values <- sort(x, decreasing = TRUE)
    topk <- sort_values[k]
    btmk <- sort_values[length(x) - k + 1]
    ifelse(x >= topk | x <= btmk, 1, 0)
  })
  
  return(bin_mat)
}




# TODO:

colwise_topk_intersect <- function(mat, k = 1000) {
  
  col_list <- asplit(mat, 2)
  topk_list <- lapply(col_list, topk_sort, k)
  topk_mat <- outer(topk_list, topk_list, Vectorize(topk_intersect))
  return(topk_mat)
}


# TODO:

colwise_cor <- function(mat, cor_method = "spearman", ncores = 1) {
  cor_mat <- WGCNA::cor(x = mat, method = cor_method, nThreads = ncores)
  return(cor_mat)
}



# TODO:

colwise_jaccard <- function(mat) {
  
  jaccard <- function(vec1, vec2) {
    sum(vec1 & vec2, na.rm = TRUE) / sum(vec1 | vec2, na.rm = TRUE)
  }
  
  col_list <- asplit(mat, 2)
  jacc_mat <- outer(col_list, col_list, Vectorize(jaccard))
  return(jacc_mat)
}



# Nested loop faster than outer(): skipping diagonal outweighs vectorized

colwise_topk_auprc <- function(mat, k = 1000) {
  
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



# TODO: smarter

# get_rank_df <- function(gene_vec, 
#                         evidence_df, 
#                         evidence_col,
#                         k = 1000) {
#   
#   
#   rank_df <- if (evidence_col == "Curated_target") {
#     
#     data.frame(Score = gene_vec) %>%
#       rownames_to_column(var = "Symbol") %>%
#       left_join(evidence_df[, c("Symbol", evidence_col)], by = "Symbol") %>%
#       # filter(Symbol != tf) %>%
#       arrange(!!sym(evidence_col)) %>%
#       mutate(Label = Symbol %in% filter(evidence_df, Curated_target)$Symbol) %>%
#       arrange(desc(Score))
#     
#   } else {
#     
#     data.frame(Score = gene_vec) %>%
#       rownames_to_column(var = "Symbol") %>%
#       left_join(evidence_df[, c("Symbol", evidence_col)], by = "Symbol") %>%
#       # filter(Symbol != tf) %>%
#       arrange(!!sym(evidence_col)) %>%
#       mutate(Label = c(rep(1, k), rep(0, nrow(.) - k))) %>%
#       arrange(desc(Score))
#     
#   }
#   
#   return(rank_df)
# }
# 
# 
# # TODO:
# 
# get_tf_performance <- function(agg_l,
#                                evidence_df,
#                                evidence_col,
#                                tf,
#                                k = 1000) {
#   
#   # TODO:
#   sc_tf_mat <- gene_vec_to_mat(agg_l, tf)
#   sc_tf_mat <- cbind(sc_tf_mat, Aggregate_rank = rowMeans(sc_tf_mat))
#   
#   # TODO:
#   auprc_l <- lapply(colnames(sc_tf_mat), function(x) {
#     gene_vec <- sc_tf_mat[, x]
#     rank_df <- get_rank_df(gene_vec, evidence_df, evidence_col, k)
#     get_au_perf(rank_df, label_col = "Label", score_col = "Score", measure = "AUPRC")
#   })
#   names(auprc_l) <- colnames(sc_tf_mat)
#   auprc_l <- unlist(auprc_l)
#   
#   return(auprc_l)
# }
# 
# 
# 
# # TODO: 
# 
# get_all_performance <- function(agg_l,
#                                 evidence_df,
#                                 evidence_col,
#                                 genes,
#                                 k = 1000,
#                                 ncores = 1) {
#   
#   # TODO:
#   
#   perf_l <- lapply(agg_l, function(agg_mat) {
#     
#     # TODO:
#     
#     perf <- mclapply(genes, function(x) {
#       
#       gene_vec <- agg_mat[, x]
#       rank_df <- get_rank_df(gene_vec, evidence_df, evidence_col, k)
#       auprc <- get_au_perf(rank_df, label_col = "Label", score_col = "Score", measure = "AUPRC")
#       topk <- sum(rank_df$Label[1:k])
#       data.frame(AUPRC = auprc, Topk = topk)
#       
#     }, mc.cores = ncores)
#     
#     perf_df <- data.frame(Symbol = genes, do.call(rbind, perf)) %>%
#       arrange(desc(AUPRC))
#     
#     return(perf_df)
#     
#   })
#   
#   names(perf_l) <- names(agg_l)
#   return(perf_l)
# }
# 
# 
# 
# # This helper takes the output list of get_all_performance(), and returns an
# # experiment by tf matrix, where each element represents the AUPRC rank of the
# # given TF's aggregate vector at recovering evidence, relative to all other
# # gene aggregate vectors in the given dataset.
# # Eg [MKA, Ascl1] == 6153: Ascl1 had the 6153th best AUPRC rank at recovering
# # Ascl1 evidence among all genes in the MKA dataset.
# 
# which_rank_mat <- function(perf_l) {
#   
#   tfs <- names(perf_l)
#   
#   which_l <- lapply(tfs, function(tf) {
#     lapply(perf_l[[tf]], function(x) {
#       which(x$Symbol == tf)
#     })
#   })
#   
#   which_mat <- do.call(cbind, which_l)
#   colnames(which_mat) <- tfs
#   
#   return(which_mat)
# }

