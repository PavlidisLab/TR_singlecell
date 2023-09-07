## Current house of functions dedicated to calculating similarity between
## ranked vectors

# TODO: consider removal of self gene (and k+1) # rank_vec <- rank_vec[names(rank_vec) != gene]
# TODO: topk_intersect() is just length of intersect, either rename or include topk sort
# TODO: binarize top/bottom k needs to use check k arg

source("R/utils/functions.R")



# Aggregate vectors are ranked with ties which causes clumping of NA==0 values,
# and so selecting the top k elements may include uninformative NAs. This checks
# for ties at the kth position, and if found, returns the first non-tied k index
# of vec_sort
# vec: named numeric vector
# k: an integer
# returns: an integer

check_k <- function(vec_sort, k) {
  
  if (vec_sort[k] == vec_sort[k - 1]) {
    tie_start <- head(sort(table(vec_sort), decreasing = TRUE))[1]
    k <- sum(vec_sort > as.numeric(names(tie_start)), na.rm = TRUE)
  }
  
  return(k)
}



# Sort (bigger = more important) the numeric vec and return the names of the top
# k elements. 
# vec: named numeric vector
# k: an integer
# check_k_arg: logical controls whether check_k() will be used
# return: a character vector of the names of the top k elements of vec

topk_sort <- function(vec, k, check_k_arg = TRUE) {
  vec_sort <- sort(vec, decreasing = TRUE)
  if (check_k_arg) k <- check_k(vec_sort, k)
  names(vec_sort[1:k])
}



# Return the length of the intersect of vec1 and vec2
# TODO: This is poorly named or should also perform sorting

topk_intersect <- function(vec1, vec2) length(intersect(vec1, vec2))



# Return the AUPRC of vec1 as scores and vec2 as labels
# TODO doc and better arg names

vec_auc <- function(vec_scores, vec_labels, measure) {
  
  stopifnot(measure %in% c("AUROC", "AUPRC"))
  
  rank_df <- data.frame(Score = vec_scores, 
                        Label = names(vec_scores) %in% vec_labels)
  
  if (all(rank_df$Label)) return(1)
  if (all(!rank_df$Label)) return(0)
  
  get_au_perf(rank_df, 
              label_col = "Label", 
              score_col = "Score", 
              measure = measure)
}


  
# Binarize matrix such that top k and bottom k is 1, everything else 0
# TODO: doc

binarize_topk_btmk <- function(mat, k = 1000) {
  
  bin_mat <- apply(mat, 2, function(x) {
    sort_values <- sort(x, decreasing = TRUE)
    topk <- sort_values[k]
    btmk <- sort_values[length(x) - k + 1]
    ifelse(x >= topk | x <= btmk, 1, 0)
  })
  
  return(bin_mat)
}




# Return a matrix of the size of the top k intersect between all columns of mat
# mat: a named numeric matrix
# k: an integer
# check_k_arg: logical controls whether check_k() will be used
# returns: an ncol(mat) * ncol(mat) integer matrix of top k overlap of mat's columns 

colwise_topk_intersect <- function(mat, k = 1000, check_k_arg = TRUE) {
  
  col_list <- asplit(mat, 2)
  topk_list <- lapply(col_list, topk_sort, k = k, check_k_arg = check_k_arg)
  topk_mat <- outer(topk_list, topk_list, Vectorize(topk_intersect))
  
  return(topk_mat)
}


# TODO: doc

colwise_cor <- function(mat, cor_method = "spearman", ncores = 1) {
  cor_mat <- WGCNA::cor(x = mat, method = cor_method, nThreads = ncores)
  return(cor_mat)
}


# TODO: doc

colwise_jaccard <- function(mat) {
  
  jaccard <- function(vec1, vec2) {
    sum(vec1 & vec2, na.rm = TRUE) / sum(vec1 | vec2, na.rm = TRUE)
  }
  
  col_list <- asplit(mat, 2)
  jacc_mat <- outer(col_list, col_list, Vectorize(jaccard))
  return(jacc_mat)
}



# Nested loop faster than outer(): skipping diagonal outweighs vectorized
# TODO: doc

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




# Convert similarity matrices into a dataframe of unique pairs and values

get_similarity_pair_df <- function(cor_mat, topk_mat, jacc_mat) {
  
  df <-
    mat_to_df(cor_mat, symmetric = TRUE, value_name = "Scor") %>%
    cbind(
      Topk = mat_to_df(topk_mat, symmetric = TRUE)$Value,
      Jaccard = mat_to_df(jacc_mat, symmetric = TRUE)$Value
    )
  
  return(df)
}




###



# Performance of rankings in a data frame. 
# ------------------------------------------------------------------------------


# This uses the ROCR package to return a df of either precision and recall (PR) 
# or true and false positive rate (ROC) calculated for each prediction in rank_df.
# Label_col assumed to be value that can be represented as binary for 0=negative
# and 1=positive. 
# Score_col assumed to be numeric value where higher positive values carry more
# importance.

get_perf_df <- function(rank_df,
                        label_col,
                        score_col = NULL,
                        measure) {
  
  stopifnot(c(label_col, score_col) %in% colnames(rank_df),
            measure %in% c("ROC", "PR"))
  
  # Negatives as 0, positives as 1
  labels <- factor(as.numeric(rank_df[[label_col]]),
                   levels = c(0, 1),
                   ordered = TRUE)
  
  # If no scores column provided, assume row order as measure of importance
  scores <- if (is.null(score_col)) {
    message("No scores provided, using row order of rank_df as relative importance")
    1 / (1:length(labels))
  } else {
    rank_df[[score_col]]
  }
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  
  if (measure == "ROC") {
    
    perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    perf_df <- data.frame(TPR = unlist(perf@y.values),
                          FPR = unlist(perf@x.values))
  } else {
    
    perf <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
    perf_df <- data.frame(Precision = unlist(perf@y.values),
                          Recall = unlist(perf@x.values))
  }
  
  return(perf_df)
}


# This uses the ROCR package to return either the area under the PR curve (AUPRC)
# or area under ROC (AUROC) for every prediction made in rank_df.
# Label_col assumed to be value that can be represented as binary for 0=negative
# and 1=positive. 
# Score_col assumed to be numeric value where higher positive values carry more
# importance.

get_au_perf <- function(rank_df,
                        label_col,
                        score_col = NULL,
                        measure) {
  
  stopifnot(c(label_col, score_col) %in% colnames(rank_df),
            measure %in% c("AUROC", "AUPRC"))
  
  # Negatives as 0, positives as 1
  labels <- factor(as.numeric(rank_df[[label_col]]),
                   levels = c(0, 1),
                   ordered = TRUE)
  
  # If no scores column provided, assume row order as measure of importance
  scores <- if (is.null(score_col)) {
    message("No scores provided, using row order of rank_df as relative importance")
    1 / (1:length(labels))
  } else {
    rank_df[[score_col]]
  }
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  
  if (measure == "AUROC") {
    perf <- ROCR::performance(pred, measure = "auc")@y.values[[1]]
  } else {
    perf <- ROCR::performance(pred, measure = "aucpr")@y.values[[1]]
  }
  
  return(perf)
}





# TODO: check k
# TODO: name structure more explicit
# TODO: more explicit about removal of TF


get_rank_df <- function(gene_vec,
                        tf,
                        evidence_df,
                        evidence_col,
                        k = NULL) {

  gene_df <- data.frame(Symbol = names(gene_vec),
                        Score = unname(gene_vec))
  
  rank_df <- evidence_df %>%
    dplyr::select(Symbol, !!sym(evidence_col)) %>% 
    left_join(., gene_df, by = "Symbol") %>%
    filter(Symbol != tf) %>%
    arrange(!!sym(evidence_col))
  
  
  if (evidence_col == "Curated_target") {
    
    labels <- filter(evidence_df, Curated_target)[["Symbol"]]
    
    rank_df <- rank_df %>%
      mutate(Label = Symbol %in% labels) %>%
      arrange(desc(Score))
    
  } else {
    
    labels <- c(rep(1, k), rep(0, nrow(rank_df) - k))
    
    rank_df <- rank_df %>%
      mutate(Label = labels) %>%
      arrange(desc(Score))
  }

  return(rank_df)
}



# TODO: better name
# TODO: measure check

get_tf_performance <- function(agg_l,
                               msr_mat,
                               evidence_df,
                               evidence_col,
                               tf,
                               measure,
                               k = NULL) {
  
  # TODO:
  gene_mat <- gene_vec_to_mat(agg_l, tf)
  gene_mat <- gene_mat[, which(msr_mat[tf, ] == 1), drop = FALSE]
  gene_mat <- cbind(gene_mat, Aggregate = rowMeans(gene_mat))

  # TODO:
  auprc_l <- lapply(colnames(gene_mat), function(x) {
    
    message(x, Sys.time())
    
    gene_vec <- gene_mat[, x]
    rank_df <- get_rank_df(gene_vec, tf, evidence_df, evidence_col, k)
    get_au_perf(rank_df, label_col = "Label", score_col = "Score", measure = measure)
    
  })
  
  names(auprc_l) <- colnames(gene_mat)
  auprc_l <- unlist(auprc_l)
  return(auprc_l)
}




# TODO:

get_all_perf_df <- function(agg_l,
                            msr_mat,
                            evidence_df,
                            evidence_col,
                            tf,
                            measure,
                            k = NULL) {
  
  # TODO:
  gene_mat <- gene_vec_to_mat(agg_l, tf)
  gene_mat <- gene_mat[, which(msr_mat[tf, ] == 1), drop = FALSE]
  gene_mat <- cbind(gene_mat, Aggregate = rowMeans(gene_mat))
  
  # TODO:
  df_l <- lapply(colnames(gene_mat), function(x) {
    
    message(x, Sys.time())
    
    gene_vec <- gene_mat[, x]
    rank_df <- get_rank_df(gene_vec, tf, evidence_df, evidence_col, k)
    
    perf_df <- get_perf_df(rank_df, 
                           label_col = "Label", 
                           score_col = "Score",
                           measure = measure)
    perf_df$ID <- x
    perf_df
    
  })
  
  all_perf_df <- do.call(rbind, df_l)
  
  return(all_perf_df)
}








# TODO:
# The idea is to get performance for every gene vector in a network at 
# recovering a given TF's evidence, and seeing how well the TF vector itself
# compares

get_all_performance <- function(agg_l,
                                evidence_df,
                                evidence_col,
                                genes,
                                k = 1000,
                                ncores = 1) {

  # TODO:

  perf_l <- lapply(agg_l, function(agg_mat) {

    # TODO:

    perf <- mclapply(genes, function(x) {

      gene_vec <- agg_mat[, x]
      rank_df <- get_rank_df(gene_vec, evidence_df, evidence_col, k)
      auprc <- get_au_perf(rank_df, label_col = "Label", score_col = "Score", measure = "AUPRC")
      topk <- sum(rank_df$Label[1:k])
      data.frame(AUPRC = auprc, Topk = topk)

    }, mc.cores = ncores)

    perf_df <- data.frame(Symbol = genes, do.call(rbind, perf)) %>%
      arrange(desc(AUPRC))

    return(perf_df)

  })

  names(perf_l) <- names(agg_l)
  return(perf_l)
}



# This helper takes the output list of get_all_performance(), and returns an
# experiment by tf matrix, where each element represents the AUPRC rank of the
# given TF's aggregate vector at recovering evidence, relative to all other
# gene aggregate vectors in the given dataset.
# Eg [MKA, Ascl1] == 6153: Ascl1 had the 6153th best AUPRC rank at recovering
# Ascl1 evidence among all genes in the MKA dataset.

which_rank_mat <- function(perf_l) {

  tfs <- names(perf_l)

  which_l <- lapply(tfs, function(tf) {
    lapply(perf_l[[tf]], function(x) {
      which(x$Symbol == tf)
    })
  })

  which_mat <- do.call(cbind, which_l)
  colnames(which_mat) <- tfs

  return(which_mat)
}

