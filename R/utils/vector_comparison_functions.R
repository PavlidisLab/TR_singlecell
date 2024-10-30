## Current house of functions dedicated to calculating similarity between
## ranked vectors

# TODO: consider removal of self gene (and k+1) # rank_vec <- rank_vec[names(rank_vec) != gene]
# TODO: topk_intersect() is just length of intersect, either rename or include topk sort
# TODO: binarize top/bottom k needs to use check k arg

source("R/utils/functions.R")
library(ROCR)



# Aggregate vectors are ranked with ties which causes clumping of NA==0 values,
# and so selecting the top k elements may include uninformative NAs. This checks
# for ties at the kth position, and if found, returns the first non-tied k index
# of vec_sort
# vec_sort: named numeric vector assumed to be sorted
# k: an integer
# returns: an integer

check_k <- function(vec_sort, k, decreasing = TRUE) {
  
  if (!decreasing) vec_sort <- -vec_sort
  
  vec_rank <- rank(-vec_sort, ties.method = "min")
  tally_ranks <- sort(table(vec_rank), decreasing = TRUE)
  
  if (tally_ranks[1] > 1) {  # ties were found, check position relative to k
    tie_start <- as.integer(names(head(tally_ranks)[1]))
    k <- min(k, tie_start - 1)
  }
  
  return(k)
}



# Sort (bigger = more important) the numeric vec and return the names of the top
# k elements. 
# vec: named numeric vector
# k: an integer
# check_k_arg: logical controls whether check_k() will be used
# return: a character vector of the names of the top k elements of vec

topk_sort <- function(vec, k, check_k_arg = TRUE, decreasing = TRUE) {
  
  vec_sort <- sort(vec, decreasing = decreasing)
  if (check_k_arg) k <- check_k(vec_sort, k = k, decreasing = decreasing)
  topk <- names(vec_sort[1:k])
  
  return(topk)
}




# Return the length of the intersect of vec1 and vec2
# TODO: This is poorly named or should also perform sorting

topk_intersect <- function(vec1, vec2) length(intersect(vec1, vec2))


  
# Binarize matrix such that top k and bottom k is 1, everything else 0
# Speed of name overlap versus logical gt/lt scores was basically identical
# TODO: 

binarize_topk_btmk <- function(mat, k, check_k_arg = TRUE) {
  
  bin_mat <- apply(mat, 2, function(vec) {
    
    vec_sort <- sort(vec, decreasing = TRUE)
  
    # Check if ties are encountered before k in either direction
    if (check_k_arg) {
      k_upper <- check_k(vec_sort, k = k)
      k_lower <- check_k(vec_sort, k = k, decreasing = FALSE)
    } else {
      k_upper <- k_lower <- k
    }
    
    # k_lower is relative to lowest value as k=1, flip so relative to dec. sort
    k_lower <- length(vec_sort) - (k_lower - 1)
    
    # Get the top and bottom genes and return binary vector of their presence
    topk <- names(vec_sort[1:k_upper])
    btmk <- names(vec_sort[k_lower:length(vec_sort)])
    bin_vec <- ifelse(names(vec) %in% c(topk, btmk), 1, 0)
    names(bin_vec) <- names(vec)
    
    return(bin_vec)
  })
  
  return(bin_mat)
}




# Return a matrix of the size of the top k intersect between all columns of mat
# mat: a named numeric matrix
# k: an integer
# check_k_arg: logical controls whether check_k() will be used
# returns: an ncol(mat) * ncol(mat) integer matrix of top k overlap of mat's columns 

colwise_topk_intersect <- function(mat, 
                                   k, 
                                   check_k_arg = TRUE,
                                   decreasing = TRUE) {
  
  col_list <- asplit(mat, 2)
  
  topk_list <- lapply(col_list, function(x) {
    topk_sort(x, k = k, check_k_arg = check_k_arg, decreasing = decreasing)
  })
  
  topk_mat <- outer(topk_list, topk_list, Vectorize(topk_intersect))
  
  return(topk_mat)
}



# TODO:

colwise_cor <- function(mat, cor_method = "spearman", ncores = 1) {
  cor_mat <- WGCNA::cor(x = mat, method = cor_method, nThreads = ncores)
  return(cor_mat)
}


# TODO:
# https://stackoverflow.com/a/66594545 Jaccard; outer faster than nested loop 

colwise_jaccard <- function(mat, k, check_k_arg = TRUE) {
  
  jaccard <- function(vec1, vec2) {
    sum(vec1 & vec2, na.rm = TRUE) / sum(vec1 | vec2, na.rm = TRUE)
  }
  
  bin_mat <- binarize_topk_btmk(mat, k = k, check_k_arg = check_k_arg)
  col_list <- asplit(bin_mat, 2)
  jacc_mat <- outer(col_list, col_list, Vectorize(jaccard))
  
  return(jacc_mat)
}



# Nested loop faster than outer(): skipping diagonal outweighs vectorized
# TODO: doc

colwise_topk_auprc <- function(mat, k) {
  
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




# TODO:

pair_colwise_cor <- function(mat1, mat2, cor_method = "spearman", ncores = 1) {
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  
  cor_l <- mclapply(1:ncol(mat1), function(x) {
    cor(mat1[, x], mat2[, x], method = cor_method, use = "pairwise.complete.obs")
  }, mc.cores = ncores)
  
  names(cor_l) <- colnames(mat1)
  return(unlist(cor_l))
}



# TODO:

pair_colwise_topk <- function(mat1, mat2, k = 200, ncores = 1) {
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  
  topk_l <- mclapply(1:ncol(mat1), function(x) {
    topk_intersect(topk_sort(vec = mat1[, x], k = k),
                   topk_sort(vec = mat2[, x], k = k))
  }, mc.cores = ncores)
  
  names(topk_l) <- colnames(mat1)
  return(unlist(topk_l))
}



# TODO:

pair_shuffle_cor <- function(mat1, mat2, cor_method = "spearman", ncores = 1) {
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  sample1 <- sample(1:ncol(mat1), ncol(mat1), replace = TRUE)
  sample2 <- sample(1:ncol(mat2), ncol(mat2), replace = TRUE)
  
  cor_l <- mclapply(1:ncol(mat1), function(x) {
    cor(mat1[, sample1[x]], mat2[, sample2[x]], 
        method = cor_method, use = "pairwise.complete.obs")
  }, mc.cores = ncores)
  
  return(unlist(cor_l))
}



# TODO:

pair_shuffle_topk <- function(mat1, mat2, k = 200, ncores = 1) {
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  sample1 <- sample(1:ncol(mat1), ncol(mat1), replace = TRUE)
  sample2 <- sample(1:ncol(mat2), ncol(mat2), replace = TRUE)
  
  
  topk_l <- mclapply(1:ncol(mat1), function(x) {
    topk_intersect(topk_sort(vec = mat1[, sample1[x]], k = k),
                   topk_sort(vec = mat2[, sample2[x]], k = k))
  }, mc.cores = ncores)
  
  return(unlist(topk_l))
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





# Functions using ROCR:: for measuring ranked performance
# TODO: consider just returning both; null as df instead of nested list
# ------------------------------------------------------------------------------



# This uses the ROCR package to a data.frame of precision and recall (PR), TPR
# and FPR (ROC), or both, calculated at each step of score/label vec
# score_vec is a numeric vector of scores where higher values == more important
# label_vec is a binary vector of equal length to score_vec where 1/TRUE == pos
# measure is one of "AUROC", "AUPRC", "both"
# returns a list

get_performance_df <- function(score_vec,
                               label_vec,
                               measure) {
  
  stopifnot(identical(length(score_vec), length(label_vec)),
            is.numeric(score_vec),
            measure %in% c("ROC", "PR", "both"))
  
  # Negatives as 0, positives as 1
  label_vec <- factor(as.integer(label_vec), levels = c(0, 1), ordered = TRUE)
  
  pred <- ROCR::prediction(predictions = score_vec, labels = label_vec)
  
  if (measure == "ROC") {
    
    roc <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    
    perf <- data.frame(TPR = unlist(roc@y.values), 
                       FPR = unlist(roc@x.values))
    
  } else if (measure == "PR") {
    
    pr <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
    
    perf <- data.frame(Precision = unlist(perf@y.values),
                       Recall = unlist(perf@x.values))
    
  } else {
    
    roc <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    pr <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
    
    perf <- data.frame(TPR = unlist(roc@y.values),
                       FPR = unlist(roc@x.values),
                       Precision = unlist(pr@y.values),
                       Recall = unlist(pr@x.values))
  }
  
  return(perf)
}



# This uses the ROCR package to return the area under the PR curve (AUPRC),
# the area under ROC (AUROC), or both.
# score_vec is a numeric vector of scores where higher values == more important
# label_vec is a binary vector of equal length to score_vec where 1/TRUE == pos
# measure is one of "AUROC", "AUPRC", "both"
# returns a list

get_auc <- function(score_vec,
                    label_vec,
                    measure) {
  
  stopifnot(identical(length(score_vec), length(label_vec)),
            is.numeric(score_vec),
            measure %in% c("AUROC", "AUPRC", "both"))
  
  # Negatives as 0, positives as 1
  label_vec <- factor(as.integer(label_vec), levels = c(0, 1), ordered = TRUE)
  
  pred <- ROCR::prediction(predictions = score_vec, labels = label_vec)
  
  if (measure == "AUROC") {
    perf <- list(
      AUROC = ROCR::performance(pred, measure = "auc")@y.values[[1]]
    )
  } else if (measure == "AUPRC") {
    perf <- list(
      AUPRC = ROCR::performance(pred, measure = "aucpr")@y.values[[1]]
    )
  } else {
    perf <- list(
      AUROC = ROCR::performance(pred, measure = "auc")@y.values[[1]],
      AUPRC = ROCR::performance(pred, measure = "aucpr")@y.values[[1]]
    )
  }
  
  return(perf)
}



# Retrieve curated targets for the given TF, casting wide enough to capture
# genes with 1:1 orthologous matches between mouse and human

get_curated_labels <- function(tf,
                               curated_df,
                               ortho_df,
                               pc_df,
                               species,
                               remove_self = TRUE) {
  
  stopifnot(species %in% c("Human", "Mouse"))
  
  # Get all valid orthologous symbols for input TF
  ortho_tf <- filter(ortho_df, 
                     Symbol_hg == str_to_upper(tf) | 
                     Symbol_mm == str_to_title(tf))
  
  tf <- union(
    filter(pc_df, Symbol == tf)$Symbol,
    c(ortho_tf$Symbol_hg, ortho_tf$Symbol_mm)
  )
  
  # Extract all target genes matching any of the ortho TF symbols
  labels <- curated_df %>%
    filter(str_to_upper(TF_Symbol) %in% str_to_upper(tf)) %>%
    distinct(Target_Symbol) %>%
    pull(Target_Symbol)
  
  # Remove the TF if it is also its own target
  if (remove_self) labels <- setdiff(str_to_upper(labels), str_to_upper(tf))

  # For labels, get the correct species symbol if it exists
  if (species == "Human") {
    
    ortho_labels <- 
      filter(pc_ortho, Symbol_hg %in% str_to_upper(labels))$Symbol_hg
    
    labels <- union(
      ortho_labels,
      filter(pc_df, Symbol %in% str_to_upper(labels))$Symbol)

  } else {
    
    ortho_labels <- 
      filter(pc_ortho, Symbol_mm %in% str_to_title(labels))$Symbol_mm
    
    labels <- union(
      ortho_labels,
      filter(pc_df, Symbol %in% str_to_title(labels))$Symbol)
  }
  
  return(labels)
}



# TODO:

get_null_performance <- function(score_vec,
                                 label_all,
                                 measure,
                                 n_target, 
                                 n_samps, 
                                 ncores = 1) {
  
  stopifnot(is.numeric(score_vec),
            length(label_all %in% names(score_vec)) > 0,
            measure %in% c("AUROC", "AUPRC", "both"))
  
  null_perf <- mclapply(1:n_samps, function(x) {
    
    null_labels <- sample(label_all, size = n_target, replace = FALSE)
    null_vec <- names(score_vec) %in% null_labels
    get_auc(score_vec = score_vec, label_vec = null_vec, measure = measure)
    
  }, mc.cores = ncores)
  
  return(null_perf)
}




# TODO:

summarize_obs_and_null_auc <- function(tf,
                                       score_vec,
                                       label_vec,
                                       label_all,
                                       n_samps = 1000,
                                       ncores = 1) {
  
  n_target <- sum(label_vec)
  
  auc <- get_auc(score_vec, label_vec = label_vec, measure = "both")
  
  null <- get_null_performance(score_vec = score_vec,
                               label_all = label_all,
                               measure = "both",
                               n_target = n_target, 
                               n_samps = n_samps, 
                               ncores = ncores)
  
  null_auprc <- unlist(lapply(null, `[[`, "AUPRC"))
  null_auroc <- unlist(lapply(null, `[[`, "AUROC"))
  
  df <- data.frame(
    Symbol = tf,
    N_targets = n_target,
    AUPRC = auc$AUPRC,
    AUPRC_quantile = ecdf(null_auprc)(auc$AUPRC),
    AUPRC_diff = auc$AUPRC - median(null_auprc),
    AUROC = auc$AUROC,
    AUROC_quantile = ecdf(null_auroc)(auc$AUROC),
    AUROC_diff = auc$AUROC - median(null_auroc)
  )
  
  return(list(Perf_df = df, Null = null)) 
}




# TODO:

curated_obs_and_null_auc <- function(tf,
                                     rank_df,
                                     score_col,
                                     curated_df,
                                     label_all,
                                     ortho_df,
                                     pc_df,
                                     species,
                                     n_samps = 1000,
                                     ncores = 1) {
  
  stopifnot(species %in% c("Mouse", "Human"))
  
  # Extract the ranking/scores for the given TF, removing the TF itself
  
  rank <- rank_df %>%
    filter(Symbol != tf) %>%
    arrange(desc(!!sym(score_col)))
    
  score_vec <- setNames(rank[[score_col]], rank$Symbol)
    
  # Extract the labels for the given TF, removing the TF itself as a target
    
  labels <- get_curated_labels(tf = tf,
                               curated_df = curated_df,
                               ortho_df = ortho_df,
                               pc_df = pc_df,
                               species = species,
                               remove_self = TRUE)
  
  label_all <- label_all[label_all != tf]
  
  # No labels may result if the TF itself was the only label
  if (length(labels) == 0) {
    paste("No labels retrieved for", tf)
    return(NA)
  }
  
  label_vec <- names(score_vec) %in% labels

  auc_df <- summarize_obs_and_null_auc(tf = tf,
                                       score_vec = score_vec,
                                       label_vec = label_vec,
                                       label_all = label_all,
                                       n_samps = n_samps,
                                       ncores = ncores)
  
  return(auc_df)
}
  
  


# TODO:

curated_obs_and_null_auc_list <- function(tfs,
                                          rank_l,
                                          score_col,
                                          curated_df,
                                          label_all,
                                          ortho_df,
                                          pc_df,
                                          species,
                                          n_samps = 1000,
                                          ncores = 1,
                                          verbose = TRUE) {
  
  stopifnot(all(tfs %in% names(rank_l)))
  
  auc_l <- lapply(tfs, function(tf) {
    
    if (verbose) message(paste(tf, Sys.time()))
    
    curated_obs_and_null_auc(
      tf = tf,
      rank_df = rank_l[[tf]],
      score_col = score_col,
      curated_df = curated_df,
      label_all = label_all,
      ortho_df = ortho_df,
      pc_df = pc_df,
      species = species,
      n_samps = n_samps,
      ncores = ncores
    )
    
  })
  
  names(auc_l) <- tfs
  auc_l <- auc_l[!is.na(auc_l)]
  
  return(auc_l)
}
  



# TODO: 

save_curated_auc_list <- function(path,
                                  tfs,
                                  rank_l,
                                  score_col,
                                  curated_df,
                                  label_all,
                                  ortho_df,
                                  pc_df,
                                  species,
                                  n_samps = 1000,
                                  ncores = 1,
                                  verbose = TRUE,
                                  force_resave = FALSE)  {
  
  if (!file.exists(path) || force_resave) {
    
    auc_l <- curated_obs_and_null_auc_list(
      tfs = tfs,
      rank_l = rank_l,
      score_col = score_col,
      curated_df = curated_df,
      label_all = label_all,
      ortho_df = ortho_df,
      pc_df = pc_df,
      species = species,
      n_samps = n_samps,
      ncores = ncores,
      verbose = verbose)
    
    saveRDS(auc_l, path)
  
  }
    
  return(invisible(NULL))
}




# TODO:

get_colwise_auc <- function(score_mat,
                            labels,
                            ncores = 1) {
  
  auc_l <- mclapply(colnames(score_mat), function(x) {
    
    score_vec <- sort(score_mat[, x], decreasing = TRUE)
    label_vec <- names(score_vec) %in% labels
    get_auc(score_vec = score_vec, label_vec = label_vec, measure = "both")
    
  }, mc.cores = ncores)
  
  auc_df <- data.frame(
    ID = colnames(score_mat),
    AUROC = vapply(auc_l, `[[`, "AUROC", FUN.VALUE = numeric(1)),
    AUPRC = vapply(auc_l, `[[`, "AUPRC", FUN.VALUE = numeric(1)),
    row.names = NULL)
  
  return(auc_df)
}




# TODO: The names of this and the actual implementation seem to swap intent
# TODO:

summarize_avg_and_individual_auc <- function(auc_df, labels) {
  
  auroc_avg <- filter(auc_df, ID == "Average")$AUROC
  auprc_avg <- filter(auc_df, ID == "Average")$AUPRC
  auc_df_no_avg <- filter(auc_df, ID != "Average")
  
  summary_df <- data.frame(
    N_targets = length(labels),
    N_datasets = nrow(auc_df_no_avg),
    AUROC_quantile = ecdf(auc_df_no_avg$AUROC)(auroc_avg),
    AUPRC_quantile = ecdf(auc_df_no_avg$AUPRC)(auprc_avg)
  )
  
  return(list(AUC_df = auc_df, Summary_df = summary_df))
}





# TODO:

get_colwise_curated_auc_list <- function(tfs,
                                         agg_l,
                                         msr_mat,
                                         curated_df,
                                         ortho_df,
                                         pc_df,
                                         species,
                                         ncores = 1,
                                         verbose = TRUE) {
  
  tf_auc_l <- lapply(tfs, function(tf) {
    
    if (verbose) message(paste(tf, Sys.time()))
    
    # Prepare curated labels, removing the TF itself if it is a target
    labels <- get_curated_labels(tf = tf,
                                 curated_df = curated_df,
                                 ortho_df,
                                 pc_df = pc_df,
                                 species = species,
                                 remove_self = TRUE)
    
    # No labels may result if the TF itself was the only label
    if (length(labels) == 0) {
      return(NA)
    }
    
    # Prepare matrix of aggregate coexpr vectors and their average score
    score_mat <- gene_vec_to_mat(agg_l, gene = tf, msr_mat = msr_mat)
    score_mat <- cbind(score_mat, Average = rowMeans(score_mat))
    
    # Calculating the AUC by using each column as a score
    auc_df <- get_colwise_auc(score_mat, labels = labels, ncores = ncores)
    
    # Summarize
    summarize_avg_and_individual_auc(auc_df, labels)
    
  })
  
  names(tf_auc_l) <- tfs
  tf_auc_l <- tf_auc_l[!is.na(tf_auc_l)]
  
  return(tf_auc_l)
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
