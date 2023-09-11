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
    
    old_k <- k
    k <- 1
    
    for (i in 2:length(vec_sort)) {
      if (vec_sort[i] == vec_sort[i - 1]) {
        k <- i - 2  # -2 b/c trailing index and want the position before tie
        break()
      }
    }
    
    if (k > old_k) stop("New index position k is larger than input k")
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




# Functions using ROCR:: for measuring ranked performance
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



# TODO:

get_curated_labels <- function(tf,
                               curated_df,
                               pc_df,
                               species,
                               remove_self = TRUE) {
  
  labels <- curated_df %>%
    filter(str_to_upper(TF_Symbol) == str_to_upper(tf)) %>% 
    distinct(Target_Symbol) %>%
    pull(Target_Symbol)
  
  if (remove_self) labels <- setdiff(str_to_upper(labels), str_to_upper(tf))

  labels <- if (species == "Human") {
    intersect(str_to_upper(labels), pc_df$Symbol)
  } else {
    intersect(str_to_title(labels), pc_df$Symbol)
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
    AUPRC_percentile_observed = ecdf(null_auprc)(auc$AUPRC),
    AUPRC_ratio_observed = auc$AUPRC / median(null_auprc),
    AUPRC_diff_observed = auc$AUPRC - median(null_auprc),
    AUROC = auc$AUROC,
    AUROC_percentile_observed = ecdf(null_auroc)(auc$AUROC),
    AUROC_ratio_observed = auc$AUROC / median(null_auroc),
    AUROC_diff_observed = auc$AUROC - median(null_auroc)
  )
  
  return(list(Perf_df = df, Null = null)) 
}




# TODO:

curated_obs_and_null_auc <- function(tf,
                                     rank_df,
                                     score_col,
                                     curated_df,
                                     label_all,
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
                               pc_df = pc_df,
                               species = species,
                               remove_self = TRUE)
  
  if (length(labels) == 0) stop(paste("No labels retrieved for", tf))
                   
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
                                          pc_df,
                                          species,
                                          n_samps = 1000,
                                          ncores = 1,
                                          verbose = TRUE) {
  
  stopifnot(all(tfs %in% names(rank_l)))
  
  auc_l <- lapply(tfs, function(tf) {
    
    if (verbose) message(paste(tf, Sys.time()))
    
    try(curated_obs_and_null_auc(
      tf = tf,
      rank_df = rank_l[[tf]],
      score_col = score_col,
      curated_df = curated_df,
      label_all = label_all,
      pc_df = pc_df,
      species = species,
      n_samps = n_samps,
      ncores = ncores))
    
  })
  
  # auc_df <- do.call(rbind, lapply(auc_l, `[[`, "Perf_df"))
  names(auc_l) <- tfs
  return(auc_l)
}
  



# TODO: 

save_curated_auc_list <- function(path,
                                  tfs,
                                  rank_l,
                                  score_col,
                                  curated_df,
                                  label_all,
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
  
  names(auc_l) <- colnames(score_mat)
  return(auc_l)
}



# TODO:

get_colwise_curated_auc_list <- function(tfs,
                                         agg_l,
                                         msr_mat,
                                         curated_df,
                                         species,
                                         ncores = 1,
                                         verbose = TRUE) {
  
  tf_auc_l <- lapply(tfs, function(tf) {
    
    if (verbose) message(paste(tf, Sys.time()))
    
    # Prepare matrix of aggregate coexpr vectors and their average score
    score_mat <- gene_vec_to_mat(agg_l, tf)
    score_mat <- score_mat[, which(msr_mat[tf, ] == 1), drop = FALSE]
    score_mat <- cbind(score_mat, Aggregate = rowMeans(score_mat))
    
    # Prepare curated labels, removing the TF itself if it is a target
    labels <- get_curated_labels(tf = tf,
                                 curated_df = curated_df,
                                 pc_df = pc_df,
                                 species = species,
                                 remove_self = TRUE)
    
    get_colwise_auc(score_mat = score_mat, labels = labels, ncores = ncores)
    
  })
  
  names(tf_auc_l) <- tfs
  return(tf_auc_l)
}
