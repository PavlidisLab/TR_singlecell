## House functions dedicated to calculating similarity and AUC performance 
## between vectors and matrices

source("R/utils/functions.R")
library(ROCR)



# Aggregate vectors have NA values imputed to 0, and are ranked with ties 
# causing clumping of these uninformative values. Selecting the TopK elements
# may include these tied values. This checks for ties at the kth position, and 
# if found, returns a new smaller k value that excludes these ties.
# vec_sort: named numeric vector assumed to be sorted
# k: an integer specifying how many top elements to select
# decreasing: logical; if TRUE, selects the highest k values, else lowest
# returns: an integer

check_k <- function(vec_sort, k, decreasing = TRUE) {
  
  stopifnot(k <= length(vec_sort))
  
  if (!decreasing) vec_sort <- -vec_sort
  
  vec_rank <- rank(-vec_sort, ties.method = "min")[1:k]
  rank_at_k <- vec_rank[k]
  
  if (rank_at_k < k) {  # min ties results in tied value at rank value under k
     k <- rank_at_k - 1 # the value before the block of ties
  }

  return(k)
}




# Sort vec and return the names of the Top k elements. 
# vec: named numeric vector
# k: an integer specifying how many top elements to select
# check_k_arg: logical controls whether check_k() will be used
# decreasing: logical; if TRUE, selects the highest k values, else lowest
# return: a character vector of the names of the top k elements of vec

topk_sort <- function(vec, k, check_k_arg = TRUE, decreasing = TRUE) {
  
  vec_sort <- sort(vec, decreasing = decreasing)
  if (check_k_arg) k <- check_k(vec_sort, k = k, decreasing = decreasing)
  topk <- names(vec_sort[1:k])
  
  return(topk)
}




# Return the size of the intersect of the top k named elements of vec1 and vec2
# vec1: named numeric vector
# vec2: named numeric vector
# k: an integer specifying how many top elements to select
# check_k_arg: logical; if TRUE, adjust k to avoid ties at the cutoff
# decreasing: logical; if TRUE, selects the highest k values, else lowest
# return: an integer of the size of the top K overlap

topk_intersect <- function(vec1, 
                           vec2, 
                           k, 
                           check_k_arg = TRUE, 
                           decreasing = TRUE) {
  
  vec_sort1 <- topk_sort(vec1, k, check_k_arg, decreasing)
  vec_sort2 <- topk_sort(vec2, k, check_k_arg, decreasing)
  topk <- length(intersect(vec_sort1, vec_sort2))
  
  return(topk)
}




# This helper is used to get length of vectors already subset for top k, 
# to avoid redundant sorting when calculated over all pairs of columns

length_intersect <- function(topk_vec1, topk_vec2) {
  length(intersect(topk_vec1, topk_vec2))
}




# Binarize matrix such that top k and bottom k *of each column* is 1, else 0.
# Used for downstream Jaccard similarity calculations.
# mat: numeric matrix with named elements
# k: an integer specifying how many top and bottom elements to select
# check_k_arg: logical; if TRUE, adjust k to avoid ties at the cutoff
# return: a binary matrix where the top k and bottom k elements in each column 
# are 1, all others are 0.

binarize_topk_btmk <- function(mat, k, check_k_arg = TRUE) {
  
  bin_mat <- apply(mat, 2, function(vec) {
    
    # sort values in decreasing order, preserving names
    vec_sort <- sort(vec, decreasing = TRUE)
  
    # check if k needs adjustment due to ties at the cutoff
    if (check_k_arg) {
      k_upper <- check_k(vec_sort, k = k)
      k_lower <- check_k(vec_sort, k = k, decreasing = FALSE)
    } else {
      k_upper <- k_lower <- k
    }
    
    # Convert k_lower from "relative to ascending order" to "relative to descending order"
    # Since vec_sort is in decreasing order, we need to count from the end
    k_lower <- length(vec_sort) - (k_lower - 1)
    
    # Select the top-k and bottom-k elements based on adjusted indices
    topk <- names(vec_sort[1:k_upper])
    btmk <- names(vec_sort[k_lower:length(vec_sort)])
    
    # Generate a binary vector indicating presence in top k or bottom k
    bin_vec <- ifelse(names(vec) %in% c(topk, btmk), 1, 0)
    names(bin_vec) <- names(vec)
    bin_vec
    
  })
  
  return(bin_mat)
}



# The following functions are used to generate measures of similarity between
# columns of a matrix
# ------------------------------------------------------------------------------



# Compute the top-k intersection matrix for all pairs of columns in a matrix.
# mat: a named numeric matrix
# k: an integer specifying how many top elements to select.
# check_k_arg: logical; if TRUE, adjust k to avoid ties at the cutoff
# decreasing: logical; if TRUE, selects the highest k values, else lowest
# returns: an ncol(mat) * ncol(mat) integer matrix of top k overlap of mat's columns 

colwise_topk_intersect <- function(mat, 
                                   k, 
                                   check_k_arg = TRUE,
                                   decreasing = TRUE) {
  
  # Convert matrix into a list of numeric vectors
  col_list <- asplit(mat, 2)
  
  # Pre-subset the top k elements for each column
  topk_list <- lapply(col_list, function(x) {
    topk_sort(x, k = k, check_k_arg = check_k_arg, decreasing = decreasing)
  })
  
  # Compute pairwise intersection sizes of top k elements
  topk_mat <- outer(topk_list, topk_list, Vectorize(length_intersect))
  
  return(topk_mat)
}




# Wrapper to WGCNA's correlation -- defaults to calculating Spearman's cor
# between each pair of columns in mat

colwise_cor <- function(mat, cor_method = "spearman", ncores = 1) {
  cor_mat <- WGCNA::cor(x = mat, method = cor_method, nThreads = ncores)
  return(cor_mat)
}




# Compute the pairwise Jaccard similarity between all columns of a matrix. This
# measures the proportion of shared top/bottom k elements between columns. 
# mat: numeric matrix with named elements.
# k: an integer specifying how many top and bottom elements to consider.
# check_k_arg: logical; if TRUE, adjust k to avoid ties at the cutoff.
# return: an ncol(mat) x ncol(mat) symmetric matrix of Jaccard similarities.
# Uses `outer()` instead of a nested loop for efficiency.
# https://stackoverflow.com/a/66594545 

colwise_jaccard <- function(mat, k, check_k_arg = TRUE) {
  
  # Jaccard similarity function for binary vectors
  jaccard <- function(vec1, vec2) {
    sum(vec1 & vec2, na.rm = TRUE) / sum(vec1 | vec2, na.rm = TRUE)
  }
  
  # Convert matrix to a binary matrix where top/bottom k elements are 1
  bin_mat <- binarize_topk_btmk(mat, k = k, check_k_arg = check_k_arg)
  
  # Convert columns to a list of binary vectors
  col_list <- asplit(bin_mat, 2)
  
  # Compute pairwise Jaccard similarities
  jacc_mat <- outer(col_list, col_list, Vectorize(jaccard))
  
  return(jacc_mat)
}




# Compute AUROC or AUPRC scores for all column pairs in a matrix.
# Uses one column as scores and the top k elements of another column as labels.
# mat: numeric matrix with named elements.
# k: an integer specifying how many top elements to select as labels.
# measure: one of "AUROC" or "AUPRC".
# check_k_arg: logical; if TRUE, adjust k to avoid ties at the cutoff
# decreasing: logical; if TRUE, selects the highest k values, else lowest
# return: a symmetric matrix where entry [i, j] contains the AUC score when using 
# column i as scores and column j as the binary label vector.
# Performance note:
# Nested loop faster than outer(): skipping diagonal outweighs vectorized

colwise_topk_auc <- function(mat, 
                             k, 
                             measure, 
                             check_k_arg = TRUE, 
                             decreasing = TRUE) {
  
  stopifnot(measure %in% c("AUROC", "AUPRC"))
  
  # Init result matrix
  auc_mat <- matrix(1, nrow = ncol(mat), ncol = ncol(mat))
  colnames(auc_mat) <- rownames(auc_mat) <- colnames(mat)
  
  # Pre-sort scores and top k labels
  sorted_scores <- lapply(seq_len(ncol(mat)), function(i) {
     sort(mat[, i], decreasing = decreasing)
  })
  
  topk_labels <- lapply(seq_len(ncol(mat)), function(j) {
    topk_sort(mat[, j], k, check_k_arg, decreasing)
  })
  
  # Nested loop for AUC comparison
  for (i in 1:nrow(auc_mat)) {
    for (j in 1:ncol(auc_mat)) {
      if (i == j) next  # Skip diagonal (perfect retrieval = 1)
      
      score_vec <- sorted_scores[[i]]
      label_vec <- (names(score_vec) %in% topk_labels[[j]]) * 1  # Binary vector
      auc_mat[i, j] <- get_auc(score_vec, label_vec, measure)[[1]]
    }
  }
  
  return(auc_mat)
}




# The following functions are used to generate measures of similarity between
# index-matched columns between pairs of matrices. This means that the result 
# is just a vector of similarities, instead of a matrix from taking every pair
# of columns between the two matrices. 
# ------------------------------------------------------------------------------



# Compute the correlation between matched columns of two matrices.
# mat1: numeric matrix with named elements.
# mat2: numeric matrix with named elements.
# cor_method: correlation method fed into WGCNA::cor
# ncores: integer specifying the number of cores for parallel execution.
# return: a named numeric vector where each element corresponds to the number 
# of overlapping top k elements between the two columns.

pair_colwise_cor <- function(mat1, mat2, cor_method = "spearman", ncores = 1) {
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  
  cor_l <- mclapply(1:ncol(mat1), function(x) {
    
    WGCNA::cor(mat1[, x], 
               mat2[, x], 
               method = cor_method, 
               use = "pairwise.complete.obs")
    
  }, mc.cores = ncores)
  
  names(cor_l) <- colnames(mat1)
  return(unlist(cor_l))
}




# Compute the top k intersection between matched columns of two matrices
# mat1: numeric matrix with named elements.
# mat2: numeric matrix with named elements.
# k: an integer specifying how many top elements to select
# check_k_arg: logical; if TRUE, adjust k to avoid ties at the cutoff
# decreasing: logical; if TRUE, selects the highest k values, else lowest
# ncores: integer specifying the number of cores for parallel execution.
# return: a named numeric vector where each element corresponds to the number 
# of overlapping top k elements between the two columns.

pair_colwise_topk <- function(mat1, 
                              mat2, 
                              k, 
                              check_k_arg = TRUE,
                              decreasing = TRUE,
                              ncores = 1) {
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  
  topk_l <- mclapply(1:ncol(mat1), function(x) {
    
    topk_intersect(vec1 = mat1[, x],
                   vec2 = mat2[, x],
                   k = k,
                   check_k_arg = check_k_arg,
                   decreasing = decreasing)
    
  }, mc.cores = ncores)
  
  names(topk_l) <- colnames(mat1)
  return(unlist(topk_l))
}




# Calculate the correlation between shuffled columns of two matrices (for null)
# mat1: numeric matrix with named elements.
# mat2: numeric matrix with named elements.
# cor_method: correlation method fed into WGCNA::cor
# ncores: integer specifying the number of cores for parallel execution.
# return: a named numeric vector where each element corresponds to the number 
# of overlapping top k elements between the two columns.

pair_shuffle_cor <- function(mat1, mat2, cor_method = "spearman", ncores = 1) {
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  sample1 <- sample(ncol(mat1))
  sample2 <- sample(ncol(mat2))
  
  cor_l <- mclapply(1:ncol(mat1), function(x) {
    
    WGCNA::cor(mat1[, sample1[x]], 
               mat2[, sample2[x]], 
               method = cor_method, 
               use = "pairwise.complete.obs")
    
  }, mc.cores = ncores)
  
  return(unlist(cor_l))
}




# Compute the top k intersection between shuffled columns of two matrices (for null)
# mat1: numeric matrix with named elements.
# mat2: numeric matrix with named elements.
# k: an integer specifying how many top elements to select
# check_k_arg: logical; if TRUE, adjust k to avoid ties at the cutoff
# decreasing: logical; if TRUE, selects the highest k values, else lowest
# ncores: integer specifying the number of cores for parallel execution.
# return: a named numeric vector where each element corresponds to the number 
# of overlapping top k elements between the two columns.

pair_shuffle_topk <- function(mat1, 
                              mat2, 
                              k, 
                              check_k_arg = TRUE,
                              decreasing = TRUE,
                              ncores = 1) {
  
  stopifnot(identical(colnames(mat1), colnames(mat2)))
  sample1 <- sample(ncol(mat1))
  sample2 <- sample(ncol(mat2))
  
  topk_l <- mclapply(1:ncol(mat1), function(x) {
    
    topk_intersect(vec1 = mat1[, sample1[x]],
                   vec2 = mat2[, sample2[x]],
                   k = k,
                   check_k_arg = check_k_arg,
                   decreasing = decreasing)
    
  }, mc.cores = ncores)
  
  return(unlist(topk_l))
}





# Functions using ROCR:: for generating AUCs for ranked/scored vectors using
# binary gene vectors corresponding to curated targets as labels.
# NOTE: It was silly to specify which measure (ROC/AUROC or PR/AUPRC) to use, 
# since I always called both.
# ------------------------------------------------------------------------------



# Compute and return a data frame of performance metrics (ROC/TPR+FPR or PR) at
# each step of the provided score vector using the ROCR package.
# score_vec: numeric vector of scores where higher values indicate higher importance.
# label_vec: binary numeric vector (0/1) of the same length as score_vec, where 1 = positive.
# measure: One of "ROC", "PR", or "both", specifying the desired performance metric.
# return: A data frame containing the requested performance metrics at each score threshold.
#          - If measure == "ROC", returns columns TPR (True Positive Rate) and FPR (False Positive Rate).
#          - If measure == "PR", returns columns Precision and Recall.
#          - If measure == "both", returns all four columns.

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
    
    perf <- data.frame(Precision = unlist(pr@y.values),
                       Recall = unlist(pr@x.values))
    
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




# Compute the Area Under the Curve (AUC) for ROC (AUROC) and/or Precision-Recall 
# (AUPRC) curves using ROCR package.
# score_vec: numeric vector of scores where higher values indicate higher importance.
# label_vec: binary numeric vector (0/1) of the same length as score_vec, where 1 = positive.
# measure: one of "AUROC", "AUPRC", or "both", specifying the desired metric.
# return: named list containing:
#          - AUROC (if measure == "AUROC" or "both")
#          - AUPRC (if measure == "AUPRC" or "both")

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
  
  tf_auc_l <- mclapply(tfs, function(tf) {
    
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
    
    # Matrix of individual experiment coexpr profiles for given TR
    score_mat <- gene_vec_to_mat(agg_l, gene = tf, msr_mat = msr_mat)
    score_mat[tf, ] <- NA  # prevent self cor from being #1 rank
    
    # When a TF-gene was not co-measured, impute to the median NA
    msr_mat <- msr_mat[, colnames(score_mat)]
    med <- median(score_mat, na.rm = TRUE)
    score_mat[msr_mat == 0] <- med
    
    # Add average (same as global rankings) to matrix for comparison
    score_mat <- cbind(score_mat, Average = rowMeans(score_mat))
    
    # Calculating the AUC by using each column as a score
    auc_df <- get_colwise_auc(score_mat, labels = labels, ncores = ncores)
    
    # Summarize
    summarize_avg_and_individual_auc(auc_df, labels)
    
  }, mc.cores = ncores)
  
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
