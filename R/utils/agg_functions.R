library(WGCNA)
library(tidyverse)
source("R/utils/functions.R")




# TODO: ranksumrank function
# TODO: need to distinguish core agg function from apply over all
# TODO: binthresh function
# TODO: minrank function (should return best cell type and value as well)
# TODO: does `...` with col_subset declared work?
# TODO: aggregate mat symmetric?
# TODO: col subset is clumsy, actually needed?
# TODO: be explicit about shape of input mat (is t() burried bad?)
# TODO: col/rowrank sign confusing?



# Rank matrix columns such that 1=best. Return as same dimension as input.

colrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {
  rank_mat <- apply(-mat, 2, rank, ties.method = ties_arg, na.last = na_arg)
  return(rank_mat)
}


# Rank matrix rows such that 1=best. Return as same dimension as input.

rowrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {
  rank_mat <- apply(-mat, 1, rank, ties.method = ties_arg, na.last = na_arg)
  return(t(rank_mat))
}


# TODO: 

allrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {
  
  rmat <- rank(-mat, ties.method = ties_arg, na.last = na_arg)
  rmat <- matrix(rmat, nrow = nrow(mat), ncol = ncol(mat))
  rownames(rmat) <- colnames(rmat) <- rownames(mat)
  return(rmat)
}


# TODO:

get_cor_mat <- function(mat, 
                        cor_method = "pearson",
                        lower_tri = TRUE,
                        ncores = 1) {
  
  cmat <- WGCNA::cor(
    mat, method = cor_method, use = "pairwise.complete.obs", nThreads = ncores)
  
  if (lower_tri) {
    diag(cmat) <- NA
    cmat[upper.tri(cmat)] <- NA
  }
  
  return(cmat)
}


# If a col of mat has fewer than min_count elements that are non-zero, set that
# col to NAs. This is done to produce an NA during correlation, instead of 
# allowing values resulting from fewer observations.

under_min_count_to_na <- function(mat, min_count = 20) {
  
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) < min_count
  mat[, na_genes] <- NA
  return(mat)
}


# Return a gene x gene matrix of 0s used to track coexpression aggregation
# across cell types. Rows are the full set of provided genes, while columns
# are an optional subset of genes.

init_agg_mat <- function(row_genes, col_genes = NULL) {
  
  if (!is.null(col_genes)) {
    
    stopifnot(all(col_genes %in% row_genes))
    
    amat <- matrix(0, nrow = length(row_genes), ncol = length(col_genes))
    rownames(amat) <- row_genes
    colnames(amat) <- col_genes
    
  } else {
      
    amat <- matrix(0, nrow = length(row_genes), ncol = length(row_genes))
    rownames(amat) <- colnames(amat) <- row_genes
  
  }
  
  return(amat)
}



# 

cor_and_rank <- function(mat,
                         cor_method = "pearson",
                         na_arg = "last",
                         ties_arg = "min") {

  stopifnot(na_arg %in% c("last", "mean"))

  cmat <- WGCNA::cor(
    x = mat, method = cor_method, use = "pairwise.complete.obs")

  if (na_arg == "last") {
    rmat <- colrank_mat(cmat, na_arg = TRUE, ties_arg = ties_arg)
  } else {
    rmat <- colrank_mat(cmat, na_arg = "keep", ties_arg = ties_arg)
    rmat <- na_to_mean(rmat)
  }

  return(rmat)
}



# Rank sum rank from Harris et al., 2021 (Jesse Gillis) 
# https://pubmed.ncbi.nlm.nih.gov/34015329/
# Rank coexpression (1=best) across cell types. Set NAs to network mean. Sum 
# the cell-type ranks, and then rank order these sums (1=best).


# 1: Column rank where NAs cors are set to 0.
all_RSR_aggregate1 <- function(mat,
                              meta, 
                              min_cell = 20, 
                              ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    cmat[is.na(cmat)] <- 0
    # diag(cmat) <- NA
    rmat <- colrank_mat(cmat)
    amat <- amat + rmat
    
  }
  
  amat <- colrank_mat(-amat, na_arg = "last", ties_arg = "min")
  
  return(amat)
}




# 2: All rank where NAs cors are set to 0.
all_RSR_aggregate2 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    # Get full symmetric corr mat. Then replace NAs with 0 and make triangular
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    cmat[is.na(cmat)] <- 0
    # diag(cmat) <- NA
    cmat[upper.tri(cmat)] <- NA
    
    rmat <- allrank_mat(cmat)
    amat <- amat + rmat
    
  }
  
  amat <- allrank_mat(-amat)
  
  return(amat)
}



# 3: Col rank where NA rank sent to mean rank
all_RSR_aggregate3 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    # diag(cmat) <- NA
    rmat <- colrank_mat(cmat, na_arg = "keep")
    rmat <- na_to_mean(rmat)
    amat <- amat + rmat
    
  }
  
  amat <- colrank_mat(-amat, na_arg = "last", ties_arg = "min")
  
  return(amat)
}



# 4: All rank where NA rank sent to mean rank
all_RSR_aggregate4 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    #
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    cmat[upper.tri(cmat)] <- NA
    
    rmat <- allrank_mat(cmat, na_arg = "keep")
    rmat[lower.tri(rmat) & is.na(rmat)] <- mean(rmat[lower.tri(rmat)], na.rm = TRUE)
  
    amat <- amat + rmat
    
  }
  
  amat <- allrank_mat(-amat)
  
  return(amat)
}



# 5: Col rank where NA rank sent to last
all_RSR_aggregate5 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    #
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    # diag(cmat) <- NA
    
    # na_arg = TRUE results in placing NAs in random order after last non-NA rank
    # rmat <- colrank_mat(cmat, na_arg = TRUE)
    
    # This sets all NA ranks to length of genes
    rmat <- colrank_mat(cmat, na_arg = "keep")
    rmat[is.na(rmat)] <- length(genes)
    
    amat <- amat + rmat
    
  }
  
  amat <- colrank_mat(-amat, na_arg = "last", ties_arg = "min")
  
  return(amat)
}



# 6: All rank where NA rank sent to last
all_RSR_aggregate6 <- function(mat,
                               meta,
                               min_cell = 20,
                               ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    cmat[upper.tri(cmat)] <- NA
    
    rmat <- allrank_mat(cmat, na_arg = "keep")
    
    n_comparisons <- length(genes) * (length(genes) - 1) / 2
    max_nonNA_rank <- max(rmat, na.rm = TRUE)
    n_na <- length(rmat[lower.tri(rmat) & is.na(rmat)])
    
    # This sets NA ranks to +1 increasing values after last non-NA rank
    # rmat[lower.tri(rmat) & is.na(rmat)] <- (max_nonNA_rank + 1) : (max_nonNA_rank + n_na)
    
    # This sets all NA ranks to the fixed last non-NA rank+1
    rmat[lower.tri(rmat) & is.na(rmat)] <- max_nonNA_rank + 1
    
    amat <- amat + rmat
    
  }
  
  amat <- allrank_mat(-amat)
  
  return(amat)
}


# Average Z-score
all_zscore_aggregate <- function(mat,
                                 meta,
                                 min_cell = 20,
                                 ...) {
  
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  na_mat <- amat
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Subset cell type counts and transpose matrix (genes as columns) for cor
    ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])
    ct_mat <- under_min_count_to_na(ct_mat, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE, ncores = 8)
    
    na_ix <- which(is.na(cmat), arr.ind = TRUE)
    na_mat[na_ix] <- na_mat[na_ix] + 1
    
    cmat[is.na(cmat)] <- 0
    # diag(cmat) <- NA
  
    zmat <- scale(cmat)
    zmat[is.na(zmat)] <- 0
    
    amat <- amat + zmat

  }
  
  agg_list <- list(
    Amat = amat,
    Avg_all = (amat / length(cts)),
    Avg_nonNA = (amat / (length(cts) - na_mat)),
    NA_mat = na_mat)
  
  return(agg_list)
}



# Get the cell type cors for given genes

all_celltype_cor <- function(mat,
                             meta,
                             min_cell = 20,
                             gene1,
                             gene2) {
  
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)

  cor_l <- lapply(cts, function(ct) {
    
    ct_mat <- t(mat[c(gene1, gene2), filter(meta, Cell_type == ct)$ID])
    
    if (sum(ct_mat[, 1] != 0) < min_cell || sum(ct_mat[, 2] != 0) < min_cell) {
      return(NA)
    }
    
    WGCNA::cor(ct_mat[, gene1], ct_mat[, gene2])

  })
  names(cor_l) <- cts
  
  cor_vec <- sort(unlist(cor_l), decreasing = TRUE)

  return(cor_vec)
}





binary_threshold <- function(mat, threshold) {
  
  # Threshold by pval
  
  binmat <- mclapply(1:ncol(cmat$p), function(x) {
    
    vec <- cmat$p[, x]
    passes <- names(sort(vec)[1:thresh])
    vec <- ifelse(names(vec) %in% passes, 1, 0)
    return(vec)
    
  }, mc.cores = ncore)
  
  binmat <- do.call(cbind, binmat)
  
}



all_threshold_aggregate <- function(mat,
                                    meta, 
                                    min_cell = 20, 
                                    cor_method = "pearson",
                                    threshold = 0,
                                    ...) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    cells <- filter(meta, Cell_type == ct)$ID
    ct_mat <- t(mat[, cells])
    na_genes <- apply(ct_mat, 2, function(x) sum(x != 0)) <= min_cell
    ct_mat[, na_genes] <- NA
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "all NAs"))
      next()
    }
    
    cmat <- WGCNA::corAndPvalue(
      x = mat, 
      method = cor_method, 
      use = "pairwise.complete.obs", 
      alternative = "greater")
    
    # threshold
    
    # amat <- amat + rmat
    
  }
  
  # amat <- colrank_mat(-amat)
  
  return(amat)
}
