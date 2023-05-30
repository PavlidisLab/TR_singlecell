## Project functions
## TODO: init_agg_mat() reconsider mat instead of row_genes as input
## -----------------------------------------------------------------------------

library(parallel)
library(DescTools)
library(tidyverse)
library(Seurat)
library(WGCNA)
library(Matrix)


# Single cell coexpression aggregation
# ------------------------------------------------------------------------------



# Rank matrix columns such that 1=best. Return as same dimension as input.

colrank_mat <- function(mat, ties_arg = "random", na_arg = "keep") {
  rank_mat <- apply(-mat, 2, rank, ties.method = ties_arg, na.last = na_arg)
  return(rank_mat)
}



# Rank matrix rows such that 1=best. Return as same dimension as input.

rowrank_mat <- function(mat, ties_arg = "random", na_arg = "keep") {
  rank_mat <- apply(-mat, 1, rank, ties.method = ties_arg, na.last = na_arg)
  return(t(rank_mat))
}



# Rank the entire matrix jointly such that 1=best. Return as same dimension as input.

allrank_mat <- function(mat, ties_arg = "random", na_arg = "keep") {
  
  rank_mat <- rank(-mat, ties.method = ties_arg, na.last = na_arg)
  rank_mat <- matrix(rank_mat, nrow = nrow(mat), ncol = ncol(mat))
  rownames(rank_mat) <- colnames(rank_mat) <- rownames(mat)
  
  return(rank_mat)
}



# Generate column-wise correlation of mat, with arguments for whether to return
# the full symmetric matrix or a lower triangular with NAs on diagonal and upper.

get_cor_mat <- function(mat, 
                        cor_method = "pearson",
                        lower_tri = FALSE,
                        ncores = 1) {
  
  if (cor_method == "bicor") {
    cmat <- WGCNA::bicor(
      mat, use = "pairwise.complete.obs", nThreads = ncores)
    
  } else {
    cmat <- WGCNA::cor(
      mat, method = cor_method, use = "pairwise.complete.obs", nThreads = ncores)
  }
  
  if (lower_tri) {
    diag(cmat) <- NA
    cmat[upper.tri(cmat)] <- NA
  }
  
  return(cmat)
}



# Replace upper tri of mat with lower tri.

lowertri_to_symm <- function(mat, na_diag = FALSE) {
  
  mat[upper.tri(mat)] <-  t(mat)[upper.tri(mat)]
  if (na_diag) diag(mat) <- NA
  return(mat)
}



# If a col of mat has fewer than min_count non-zero elements, set that col to
# all NAs. This is done to produce an NA during correlation, instead of 
# allowing values resulting from fewer observations.

under_min_count_to_na <- function(mat, min_count = 20) {
  
  na_genes <- apply(mat, 2, function(x) sum(x != 0)) < min_count
  mat[, na_genes] <- NA
  return(mat)
}



# Subset mat to cell_type, set under min count genes to NA, and transpose.
# Expects mat is genes x cells, and will return a cells x genes mat for corr.

subset_and_filter <- function(mat, meta, cell_type, min_count = 20) {
  
  ids <- dplyr::filter(meta, Cell_type == cell_type)$ID
  ct_mat <- t(mat[, ids])
  ct_mat <- under_min_count_to_na(ct_mat, min_count)
  stopifnot(all(rownames(ct_mat) %in% meta$ID))
  return(ct_mat)
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



# Set NAs in matrix to 0

na_to_zero <- function(mat) {
  mat[is.na(mat)] <- 0
  return(mat)
}


# Set diag in matrix to 1

diag_to_one <- function(mat) {
  diag(mat) <- 1
  return(mat)
}


# Set upper triangle to NA

upper_to_na <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
}



# Rank sum rank aggregation from Jesse Gillis's group.
# https://pubmed.ncbi.nlm.nih.gov/27165153/
# https://pubmed.ncbi.nlm.nih.gov/25717192/
# https://pubmed.ncbi.nlm.nih.gov/34015329/
# A gene-gene correlation matrix is generated for each cell type (defined in
# meta), and these correlations are ranked, summed, and then re-ranked such
# that 1=best. Genes with counts under the min_count threshold for a given cell
# type are set to NA, to produce NA cor. Then, NA cors are set to 0 before 
# ranking. The output is a gene-gene matrix representing the rank of correlation
# aggregated across cell types.


# Here, ranking is done column-wise such that each column/gene vector is ranked
# only within that vector. Standardize controls whether the final ranks should
# be kept in their original integer order, or standardized to be [0, 1]

RSR_colrank <- function(mat,
                        meta,
                        min_cell = 20,
                        standardize = TRUE) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Get count matrix for current cell type, coercing low count genes to NAs
    
    ct_mat <- subset_and_filter(mat, meta, ct, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    # Get cell-type cor matrix: full symmetric for ranking, NA cors to 0, and
    # diag (self-cor) coerced to 1
    
    cmat <- ct_mat %>% 
      get_cor_mat(lower_tri = FALSE) %>% 
      na_to_zero() %>% 
      diag_to_one()
    
    # Column-wise rank the correlations and add to aggregate matrix
    
    rmat <- colrank_mat(cmat)
    amat <- amat + rmat
    
  }
  
  if (standardize) {
    amat <- colrank_mat(amat)
    ngene <- length(genes)
    amat <- apply(amat, 2, function(x) x/ngene)
  } else {
    amat <- colrank_mat(-amat)
  }
  
  return(amat)
}



# Here, ranking is done for the entire matrix jointly. Standardize controls 
# whether the final ranks should be kept in their original integer order, or 
# standardized to be [0, 1]

RSR_allrank <- function(mat,
                        meta,
                        min_cell = 20,
                        standardize = TRUE) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Get count matrix for current cell type, coercing low count genes to NAs
    
    ct_mat <- subset_and_filter(mat, meta, ct, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }

    # Get cell-type cor matrix: full symmetric, NA cors to 0, diag (self-cor)
    # coerced to 1, then to triangular to prevent double ranking symmetric matrix

    cmat <- ct_mat %>%
      get_cor_mat(lower_tri = FALSE) %>%
      na_to_zero() %>%
      diag_to_one() %>%
      upper_to_na()

    # Rank the tri matrix and add to aggregate matrix

    rmat <- allrank_mat(cmat)
    amat <- amat + rmat
    
  }
  
  if (standardize) {
    amat <- allrank_mat(amat) / sum(!is.na(amat))
  } else {
    amat <- allrank_mat(-amat)
  }
  
  return(amat)
}



# Here, correlations are transformed using Fisher' Z transformation, and these 
# are then averaged across cell types. A list of three matrices are returned: 
# the average calculated across all cell types, the average using only non-NA 
# observations, and a matrix tracking the count of NAs for each gene-gene pair 
# across cell types.
# https://bookdown.org/mwheymans/bookmi/pooling-correlation-coefficients-1.html

fishersZ_aggregate <- function(mat,
                               meta,
                               min_cell = 20) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta))
  
  cts <- unique(meta$Cell_type)
  genes <- rownames(mat)
  
  amat <- init_agg_mat(row_genes = genes)
  na_mat <- amat  # matrix of 0s for tracking NA gene-gene cor pairs
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    ct_mat <- subset_and_filter(mat, meta, ct, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      next()
    }
    
    # Generate cor matrix, and count NAs before setting to 0
    
    cmat <- get_cor_mat(ct_mat, lower_tri = FALSE)
    
    na_ix <- which(is.na(cmat), arr.ind = TRUE)
    na_mat[na_ix] <- na_mat[na_ix] + 1
    
    cmat <- cmat %>% 
      na_to_zero() %>% 
      diag_to_one()
    
    # Fisher's z transformation and add result to aggregate matrix
    
    zmat <- DescTools::FisherZ(cmat)
    amat <- amat + zmat
    
  }
  
  agg_list <- list(
    Avg_all = (amat / length(cts)),
    Avg_nonNA = (amat / (length(cts) - na_mat)),
    NA_mat = na_mat)
  
  return(agg_list)
}



# Get correlation of gene1 and gene2 across cell types represented in meta

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



# Misc helpers
# ------------------------------------------------------------------------------


# Convert a matrix into a long and skinny df. If symmetric, only return the
# unique values.

mat_to_df <- function(mat, symmetric = TRUE) {
  
  if (symmetric) {
    df <- data.frame(
      Row = rownames(mat)[row(mat)[lower.tri(mat)]],
      Col = colnames(mat)[col(mat)[lower.tri(mat)]],
      Value = mat[lower.tri(mat)],
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Row = rownames(mat)[row(mat)],
      Col = colnames(mat)[col(mat)],
      Value = c(mat),
      stringsAsFactors = FALSE
    )
  }
  return(df)
}


# Return a df of the max cor for each cell type for the given gene, removing NAs
# and the gene itself to prevent cor=1

max_cor_df <- function(cmat_list, tf) {
  
  cor_max <- lapply(names(cmat_list), function(x) {
    
    cor_mat <- cmat_list[[x]]
    
    vec <- cor_mat[tf, setdiff(colnames(cor_mat), tf)]
    
    if (all(is.na(vec))) {
      return(NA)
    }
    
    data.frame(
      Cell_type = x,
      Value = max(vec, na.rm = TRUE),
      Symbol = names(vec)[which.max(vec)])
  })
  
  cor_max <- cor_max[!is.na(cor_max)]
  
  cor_max <- data.frame(do.call(rbind, cor_max)) %>% 
    arrange(desc(Value))
  
  return(cor_max)
}


# Working with Seurat object
# TODO: sample_expression_level() impl requires subset (costly) - try pre-format input mat
# ------------------------------------------------------------------------------


# TODO:

top_expr_quantile <- function(sdat, gene, qtl = 0.9) {
  
  expr_cutoff <- quantile(sdat@assays$RNA@data[gene, ], qtl)
  expr_cutoff <- ifelse(expr_cutoff == 0L, 1, expr_cutoff)
  sdat$Top_expr_quantile <- sdat@assays$RNA@data[gene, ] >= expr_cutoff
  
  return(sdat)
}



# TODO:

top_expr_celltype <- function(sdat, avg_mat, gene) {
  
  stopifnot("Cell_type" %in% colnames(sdat@meta.data))
  
  top_ct <- names(which.max(avg_mat[gene , ]))
  sdat$Top_expr_celltype <- sdat$Cell_type == top_ct
  
  return(sdat)
}



# TODO: doc and use slot get method

get_ct_avg <- function(sdat, assay = "RNA", scale = FALSE, ncores = 8) {
  
  stopifnot("Cell_type" %in% colnames(sdat@meta.data))
  
  cts <- unique(sdat$Cell_type)
  
  if (scale) {
    mat <- sdat@assays[[assay]]@scale.data
  } else {
    mat <- sdat@assays[[assay]]@data
  }
  
  ct_avg <- mclapply(cts, function(x) {
    rowMeans(mat[, sdat$Cell_type == x])
  }, mc.cores = ncores)
  
  ct_avg <- do.call(cbind, ct_avg)
  rownames(ct_avg) <- rownames(mat)
  colnames(ct_avg) <- cts
  
  return(ct_avg)
}



# TODO: doc and use slot get method


sample_expression_level <- function(sdat, targets, rank_window = 200) {
  
  stopifnot(all(targets %in% rownames(sdat)))
  
  # get the ordered average expression of genes across all cells
  avg_all <- sort(rowMeans(sdat@assays$RNA@data))
  
  # for each target sample a gene whose average expression is within rank window
  
  sample_ix <- vapply(targets, function(x) {
    
    ix_orig <- which(names(avg_all) == x)
    ix_new <- ix_orig
    
    while (ix_new == ix_orig) {
      ix_new <- ix_orig + sample(-rank_window:rank_window, 1) 
    }
    
    return(ix_new)
    
  }, FUN.VALUE = integer(1), USE.NAMES = FALSE)
  
  return(names(avg_all[sample_ix]))
}



# Performance of rankings in a data frame. 
# ------------------------------------------------------------------------------


# TODO:

get_perf_df <- function(rank_df, label_col, measure = NULL) {
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("ROC", "PR"))
  
  # positives/has evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(rank_df[[label_col]]))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
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


# TODO:

get_au_perf <- function(rank_df, label_col, measure = NULL) {
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("AUROC", "AUPRC"))
  
  # positives/has evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(rank_df[[label_col]]))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  
  if (measure == "AUROC") {
    perf <- ROCR::performance(pred, measure = "auc")@y.values[[1]]
  } else {
    perf <- ROCR::performance(pred, measure = "aucpr")@y.values[[1]]
  }
  
  return(perf)
}



# TODO:

all_perf_df <- function(rank_df,
                        keep_cols,
                        label_col,
                        measure = NULL) {
  
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("ROC", "PR"))
 
  perf_l <- lapply(keep_cols, function(x) {
    
    get_perf_df(
      rank_df = dplyr::arrange(rank_df, !!sym(x)),
      label_col,
      measure) %>%
      mutate(Group = x)
  })
  
  perf_df <- do.call(rbind, perf_l) %>% 
    mutate(Group = factor(Group, levels = keep_cols))
   
  return(perf_df)
}


# TODO:

all_au_perf <- function(rank_df,
                        keep_cols,
                        label_col,
                        measure = NULL) {
  
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("AUROC", "AUPRC"))
  
  perf_l <- lapply(keep_cols, function(x) {
    
    get_au_perf(
      rank_df = dplyr::arrange(rank_df, !!sym(x)),
      label_col,
      measure)
  })
  
  names(perf_l) <- keep_cols
  
  return(perf_l)
}




# Functions for QC/preprocessing count matrices and metadata
# ------------------------------------------------------------------------------



# Gives a matrix with ENSEMBL IDs as rownames, return the matrix with the 
# corresponding gene symbols as rownames. Blank gene symbols are removed

ensembl_to_symbol <- function(mat, ensembl_df) {
  
  stopifnot(c("Gene_ID", "Symbol") %in% colnames(ensembl_df))
  
  ids <- intersect(pc$Gene_ID, rownames(mat))
  
  if (length(ids) == 0) stop("No common ENSEMBL IDs in rownames of matrix")
  
  common_genes <- data.frame(
    ID = ids,
    Symbol = pc$Symbol[match(ids, pc$Gene_ID)]) %>% 
    filter(Symbol != "")
  
  # dupl_genes <- common_genes$Symbol[which(duplicated(common_genes$Symbol))]
  
  mat <- mat[common_genes$ID, ]
  rownames(mat) <- common_genes$Symbol
  
  return(mat)
}



# Assumes that mat is sparse gene x cell count matrix. Filters the matrix for 
# unique gene symbols in pcoding_df, and fills the missing genes as 0s. 

get_pcoding_only <- function(mat, pcoding_df) {
  
  stopifnot("Symbol" %in% colnames(pcoding_df))
  
  genes <- unique(pcoding_df$Symbol)
  common <- intersect(rownames(mat), genes)
  missing <- setdiff(genes, rownames(mat))
  
  if (length(common) == 0) stop("No common symbols in rownames of mat")
  
  pc_mat <- mat[common, ]
  pc_mat <- rbind(pc_mat, Matrix(0, nrow = length(missing), ncol = ncol(mat)))
  rownames(pc_mat) <- c(common, missing)
  pc_mat <- pc_mat[genes, ]
  
  return(pc_mat)
}



# Given a vector of genes that have either common gene symbols or ensembl
# ids, return a subst of gene_vec only containing the mitochondrial genes. 
# Assumes gene_vec has only mouse or human symbols/ensembl IDs.

get_mt_genes <- function(gene_vec,
                         mt_path = "/home/amorin/Data/Metadata/mitochondrial_genes_all.tsv") {
  
  mt_table <- read.delim(mt_path, stringsAsFactors = FALSE)
  mt_genes <- gene_vec[gene_vec %in% c(mt_table$Gene_stable_ID, mt_table$Gene_name)]
  
  return(mt_genes)
}



# This adds columns to metadata: the number of total UMI counts for each 
# cell/column of mat, the number of non-zero expressing genes, and the RNA
# novelty/compexity, which is the ratio of the log10 gene counts to log10 umi 
# counts. It additionally adds ratio of mitochondrial if available. 
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

add_count_info <- function(mat, meta) {
  
  mt_genes <- get_mt_genes(rownames(mat))
  
  meta <- meta %>% 
    mutate(
      UMI_counts = colSums(mat),
      Gene_counts = colSums(mat > 0),
      RNA_novelty = log10(Gene_counts) / log10(UMI_counts)
    )
  
  if (length(mt_genes) > 0) {
    mt_ratio <- colSums(mat[mt_genes, , drop = FALSE]) / meta$UMI_counts
    meta$MT_ratio = mt_ratio
  }
  
  # Remove cell x gene features of this type, if present
  meta[, c("nFeature_RNA", "nCount_RNA")] <- NULL
  
  return(meta)
}



# This removes from mat cells that fail any of the filters as laid out in 
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

rm_low_qc_cells <- function(mat, 
                            meta,
                            min_counts = 500,
                            min_genes = 250,
                            min_novelty = 0.8,
                            max_mt_ratio = 0.2) {
  
  keep_id <- meta %>% 
    filter(
      UMI_counts >= min_counts,
      Gene_counts >= min_genes,
      RNA_novelty > min_novelty) %>% 
    pull(ID)
  
  if ("MT_ratio" %in% colnames(meta)) {
    keep_id <- intersect(keep_id, filter(meta, MT_ratio < max_mt_ratio)$ID)
  }
  
  if (length(keep_id) == 0) stop("No remaining cells after filtering")
  
  return(mat[, keep_id])
}



# Assumes that mat is a gene x cell count matrix. Filters the matrix of genes 
# that are not detected (nonzero) in at least min_count cells

filter_low_expressed_genes <- function(mat, min_count = 20) {
  
  mat <- mat > 0
  counts <- rowSums(mat)
  keep <- which(counts >= min_count)
  
  return(mat[keep, ])
}




# Old/unused but kept for reference
# ------------------------------------------------------------------------------


# Count NAs for each element of a list of matrices
# No longer used as correlation is generated iteratively, not in a full list

# count_nas <- function(cmat_list) {
#   mat <- apply(simplify2array(cmat_list), 1:2, function(x) sum(is.na(x)))
#   return(mat)
# }


# Set NAs in matrix to the average value of the matrix
# Sticking with setting NA cors to 0 - extremely similar result

# na_to_mean <- function(mat) {
#   mat[is.na(mat)] <- mean(mat, na.rm = TRUE)
#   return(mat)
# }

