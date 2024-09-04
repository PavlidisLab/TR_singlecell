## Project functions
## -----------------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(parallel)
library(WGCNA)
library(Matrix)
library(qlcMatrix)
library(DescTools)



# Not in

'%!in%' <- function(x, y) !('%in%'(x, y))



# Loading and saving data objects
# ------------------------------------------------------------------------------


# Execute a function with provided args and save an RDS to path

save_function_results <- function(path, 
                                  fun,
                                  args,
                                  force_resave = FALSE) {
  
  if (!file.exists(path) || force_resave) {
    
    result <- do.call(fun, args)
    
    if (!is.null(result)) {
      saveRDS(result, path)
    }
    
    return(invisible(NULL))
  }
}



# Wrapper around data.table::fwrite() to save a matrix as a tab delim data.frame
# so that matrix rownames are saved as the first column

fwrite_mat <- function(mat, path) {
  
  fwrite(
    data.frame(mat, check.names = FALSE),
    sep = "\t",
    row.names = TRUE,
    quote = FALSE,
    verbose = FALSE,
    showProgress = FALSE,
    file = path
  )
}



# Calls data.table::fread() to read in a table where the first column corresponds
# to gene names, and convert this table to a gene x gene matrix. 
# sub_genes controls if only a subset of column genes should be loaded 

fread_to_mat <- function(path, genes, sub_genes = NULL) {
  
  if (!is.null(sub_genes)) {
    
    dat <- fread(path, sep = "\t", select = c("V1", sub_genes))
    mat <- as.matrix(dat[, -1, drop = FALSE])
    rownames(mat) <- dat$V1
    colnames(mat) <- sub_genes
    mat <- mat[genes, sub_genes, drop = FALSE]
    
  } else {
    
    dat <- fread(path, sep = "\t")
    mat <- as.matrix(dat[, -1, drop = FALSE])
    rownames(mat) <- colnames(mat) <- dat$V1
    mat <- mat[genes, genes, drop = FALSE]
  }
  
  return(mat)
}



# Loads aggregate correlation matrices or NA count matrices for the given dataset
# ids into a list. 
# Pattern: "_NA_mat.tsv" for NA counts

load_agg_mat_list  <- function(ids,
                               dir = "/space/scratch/amorin/TR_singlecell/",
                               pattern = "_RSR_allrank_CPM.tsv",  
                               genes,
                               sub_genes = NULL,
                               verbose = TRUE) {
  
  mat_l <- lapply(ids, function(x) {
    if (verbose) message(paste(x, Sys.time()))
    path <- file.path(dir, x, paste0(x, pattern))
    fread_to_mat(path, genes, sub_genes)
  })
  
  names(mat_l) <- ids
  gc(verbose = FALSE)
  
  return(mat_l)
}



# Call load_agg_mat_list to load/save a local copy

load_or_generate_agg <- function(path, 
                                 ids, 
                                 dir = "/space/scratch/amorin/TR_singlecell/",
                                 pattern = "_RSR_allrank_CPM.tsv",
                                 genes, 
                                 sub_genes = NULL) {
  
  if (!file.exists(path)) {
    
    agg_l <- load_agg_mat_list(ids = ids, 
                               dir = dir,
                               genes = genes, 
                               pattern = pattern,
                               sub_genes = sub_genes)
    
    saveRDS(agg_l, path)
  } else {
    agg_l <- readRDS(path)
  }
  
  return(agg_l)
}



# Given a vector of ids, will load the associated list of normalized matrices
# and metadata into a list. 

load_dat_list <- function(ids,
                          sc_dir = "/space/scratch/amorin/TR_singlecell/",
                          pattern = "_clean_mat_and_meta_CPM.RDS") {
  
  dat_l <- lapply(ids, function(x) {
    path <- file.path(sc_dir, x, paste0(x, pattern))
    dat <- readRDS(path)
  })
  names(dat_l) <- ids
  
  return(dat_l)
}


# Uses fread() to read from path, assuming that the introduced V1 (rownames)
# column is for genes. Coerces to matrix and assigns gene rownames. Columns
# are assumed to be cell IDs (thus distinct from fread_to_mat())

read_count_mat <- function(dat_path) {
  
  dat <- suppressWarnings(fread(dat_path))
  genes <- dat$V1
  dat$V1 <- NULL
  mat <- as.matrix(dat)
  rownames(mat) <- genes
  
  return(mat)
}



# Single cell coexpression aggregation
# ------------------------------------------------------------------------------

# All of these assume that bigger values are more important, and that the
# rank importance is the inverse such that 1=best. Minimum ranking is used for
# ties, as this is used to infer where ties pile up in the overall ranking.


# Rank matrix columns such that 1=best. Return as same dimension as input.

colrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {
  
  rank_mat <- apply(-mat, 2, rank, ties.method = ties_arg, na.last = na_arg)
  return(rank_mat)
}



# Rank the entire matrix jointly such that 1=best. Return as same dimension as input.

allrank_mat <- function(mat, ties_arg = "min", na_arg = "keep") {
  
  rank_mat <- rank(-mat, ties.method = ties_arg, na.last = na_arg)
  rank_mat <- matrix(rank_mat, nrow = nrow(mat), ncol = ncol(mat))
  rownames(rank_mat) <- colnames(rank_mat) <- rownames(mat)
  return(rank_mat)
}



# TODO: use instead of external?
# fz <- function(r) 0.5 * log((1 + r) / (1 - r))
# identical(fz(0.5), DescTools::FisherZ(0.5))



# Column-wise Pearson's correlation for sparse matrices. 
# Wrapper for qlcMatrix::sparseCor() to format resulting matrix

sparse_pcor <- function(mat) {
  
  stopifnot(inherits(mat, "dgCMatrix"))
  
  cmat <- qlcMatrix::corSparse(mat)
  colnames(cmat) <- rownames(cmat) <- colnames(mat)
  cmat <- diag_to_one(cmat)
  
  return(cmat)
}



# Column rank matrix while preserving sparsity, which is then fed into sparse
# Pearson's cor to get a sparse Spearman's impl.
# https://saket-choudhary.me/blog/2022/03/09/sparsespearman/
# https://github.com/saketkc/blog/blob/main/2022-03-10/SparseSpearmanCorrelation2.ipynb

sparse_rank <- function(mat) {
  
  stopifnot(inherits(mat, "dgCMatrix"))
  
  nonzero_per_col <- diff(mat@p)
  zero_per_col <- nrow(mat) - nonzero_per_col
  offsets <- (zero_per_col - 1) / 2
  x <- mat@x
  
  ## split entries to columns
  col_l <- split(x, f = rep.int(1:ncol(mat), nonzero_per_col))
  
  # Rank and offset
  sparse_ranks <- lapply(seq_along(col_l), function(i) {
    rank(col_l[[i]]) + offsets[i]
  })
  
  mat@x <- unlist(sparse_ranks)
  
  return(mat)
}



# Column-wise Spearman's correlation for sparse matrices. 
# Note: This is only calculated over nonzero columns, as all zero columns affect
# the offset logic employed by sparse_rank(). All zero columns are maintained as
# NAs, with 1 on the diagonal

sparse_scor <- function(mat) {
  
  stopifnot(inherits(mat, "dgCMatrix"))
  
  cmat_all <- matrix(NA, ncol(mat), ncol(mat))
  colnames(cmat_all) <- rownames(cmat_all) <- colnames(mat)
  which_nonzero <- colSums(mat > 0) > 0
  mat <- mat[, which_nonzero]
  
  rmat <- sparse_rank(mat)
  cmat <- sparse_pcor(rmat)
  cmat_all[which_nonzero, which_nonzero] <- cmat
  cmat_all <- diag_to_one(cmat_all)
  
  return(cmat_all)
}



# Wrapper for calling sparse correlation
# cor_method must be one of [pearson, spearman]

calc_sparse_correlation <- function(mat, cor_method) {
  
  stopifnot(inherits(mat, "dgCMatrix"))
  stopifnot(cor_method %in% c("pearson", "spearman"))
  
  cmat <- if (cor_method == "pearson") {
    sparse_pcor(mat) 
  } else {
    sparse_scor(mat)
  }

  return(cmat)
}



# Set NAs in matrix to 0

na_to_zero <- function(mat) {
  
  stopifnot(is.matrix(mat))
  mat[is.na(mat)] <- 0
  return(mat)
}



# Set diag in matrix to 1

diag_to_one <- function(mat) {
  
  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  diag(mat) <- 1
  return(mat)
}



# Set upper triangle to NA

uppertri_to_na <- function(mat) {
  
  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  mat[upper.tri(mat)] <- NA
  return(mat)
}



# Replace upper tri of mat with lower tri.

lowertri_to_symm <- function(mat) {
  
  stopifnot(is.matrix(mat), ncol(mat) == nrow(mat))
  mat[upper.tri(mat)] <-  t(mat)[upper.tri(mat)]
  return(mat)
}



# If a col of mat has fewer non-zero elements than min_count, set that col to
# 0. This is done to produce an NA during correlation, instead of allowing cors 
# derived from small n. 

zero_sparse_cols <- function(mat, min_count = 20) {
  
  stopifnot(inherits(mat, "dgCMatrix"))
  stopifnot(is.numeric(min_count), min_count >= 0 & min_count <= nrow(mat))
  
  nonzero_cells <- colSums(mat != 0)
  filt_genes <- nonzero_cells < min_count
  if (any(filt_genes)) mat[, filt_genes] <- 0

  return(mat)
}



# Subset mat to cell_type, set under min count genes to 0 and transpose.
# Expects mat is genes x cells, and will return a cells x genes mat for corr.

prepare_celltype_mat <- function(mat, meta, cell_type, min_count = 20) {
  
  stopifnot(inherits(mat, "dgCMatrix"))
  stopifnot(c("ID", "Cell_type") %in% colnames(meta), cell_type %in% meta$Cell_type)
  
  ids <- dplyr::filter(meta, Cell_type %in% cell_type)[["ID"]]
  stopifnot(all(ids %in% colnames(mat)))
  
  ct_mat <- t(mat[, ids])
  ct_mat <- zero_sparse_cols(ct_mat, min_count)
  stopifnot(all(rownames(ct_mat) %in% meta$ID))
  
  return(ct_mat)
}



# Initate a matrix of 0s for holding aggregate correlations. Assumes mat is a 
# gene x cell count matrix and returns a gene x gene matrix of 0s with the names
# and dimension of the gene rows of mat

init_agg_mat <- function(mat) {
  
  stopifnot(is.matrix(mat) || inherits(mat, "dgCMatrix")) 
  stopifnot(!is.null(rownames(mat)))
  amat <- matrix(0, nrow = nrow(mat), ncol = nrow(mat))
  rownames(amat) <- colnames(amat) <- rownames(mat)
  return(amat)
}



# Count the NAs in cmat and increment the corresponding indices of na_mat

increment_na_mat <- function(cmat, na_mat) {
  
  na_ix <- which(is.na(cmat), arr.ind = TRUE)
  na_mat[na_ix] <- na_mat[na_ix] + 1
  return(na_mat)
}



# Take the input correlation matrix and transform for ranking.
# all methods set NA cors to 0, and ensure the diagonal equals 1.
# allrank: make the matrix tri. to prevent double ranking symmetric elements and
#          then jointly rank the matrix (lower rank = positive cor).  
# colrank: rank each column seperately (lower rank = positive cor).
# FZ:      perform Fisher's Z transform on the raw correlations

transform_correlation_mat <- function(cmat, agg_method) {
  
  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"))
  stopifnot(is.matrix(cmat), identical(rownames(cmat), colnames(cmat)))
  
  cmat <- cmat %>%
    na_to_zero() %>%
    diag_to_one()
  
  if (agg_method == "allrank") {
    
    cmat <- cmat %>% 
      uppertri_to_na() %>% 
      allrank_mat()
    
  } else if (agg_method == "colrank") {
    
    cmat <- colrank_mat(cmat)
    
  } else if (agg_method == "FZ") {
    
    cmat <- DescTools::FisherZ(cmat) 
  }
  
  return(cmat)
}



# Format the summed aggregate correlations.
# all_rank: jointly re-rank the summed ranks and standardize into [0, 1], then 
#           covert back to a symmetric matrix for ease of downstream operations
# col_rank: re-rank the summed ranks for each column separately and standardize
#           into [0, 1]
# FZ:       divide each element by its count of non-NA observations

finalize_agg_mat <- function(amat, agg_method, n_celltypes, na_mat) {
  
  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"))
  stopifnot(is.matrix(amat), identical(rownames(amat), colnames(amat)))
  
  if (agg_method == "allrank") {
    
    amat <- allrank_mat(amat) / sum(!is.na(amat))
    amat <- diag_to_one(amat)
    amat <- lowertri_to_symm(amat)
    
  } else if (agg_method == "colrank") {
    
    amat <- colrank_mat(amat)
    ngene <- nrow(amat)
    amat <- apply(amat, 2, function(x) x/ngene)
    
  } else if (agg_method == "FZ") {
    
    amat <- amat / (n_celltypes - na_mat)
    
  }
  
  return(amat)
}




# Main function to aggregate gene correlation across all cell types of dataset.
# Returns a list of two matrices:
# -- Agg_mat: symmetric gene x gene matrix of aggregate correlation scores
# -- NA_mat: symmetric gene x gene matrix tracking count of NA gene cor pairs

# mat: cell by gene count matrix
# meta: data frame with cell IDs and cell types
# cor_method:
# agg_method:
# min_cell: count of cells in a cell type that must have non-zero counts for a 
#           gene considered to be measured and thus viable for calculating cor
# ncores: number of cores to use

aggregate_celltype_correlation <- function(mat,
                                           meta,
                                           cor_method,
                                           agg_method,
                                           min_cell = 20,
                                           verbose = TRUE) {
  
  stopifnot(cor_method %in% c("pearson", "spearman"))
  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"))
  
  cts <- unique(meta$Cell_type)
  n_cts <- length(cts)
  
  # Matrices of 0s for tracking aggregate correlation and count of NAs
  
  amat <- init_agg_mat(mat)
  na_mat <- amat  
  
  for (ct in cts) {
    
    if (verbose) message(paste(ct, Sys.time()))
    
    # Get count matrix for current cell type, coercing low count genes to NA/0
    
    ct_mat <- prepare_celltype_mat(mat, meta, ct, min_cell)
    n_nonmsr <- max(sum(ct_mat == 0), sum(is.na(ct_mat)))
    
    if (identical(n_nonmsr, length(ct_mat))) {
      message(paste(ct, "skipped due to insufficient counts"))
      na_mat <- na_mat + 1
      next()
    }
    
    # Get cell-type cor matrix and increment count of NAs before imputing to 0
    
    cmat <- calc_sparse_correlation(ct_mat, cor_method)
    na_mat <- increment_na_mat(cmat, na_mat)
    
    # Transform raw correlation matrix, add to aggregate and clean up
    
    cmat <- transform_correlation_mat(cmat, agg_method)
    amat <- amat + cmat
    rm(cmat, ct_mat)
    gc(verbose = FALSE)
    
  }

  # Final format of aggregate matrix and return along with the tracked NA mat
  
  amat <- finalize_agg_mat(amat, agg_method, n_cts, na_mat)
  return(list(Agg_mat = amat, NA_mat = na_mat))
  
}




# Rank sum rank aggregation from Jesse Gillis's group.
# Returns a gene-gene matrix representing aggregated gene-gene correlation across
# cell types.
# https://pubmed.ncbi.nlm.nih.gov/27165153/
# https://pubmed.ncbi.nlm.nih.gov/25717192/
# https://pubmed.ncbi.nlm.nih.gov/34015329/

# A gene-gene correlation matrix is generated for each cell type group in a 
# gene x cell count matrix. Each correlation matrix is ranked (1=best), summed,
# and then re-ranked (1=best). Standardized ranks results in gene pair 
# correlations scored as [0,1] (1=best). Genes with counts under the min_cell 
# threshold for a given cell type are set to NA, to produce NA cors instead of
# real but small n values. These NAs are then imputed to 0 for ranking. This
# introduces many ties, and there will be variable gene coverage across cell
# types and datasets. For this reason, ties are set to min, as randomly breaking
# ties varies from run to run. Min ranking also results in "clumping" of ranks
# which allows inference between real positive correlations, imputed 0s, and
# real negative correlations. However, the count of imputed gene pairs are also
# directly tracked in a separate matrix.



# Here, ranking is done for the entire matrix jointly. This returns a list of
# two matrices: the aggregate correlation matrix and the NA counts matrix





# Here, ranking is done column-wise such that each column/gene vector is ranked
# only within that vector.  This does not track NA counts (assumed to already be
# done with allrank) and returns the aggregate matrix directly. 




# Here, correlations are transformed using Fisher' Z transformation, and these 
# are then averaged across cell types. A list of two matrices are returned: 
# the average using only non-NA observations, and a matrix tracking the count of NAs for each gene-gene pair 
# across cell types.
# https://bookdown.org/mwheymans/bookmi/pooling-correlation-coefficients-1.html
# https://rdrr.io/cran/DescTools/man/FisherZ.html
