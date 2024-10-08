## Project functions
## -----------------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(data.table, quietly = TRUE)
library(parallel)
library(WGCNA)
library(Matrix)



# Not in

'%!in%' <- function(x, y) !('%in%'(x, y))



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



# Loading data objects
# ------------------------------------------------------------------------------


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

# NOTE: All of these assume that bigger values are more important, and that the
# rank importance is the inverse such that 1=best. Seed is set for default
# random ties.


# Rank matrix columns such that 1=best. Return as same dimension as input.

colrank_mat <- function(mat, ties_arg = "random", na_arg = "keep", seed = 3) {
  
  set.seed(seed)
  rank_mat <- apply(-mat, 2, rank, ties.method = ties_arg, na.last = na_arg)
  return(rank_mat)
}



# Rank matrix rows such that 1=best. Return as same dimension as input.

rowrank_mat <- function(mat, ties_arg = "random", na_arg = "keep", seed = 3) {
  
  set.seed(seed)
  rank_mat <- apply(-mat, 1, rank, ties.method = ties_arg, na.last = na_arg)
  return(t(rank_mat))
}



# Rank the entire matrix jointly such that 1=best. Return as same dimension as input.

allrank_mat <- function(mat, ties_arg = "random", na_arg = "keep", seed = 3) {
  
  set.seed(seed)
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


  } else if (cor_method == "spearman") {  # pairwise prohib. slow for scor
    cmat <- WGCNA::cor(
      mat, method = "spearman", use = "everything", nThreads = ncores)

  } else {
    cmat <- WGCNA::cor(
      mat, method = cor_method, use = "everything", nThreads = ncores)
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
  
  binary_counts <- mat > 0
  nonzero_cells <- colSums(binary_counts)
  filt_genes <- nonzero_cells < min_count
  mat[, filt_genes] <- NA
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



# Assumes mat is a gene x cell count matrix and returns a gene x gene matrix of 
# 0s with the names and dimension of the gene rows of mat

init_agg_mat <- function(mat) {
  
  stopifnot(length(rownames(mat)) > 1)
  amat <- matrix(0, nrow = nrow(mat), ncol = nrow(mat))
  rownames(amat) <- colnames(amat) <- rownames(mat)
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



# Set diag in matrix to NA

diag_to_na <- function(mat) {
  diag(mat) <- NA
  return(mat)
}



# Set upper triangle to NA

upper_to_na <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
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

RSR_allrank <- function(mat,
                        meta,
                        cor_method = "pearson",
                        min_cell = 20,
                        standardize = TRUE,
                        ncores = 1) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta), length(rownames(mat)) > 0)
  
  cts <- unique(meta$Cell_type)

  # Matrices of 0s for tracking aggregate correlation and count of NAs
  
  amat <- init_agg_mat(mat)
  na_mat <- amat  
  
  for (ct in cts) {
    
    message(paste(ct, Sys.time()))
    
    # Get count matrix for current cell type, coercing low count genes to NAs
    
    ct_mat <- subset_and_filter(mat, meta, ct, min_cell)
    
    if (sum(is.na(ct_mat)) == length(ct_mat)) {
      message(paste(ct, "skipped due to insufficient counts"))
      na_mat <- na_mat + 1
      next()
    }
    
    # Get cell-type cor matrix and increment count of NAs before imputing to 0,
    # set diag (self-cor) to 1, and set to triangular to prevent symmetric ranking
    
    cmat <- get_cor_mat(ct_mat, 
                        lower_tri = FALSE, 
                        cor_method = cor_method,
                        ncores = ncores)
    
    na_ix <- which(is.na(cmat), arr.ind = TRUE)
    na_mat[na_ix] <- na_mat[na_ix] + 1
    
    cmat <- cmat %>%
      na_to_zero() %>%
      diag_to_one() %>%
      upper_to_na()
    
    # Rank the tri matrix and add to aggregate matrix
    
    rmat <- allrank_mat(cmat, ties_arg = "min")
    amat <- amat + rmat
    rm(cmat, ct_mat, rmat)
    gc(verbose = FALSE)
    
  }
  
  # Final rank of aggregate and convert triangular matrix back to symmetric
  
  if (standardize) {
    amat <- allrank_mat(amat, ties_arg = "min") / sum(!is.na(amat))
  } else {
    amat <- allrank_mat(amat, ties_arg = "min")
  }
  
  amat <- lowertri_to_symm(amat)
  
  return(list(Agg_mat = amat, NA_mat = na_mat))
}



# Get correlation of gene1 and gene2 across cell types represented in meta

all_celltype_cor <- function(mat,
                             meta,
                             gene1,
                             gene2,
                             min_cell = 20,
                             cor_method = "pearson") {
  
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta),
            c(gene1, gene2) %in% rownames(mat))
  
  cts <- unique(meta$Cell_type)
  
  cor_l <- lapply(cts, function(ct) {
    
    ct_mat <- t(mat[c(gene1, gene2), filter(meta, Cell_type == ct)$ID])
    
    if (sum(ct_mat[, 1] != 0) < min_cell || sum(ct_mat[, 2] != 0) < min_cell) {
      return(NA)
    }
    
    WGCNA::cor(ct_mat[, gene1], ct_mat[, gene2], method = cor_method)
    
  })
  names(cor_l) <- cts
  
  cor_vec <- round(sort(unlist(cor_l), decreasing = TRUE), 5)
  
  return(cor_vec)
}



# Get a list of correlations of gene1 and gene2 across cell types for every
# experiment in ids.

get_all_cor_l <- function(ids, gene1, gene2, cor_method = "pearson") {
  
  cor_l <- lapply(ids, function(x) {
    dat <- load_dat_list(x)[[1]]
    mat <- dat$Mat
    meta <- dat$Meta
    all_celltype_cor(mat, meta, gene1, gene2, cor_method = cor_method)
  })
  
  names(cor_l) <- ids
  keep <- vapply(cor_l, function(x) length(x) > 0, logical(1))
  cor_l <- cor_l[keep]
  
  return(cor_l)
}



# Misc helpers
# ------------------------------------------------------------------------------


# Convert a matrix into a long and skinny df. If symmetric, only return the
# unique values.

mat_to_df <- function(mat, symmetric = TRUE, value_name = NULL) {
  
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
  
  if (!is.null(value_name)) colnames(df)[colnames(df) == "Value"] <- value_name
  
  return(df)
}



# Mat assumed to be a gene by experiment matrix of the aggregate profiles for 
# the given gene across experiments, and msr_mat a binary gene by experiment 
# matrix tracking whether or not genes were measured in a given dataset. Subset
# the input mat to only include the experiments/columns in which the given gene
# is measured.

subset_to_measured <- function(mat, gene, msr_mat) {
  
  stopifnot(gene %in% rownames(msr_mat), 
            length(intersect(rownames(mat), rownames(msr_mat))) > 0,
            length(intersect(colnames(mat), colnames(msr_mat))) > 0)
  
  msr_exps <- names(which(msr_mat[gene, ] == 1))
  mat <- mat[, msr_exps, drop = FALSE]
  
  return(mat)
}



# For the given gene, bind its coexpression profiles for all datasets in agg_l
# into a matrix, keeping only the profiles that are measured.

gene_vec_to_mat <- function(agg_l, gene, msr_mat) {
  
  genes <- rownames(agg_l[[1]])
  stopifnot(gene %in% genes)
  
  gene_mat <- do.call(cbind, lapply(agg_l, function(x) x[genes, gene]))
  gene_mat <- subset_to_measured(gene_mat, msr_mat = msr_mat, gene = gene)
  
  return(gene_mat)
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
  meta <- meta[, !(colnames(meta) %in% c("nFeature_RNA", "nCount_RNA"))]
  
  return(meta)
}



# This subsets mat to remove cells that fail any of the filters laid out in: 
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html
# The RNA novelty filter is relaxed for Smart-seq runs as more reads can go to
# genes like mitochondrial https://pubmed.ncbi.nlm.nih.gov/33662621/

rm_low_qc_cells <- function(mat, 
                            meta,
                            min_counts = 500,
                            min_genes = 250,
                            min_novelty = NULL,
                            max_mt_ratio = 0.2) {
  
  keep_id <- meta %>%
    mutate(Is_smartseq = str_detect(str_to_lower(assay), "smart-seq")) %>%
    filter(
      UMI_counts >= min_counts,
      Gene_counts >= min_genes,
      ifelse(Is_smartseq, RNA_novelty > 0.5, RNA_novelty > 0.8)) %>%
    pull(ID)
  
  if ("MT_ratio" %in% colnames(meta)) {
    keep_id <- intersect(keep_id, filter(meta, MT_ratio < max_mt_ratio)$ID)
  }
  
  if (length(keep_id) == 0) stop("No remaining cells after filtering")
  
  return(mat[, keep_id])
}
