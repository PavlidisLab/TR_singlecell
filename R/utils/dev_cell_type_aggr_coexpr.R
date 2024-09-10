source("R/utils/dev_functions.R")


# TODO: check before merge

init_agg_mat <- function(pc_df) {
  
  amat <- matrix(0, nrow = nrow(pc_df), ncol = nrow(pc_df))
  rownames(amat) <- colnames(amat) <- pc_df$Symbol
  return(amat)
}



# Assumes path leads to an RDS of a list with 2 elements: the sparse count
# matrix and the metadata mapping cell IDs to cell type

load_dataset <- function(path) {

  dat <- readRDS(path)
  meta <- dat$Meta
  mat <- dat$Mat
  stopifnot(identical(colnames(mat), meta$ID))
  
  return(list(Mat = mat, Meta = meta))
}



# TODO:

aggr_coexpr_across_datasets <- function(ct_df,
                                        pc_df,
                                        cor_method = "pearson",
                                        agg_method = "FZ",
                                        min_cell = 20,
                                        verbose = TRUE) {
  
  stopifnot(cor_method %in% c("pearson", "spearman"))
  stopifnot(agg_method %in% c("allrank", "colrank", "FZ"))
  
  data_ids <- unique(ct_df[["ID"]])
  n_cts <- length(data_ids) # All cell types for a dataset are collapsed
  
  # Matrices of 0s for tracking aggregate correlation and count of NAs
  
  amat <- init_agg_mat(pc_df)
  na_mat <- amat  
  
  for (id in data_ids) {
    
    if (verbose) message(paste(id, Sys.time()))
    
    ct <- filter(ct_df, ID == id)[["Cell_type"]]
    dat_path <- unique(filter(ct_df, ID == id)[["Path"]])
    
    # Load dataset and get count matrix for current cell type
      
    dat = load_dataset(dat_path)
    
    ct_mat <- prepare_celltype_mat(mat = dat$Mat, 
                                   meta = dat$Meta, 
                                   cell_type = ct, 
                                   min_count = min_cell)

    # Check if filtering removed all genes
    
    no_msr <- all(ct_mat == 0)
    
    if (no_msr) {
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
