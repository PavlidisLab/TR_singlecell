# Generating aggregate correlation for a specific cell type
# ------------------------------------------------------------------------------


# TODO: note that this allows multiple matches; check before merge

subset_and_filter <- function(mat, meta, cell_type, min_count = 20) {
  
  ids <- dplyr::filter(meta, Cell_type %in% cell_type)$ID
  ct_mat <- t(mat[, ids])
  ct_mat <- under_min_count_to_na(ct_mat, min_count)
  stopifnot(all(rownames(ct_mat) %in% meta$ID))
  
  return(ct_mat)
}



# TODO: check before merge

init_agg_mat <- function(pc_df) {
  
  amat <- matrix(0, nrow = nrow(pc_df), ncol = nrow(pc_df))
  rownames(amat) <- colnames(amat) <- pc_df$Symbol
  return(amat)
}



# TODO:

aggr_celltype_coexpr <- function(ct_df, 
                                 pc_df, 
                                 cor_method = "pearson",
                                 min_cell = 20,
                                 ncores = 1,
                                 in_dir,
                                 out_dir) {
  
  ids <- unique(ct_df[["ID"]])
  
  # Matrices of 0s for tracking aggregate correlation and count of NAs
  
  amat <- init_agg_mat(pc_df)
  na_mat <- amat  
  
  # Loop through every dataset ID, generating ranked cor and adding to aggregate
  
  for (id in ids) {
    
    message(paste(id, Sys.time()))
    cmat_path <- file.path(out_dir, paste0(id, "_cormat.tsv")) 
    
    # Load or generate raw correlation matrix
    
    if (!file.exists(cmat_path)) {
      
      # Load dataset
      
      processed_path <- file.path(in_dir, id, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
      dat <- readRDS(processed_path)
      meta <- dat$Meta
      mat <- dat$Mat
      stopifnot(identical(colnames(mat), meta$ID))
      mat <- as.matrix(mat)
      
      # Subset to given cell type and set low count genes to NA
      
      ct <- filter(ct_df, ID == id)[["Cell_type"]]
      ct_mat <- subset_and_filter(mat, meta, ct, min_cell) 
      
      if (sum(is.na(ct_mat)) == length(ct_mat)) {
        message(paste(ct, "skipped due to insufficient counts"))
        na_mat <- na_mat + 1
        next()
      }
      
      # Generate cor matrix and save
      
      cmat <- get_cor_mat(ct_mat,
                          cor_method = cor_method,
                          lower_tri = FALSE, 
                          ncores = ncores)
      
      fwrite(
        data.frame(cmat, check.names = FALSE),
        sep = "\t",
        row.names = TRUE,
        quote = FALSE,
        verbose = FALSE,
        showProgress = FALSE,
        file = cmat_path
      )
      
      rm(ct_mat)
      
    } else {
      
      cmat <- fread_to_mat(path = cmat_path, genes = pc_df$Symbol)
      
    }
    
    # Track NAs from non-measured gene pairs
    
    na_ix <- which(is.na(cmat), arr.ind = TRUE)
    na_mat[na_ix] <- na_mat[na_ix] + 1
    
    # Formatting cor matrix: NAs cors to 0 (rank between pos/neg), ensure
    # self cor is 1, and make triangular to prevent double ranking symmetry
    
    cmat <- cmat %>%
      na_to_zero() %>%
      diag_to_one() %>%
      upper_to_na()
    
    # Rank the whole tri. matrix jointly and add to aggregate matrix
    
    rmat <- allrank_mat(cmat, ties_arg = "min")
    amat <- amat + rmat
    rm(cmat, rmat)
    gc(verbose = FALSE)
    
  }
  
  # Standardize rank of aggregate and convert tri. matrix back to symmetric
  
  amat <- allrank_mat(amat, ties_arg = "min") / sum(!is.na(amat))
  amat <- lowertri_to_symm(amat)
  
  return(list(Agg_mat = amat, NA_mat = na_mat))
}





aggr_celltype_coexpr_colrank <- function(ct_df,
                                         pc_df,
                                         cor_method = "pearson",
                                         min_cell = 20,
                                         ncores = 1,
                                         in_dir,
                                         out_dir) {
  
  
  ids <- unique(ct_df[["ID"]])
  
  # Matrices of 0s for tracking aggregate correlation and count of NAs
  
  amat <- init_agg_mat(pc_df)
  na_mat <- amat  
  
  # Loop through every dataset ID, generating ranked cor and adding to aggregate
  
  for (id in ids) {
    
    message(paste(id, Sys.time()))
    cmat_path <- file.path(out_dir, paste0(id, "_cormat.tsv")) 
    
    # Load or generate raw correlation matrix
    
    if (!file.exists(cmat_path)) {
      
      # Load dataset
      
      processed_path <- file.path(in_dir, id, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
      dat <- readRDS(processed_path)
      meta <- dat$Meta
      mat <- dat$Mat
      stopifnot(identical(colnames(mat), meta$ID))
      mat <- as.matrix(mat)
      
      # Subset to given cell type and set low count genes to NA
      
      ct <- filter(ct_df, ID == id)[["Cell_type"]]
      ct_mat <- subset_and_filter(mat, meta, ct, min_cell) 
      
      if (sum(is.na(ct_mat)) == length(ct_mat)) {
        message(paste(ct, "skipped due to insufficient counts"))
        na_mat <- na_mat + 1
        next()
      }
      
      # Generate cor matrix and save
      
      cmat <- get_cor_mat(ct_mat,
                          cor_method = cor_method,
                          lower_tri = FALSE, 
                          ncores = ncores)
      
      fwrite(
        data.frame(cmat, check.names = FALSE),
        sep = "\t",
        row.names = TRUE,
        quote = FALSE,
        verbose = FALSE,
        showProgress = FALSE,
        file = cmat_path
      )
      
      rm(ct_mat)
      
    } else {
      
      cmat <- fread_to_mat(path = cmat_path, genes = pc_df$Symbol)
      
    }
    
    # Track NAs from non-measured gene pairs
    
    na_ix <- which(is.na(cmat), arr.ind = TRUE)
    na_mat[na_ix] <- na_mat[na_ix] + 1
    
    # Formatting cor matrix: NAs cors to 0 (rank between pos/neg), ensure
    # self cor is 1, and make triangular to prevent double ranking symmetry
    
    cmat <- cmat %>%
      na_to_zero() %>%
      diag_to_one()
    
    # Rank the whole tri. matrix jointly and add to aggregate matrix
    
    rmat <- colrank_mat(cmat, ties_arg = "min")
    amat <- amat + rmat
    rm(cmat, rmat)
    gc(verbose = FALSE)
    
  }
  
  # Standardize rank of aggregate and convert tri. matrix back to symmetric
  
  amat <- colrank_mat(amat, ties_arg = "min")
  ngene <- nrow(amat)
  amat <- apply(amat, 2, function(x) x/ngene)

  return(list(Agg_mat = amat, NA_mat = na_mat))
}




aggr_celltype_coexpr_fishersZ <- function(ct_df,
                                          pc_df,
                                          cor_method = "pearson",
                                          min_cell = 20,
                                          ncores = 1,
                                          in_dir,
                                          out_dir) {
  
  
  ids <- unique(ct_df[["ID"]])
  
  # Matrices of 0s for tracking aggregate correlation and count of NAs
  
  amat <- init_agg_mat(pc_df)
  na_mat <- amat  
  
  # Loop through every dataset ID, generating ranked cor and adding to aggregate
  
  for (id in ids) {
    
    message(paste(id, Sys.time()))
    cmat_path <- file.path(out_dir, paste0(id, "_cormat.tsv")) 
    
    # Load or generate raw correlation matrix
    
    if (!file.exists(cmat_path)) {
      
      # Load dataset
      
      processed_path <- file.path(in_dir, id, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
      dat <- readRDS(processed_path)
      meta <- dat$Meta
      mat <- dat$Mat
      stopifnot(identical(colnames(mat), meta$ID))
      mat <- as.matrix(mat)
      
      # Subset to given cell type and set low count genes to NA
      
      ct <- filter(ct_df, ID == id)[["Cell_type"]]
      ct_mat <- subset_and_filter(mat, meta, ct, min_cell) 
      
      if (sum(is.na(ct_mat)) == length(ct_mat)) {
        message(paste(ct, "skipped due to insufficient counts"))
        na_mat <- na_mat + 1
        next()
      }
      
      # Generate cor matrix and save
      
      cmat <- get_cor_mat(ct_mat,
                          cor_method = cor_method,
                          lower_tri = FALSE, 
                          ncores = ncores)
      
      fwrite(
        data.frame(cmat, check.names = FALSE),
        sep = "\t",
        row.names = TRUE,
        quote = FALSE,
        verbose = FALSE,
        showProgress = FALSE,
        file = cmat_path
      )
      
      rm(ct_mat)
      
    } else {
      
      cmat <- fread_to_mat(path = cmat_path, genes = pc_df$Symbol)
      
    }
    
    # Track NAs from non-measured gene pairs
    
    na_ix <- which(is.na(cmat), arr.ind = TRUE)
    na_mat[na_ix] <- na_mat[na_ix] + 1
    
    # Formatting cor matrix: NAs cors to 0 (rank between pos/neg), ensure
    # self cor is 1, and make triangular to prevent double ranking symmetry
    
    cmat <- cmat %>%
      na_to_zero() %>%
      diag_to_one()
    
    # Fisher's z transformation and add result to aggregate matrix
    
    zmat <- DescTools::FisherZ(cmat)
    amat <- amat + zmat
    rm(cmat, zmat)
    gc(verbose = FALSE)
    
  }
  
  agg_list <- list(
    Avg_all = (amat / length(ids)),
    Avg_nonNA = (amat / (length(ids) - na_mat)),
    NA_mat = na_mat)
  
  return(agg_list)
}


