## Process count matrix and get aggregate correlation for GSE180928
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE180928"
sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata", id)
dat_path <- file.path(sc_dir, paste0(id, "_filtered_cell_counts.csv"))
meta_path <- file.path(sc_dir, paste0(id, "_metadata.csv"))
out_dir <- file.path("/space/scratch/amorin/TR_singlecell", id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank_minrank_nostd.RDS"))
colrank_path <- file.path(out_dir, paste0(id, "_RSR_colrank_minrank_nostd.RDS"))
zcor_path <- file.path(out_dir, paste0(id, "_fishersZ.RDS"))
pc <- read.delim(ref_hg_path, stringsAsFactors = FALSE)



RSR_allrank <- function(mat,
                        meta,
                        min_cell = 20,
                        standardize = FALSE) {
  
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
    
    rmat <- allrank_mat(cmat, ties_arg = "min")
    amat <- amat + rmat
    rm(cmat, ct_mat)
    gc(verbose = FALSE)
    
  }
  
  if (standardize) {
    amat <- allrank_mat(amat, ties_arg = "min") / sum(!is.na(amat))
  } else {
    amat <- allrank_mat(amat, ties_arg = "min")
  }
  
  return(amat)
}




RSR_colrank <- function(mat,
                        meta,
                        min_cell = 20,
                        standardize = FALSE) {
  
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
    
    rmat <- colrank_mat(cmat, ties_arg = "min")
    amat <- amat + rmat
    rm(cmat, ct_mat)
    gc(verbose = FALSE)
    
  }
  
  if (standardize) {
    amat <- colrank_mat(amat, ties_arg = "min")
    ngene <- length(genes)
    amat <- apply(amat, 2, function(x) x/ngene)
  } else {
    amat <- colrank_mat(amat, ties_arg = "min")
  }
  
  return(amat)
}




if (!file.exists(processed_path)) {
  
  dat <- read.delim(dat_path, sep = ",")
  meta <- read.delim(meta_path, sep = ",")
  
  # Preparing count matrix
  
  rownames(dat) <- dat$X
  dat$X <- NULL
  colnames(dat) <- str_replace_all(colnames(dat), "\\.", "_")
  mat <- Matrix(as.matrix(dat), sparse = TRUE)
  
  # Ready metadata
  
  meta <- meta %>% 
    dplyr::rename(ID = X, Cell_type = Lineage) %>% 
    mutate(ID = str_replace_all(ID, "-", "_")) %>% 
    add_count_info(mat = mat)
  
  # QC plots
  
  p1 <- all_hist(meta)
  p2 <- qc_scatter(meta)
  
  ggsave(p1, device = "png", dpi = 300, height = 12, width = 16, bg = "white",
         filename = file.path(out_dir, paste0(id, "_QC_histograms.png")))
  
  ggsave(p2, device = "png", dpi = 300, height = 8, width = 8,
         filename = file.path(out_dir, paste0(id, "_QC_scatter.png")))
  
  # Remove cells failing QC, keep only protein coding genes, and normalize
  
  mat <- rm_low_qc_cells(mat, meta) %>%
    get_pcoding_only(pcoding_df = pc) %>% 
    Seurat::LogNormalize(., verbose = FALSE)
    
  meta <- filter(meta, ID %in% colnames(mat))
  
  stopifnot(all(colnames(mat) %in% meta$ID), length(meta$ID) > 0)
  
  saveRDS(list(Mat = mat, Meta = meta), file = processed_path)
  
  rm(dat)
  gc()
  
} else {
  
  dat <- readRDS(processed_path)
  meta <- dat[[2]]
  mat <- dat[[1]]
  
}


stopifnot(identical(colnames(mat), meta$ID))


mat <- as.matrix(mat)


if (!file.exists(allrank_path)) {
  rsr_all <- RSR_allrank(mat, meta)
  saveRDS(rsr_all, allrank_path)
}


if (!file.exists(colrank_path)) {
  rsr_col <- RSR_colrank(mat, meta)
  saveRDS(rsr_col, colrank_path)
}
