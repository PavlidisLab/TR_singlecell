## Generate cell type correlation into a list for inspection of intermediates of
## aggregation.
## -----------------------------------------------------------------------------

source("R/00_config.R")
source("R/utils/functions.R")

id <- "GSE180928"
species <- "Human" 

sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets", id)
dat_path <- file.path(sc_dir, paste0(id, "_cellxgene_seurat.RDS"))
out_dir <- file.path("/space/scratch/amorin/TR_singlecell", id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank.RDS"))
corlist_path <- file.path(out_dir, paste0(id, "_celltype_corlist.RDS"))


dat <- readRDS(processed_path)
meta <- dat[[2]]
mat <- as.matrix(dat[[1]])
rsr <- readRDS(allrank_path)

genes <- rownames(mat)


get_cor_list <- function(mat, meta, ncores = 1) {
  
  cts <- unique(meta$Cell_type)
  
  # Split mat into list of cell type matrices, setting low count genes to NA,
  # calculating gene-cor cor for each cell type, and setting NA cors to 0
  
  ct_l <- mclapply(cts, function(x) {
    
    ct_mat <- subset_and_filter(mat = mat, 
                                meta = meta, 
                                cell_type = x, 
                                min_count = 20)
    
    cor_mat <- ct_mat %>%
      get_cor_mat(lower_tri = FALSE) %>%
      na_to_zero() %>%
      diag_to_one()
    
    return(cor_mat)
    
  }, mc.cores = ncores)
  
  names(ct_l) <- cts
  
  return(ct_l)
}



if (!file.exists(corlist_path)) {
  cor_l <- get_cor_list(mat = mat, meta = meta, ncores = ncore)
  saveRDS(cor_l, corlist_path)
} else {
  cor_l <- readRDS(corlist_path)
}

  
# Get the average cor across cell types

# cor_avg <- apply(simplify2array(cor_l), 1:2, mean)
