## Load a google sheets metadata table tracking scRNA-seq experiments, calculate
## additional details from the finished experiments, format and save
## 1) X
## 2) X
## 3) X
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
library(parallel)
source("R/00_config.R")
source("R/utils/functions.R")

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)


# Load and save a local version of the raw metadata on Gsheets

if (!file.exists(ghseets_meta_raw_path)) {
  
  meta <- read_sheet(gsheets_id, sheet = "Main", trim_ws = TRUE)
  
  write.table(meta,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE,
              file = ghseets_meta_raw_path)
} else {
  meta <- read.delim(ghseets_meta_raw_path, stringsAsFactors = FALSE)
}




# Load a dataset's count matrix and meta and return a list of two data.frames: 
# 1) 1 x 3 of the experiment ID and the count of unique cell types and total number of cells 
# 2) 1 x 2 of the count of cells for each unique cell type 
# ------------------------------------------------------------------------------


add_meta <- function(id) {
  
  dat <- load_dat_list(id)[[1]]
  
  ct_count <- dplyr::count(dat$Meta, Cell_type, name = "N_cells")
  
  dat_count <- data.frame(
    ID = id,
    N_cells = ncol(dat$Mat),
    N_celltypes = nrow(ct_count)
  )
  
  return(list(Dat_count = dat_count, Ct_count = ct_count))
  
}


# Create a binary gene by experiment matrix that tracks whether a given gene
# was measured (at least one count in at least 20 cells in at least one cell
# type) for each dataset/id.
# ------------------------------------------------------------------------------


get_gene_msr_mat <- function(ids, meta, genes) {
  
  msr_mat <- matrix(0, nrow = length(genes), ncol = length(ids))
  colnames(msr_mat) <- ids
  rownames(msr_mat) <- genes
  
  # For each dataset, load the matrix tracking NA counts of gene msrmt. pairs, 
  # using the diagonal (identity) to binarize if a gene was measured at least
  # once (1) or not (0)
  
  for (id in ids) {
    
    na_mat <- load_agg_mat_list(
      id, genes = genes, pattern = "_NA_mat_CPM.tsv")[[1]]
    
    n_celltype <- filter(meta, ID == id)$N_celltype
    
    binary_msr <- as.integer(diag(na_mat) != n_celltype)
    msr_mat[, id] <- binary_msr
    
    rm(na_mat)
    gc(verbose = FALSE)
    
  }
  
  return(msr_mat)
}




# Check whether the given IDs have an aggregate coexpression matrix, then load,
# inspect, and add count info
# ------------------------------------------------------------------------------


loaded <- lapply(meta$ID, function(x) {
  file.exists(file.path(amat_dir, x, paste0(x, "_RSR_allrank_CPM.tsv")))
})


meta_loaded <- meta[unlist(loaded), ]


missing <- setdiff(meta$ID, meta_loaded$ID)


stopifnot(all(meta_loaded$Species %in% c("Human", "Mouse")))


# Inspect assay/technology for consistency
n_platform <- sort(table(meta$Platform))



# Iteratively load each experiment and get cell/cell type counts
meta_counts <- lapply(meta_loaded$ID, add_meta)



# Extract the list of cell counts per cell type for each dataset for later inspection
counts_l <- lapply(meta_counts, `[`, "Ct_count")
names(counts_l) <- meta_loaded$ID


# Bind cell counts and enter into meta
meta_cols <- do.call(rbind, lapply(meta_counts, `[[`, "Dat_count"))
stopifnot(identical(meta_cols$ID, meta_loaded$ID))
meta_loaded$N_cells <- meta_cols$N_cells
meta_loaded$N_celltypes <- meta_cols$N_celltypes



# Save meta and list of cell types

write.table(meta_loaded,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = sc_meta_path)


saveRDS(counts_l, file = celltype_list_path)




# 
if (!file.exists(msr_mat_hg_path) || force_resave) {
  msr_mat_hg <- get_gene_msr_mat(ids_hg, sc_meta, pc_hg$Symbol)
  saveRDS(msr_mat_hg, msr_mat_hg_path)
} else {
  msr_mat_hg <- readRDS(msr_mat_hg_path)
}


if (!file.exists(msr_mat_mm_path) || force_resave) {
  msr_mat_mm <- get_gene_msr_mat(ids_mm, sc_meta, pc_mm$Symbol)
  saveRDS(msr_mat_mm, msr_mat_mm_path)
} else {
  msr_mat_mm <- readRDS(msr_mat_mm_path)
}
