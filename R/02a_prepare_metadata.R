## Load a google sheets metadata table tracking scRNA-seq experiments, calculate
## details about the count of cell/celltypes and gene measurement, and export
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
library(data.table)
source("R/00_config.R")
source("R/utils/functions.R")

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)


# Load and save a local version of the raw metadata on Gsheets

if (!file.exists(gsheets_meta_raw_path)) {
  
  meta <- read_sheet(gsheets_id, sheet = "Main", trim_ws = TRUE)
  
  write.table(meta,
              sep = "\t",
              row.names = FALSE,
              quote = FALSE,
              file = ghseets_meta_raw_path)
} else {
  meta <- read.delim(gsheets_meta_raw_path, stringsAsFactors = FALSE)
}



# Functions
# ------------------------------------------------------------------------------


# Load a dataset's count matrix and meta and return a list of two data.frames: 
# 1) 1 x 3 of the experiment ID and the count of cells and unique cell types 
# 2) 1 x 2 of the count of cells for each unique cell type 
# This information is coupled due to the cost of loading data.

calc_cell_counts <- function(id) {
  
  dat <- load_dat_list(id)[[1]]
  
  # Count of cells per unique cell type
  ct_count <- dplyr::count(dat$Meta, Cell_type, name = "N_cells")
  
  # Count of total cells and count of unique cell types for the ID
  dat_count <- data.frame(
    ID = id,
    N_cells = ncol(dat$Mat),
    N_celltypes = nrow(ct_count)
  )
  
  return(list(Dat_count = dat_count, Ct_count = ct_count))
  
}



# For the given ID, load the matrix tracking NA counts of gene msrmt pairs, 
# using the diagonal (identity) to binarize if a gene was measured at least
# once (1) or not (0)

calc_gene_msr_vec <- function(id, meta, genes) {
  
  na_mat <- load_agg_mat_list(id, genes = genes, pattern = "_NA_mat_CPM.tsv")[[1]]
  n_celltype <- filter(meta, ID == id)$N_celltype
  binary_msr <- as.integer(diag(na_mat) != n_celltype)
  rm(na_mat)
  gc(verbose = FALSE)
    
  return(binary_msr)
}



# Create a binary gene by experiment matrix that tracks whether a given gene
# was measured (at least one count in at least 20 cells in at least one cell
# type) for each dataset/id.

calc_gene_msr_mat <- function(ids, meta, genes) {
  
  msr_l <- lapply(ids, calc_gene_msr_vec, meta = meta, genes = genes)
  msr_mat <- do.call(cbind, msr_l)
  colnames(msr_mat) <- ids
  rownames(msr_mat) <- genes
  
  return(msr_mat)
}



# Check metadata before adding to metadata 
# ------------------------------------------------------------------------------


# Ensure all IDs in meta have a generated aggregated coexpr matrix
ids <- meta$ID


is_loaded <- unlist(lapply(ids, function(x) {
  file.exists(file.path(amat_dir, x, paste0(x, "_RSR_allrank_CPM.tsv")))
}))


not_loaded <- meta[!is_loaded, "ID"]


if (length(not_loaded) > 0) {
  stop("IDs missing data: ", paste(toString(not_loaded)))
}



# Ensure mouse or human and isolate IDs for gene measurement tracking
stopifnot(all(meta$Species %in% c("Human", "Mouse")))
ids_hg <- filter(meta, Species == "Human")$ID
ids_mm <- filter(meta, Species == "Mouse")$ID



# Inspect assay/technology for consistency
n_platform <- sort(table(meta$Platform))



# Iteratively load each experiment and get count of cells/cell types
meta_counts <- lapply(ids, calc_cell_counts)



# Extract the list of cell counts per cell type for each dataset for later inspection
counts_l <- lapply(meta_counts, `[`, "Ct_count")
names(counts_l) <- ids



# Bind cell counts into meta
meta_cols <- do.call(rbind, lapply(meta_counts, `[[`, "Dat_count"))
stopifnot(identical(meta_cols$ID, meta$ID))
meta$N_cells <- meta_cols$N_cells
meta$N_celltypes <- meta_cols$N_celltypes
stopifnot(is.integer(meta$N_cells), is.integer(meta$N_celltypes))



# Generate binary measurement matrices
msr_mat_hg <- calc_gene_msr_mat(ids_hg, meta, pc_hg$Symbol)
msr_mat_mm <- calc_gene_msr_mat(ids_mm, meta, pc_mm$Symbol)



# Bind gene measurement counts into meta
n_msr_hg <- colSums(msr_mat_hg)
n_msr_mm <- colSums(msr_mat_mm)
n_msr <- c(n_msr_hg, n_msr_mm)[ids]
stopifnot(is.numeric(n_msr))
meta$N_genes <- n_msr



# Save out
# ------------------------------------------------------------------------------



write.table(meta,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = sc_meta_path)


saveRDS(counts_l, file = celltype_list_path)


saveRDS(msr_mat_hg, msr_mat_hg_path)


saveRDS(msr_mat_mm, msr_mat_mm_path)
