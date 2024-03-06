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



# Generates a list of 2 matrices tracking gene measurement, bundled as they
# load and use the same experiment-specific NA tracking matrices.

# msr_mat is a binary gene x experiment matrix tracking gene measurement. this
# is used to describe gene coverage for each experiment.

# comsr_mat is a gene x gene matrix tallying the number of experiments in
# which pairs of genes were co-measured in at least on cell type. this is used
# to describe the count of contributing experiments for the TR rankings


gen_msr_mat_list <- function(ids, meta, genes) {
  
  # Init tracking matrices
  msr_mat <- matrix(0, nrow = length(genes), ncol = length(ids))
  rownames(msr_mat) <- genes
  colnames(msr_mat) <- ids
  
  comsr_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
  rownames(comsr_mat) <- colnames(comsr_mat) <- genes
  
  for (id in ids) {
    
    # NA mat counts the NA gene pairs; max value is the count of cell types
    na_mat <- load_agg_mat_list(id, genes = genes, pattern = "_NA_mat_CPM.tsv")[[1]]
    n_celltype <- filter(meta, ID == id)$N_celltype
    
    # Yes/no if gene pairs were co-measured, and add to tracking mat
    comsr <- (na_mat < n_celltype) & (na_mat < n_celltype)
    comsr_mat <- comsr_mat + comsr
    
    # Diag (self gene) equal to count of cell types means all cell types were NA
    msr_mat[, id] <- as.integer(diag(na_mat) != n_celltype)
    
    rm(na_mat)
    gc(verbose = FALSE)
  }
  
  return(list(Msr_mat = msr_mat, Comsr_mat = comsr_mat))
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



# Generate measurement matrices
msr_l_hg <- gen_msr_mat_list(ids_hg, meta, pc_hg$Symbol)
msr_l_mm <- gen_msr_mat_list(ids_mm, meta, pc_mm$Symbol)



# Bind gene measurement counts into meta
n_msr_hg <- colSums(msr_l_hg$Msr_mat)
n_msr_mm <- colSums(msr_l_mm$Msr_mat)
n_msr <- c(n_msr_hg, n_msr_mm)[ids]
stopifnot(is.numeric(n_msr))
meta$N_genes <- n_msr



# Save out
# ------------------------------------------------------------------------------


# Meta
write.table(meta,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = sc_meta_path)


# Cell type counts
saveRDS(counts_l, file = celltype_list_path)


# Measurement matrices
saveRDS(msr_l_hg$Msr_mat, msr_mat_hg_path)
saveRDS(msr_l_mm$Msr_mat, msr_mat_mm_path)


# Co-measurement matrices
saveRDS(msr_l_hg$Comsr_mat, comsr_mat_hg_path)
saveRDS(msr_l_mm$Comsr_mat, comsr_mat_mm_path)
