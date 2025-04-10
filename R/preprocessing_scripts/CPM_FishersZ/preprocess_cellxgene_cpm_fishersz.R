## This script takes in arguments for an experiment ID and species, loading the
## associated ID's seurat object and pre-processing the data before generating
## the aggregate correlation matrix and a matrix tracking NAs
## -----------------------------------------------------------------------------

library(aggtools)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
species <- args[2]
cor_method <- "pearson"
agg_method <- "FZ"
force_resave <- FALSE

dat_dir <- file.path(sc_dir, id)
dat_path <- list.files(dat_dir, pattern = "^.*_seurat.*.RDS", full.names = TRUE)
out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
amat_path <- file.path(out_dir, paste0(id, "_FZ_mat_CPM.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat_CPM.tsv"))


# Protein coding gene table (all CxG are ENSEMBL)
pc <- if (str_to_lower(species) %in% c("human", "hg")) {
  read.delim(ens_hg_path, stringsAsFactors = FALSE)
} else if (str_to_lower(species) %in% c("mouse", "mm")) {
  read.delim(ens_mm_path, stringsAsFactors = FALSE)
} else {
  stop("Species not recognized")
}



if (!file.exists(processed_path) || force_resave) {
  
  # Load Seurat object -- if multiple, merge into one
  dat <- lapply(dat_path, readRDS)
  dat <- reduce(dat, merge)
  
  # Extract count matrix: default counts slot, but use data slot if counts empty
  mat <- GetAssayData(dat, slot = "counts")
  
  if (length(mat) == 0 || all(rowSums(mat) == 0)) {
    mat <- GetAssayData(dat, slot = "data")
  }
  
  stopifnot(typeof(mat@i) == "integer")
  
  
  # Ready metadata
  change_colnames <- c(Cell_type = "cell_type", Old_ID = "ID")
  
  meta <- dat[[]] %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    rownames_to_column(var = "ID") %>% 
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
    ensembl_to_symbol(ensembl_df = pc) %>% 
    get_pcoding_only(pcoding_df = pc) %>% 
    Seurat::NormalizeData(., normalization.method = "RC", scale.factor = 1e6, verbose = FALSE)
  
  meta <- filter(meta, ID %in% colnames(mat))
  mat <- mat[, meta$ID]
  
  stopifnot(identical(colnames(mat), meta$ID), length(meta$ID) > 0)
  
  message(paste("Count of cells:", ncol(mat),
                "Count unique cell types: ", n_distinct(meta$Cell_type)))
  
  saveRDS(list(Mat = mat, Meta = meta), file = processed_path)
  
  rm(dat)
  gc()
  
} 



if (!file.exists(amat_path)) {
  
  # Load processed data
  dat <- readRDS(processed_path)
  meta <- dat$Meta
  mat <- dat$Mat
  stopifnot(identical(colnames(mat), meta$ID))
  
  # Run coexpression aggregate
  agg_l <- aggr_coexpr_single_dataset(
    mat = mat, 
    meta = meta,
    pc_df = pc,
    cor_method = cor_method,
    agg_method = agg_method
  )

  
  # Save output matrices as .tsv
  fwrite_mat(agg_l$Agg_mat, amat_path)
  fwrite_mat(agg_l$NA_mat, namat_path)
  
}
