## GSE173096
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE173096"
species <- "Mouse"

dat_dir <- file.path(sc_dir, id)

out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat.tsv"))

pc <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

# Downloaded from cxg but need to remove duplicated cells from GSE145929
dat_path <- file.path(dat_dir, paste0(id, "_cellxgene_seurat.RDS"))
meta_path <- file.path(dat_dir, paste0(id, "_metadata.txt"))
dupl_path <- file.path(dat_dir, paste0(id, "_dupl_ids_to_remove.tsv"))



if (!file.exists(processed_path)) {
  
  # Load Seurat object and duplicated IDs
  
  dat <- readRDS(dat_path)
  dupl_ids <- read.delim(dupl_path, stringsAsFactors = FALSE)
  
  
  # Extract count matrix: default counts slot, but use data slot if counts empty
  
  mat <- GetAssayData(dat, slot = "counts")
  
  if (length(mat) == 0 || all(rowSums(mat) == 0)) {
    mat <- GetAssayData(dat, slot = "data")
  }
  
  
  # Ready metadata
  # "GSE173096" Remove duplicated ID
  
  change_colnames <- c(Cell_type = "cell_type", ID = "V1")
  
  meta <- dat[[]] %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(assay = "10x 3' v2/v3") %>% 
    rownames_to_column(var = "ID") %>% 
    filter(!(ID %in% dupl_ids[[1]]))
  
  mat <- mat[, meta$ID]
  meta <- add_count_info(mat = mat, meta = meta)
  
  
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
    Seurat::LogNormalize(., verbose = FALSE)
  
  meta <- filter(meta, ID %in% colnames(mat))
  mat <- mat[, meta$ID]
  
  stopifnot(identical(colnames(mat), meta$ID), length(meta$ID) > 0)
  
  message(paste("Count of cells:", ncol(mat),
                "Count unique cell types: ", n_distinct(meta$Cell_type)))
  
  saveRDS(list(Mat = mat, Meta = meta), file = processed_path)
  gc()
  
} else {
  
  dat <- readRDS(processed_path)
  meta <- dat$Meta
  mat <- dat$Mat
  
}


stopifnot(identical(colnames(mat), meta$ID))


mat <- as.matrix(mat)


if (!file.exists(allrank_path)) {
  
  rsr_all <- RSR_allrank(mat, meta)
  
  # Write as data.frames (preserve rownames) with data.table fwrite (fast)
  
  fwrite(
    data.frame(rsr_all$Agg_mat, check.names = FALSE),
    sep = "\t",
    row.names = TRUE,
    quote = FALSE,
    verbose = FALSE,
    showProgress = FALSE,
    file = allrank_path
  )
  
  
  fwrite(
    data.frame(rsr_all$NA_mat, check.names = FALSE),
    sep = "\t",
    row.names = TRUE,
    quote = FALSE,
    verbose = FALSE,
    showProgress = FALSE,
    file = namat_path
  )
}
