## Process count matrix and get aggregate correlation for HypoMap
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "Tabula_Sapiens"
sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets", id)
dat_path <- file.path(sc_dir, paste0(id, "_cellxgene_seurat.RDS"))
out_dir <- file.path("/space/scratch/amorin/TR_singlecell", id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank.RDS"))
colrank_path <- file.path(out_dir, paste0(id, "_RSR_colrank.RDS"))
zcor_path <- file.path(out_dir, paste0(id, "_fishersZ.RDS"))
pc <- read.delim(ens_hg_path, stringsAsFactors = FALSE)


if (!file.exists(processed_path)) {
  
  dat <- readRDS(dat_path)
  
  # Extract count matrix
  
  mat <- GetAssayData(dat, slot = "counts")
  
  # Ready metadata
  
  meta <- dat@meta.data %>% 
    dplyr::rename(Cell_type = cell_type) %>% 
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


# if (!file.exists(colrank_path)) {
#   rsr_col <- RSR_colrank(mat, meta)
#   saveRDS(rsr_col, colrank_path)
# }
# 
# 
# if (!file.exists(zcor_path)) {
#   zcor <- fishersZ_aggregate(mat, meta)
#   saveRDS(zcor, zcor_path)
# }
