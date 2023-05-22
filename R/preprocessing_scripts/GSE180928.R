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
pc <- read.delim(ref_hg_path, stringsAsFactors = FALSE)


if (!file.exists(processed_path)) {
  
  dat <- read.delim(dat_path, sep = ",")
  meta <- read.delim(meta_path, sep = ",")
  
  # Preparing count matrix
  
  rownames(dat) <- dat$X
  dat$X <- NULL
  colnames(dat) <- str_replace_all(colnames(dat), "\\.", "_")
  mat <- as.matrix(dat)
  
  # Ready metadata
  
  meta <- meta %>% 
    dplyr::rename(ID = X, Cell_type = Cluster) %>% 
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
  
  saveRDS(list(mat, meta), file = processed_path)
  
  rm(dat)
  gc()
  
} else {
  
  dat <- readRDS(processed_path)
  meta <- dat[[2]]
  mat <- dat[[1]]
  
}


stopifnot(identical(colnames(mat), meta$ID))


rsr1 <- all_RSR_aggregate1(mat, meta)
saveRDS(rsr1, file = "/space/scratch/amorin/R_objects/GSE180928_RSR1.RDS")

z1 <- all_zscore_aggregate(mat, meta)
saveRDS(z1, file = "/space/scratch/amorin/R_objects/GSE180928_Z1.RDS")

rsr2 <- all_RSR_aggregate2(mat, meta)
saveRDS(rsr2, file = "/space/scratch/amorin/R_objects/GSE180928_RSR2.RDS")

rsr3 <- all_RSR_aggregate3(mat, meta)
saveRDS(rsr3, file = "/space/scratch/amorin/R_objects/GSE180928_RSR3.RDS")

rsr4 <- all_RSR_aggregate4(mat, meta)
saveRDS(rsr4, file = "/space/scratch/amorin/R_objects/GSE180928_RSR4.RDS")

rsr5 <- all_RSR_aggregate5(mat, meta)
saveRDS(rsr5, file = "/space/scratch/amorin/R_objects/GSE180928_RSR5.RDS")

rsr6 <- all_RSR_aggregate6(mat, meta)
saveRDS(rsr6, file = "/space/scratch/amorin/R_objects/GSE180928_RSR6.RDS")
