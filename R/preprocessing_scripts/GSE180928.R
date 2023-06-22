## Process count matrix and get aggregate correlation for GSE180928
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(data.table)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE180928"

sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata", id)
dat_path <- file.path(sc_dir, paste0(id, "_filtered_cell_counts.csv"))
meta_path <- file.path(sc_dir, paste0(id, "_metadata.csv"))
out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat.tsv"))

# Used for generating/comparing column ranked and avg Fisher's z matrices 
colrank_path <- file.path(out_dir, paste0(id, "_RSR_colrank.tsv"))
zcor_path <- file.path(out_dir, paste0(id, "_fishersZ.RDS"))

pc <- read.delim(ref_hg_path, stringsAsFactors = FALSE)


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
  meta <- dat$Meta
  mat <- dat$Mat
  
}


stopifnot(identical(colnames(mat), meta$ID))


mat <- as.matrix(mat)


if (!file.exists(allrank_path)) {
  
  rsr_all <- RSR_allrank(mat, meta)
  
  suppressMessages(fwrite(
    rsr_all$Agg_mat,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    file = allrank_path
  ))
  
  suppressMessages(fwrite(
    rsr_all$NA_mat,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    file = namat_path
  ))
}



if (!file.exists(colrank_path)) {
  
  rsr_col <- RSR_colrank(mat, meta)
  
  suppressMessages(fwrite(
    rsr_col,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    file = colrank_path
  ))
}


if (!file.exists(zcor_path)) {
  zcor <- fishersZ_aggregate(mat, meta)
  saveRDS(zcor, zcor_path)
}
