## Process count matrix and get aggregate correlation for GSE180928
## -----------------------------------------------------------------------------

library(aggtools)
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
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
amat_path <- file.path(out_dir, paste0(id, "_FZ_mat_CPM.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat.tsv"))


pc <- read.delim(ref_hg_path, stringsAsFactors = FALSE)



if (!file.exists(processed_path)) {
  
  meta <- read.delim(meta_path, sep = ",")
  
  mat <- read_count_mat(dat_path)
  colnames(mat) <- str_replace_all(colnames(mat), "\\.", "-")
  mat <- mat[, meta$X]
  
  stopifnot(identical(colnames(mat), meta$X))
  
  
  # Ready metadata
  
  change_colnames <- c(Cell_type = "Lineage", ID = "X")
  
  meta <- meta %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(assay = "10x 3' v2/v3") %>% 
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
    Seurat::NormalizeData(., normalization.method = "RC", scale.factor = 1e6, verbose = FALSE)
    
  meta <- filter(meta, ID %in% colnames(mat))
  mat <- mat[, meta$ID]
  
  stopifnot(identical(colnames(mat), meta$ID), length(meta$ID) > 0)
  
  saveRDS(list(Mat = mat, Meta = meta), file = processed_path)
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
    cor_method = "pearson",
    agg_method = "FZ"
  )

  
  # Save output matrices as .tsv
  fwrite_mat(agg_l$Agg_mat, amat_path)
  fwrite_mat(agg_l$NA_mat, namat_path)
  
}
