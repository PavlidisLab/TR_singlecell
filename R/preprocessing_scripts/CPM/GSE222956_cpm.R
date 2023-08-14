## GSE222956
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE222956"
species <- "Mouse"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)

out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank_CPM.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat_CPM.tsv"))


pc <- read.delim(ref_mm_path, stringsAsFactors = FALSE)


# Files were directly downloaded from GEO, see GSE222956_download.sh in dat_dir
meta_path <- file.path(dat_dir, paste0(id, "_metadata.tsv"))
mat_path <- list.files(dat_dir, pattern = "mtx", full.names = TRUE)
features_path <- list.files(dat_dir, pattern = "features", full.names = TRUE)
barcodes_path <- list.files(dat_dir, pattern = "barcodes", full.names = TRUE)



if (!file.exists(processed_path)) {
  
  # Load metadata and the count matrix
  
  meta <- read.delim(meta_path)
  meta <- data.frame(ID = str_replace(meta$cell.ID, "^cTG_|^ctr_", ""), meta)
  
  features_l <- lapply(features_path, read.delim, header = FALSE)
  barcodes_l <- lapply(barcodes_path, read.delim, header = FALSE)
  
  mat_l <- lapply(1:length(mat_path), function(x) {
    mat <- Matrix::readMM(mat_path[[x]])
    rownames(mat) <- features_l[[x]]$V2
    colnames(mat) <- barcodes_l[[x]]$V1
    return(mat)
  })
  
  stopifnot(identical(rownames(mat_l[[1]]), rownames(mat_l[[length(mat_l)]])))
  
  mat <- do.call(cbind, mat_l)
  mat <- mat[, meta$ID]
  
  stopifnot(identical(colnames(mat), meta$ID))
  
  
  # Ready metadata
  
  change_colnames <- c(Cell_type = "cell.type")
  
  meta <- meta %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(assay = "10x 3' v3") %>% 
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
