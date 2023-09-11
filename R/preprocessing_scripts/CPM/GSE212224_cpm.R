## GSE212224
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE212224"
species <- "Mouse"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)

out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank_CPM.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat_CPM.tsv"))


pc <- read.delim(ref_mm_path, stringsAsFactors = FALSE)


# Files were directly downloaded from GEO, see GSE212224_download.sh in dat_dir
meta_path <- list.files(dat_dir, pattern = "metadata", full.names = TRUE)
mat_path <- list.files(dat_dir, pattern = "counts", full.names = TRUE)
features_path <- list.files(dat_dir, pattern = "features", full.names = TRUE)
barcodes_path <- list.files(dat_dir, pattern = "barcodes", full.names = TRUE)



if (!file.exists(processed_path)) {
  
  # Load metadata and the count matrix
  # "GSE212224" extensive duplicates
  
  meta_l <- lapply(meta_path, fread)
  
  features_l <- lapply(features_path, read.delim, header = FALSE)
  barcodes_l <- lapply(barcodes_path, read.delim, header = FALSE)
  
  mat_l <- lapply(1:length(mat_path), function(x) {
    mat <- Matrix::readMM(mat_path[[x]])
    rownames(mat) <- features_l[[x]]$V1
    colnames(mat) <- barcodes_l[[x]]$V1
    return(mat)
  })
  
  
  dupl <- intersect(colnames(mat_l[[1]]), colnames(mat_l[[2]]))
  set1 <- setdiff(colnames(mat_l[[1]]), colnames(mat_l[[2]]))
  set2 <- setdiff(colnames(mat_l[[2]]), colnames(mat_l[[1]]))
  common_genes <- intersect(rownames(mat_l[[1]]), rownames(mat_l[[2]]))
  
  
  mat <- cbind(mat_l[[1]][common_genes, c(dupl, set1)], 
               mat_l[[2]][common_genes, set2])
  
  meta <- rbind(
    filter(meta_l[[1]], V1 %in% c(dupl, set1)) %>% dplyr::select(c(V1, CellType)), 
    filter(meta_l[[2]], V1 %in% set2) %>% dplyr::select(c(V1, CellType))
  )
  
  mat <- as.matrix(mat[, meta$V1])

  stopifnot(identical(colnames(mat), meta$V1))
  
  
  # Ready metadata
  # "GSE212224" collapse and clean cell type
  
  change_colnames <- c(Cell_type = "CellType", ID = "V1")
  
  meta <- as.data.frame(meta) %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(
      assay = "Seq-well",
      Cell_type = str_replace(Cell_type, "_[:digit:]$", ""),
      Cell_type = str_replace(Cell_type, "[:digit:]$", ""),  
      Cell_type = str_replace(Cell_type, "Lineage", "lineage"),
      Cell_type = str_replace(Cell_type, "Neutropils", "Neutrophils")
    ) %>% 
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
