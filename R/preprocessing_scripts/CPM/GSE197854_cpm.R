## GSE197854
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "GSE197854"
species <- "Mouse"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)

out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank_CPM.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat_CPM.tsv"))


pc <- read.delim(ref_mm_path, stringsAsFactors = FALSE)


# Files were directly downloaded from GEO, see GSE197854_download.sh in dat_dir
dat_path <- list.files(dat_dir, pattern = "counts", full.names = TRUE)
meta_path <- list.files(dat_dir, pattern = "metadata", full.names = TRUE)


read_data <- function(dat_path) {
  
  mat <- fread(dat_path)
  dupl <- which(duplicated(mat$gene))
  mat <- mat[-dupl, ]
  genes <- mat$gene
  mat$gene <- NULL
  mat <- as.matrix(mat)
  rownames(mat) <- genes
  colnames(mat) <- str_replace(colnames(mat), "\\.", "-")
  
  return(mat)
}




if (!file.exists(processed_path)) {
  
  # Load metadata and the count matrix
  
  meta_l <- lapply(meta_path, function(x) as.data.frame(fread(x)))
  meta <- do.call(rbind, lapply(meta_l, `[`, c("barcodes", "main.pop")))
  
  
  mat_l <- lapply(dat_path, read_data)
  stopifnot(identical(rownames(mat_l[[1]]), rownames(mat_l[[2]])))
  mat <- do.call(cbind, mat_l)
  mat <- mat[, meta$barcodes]
  
  stopifnot(identical(colnames(mat), meta$barcodes))

  
  # Ready metadata
  # "GSE197854" collapse cell types with integer delim and remove unlabeled
  
  change_colnames <- c(Cell_type = "main.pop", ID = "barcodes")
  
  meta <- meta %>% 
    dplyr::rename(any_of(change_colnames)) %>% 
    mutate(
      assay = "10x 5' v2",
      Cell_type = str_replace(Cell_type, "_[:digit:]$", "")
      ) %>% 
    add_count_info(mat = mat)
  
  
  meta <- filter(meta, !is.na(Cell_type) & Cell_type != "-")
  mat <- mat[, meta$ID]
  
  
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
  
  stopifnot(all(colnames(mat) %in% meta$ID), length(meta$ID) > 0)
  
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
