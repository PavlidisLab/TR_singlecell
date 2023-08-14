## This script downloads, preprocesses, and generates the aggregate correlation
## network for the Tabula Muris dataset
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6642641/
## https://figshare.com/articles/dataset/Robject_files_for_tissues_processed_by_Seurat/5821263
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(Seurat)
library(data.table)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "TabulaMuris"
species <- "Mouse"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)

dl_table <- read.delim(file.path(dat_dir, "TabulaMuris_download.tsv"), stringsAsFactors = FALSE, header = FALSE)


# Download each tissue as an Robj file

for (i in 1:nrow(dl_table)) {
  
  dl_path <- file.path(dat_dir, dl_table[i, 1])
  dl_url <- dl_table[i, 2]
  
  if (!file.exists(dl_path)) {
    curl::curl_download(url = dl_url, destfile = dl_path, quiet = FALSE)
  }
}


dat_path <- list.files(dat_dir, pattern = ".Robj", full.names = TRUE)
out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta_CPM.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank_CPM.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat_CPM.tsv"))


pc <- read.delim(ref_mm_path, stringsAsFactors = FALSE)



if (!file.exists(processed_path)) {
  
  # Load and update Seurat objects, each of which is named "tiss", then merge
  
  dat_l <- lapply(dat_path, function(x) {
    load(x)
    UpdateSeuratObject(tiss)
  })
  
  dat <- reduce(dat_l, merge)
  gc(verbose = FALSE)

  
  # Extract count matrix: default counts slot, but use data slot if counts empty
  
  mat <- GetAssayData(dat, slot = "counts")
  
  if (length(mat) == 0 || all(rowSums(mat) == 0)) {
    mat <- GetAssayData(dat, slot = "data")
  }
  
  
  # Ready metadata
  
  meta <- dat[[]] %>% 
    dplyr::rename(Cell_type = cell_ontology_class) %>% 
    rownames_to_column(var = "ID") %>% 
    add_count_info(mat = mat)
  
  
  # QC plots
  
  p1 <- all_hist(meta)
  p2 <- qc_scatter(meta)
  
  ggsave(p1, device = "png", dpi = 300, height = 12, width = 16, bg = "white",
         filename = file.path(out_dir, paste0(id, "_QC_histograms.png")))
  
  ggsave(p2, device = "png", dpi = 300, height = 8, width = 8,
         filename = file.path(out_dir, paste0(id, "_QC_scatter.png")))
  
  
  # This hacky addition is to allow a more relaxed RNA novelty filter, given the
  # sequencing is deeper and shows an RNA novelty distn similar to smart-seq
  
  meta <- mutate(meta, assay = "Smart-seq")
  
  # Remove cells failing QC, keep only protein coding genes, and normalize
  
  mat <- rm_low_qc_cells(mat, meta) %>%
    get_pcoding_only(pcoding_df = pc) %>% 
    Seurat::NormalizeData(., normalization.method = "RC", scale.factor = 1e6, verbose = FALSE)
  
  meta <- meta %>% 
    filter(ID %in% colnames(mat)) %>% 
    mutate(assay = "FACS_Nova-seq")  # update with actual assay
  
  stopifnot(all(colnames(mat) %in% meta$ID), length(meta$ID) > 0)
  
  message(paste("Count of cells:", ncol(mat),
                "Count unique cell types: ", n_distinct(meta$Cell_type)))
  
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
