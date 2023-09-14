## HPA
## -----------------------------------------------------------------------------

library(WGCNA)
library(tidyverse)
library(data.table)
library(Seurat)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "HPA"
species <- "Human"

dat_dir <- file.path(sc_dir, id)
if (!dir.exists(dat_dir)) dir.create(dat_dir)

out_dir <- file.path(amat_dir, id)
processed_path <- file.path(out_dir, paste0(id, "_clean_mat_and_meta.RDS"))
allrank_path <- file.path(out_dir, paste0(id, "_RSR_allrank.tsv"))
namat_path <- file.path(out_dir, paste0(id, "_NA_mat.tsv"))


pc <- read.delim(ens_hg_path, stringsAsFactors = FALSE)


# Files were directly downloaded from human protein atlas, see HPA_download.sh in dat_dir
tissue_dir <- list.dirs(dat_dir, full.names = TRUE, recursive = FALSE)
cluster_meta_path <- file.path(dat_dir, "rna_single_cell_cluster_description.tsv")


# Load each tissue's count matrix and cell IDs into a list, then match the 
# cell ID/cluster numbers with the respective cell type in cluster_meta.

load_and_match_data <- function(cluster_meta, tissue_dir) {
  
  
  tissues <- unique(cluster_meta$Tissue)
  
  
  dat_l <- lapply(1:length(tissues), function(x) {
    
    # Load mat and meta
    
    tissue <- tissues[x]
    cell_ids <- read.delim(file.path(tissue_dir[x], "cell_data.tsv"))
    mat <- fread(file.path(tissue_dir[x], "read_count.tsv"), header = TRUE)
    
    # Match cluster number to cell type, cleaning up names and setting unique ID
    
    ct <- cluster_meta %>% 
      filter(Tissue == tissue) %>% 
      mutate(
        Cell_type = paste0(tissue, "_", Cell_type),
        Cell_type = str_replace_all(Cell_type, " ", "_"))
    
    cell_meta <- cell_ids %>%
      mutate(cluster = paste0("c-", cluster)) %>%
      dplyr::rename(Cluster = cluster) %>%
      left_join(.,
                ct[, c("Tissue", "Cluster", "Cell_type", "Cell_type_group")],
                by = "Cluster") %>% 
      dplyr::rename(ID = cell_id) %>% 
      mutate(ID = paste0(ID, "_", Cell_type))
    
    # Format matrix
    
    genes <- mat[[1]]
    mat <- as.matrix(mat[, -1])
    rownames(mat) <- genes
    colnames(mat) <- cell_meta$ID
    
    gc(verbose = FALSE)
    list(Mat = mat, Meta = cell_meta)
    
  })
  
  names(dat_l) <- tissues
  
  return(dat_l)
  
}



if (!file.exists(processed_path)) {
  
  # For HPA must load each tissue individually before binding into one mat
  
  cluster_meta <- read.delim(cluster_meta_path)
  colnames(cluster_meta) <- str_replace_all(colnames(cluster_meta), "\\.", "_")
  
  
  # Tissue names slightly differ between paths and meta, inspect to ensure match
  
  tissues <- data.frame(
    Dir = list.dirs(dat_dir, full.names = FALSE, recursive = FALSE),
    Meta = unique(cluster_meta$Tissue)
  )
  
  stopifnot(identical(
    str_to_lower(str_replace(tissues$Dir, "_", "")),
    str_to_lower(str_replace(tissues$Meta, " ", ""))
  ))

  
  dat <- load_and_match_data(cluster_meta, tissue_dir)
  
  
  mat <- do.call(cbind, lapply(dat, `[[`, "Mat"))
  meta <- do.call(rbind, lapply(dat, `[[`, "Meta"))

  
  # Ready metadata
  
  meta <- meta %>% 
    mutate(assay = "NA") %>% 
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
    ensembl_to_symbol(., ensembl_df = pc) %>% 
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
