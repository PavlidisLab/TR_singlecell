library(WGCNA)
library(tidyverse)
library(Matrix)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")

sc_dir <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/jules_garreau_sc_datasets/"
dat_path <- file.path(sc_dir, "GSE144136/GSE144136_GeneBarcodeMatrix_Annotated.mtx.gz")
meta_path <- file.path(sc_dir, "GSE144136/GSE144136_CellNames.csv.gz")
genes_path <- file.path(sc_dir, "GSE144136/GSE144136_GeneNames.csv.gz")
out_path <- "/space/scratch/amorin/R_objects/GSE144136_mat_and_meta.RDS"
pc <- read.delim("/home/amorin/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)


if (!file.exists(out_path)) {
  
  dat <- as.matrix(Matrix::readMM(dat_path))
  meta <- read.delim(meta_path, sep = ",")
  genes <- read.delim(genes_path, sep = ",")
  
  stopifnot(nrow(meta) == ncol(dat), nrow(genes) == nrow(dat))
  
  # Meta is concat in name. Need to split. Variable delimiter names so done
  # step-wise.
  
  split_meta <- str_split(meta$x, "_", simplify = FALSE) 
    
  cell_tag <- unlist(lapply(split_meta, function(x) tail(x, 1)))
  batch <- unlist(lapply(split_meta, function(x) x[length(x) - 1]))
  condition <- unlist(lapply(split_meta, function(x) x[length(x) - 2]))
  
  # Remove batch/condition info from shorter cell type labels
  rm <- unique(c(batch, condition))
  
  # Different granularity of clusters/celltypes
  
  cluster <- lapply(split_meta, function(x) (x[1:4]))
  cluster <- lapply(cluster, function(x) paste(x[!x %in% rm], collapse = "_"))
  cluster <- lapply(cluster, function(x) str_replace(x, "(\\.[:digit:]+).*", "\\1"))
  cluster <- lapply(cluster, function(x) str_replace(x, "\\/", "_"))
  
  celltype <- lapply(cluster, function(x) str_replace(x, "\\..*", ""))
  
  
  meta_final <- data.frame(
    ID = meta$x,
    Cell_tag = unlist(cell_tag),
    Batch = unlist(batch),
    Condition = unlist(condition),
    Cluster = unlist(cluster),
    Cell_type = unlist(celltype)
  )
  
  
  mat <- dat
  colnames(mat) <- meta_final$ID
  rownames(mat) <- genes$x
  
  common_genes <- intersect(pc$Symbol, genes$x)
  mat <- mat[common_genes, ]

  stopifnot(all(colnames(mat) %in% meta_final$ID))
  
  saveRDS(list(mat, meta_final), file = out_path)
  
  rm(dat)
  gc()
  
} else {
  
  dat <- readRDS(out_path)
  meta <- dat[[2]]
  mat <- dat[[1]]
  
}


stopifnot(identical(colnames(mat), meta$ID))



rsr1 <- all_RSR_aggregate1(mat, meta)
saveRDS(rsr1, file = "/space/scratch/amorin/R_objects/GSE144136_RSR1.RDS")

z1 <- all_zscore_aggregate(mat, meta)
saveRDS(z1, file = "/space/scratch/amorin/R_objects/GSE144136_Z1.RDS")

rsr2 <- all_RSR_aggregate2(mat, meta)
saveRDS(rsr2, file = "/space/scratch/amorin/R_objects/GSE144136_RSR2.RDS")