library(WGCNA)
library(tidyverse)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")


sc_dir_hg <- "/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata/"
dat_path <- file.path(sc_dir_hg, "GSE180928/GSE180928_filtered_cell_counts.csv")
meta_path <- file.path(sc_dir_hg, "GSE180928/GSE180928_metadata.csv")
out_path <- "/space/scratch/amorin/R_objects/20_04_2023_GSE180928.RDS"


if (!file.exists(out_path)) {
  
  dat <- read.delim(file.path(sc_dir_hg, "GSE180928/GSE180928_filtered_cell_counts.csv"), sep = ",")
  meta <- read.delim(file.path(sc_dir_hg, "GSE180928/GSE180928_metadata.csv"), sep = ",")
  
  rownames(dat) <- dat$X
  dat$X <- NULL
  colnames(dat) <- str_replace_all(colnames(dat), "\\.", "_")
  dat <- as.matrix(dat)
  
  
  meta <- meta %>% 
    dplyr::rename(ID = X, Cell_type = Cluster) %>% 
    mutate(ID = str_replace_all(ID, "-", "_"))
  
  stopifnot(all(colnames(dat) %in% meta$ID))
  
  dat <- dat[, meta$ID]
  
  saveRDS(list(dat, meta), file = out_path)
  
} else {
  
  dat <- readRDS(out_path)
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



# rsr1 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR1.RDS")
# rsr2 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR2.RDS")
# rsr3 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR3.RDS")
# z1 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_Z1.RDS")

