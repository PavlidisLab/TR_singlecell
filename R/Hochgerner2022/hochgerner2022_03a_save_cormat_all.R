## Generate gene-gene correlation matrix across all cells
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(future)
library(parallel)
library(WGCNA)
source("R/utils/functions.R")
source("R/00_config.R")


# Output RDS
cor_all_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_cormat_all.RDS"

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)


message(paste("Begin", Sys.time()))


cor_all <- WGCNA::cor(x = t(as.matrix(sdat@assays$RNA@data)),
                      nThreads = ncore,
                      use = "pairwise.complete.obs")


message(paste("End", Sys.time()))


saveRDS(cor_all, cor_all_path)
