## Perform a diff expr analysis for cells binarized by whether they belong to 
## the cell type with the highest average expression for the given TR.
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(future)
library(parallel)
source("R/utils/functions.R")
source("R/00_config.R")

# Output RDS
de_ct_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_TR_DEA_celltype.RDS"

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Matrix of average expression per cell type to find max average expression
ct_avg <- get_ct_avg(sdat)

# Genes to iterate over to find in vs out grouped by expression
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")


# For each TF (j), get the cell type that has the max average expression for the
# given TF, and then for each gene (i) perform a simple diff expression test 
# across all cells for group status in that max average cell type.
# (Using Seurat default Wilcoxon test)
# ------------------------------------------------------------------------------


tf_de <- lapply(tfs, function(x) {
  
  # Add metadata column for whether cell belongs to cell type with highest expr
  sdat <- top_expr_celltype(sdat, avg_mat = ct_avg, gene = x)
  
  message(paste(x, "Top cell type DE", Sys.time()))
  
  de_top_qntl <- FindMarkers(sdat,
                             group.by = "Top_expr_celltype",
                             ident.1 = "TRUE",
                             min.cells.group = 1,
                             min.cells.feature = 1,
                             min.pct = 0,
                             logfc.threshold = 0,
                             only.pos = FALSE)
  
  return(de_top_qntl)
  
})

names(tf_de) <- tfs

# Save out
saveRDS(tf_de, de_ct_path)
