## Perform a diff expr analysis for cells binarized by whether they are in the 
## top quantile of expression for the given TRs.
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(future)
source("R/utils/functions.R")
source("R/00_config.R")

# Output RDS
de_ct_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_TR_DEA_quantile.RDS"

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Genes to iterate over to find in vs out grouped by expression
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")


# For each TF (j), group cells by whether or not they are in the top quantile of
# expression for the given TF, and then for each gene (i) perform a simple diff
# expression test across all cells for group status in the top quantile cells.
# (Using Seurat default Wilcoxon test)
# ------------------------------------------------------------------------------


tf_de <- lapply(tfs, function(x) {
  
  # Add metadata column for quantile group status
  sdat <- top_expr_quantile(sdat, gene = x, qtl = 0.9)
  
  message(paste(x, "Top quantile DE", Sys.time()))
  
  de_top_qntl <- FindMarkers(sdat,
                             group.by = "Top_expr_quantile",
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
saveRDS(tf_de, outfile)
