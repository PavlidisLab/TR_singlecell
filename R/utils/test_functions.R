## For testing utils/functions
## -----------------------------------------------------------------------------

library(tidyverse)
library(assertthat)
library(Seurat)
source("R/utils/functions.R")
source("R/00_config.R")

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)

tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

targets <- filter(rank_l$Mouse$Ascl1, Curated_target)$Symbol


# Testing sample_expression_level()
# ------------------------------------------------------------------------------


rank_window <- 10
avg_all <- sort(rowMeans(sdat@assays$RNA@data))


set.seed(12)
tt1 <- sample_expression_level(sdat, targets, rank_window)
ix1_a <- which(names(avg_all) == targets[1])
ix1_b <- which(names(avg_all) == tt1[1])
assertthat::assert_that(abs(ix1_a - ix1_b) <= rank_window)
assertthat::assert_that(ix1_a != ix1_b)


set.seed(13)
tt2 <- sample_expression_level(sdat, targets, rank_window)
ix2_a <- which(names(avg_all) == targets[length(targets)])
ix2_b <- which(names(avg_all) == tt2[length(targets)])
assertthat::assert_that(abs(ix2_a - ix2_b) <= rank_window)
assertthat::assert_that(ix2_a != ix2_b)


# plot(density(log10(avg_all)))
# abline(v = log10(avg_all[ix1_a]), col = "red", lw = 2)
# abline(v = log10(avg_all[ix1_a]), col = "black", lty = 3, lw = 3)



# Add count info to metadata
# head(meta)
# meta <- add_count_info(mat, meta)
# head(meta)
# 
# 
# mat2 <- mat
# meta2 <- meta
# mat2 <- rbind(mat2, "MT-1" = runif(ncol(mat2), min = 0, max = 300))
# mat2[19213:19214, 1:10]
# meta2 <- add_count_info(mat2, meta2)
# head(meta2)


# Remove low QC cells
# dim(mat)
# mat <- rm_low_qc_cells(mat, meta)
# dim(mat)
# 
# dim(mat2)
# mat <- rm_low_qc_cells(mat2, meta2)
# dim(mat2)


# Ensembl ID matching
# mat1 <- mat_norm
# mat2 <- mat_norm
# rownames(mat1) <- 1:nrow(mat1)
# tt1 <- ensembl_to_symbol(mat1, pc)
# tt2 <- ensembl_to_symbol(mat2, pc)
# tt2[1:5, 1:5]