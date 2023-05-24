## This script looks to compare the rankings from coexpression aggregation
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

dat <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_clean_mat_and_meta.RDS")
rsr_all <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_RSR_allrank.RDS")
rsr_col <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_RSR_colrank.RDS")
zcor <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_fishersZ.RDS")


# TODO: replace with proper output

saveRDS(
  list(
    Col_std = amat_colrank_std,
    All_std = amat_allrank_std,
    All_nostd = amat_allrank_nostd
  ),
  file = "/space/scratch/amorin/R_objects/Posner2022_subset_RSR_comparison.RDS"
)

gene <- "Rpl37a"

cor(amat_colrank_std[, gene], amat_allrank_std[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_colrank_std[, gene], amat_allrank_std[, gene])


cor(amat_allrank_std[, gene], amat_allrank_nostd[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_allrank_std[, gene], amat_allrank_nostd[, gene])


cor(amat_colrank_std[, gene], amat_allrank_nostd[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_colrank_std[, gene], amat_allrank_nostd[, gene])


