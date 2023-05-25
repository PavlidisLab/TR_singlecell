## This script looks to compare the rankings from coexpression aggregation
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

dat <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_clean_mat_and_meta.RDS")
rsr_all <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_RSR_allrank.RDS")
rsr_all <- lowertri_to_symm(rsr_all)
rsr_col <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_RSR_colrank.RDS")
zcor <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_fishersZ.RDS")
na_mat <- zcor$NA_mat

stopifnot(identical(rownames(rsr_all), rownames(rsr_col)))
stopifnot(identical(rownames(rsr_all), rownames(zcor$Avg_all)))
stopifnot(identical(rownames(rsr_all), rownames(zcor$NA_mat)))


gene <- "MECP2"


rank_df <- data.frame(
  Symbol = rownames(rsr_all),
  RSR_all = rsr_all[, gene],
  RSR_col = rsr_col[, gene],
  Zcor_all = zcor$Avg_all[, gene],
  Zcor_nonNA = zcor$Avg_nonNA[, gene],
  N_na = na_mat[, gene]
)


cor(select_if(rank_df, is.numeric))


tt <- mat_to_df(rsr_all, symmetric = TRUE)


tt2 <- rank_df %>% slice_max(RSR_all, n = 1000)
cor(select_if(tt2, is.numeric))


cor(amat_colrank_std[, gene], amat_allrank_std[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_colrank_std[, gene], amat_allrank_std[, gene])


cor(amat_allrank_std[, gene], amat_allrank_nostd[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_allrank_std[, gene], amat_allrank_nostd[, gene])


cor(amat_colrank_std[, gene], amat_allrank_nostd[, gene], method = "spearman", use = "pairwise.complete.obs")
plot(amat_colrank_std[, gene], amat_allrank_nostd[, gene])


