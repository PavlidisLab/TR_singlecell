## This script looks to compare the rankings from coexpression aggregation
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "Posner2022"
dat <- readRDS(file.path("/space/scratch/amorin/TR_singlecell", id, paste0(id, "_clean_mat_and_meta.RDS")))
rsr_all <- readRDS(file.path("/space/scratch/amorin/TR_singlecell", id, paste0(id, "_RSR_allrank.RDS")))
rsr_all <- lowertri_to_symm(rsr_all)
rsr_col <- readRDS(file.path("/space/scratch/amorin/TR_singlecell", id, paste0(id, "_RSR_colrank.RDS")))
zcor <- readRDS(file.path("/space/scratch/amorin/TR_singlecell", id, paste0(id, "_fishersZ.RDS")))
na_mat <- zcor$NA_mat

stopifnot(identical(rownames(rsr_all), rownames(rsr_col)))
stopifnot(identical(rownames(rsr_all), rownames(zcor$Avg_all)))
stopifnot(identical(rownames(rsr_all), rownames(zcor$NA_mat)))


genes <- rownames(rsr_all)

med_na <- apply(na_mat, 1, median)
min_na <- min(na_mat)
max_na <- max(na_mat)
which_min_na <- which(na_mat == min_na, arr.ind = TRUE)


cor_l <- lapply(genes, function(x) {
  cor(rsr_all[, x], rsr_col[, x], use = "pairwise.complete.obs", method = "spearman")
})
names(cor_l) <- genes


summary(unlist(cor_l))


# Example of best cor

cor_l[which.max(cor_l)]


# Example of cor that disagrees even though there are not many NAs

cor_l[which.min(cor_l)]



gene <- "Rbm25"  # DDX17


rank_df <- data.frame(
  Symbol = rownames(rsr_all),
  RSR_all = rsr_all[, gene],
  RSR_col = rsr_col[, gene],
  Zcor_all = zcor$Avg_all[, gene],
  Zcor_nonNA = zcor$Avg_nonNA[, gene],
  N_na = na_mat[, gene]
)


cor(select_if(rank_df, is.numeric), use = "pairwise.complete.obs", method = "spearman")


plot(rank_df$RSR_all, rank_df$RSR_col)

# Here, splitting the ranks into 10 bins, decreasing by column or all rank. 
# Then correlating the ranks in each bin to show how the random ranking
# of NA->0 cor appears

a1 <- arrange(rank_df, desc(RSR_all)) %>% filter(Symbol != gene)
a2 <- split(a1, cut(1:nrow(a1), 10))
a3 <- lapply(a2, function(x) cor(x$RSR_all, x$RSR_col, method = "spearman"))
plot(unlist(a3))



k <- 1920
sub_all <- slice_max(rank_df, RSR_all, n = k)
sub_col <- slice_max(rank_df, RSR_col, n = k)
cor(select_if(sub_all, is.numeric), method = "spearman")
cor(select_if(sub_col, is.numeric), method = "spearman")
plot(sub_all$RSR_all, sub_all$RSR_col)
plot(sub_col$RSR_all, sub_col$RSR_col)
