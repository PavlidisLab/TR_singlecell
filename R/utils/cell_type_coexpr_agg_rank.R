## This script tests/examines the Harris 2021 et al. coexpr aggregation approach
## TODO: fix inconsistent dims of cor->rank
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(parallel)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Load correlation matrices generated per cell type and across all cells
cor_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_celltype.RDS")
cor_all <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_all.RDS")

# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)


# Generating intermediate steps of aggregation as well as w/ and w/o impute NAs
# ------------------------------------------------------------------------------


# Count NA cors and then impute NA cors  to 0 (cor0) or mean (cormean)
# Note that the Harris 2021 et al. approach imputes NA ranks, not cors.

n_na <- count_nas(cor_ct)
cor0_ct <- lapply(cor_ct, na_to_zero)
cormean_ct <- lapply(cor_ct, na_to_mean)


ct_ix <- sample(n_distinct(sdat$Cell_type), 1)
na_ix <- which(is.na(cor_ct[[1]]), arr.ind = TRUE)
sample_na_ix <- na_ix[sample(nrow(na_ix), 1), ]
assertthat::are_equal(cor0_ct[[ct_ix]][sample_na_ix[1], sample_na_ix[2]], 0)
assertthat::are_equal(cormean_ct[[ct_ix]][sample_na_ix[1], sample_na_ix[2]], mean(cor_ct[[ct_ix]], na.rm = TRUE))


# Convert +/- NA imputed cors to ranks, and for -impute convert the resultant
# NA ranks to the average of the rank matrix (as in Harris 2021)

cor_ct_rank <- mclapply(cor_ct, colrank_mat, mc.cores = ncore)
cor0_ct_rank <- mclapply(cor0_ct, colrank_mat, mc.cores = ncore)
cormean_ct_rank <- mclapply(cormean_ct, colrank_mat, mc.cores = ncore)
cor_ct_rankmean <- mclapply(cor_ct_rank, na_to_mean, mc.cores = ncore)

# Ensure max cor equals min (best) rank

ct_ix <- sample(n_distinct(sdat$Cell_type), 1)
tf_ix <- sample(ncol(cor_ct[[ct_ix]]), 1)

assertthat::are_equal(
  names(which.max(cor_ct[[ct_ix]][ , tf_ix])),
  names(which.min(cor_ct_rank[[ct_ix]][, tf_ix])))

assertthat::are_equal(
  names(which.max(cor0_ct[[ct_ix]][ , tf_ix])),
  names(which.min(cor0_ct_rank[[ct_ix]][ , tf_ix])))

assertthat::are_equal(
  names(which.max(cormean_ct[[ct_ix]][ , tf_ix])),
  names(which.min(cormean_ct_rank[[ct_ix]][ , tf_ix])))


# NAs ranked last, want min non-NA ranking

assertthat::are_equal(
  names(which.min(cor_ct[[ct_ix]][ , tf_ix])),
  names(which.max(cor_ct_rank[[ct_ix]][ , tf_ix])))

assertthat::are_equal(
  names(which.min(cor0_ct[[ct_ix]][ , tf_ix])),
  names(which.max(cor0_ct_rank[[ct_ix]][ , tf_ix])))



# Aggregate cell type rank matrices
# https://stackoverflow.com/questions/42628385/sum-list-of-matrices-with-nas
# Mean_nona divides by the number of non-NA observations
# Mean_all divides by the total number of cell types


cor_rank_sum <- apply(simplify2array(cor_ct_rank), 1:2, sum, na.rm = TRUE)
cor0_rank_sum <- apply(simplify2array(cor0_ct_rank), 1:2, sum, na.rm = TRUE)
cormean_rank_sum <- apply(simplify2array(cormean_ct_rank), 1:2, sum, na.rm = TRUE)
cor_rankmean_sum <- apply(simplify2array(cor_ct_rankmean), 1:2, sum, na.rm = TRUE)


gene_ix1 <- sample(nrow(cor_rank_sum), 1)
gene_ix2 <- sample(ncol(cor_rank_sum), 1)


assertthat::are_equal(
  sum(sapply(cor_ct_rank, function(x) x[gene_ix1, gene_ix2]), na.rm = TRUE),
  cor_rank_sum[gene_ix1, gene_ix2])


# Then convert these aggregations to ranks.

final_keepna <- colrank_mat(-cor_rank_sum)
final_cor0 <- colrank_mat(-cor0_rank_sum)
final_cormean <- colrank_mat(-cormean_rank_sum)
final_rankmean <- colrank_mat(-cor_rankmean_sum)

harris2021 <- aggregate_cor(cor_ct, impute_na = TRUE, ncores = ncore)



# Focus on single TF and organize gene rankings
# ------------------------------------------------------------------------------


tf <- "Mef2c"


# Which TF-gene pair had max cor in each cell type
max_cor <- max_cor_df(cor_ct, tf)

# NA df
na_df <- data.frame(Count_NA = n_na[, tf]) %>% 
  rownames_to_column(var = "Symbol")


final_keepna <- colrank_mat(-cor_rank_sum)
final_cor0 <- colrank_mat(-cor0_rank_sum)
final_cormean <- colrank_mat(-cormean_rank_sum)
final_rankmean <- colrank_mat(-cor_rankmean_sum)

harris2021 <- aggregate_cor(cor_ct, impute_na = TRUE, ncores = ncore)


# Gene rankings by different aggregations in single data frame

rank_df <- data.frame(
  Symbol = rownames(harris2021),
  Sum_keepNA = cor_rank_sum[, tf],
  Sum_cor0 = cor0_rank_sum[, tf],
  Sum_cormean = cormean_rank_sum[, tf],
  Sum_rankmean = cor_rankmean_sum[, tf],
  Rank_sum_keepNA = final_keepna[, tf],
  Rank_sum_cor0 = final_cor0[, tf],
  Rank_sum_cormean = final_cormean[, tf],
  Rank_sum_rankmean = final_rankmean[, tf],
  Harris2021 = harris2021[, tf]
) %>%
  left_join(na_df, by = "Symbol") %>% 
  arrange(Harris2021)


# Correlation of the ranks and the count of NAs
rank_df_cor <- 
  round(cor(select_if(rank_df, is.numeric), use = "pairwise.complete.obs"), 3)


# Gene NA count. When equal to number of cell types, it means that those genes
# and the TF were never expressed in the same cell type. The minimum NA count
# reflects the count of cell types in which the TF was not detected. Genes
# equal to this count have a non-NA cor in every cell type with TF detection.

hist(rank_df$Count_NA, breaks = 100)
sum(rank_df$Count_NA == n_distinct(sdat$Cell_type))
sum(rank_df$Count_NA == min(rank_df$Count_NA))

filter(rank_df, Count_NA < 120) %>% head()
filter(rank_df, Count_NA == min(rank_df$Count_NA)) %>% head()


# Inspecting individual genes/cell types
# ------------------------------------------------------------------------------


gene2 <- "Lhx8"

sort(sapply(cor_ct, function(x) x[gene2, tf]), decreasing = TRUE)
sort(sapply(cor_ct_rank, function(x) x[gene2, tf]))

ct <- "GABA-50-Chat-Vip"

head(sort(cor_ct[[ct]][, tf], decreasing = TRUE))


# By cell type
sdat_sub <- subset(sdat, idents = ct)
plot_scatter(sdat_sub, tf, gene2, slot = "data", jitter = TRUE)
plot_scatter(sdat_sub, tf, gene2, slot = "counts", jitter = TRUE)
plot_scatter(sdat_sub, tf, gene2, slot = "counts", jitter = FALSE)
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "data")[c(tf, gene2), ])))
cor(t(as.matrix(GetAssayData(object = sdat_sub, slot = "counts")[c(tf, gene2), ])))

sum(sdat@assays$RNA@counts[gene2, sdat$Cell_type == ct] > 0)

sum(sdat@assays$RNA@counts[tf, sdat$Cell_type == ct] > 0)

sum(sdat@assays$RNA@counts[tf, sdat$Cell_type == ct] > 0 & 
      sdat@assays$RNA@counts[gene2, sdat$Cell_type == ct] > 0)

# Across all cells
plot_scatter(sdat, tf, gene2, slot = "data", jitter = TRUE)
plot_scatter(sdat, tf, gene2, slot = "counts", jitter = TRUE)
plot_scatter(sdat, tf, gene2, slot = "counts", jitter = FALSE)
cor_all[gene2, tf]
cor(t(as.matrix(GetAssayData(object = sdat, slot = "counts")[c(gene2, tf), ])))


###

rank_df <- rank_l$Mouse[[tf]] %>% 
  dplyr::select(Symbol, Rank_integrated, Curated_target) %>% 
  left_join(rank_df, by = "Symbol") %>% 
  filter(Symbol != tf)


keep_cols <- c("Rank_integrated", "Harris2021", "Sum_keepNA", "Sum_cor0", "Sum_cormean")


pr_df <- all_perf_df(rank_df, keep_cols, label_col = "Curated_target", measure = "PR")
auprc <- all_au_perf(rank_df, keep_cols, label_col = "Curated_target", measure = "AUPRC")

roc_df <- all_perf_df(rank_df, keep_cols, label_col = "Curated_target", measure = "ROC")
auroc <- all_au_perf(rank_df, keep_cols, label_col = "Curated_target", measure = "AUROC")

cols <- c(rep("lightgrey", length(keep_cols)))
names(cols) <- keep_cols
cols["Rank_integrated"] <- "black"


plot_perf(df = roc_df, auc_l = auroc, measure = "ROC", cols = cols, title = tf, ncol_legend = 1)
plot_perf(df = pr_df, auc_l = auprc, measure = "PR", cols = cols, title = tf, ncol_legend = 1)
