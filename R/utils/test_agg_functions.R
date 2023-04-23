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
  dat <- dat[[1]]
  
}


stopifnot(identical(colnames(dat), meta$ID))

cts <- unique(meta$Cell_type)
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")
genes <- rownames(dat)


# Min count of non-zero expressing cells to keep gene for correlating
min_cell <- 20

# Threshold topn n position for binarizing
thresh <- floor(nrow(dat) * 0.005)
thresh <- thresh + 1  # +1 for the moment to account for keeping self cor


# Subset for speed of testing
mat <- dat[5001:6000, ]
meta_sub <- filter(meta, Cell_type %in% cts[1:2])




# Check init agg mat
# ------------------------------------------------------------------------------


genes <- rownames(mat)
amat1 <- init_agg_mat(row_genes = genes)
# amat2 <- init_agg_mat(row_genes = genes, col_genes = str_to_upper(tfs))




# Subset to cell type and check min to NA. Transpose so genes are columns
# ------------------------------------------------------------------------------


ct <- cts[1]
ct_mat <- t(mat[, filter(meta, Cell_type == ct)$ID])

na_mat <- under_min_count_to_na(ct_mat)

na_gene <- genes[which(is.na(na_mat), arr.ind = TRUE)[, "col"][1]]

stopifnot(sum(ct_mat[, na_gene] != 0) < min_cell)
stopifnot(sum(is.na(na_mat)) > sum(is.na(ct_mat)))


# na_mat[1:5, 1:5]
# ct_mat[1:5, 1:5]


# Check cor. Transpose mat so that genes are cols.
# ------------------------------------------------------------------------------


cmat1 <- get_cor_mat(t(mat))
cmat2 <- get_cor_mat(t(mat), lower_tri = FALSE)
cmat3 <- get_cor_mat(na_mat)
cmat4 <- get_cor_mat(na_mat, lower_tri = FALSE)

all(is.na(cmat3[na_gene, ]))
table(!is.na(cmat2[na_gene, ]))

cmat1[1:5, 1:5]
cmat2[1:5, 1:5]
cmat3[1:5, 1:5]
cmat4[1:5, 1:5]


# Check ranking of each element of matrix
# ------------------------------------------------------------------------------


allrank1 <- allrank_mat(cmat1)
allrank2 <- allrank_mat(cmat2)
allrank3 <- allrank_mat(cmat3)
allrank4 <- allrank_mat(cmat4)

allrank1[1:5, 1:5]
allrank2[1:5, 1:5]

which(cmat1 == max(cmat1, na.rm = TRUE), arr.ind = TRUE)
which(allrank1 == min(allrank1, na.rm = TRUE), arr.ind = TRUE)

cmat1[630, 498, drop = FALSE]
allrank1[630, 498, drop = FALSE]

df <- mat_to_df(allrank1) %>% 
  dplyr::rename(Rank_cor = Value) %>% 
  mutate(Cor = mat_to_df(cmat1)[, "Value"])



# Check row rank. Default is ties as min and keep NA. mat will be ordered by
# cells, na_mat will be ordered by genes
# ------------------------------------------------------------------------------


rowrank1 <- rowrank_mat(mat)
rowrank2 <- rowrank_mat(na_mat)
rowrank3 <- rowrank_mat(mat, ties_arg = "random")
rowrank4 <- rowrank_mat(na_mat, ties_arg = "random")
rowrank5 <- rowrank_mat(mat, na_arg = "last")
rowrank6 <- rowrank_mat(na_mat, na_arg = "last")

assertthat::are_equal(dim(mat), dim(rowrank1))


head(sort(mat[1, ], decreasing = TRUE))
head(sort(rowrank1[1, ], decreasing = FALSE))
head(sort(rowrank2[1, ], decreasing = FALSE))

rowrank1[1:5, 1:5]
rowrank2[1:5, 1:5]
rowrank3[1:5, 1:5]
rowrank4[1:5, 1:5]
rowrank5[1:5, 1:5]
rowrank6[1:5, 1:5]



# Check col rank
# ------------------------------------------------------------------------------


colrank1 <- colrank_mat(t(mat)) # default is min

colrank1[1:5, 1:5]
t(mat)[1:5, 1:5]


# Check cor and rank
# ------------------------------------------------------------------------------


# cr1 <- cor_and_rank(mat, na_arg = "last")
# cr2 <- cor_and_rank(na_mat, na_arg = "last")
# cr3 <- cor_and_rank(mat, na_arg = "mean")
# cr4 <- cor_and_rank(na_mat, na_arg = "mean")
# cr5 <- cor_and_rank(mat, ties_arg = "min")
# cr6 <- cor_and_rank(na_mat, ties_arg = "random")
# 
# cr1[1:5, 1:5]
# cr2[1:5, 1:5]
# cr3[1:5, 1:5]
# cr4[1:5, 1:5]
# cr5[1:5, 1:5]
# cr6[1:5, 1:5]
# 
# identical(cr1, cr2)


# Check RSR and Zscore
# ------------------------------------------------------------------------------


rsr1 <- all_RSR_aggregate1(mat, meta)
rsr2 <- all_RSR_aggregate2(mat, meta)
rsr3 <- all_RSR_aggregate3(mat, meta)
rsr4 <- all_RSR_aggregate4(mat, meta)
rsr5 <- all_RSR_aggregate5(mat, meta)
rsr6 <- all_RSR_aggregate6(mat, meta)

z1 <- all_zscore_aggregate(mat, meta)


summary(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], rsr3[,x])))
summary(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], rsr5[,x])))
summary(sapply(1:ncol(rsr4), function(x) cor(rsr6[,x], rsr5[,x], use = "pairwise.complete.obs")))
summary(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], z1$Avg_all[,x])))
summary(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], z1$Avg_nonNA[,x])))



rank_df <- data.frame(
  Symbol = rownames(rsr1),
  RSR1 = rsr1[, "ASCL1"],
  RSR2 = rsr2[, "ASCL1"],
  RSR3 = rsr3[, "ASCL1"],
  RSR4 = rsr4[, "ASCL1"],
  RSR5 = rsr5[, "ASCL1"],
  RSR6 = rsr6[, "ASCL1"],
  Z_all = z1$Avg_all[, "ASCL1"],
  Z_nonNA = z1$Avg_nonNA[, "ASCL1"]
)


plot(rsr1[,1], z1$Avg_all[,1])
plot(rsr1[,1], z1$Avg_nonNA[,1])
cor(rsr1[,1], z1$Avg_nonNA[,1], method = "spearman")
cor(rsr1[,1], z1$Avg_nonNA[,1], method = "spearman")

summary(sapply(1:ncol(z1$Avg_all), function(x) cor(z1$NA_mat[, x], z1$Avg_all[, x], method = "spearman")))
summary(sapply(1:ncol(z1$Avg_all), function(x) cor(z1$NA_mat[, x], z1$Avg_nonNA[, x], method = "spearman")))




# z1$Avg_all[1:5, 1:5]
# z1$Avg_nonNA[1:5, 1:5]
all_celltype_cor(mat, meta, gene1 = "ASCL1", gene2 = "OLIG1")



# Need to keep track of NAs -> 0s without ranking symmetric NAs
# ------------------------------------------------------------------------------


ct_mat <- t(mat[, filter(meta, Cell_type == cts[1])$ID])
ct_mat <- under_min_count_to_na(ct_mat, min_cell)

cmat_tri <- get_cor_mat(ct_mat)
cmat_full <- get_cor_mat(ct_mat, lower_tri = FALSE)

cmat_full_nato0 <- cmat_full
cmat_full_nato0[is.na(cmat_full_nato0)] <- 0

cmat_full_nato0_tri <- cmat_full_nato0
diag(cmat_full_nato0_tri) <- NA
cmat_full_nato0_tri[upper.tri(cmat_full_nato0_tri)] <- NA

rmat_full_nato0_tri <- allrank_mat(cmat_full_nato0_tri)



ct_mat[1:5, 1:5]
cmat_tri[1:5, 1:5]
cmat_full[1:5, 1:5]
cmat_full_nato0[1:5, 1:5]
cmat_full_nato0_tri[1:5, 1:5]
rmat_full_nato0_tri[1:5, 1:5]


sum(is.na(cmat_full))
which(is.na(cmat_full), arr.ind = TRUE)[1,]
cmat_full[24, 1, drop = FALSE]

sum(is.na(cmat_full_nato0))
sum(is.na(cmat_tri))
sum(is.na(cmat_full_nato0_tri))

which(cmat_full_nato0_tri == 0, arr.ind = TRUE)[1,]
cmat_full_nato0_tri[24, 1, drop = FALSE]


which(rmat_full_nato0_tri == min(rmat_full_nato0_tri, na.rm = TRUE), arr.ind = TRUE)
cmat_full_nato0_tri[630, 498, drop = FALSE]
rmat_full_nato0_tri[630, 498, drop = FALSE]




# Speed of cor. ncore8 actually slightly slows down, likely NA related.
# ------------------------------------------------------------------------------


res <- microbenchmark::microbenchmark(
  F1 = WGCNA::cor(mat, nThreads = 1),
  F2 = WGCNA::cor(mat, nThreads = 8),
  F3 = WGCNA::corAndPvalue(mat),
  F4 = stats::cor(mat),
  times = 10
)



# Checking 0 vs NA into cor (does having 0s allow faster parallel)
# Doesn't seem so

ct <- cts[1]
ct_mat1 <- t(mat[, filter(meta, Cell_type == ct)$ID])
ct_mat2 <- under_min_count_to_na(ct_mat, min_cell)



res <- microbenchmark::microbenchmark(
  F1 = get_cor_mat(ct_mat1, lower_tri = FALSE, ncores = 1),
  F2 = get_cor_mat(ct_mat1, lower_tri = FALSE, ncores = 8),
  F3 = get_cor_mat(ct_mat2, lower_tri = FALSE, ncores = 1),
  F4 = get_cor_mat(ct_mat2, lower_tri = FALSE, ncores = 8),
  times = 10
)


# Pearson vs Spearman

res <- microbenchmark::microbenchmark(
  F1 = get_cor_mat(ct_mat1, cor_method = "pearson", lower_tri = FALSE, ncores = 1),
  F2 = get_cor_mat(ct_mat1, cor_method = "pearson", lower_tri = FALSE, ncores = 8),
  F3 = get_cor_mat(ct_mat1, cor_method = "spearman", lower_tri = FALSE, ncores = 1),
  F4 = get_cor_mat(ct_mat1, cor_method = "spearman", lower_tri = FALSE, ncores = 8),
  times = 10
)