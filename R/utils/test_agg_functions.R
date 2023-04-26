library(WGCNA)
library(tidyverse)
library(parallel)
library(cowplot)
source("R/utils/dev_aggregate_functions.R")

dat <- readRDS("/space/scratch/amorin/R_objects/GSE180928_mat_and_meta.RDS")
mat <- dat[[1]]
meta <- dat[[2]]

stopifnot(identical(colnames(mat), meta$ID))

cts <- unique(meta$Cell_type)
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")
genes <- rownames(dat)

# Min count of non-zero expressing cells to keep gene for correlating
min_cell <- 20

# Subset for speed of testing
mat <- mat[5001:6000, ]
meta_sub <- filter(meta, Cell_type %in% cts[1:2])



# Check init agg mat
# ------------------------------------------------------------------------------


genes <- rownames(mat)
amat1 <- init_agg_mat(row_genes = genes)
# amat2 <- init_agg_mat(row_genes = genes, col_genes = str_to_upper(tfs))
stopifnot(all(amat1 == 0))



# Subset to cell type and check min to NA. Transposes mat so genes are columns
# ------------------------------------------------------------------------------


ct <- cts[1]
ct_mat <- ct_mat <- subset_and_filter(mat, meta, cell_type = ct, min_count = min_cell)

na_gene <- genes[which(is.na(ct_mat), arr.ind = TRUE)[, "col"][1]]

stopifnot(
  sum(mat[na_gene, filter(meta, Cell_type == ct)$ID] != 0) < min_cell)

stopifnot(
  sum(is.na(mat[na_gene, filter(meta, Cell_type == ct)$ID])) < sum(is.na(ct_mat)))


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
# rowrank6 <- rowrank_mat(na_mat, na_arg = "last")

assertthat::are_equal(dim(mat), dim(rowrank1))


# head(sort(mat[1, ], decreasing = TRUE))
# head(sort(rowrank1[1, ], decreasing = FALSE))
# head(sort(rowrank2[1, ], decreasing = FALSE))
# 
# rowrank1[1:5, 1:5]
# rowrank2[1:5, 1:5]
# rowrank3[1:5, 1:5]
# rowrank4[1:5, 1:5]
# rowrank5[1:5, 1:5]
# rowrank6[1:5, 1:5]



# Check col rank
# ------------------------------------------------------------------------------


colrank1 <- colrank_mat(t(mat)) # default is min

# colrank1[1:5, 1:5]
# t(mat)[1:5, 1:5]


# Check RSR and Zscore
# 1: Col rank NA cors to 0
# 2: All rank NA cors to 0
# 3: Col rank NA cors to mean
# 4: All rank NA cors to mean
# 5: Col rank NA cors to last
# 6: All rank NA cors to last
# ------------------------------------------------------------------------------


rsr1 <- all_RSR_aggregate1(mat, meta)
rsr2 <- all_RSR_aggregate2(mat, meta)
rsr3 <- all_RSR_aggregate3(mat, meta)
rsr4 <- all_RSR_aggregate4(mat, meta)
rsr5 <- all_RSR_aggregate5(mat, meta)
rsr6 <- all_RSR_aggregate6(mat, meta)

z1 <- all_zscore_aggregate(mat, meta)


# rsr1 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR1.RDS")
# rsr2 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR2.RDS")
# rsr3 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR3.RDS")
# z1 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_Z1.RDS")


# The all rank matrices are lower tri (NAs in upper), fill out for comparison

rsr_full2 <- lowertri_to_symm(rsr2, na_diag = FALSE)
rsr_full4 <- lowertri_to_symm(rsr4, na_diag = FALSE)
rsr_full6 <- lowertri_to_symm(rsr6, na_diag = FALSE)


stopifnot(identical(rsr_full2[1, ], rsr_full2[, 1]))



# Retrieve the column wise correlations between two aggregation matrices


plot_scatter <- function(mat1, mat2, ix) {
  
  p <- data.frame(Mat1 = mat1[, ix], Mat2 = mat2[, ix]) %>% 
    ggplot(aes(x = Mat1, y = Mat2)) +
    geom_point(shape = 21, alpha = 0.6) + 
    ggtitle(rownames(mat1)[ix]) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 20))
  
  return(p)
}


cor_summary <- function(mat1, mat2, ncores = 8) {
  
  cor_l <- mclapply(1:ncol(mat1), function(x) {
    cor(mat1[, x], mat2[, x], method = "spearman")
  }, mc.cores = ncores)
  
  min_ix = which.min(unlist(cor_l))

  p <- plot_scatter(mat1, mat2, min_ix)
  
  return(list(
    Summary = summary(unlist(cor_l)),
    Min_ix = min_ix,
    Min_plot = p
  ))
}




# Comparing all rank with column rank when setting NA cors to 0

summ_1vs2 <- cor_summary(mat1 = rsr1, mat2 = rsr_full2)
p_1vs2_a <- plot_scatter(mat1 = rsr1, mat2 = rsr_full2, ix = which(rownames(rsr1) == "ASCL1"))
p_1vs2_a <- p_1vs2_a + xlab("Column rank") + ylab("All rank")
p_1vs2_b <- summ_1vs2$Min_plot + xlab("Column rank") + ylab("All rank")
p_1vs2 <- plot_grid(p_1vs2_a, p_1vs2_b)


# Comparing all rank with column rank when setting NA ranks to mean
summ_3vs4 <- cor_summary(mat1 = rsr3, mat2 = rsr_full4)


# Comparing all rank with column rank when setting NA ranks to last
summ_5vs6 <- cor_summary(mat1 = rsr5, mat2 = rsr_full6)


# Comparing all rank setting NA cor to 0 with NA rank to mean rank 
summ_2vs4 <- cor_summary(mat1 = rsr2, mat2 = rsr_full4)


# Comparing col rank setting NA cor to 0 with NA rank to mean rank 
summ_1vs3 <- cor_summary(mat1 = rsr1, mat2 = rsr3)




# Comparing col rank setting NA cor to 0 with NA rank to Z score
summary(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], z1$Avg_all[,x], method = "spearman")))
summary(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], z1$Avg_nonNA[,x], method = "spearman")))
which.min(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], z1$Avg_all[,x], method = "spearman")))
which.min(sapply(1:ncol(rsr1), function(x) cor(rsr1[,x], z1$Avg_nonNA[,x], method = "spearman")))
plot(rsr1[, 177], z1$Avg_all[, 177])
plot(rsr1[, 763], z1$Avg_nonNA[, 763])
head(sort(z1$Avg_all[177, ], decreasing = TRUE))  # unusually low Zscore 


# Comparing Z score average by all CTs or just by non-NA CTs
summary(sapply(1:ncol(z1$Avg_all), function(x) cor(z1$Avg_all[, x], z1$Avg_nonNA[, x], method = "spearman")))
which.min(sapply(1:ncol(z1$Avg_all), function(x) cor(rsr1[,x], z1$Avg_all[,x], method = "spearman")))
plot(z1$Avg_nonNA[, 177], z1$Avg_all[, 177])


# Relationship between Zscores and count NAs
summary(sapply(1:ncol(z1$Avg_all), function(x) cor(z1$NA_mat[, x], z1$Avg_all[, x], method = "spearman")))
summary(sapply(1:ncol(z1$Avg_nonNA), function(x) cor(z1$NA_mat[, x], z1$Avg_nonNA[, x], method = "spearman")))
which.max(sapply(1:ncol(z1$Avg_all), function(x) cor(z1$NA_mat[, x], z1$Avg_all[, x], method = "spearman")))
which.max(sapply(1:ncol(z1$Avg_nonNA), function(x) cor(z1$NA_mat[, x], z1$Avg_nonNA[, x], method = "spearman")))
plot(z1$Avg_all[, 461], z1$Avg_nonNA[, 461])



# Single dataframe of all rankings for demo gene
gene <- "ASCL1"

rank_df <- data.frame(
  Symbol = rownames(rsr1),
  RSR1 = rsr1[, gene],
  RSR2 = rsr_full2[, gene],
  # RSR3 = rsr3[, gene],
  # RSR4 = rsr_full4[, gene],
  # RSR5 = rsr5[, gene],
  # RSR6 = rsr_full6[, gene],
  Z_all = z1$Avg_all[, gene],
  Z_nonNA = z1$Avg_nonNA[, gene]
)



# Examining individual cell type cors
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