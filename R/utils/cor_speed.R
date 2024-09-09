source("R/00_config.R")
source("R/utils/dev_functions.R")
library(WGCNA)
library(microbenchmark)


pc_df <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
id <- "GSE199317"
ct <- "Microglia"
dat <- readRDS(file.path(amat_dir, id, paste0(id, "_clean_mat_and_meta_CPM.RDS")))
meta <- dat$Meta
mat <- dat$Mat


mat <- t(mat[, filter(meta, Cell_type == ct)[["ID"]]])
mat_dense <- as.matrix(mat)


res <- microbenchmark(
  F1 = calc_sparse_correlation(mat, cor_method = "pearson"),
  F2 = cor(mat_dense),
  times = 5
)


print(res)