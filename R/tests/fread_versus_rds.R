library(data.table)
library(microbenchmark)
source("R/00_config.R")

id <- "GSE180928"
rds_path <- file.path(amat_dir, id, paste0(id, "_RSR_allrank.RDS"))
tsv_path <- file.path(amat_dir, id, paste0(id, "_RSR_allrank.tsv"))


sub_genes <- c("ASCL1", "RUNX1")

# dat_rds <- readRDS(rds_path)
# dat_tsv <- as.matrix(fread(tsv_path, sep = "\t"))
# dat_tsv_sub <- as.matrix(fread(tsv_path, sep = "\t", select = sub_genes))


res <- microbenchmark(
  F1 = readRDS(rds_path)[, sub_genes],
  F2 = as.matrix(fread(tsv_path, sep = "\t"))[, sub_genes],
  F3 = as.matrix(fread(tsv_path, sep = "\t", select = sub_genes)),
  times = 5
)
