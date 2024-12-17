## 
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
library(aggtools)
library(ggrepel)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 200


# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# Protein coding genes and TFs
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)


file <- "/space/scratch/amorin/R_objects/TRsc/GSE132042_10X_smartseq_comparison.RDS"

# Measurement matrices used for filtering when a gene was never expressed
# msr_mm <- readRDS(msr_mat_mm_path)


dat_ss <- load_dat_list("GSE132042SmartSeq2")[[1]]
dat_10x <- load_dat_list("GSE132042")[[1]]

# table(dat_ss$Meta$assay)
# table(dat_10x$Meta$assay)


common_counts <- inner_join(
  count(dat_ss$Meta, Cell_type),
  count(dat_10x$Meta, Cell_type),
  by = "Cell_type", suffix = c("_smartseqV2", "_10XV2")
) %>% 
  mutate(Cell_type = as.character(Cell_type)) %>% 
  filter(n_smartseqV2 >= 100 & n_10XV2 >= 100)


common_cts <- common_counts$Cell_type


meta_ss <- filter(dat_ss$Meta, Cell_type %in% common_cts)
mat_ss <- dat_ss$Mat[, meta_ss$ID]


meta_10x <- filter(dat_10x$Meta, Cell_type %in% common_cts)
mat_10x <- dat_10x$Mat[, meta_10x$ID]


if (!file.exists(file)) {
  
  agg_ss <- aggtools::aggr_coexpr_single_dataset(
    mat = mat_ss,
    meta = meta_ss,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = "allrank"
  )
  
  agg_10x <- aggtools::aggr_coexpr_single_dataset(
    mat = mat_10x,
    meta = meta_10x,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = "allrank"
  )
  
  agg_l <- list(Agg_SS = agg_ss, Agg_10x = agg_10x)
  
  saveRDS(agg_l, file)
  
} else {
  
  agg_l <- readRDS(file)
  
}

stop()



# 
# ------------------------------------------------------------------------------


mat_ss_log <- log2(mat_ss + 1)
mat_10x_log <- log2(mat_10x + 1)


avg_expr_mat <- function(mat, meta, cts, tfs) {
  
  avg_l <- lapply(cts, function(ct) {
    
    ct_ids <- filter(meta, Cell_type == ct)$ID
    ct_mat <- mat[tfs, ct_ids]
    rowMeans(ct_mat)
    
  })
  
  avg_mat <- as.matrix(do.call(cbind, avg_l))
  colnames(avg_mat) <- common_cts
  rownames(avg_mat) <- tfs
  
  return(avg_mat)
}


avg_ss <- avg_expr_mat(mat_ss_log, meta_ss, common_cts, tfs_mm$Symbol)
avg_10x <- avg_expr_mat(mat_10x_log, meta_10x, common_cts, tfs_mm$Symbol)
gc()


all0_ss <- which(rowSums(avg_ss) == 0)
all0_10x <- which(rowSums(avg_10x) == 0)
keep <- setdiff(tfs_mm$Symbol, intersect(names(all0_ss), names(all0_10x)))


avg_ss <- avg_ss[keep, ]
avg_10x <- avg_10x[keep, ]


avg_topk <- pair_colwise_topk(avg_ss, avg_10x, k = k, ncores = ncore)
avg_scor <- pair_colwise_cor(avg_ss, avg_10x, cor_method = "spearman", ncores = ncore)


avg_topk_shuffle <- unlist(lapply(1:10, function(iter) {
  pair_shuffle_topk(avg_ss, avg_10x, k = k, ncores = ncore)
}))


plot(density(avg_topk_shuffle))
abline(v = mean(avg_topk))


colnames(avg_ss) <- paste0(common_cts, "_SS")
colnames(avg_10x) <- paste0(common_cts, "_10x")
all_topk_mat <- colwise_topk_intersect(cbind(avg_ss, avg_10x), k = k)
all_topk_df <- mat_to_df(all_topk_mat)

all_scor_mat <- cor(cbind(avg_ss, avg_10x), method = "spearman")
all_scor_df <- mat_to_df(all_scor_mat)

