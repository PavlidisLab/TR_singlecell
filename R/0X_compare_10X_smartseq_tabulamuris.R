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

