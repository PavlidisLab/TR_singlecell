## Examine differences between results generated from log norm vs CPM norm
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(Seurat)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_cpm_hg <- readRDS(msr_mat_hg_path)
msr_cpm_mm <- readRDS(msr_mat_mm_path)
msr_ln_hg <- readRDS("/space/scratch/amorin/R_objects/binary_measurement_matrix_hg_lognorm.RDS")
msr_ln_mm <- readRDS("/space/scratch/amorin/R_objects/binary_measurement_matrix_mm_lognorm.RDS")

stopifnot(identical(rownames(msr_cpm_hg), rownames(msr_ln_hg)),
          identical(colnames(msr_cpm_hg), colnames(msr_ln_hg)))

stopifnot(identical(rownames(msr_cpm_mm), rownames(msr_ln_mm)),
          identical(colnames(msr_cpm_mm), colnames(msr_ln_mm)))


# 1. Differences in binary measurement
# ------------------------------------------------------------------------------


get_msr_diff <- function(msr_cpm, msr_ln) {
  
  msr_diff <- msr_cpm - msr_ln 
  msr_diff_df <- data.frame(which(msr_diff != 0, arr.ind = TRUE))
  msr_diff_df$row <- rownames(msr_cpm)[msr_diff_df$row]
  msr_diff_df$col <- colnames(msr_cpm)[msr_diff_df$col]
  rownames(msr_diff_df) <- NULL
  colnames(msr_diff_df) <- c("Gene", "Experiment")
  msr_diff_df$Diff <- msr_diff[msr_diff != 0]
  
  return(msr_diff_df)
}



msr_diff_hg <- get_msr_diff(msr_cpm_hg, msr_ln_hg)
msr_diff_mm <- get_msr_diff(msr_cpm_mm, msr_ln_mm)


# Human: SIGLEC5 most diff at 12 datasets, and a handful of experiments have 2
# genes that differ. In all cases the genes are measured in log norm not CPM.

sort(table(msr_diff_hg$Experiment))
sort(table(msr_diff_hg$Gene))

# Mouse: Arhgap26 and Fam220a most diff at 4 datasets, and GSE207848 has 12 
# genes that differ. In all cases the genes are measured in log norm not CPM.

sort(table(msr_diff_mm$Experiment))
sort(table(msr_diff_mm$Gene))


# Load and inpsect data

id <- "GSE135922"
gene <- "ACTB"

dat_cpm <- load_dat_list(id)[[1]]
dat_ln <- load_dat_list(id, suffix = "_clean_mat_and_meta.RDS")[[1]]
stopifnot(identical(dim(dat_cpm$Mat), dim(dat_ln$Mat)))

diff_df <- data.frame(
  Symbol = rownames(dat_cpm$Mat),
  Diff = rowSums(dat_cpm$Mat) - rowSums(dat_ln$Mat)
)

plot(dat_cpm$Mat[gene, ], dat_ln$Mat[gene, ])
plot(rowSums(dat_cpm$Mat), rowSums(dat_ln$Mat))
cor(rowSums(dat_cpm$Mat), rowSums(dat_ln$Mat), method = "spearman")
view(data.frame(CPM = dat_cpm$Mat[gene, ], LN = dat_ln$Mat[gene, ]))
