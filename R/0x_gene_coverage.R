## Inspect NA pairs of TFs and target genes
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID[1:3]
ids_mm <- filter(sc_meta, Species == "Mouse")$ID[1:3]



# TODO

load_na_mat_list <- function(ids,
                             dir = "/space/scratch/amorin/TR_singlecell/",
                             sub_genes = NULL) {
  
  
  if (!is.null(sub_genes)) sub_genes <- c("V1", sub_genes)
  
  mat_l <- lapply(ids, function(x) {
    path <- file.path(amat_dir, x, paste0(x, "_NA_mat.tsv"))
    dat <- fread(path, sep = "\t", select = sub_genes)
    mat <- as.matrix(dat[, -1, drop = FALSE])
    rownames(mat) <- colnames(mat) <- dat$V1
    mat
  })
  
  names(mat_l) <- ids
  gc(verbose = FALSE)
  
  return(mat_l)
}



na_hg <- load_na_mat_list(ids_hg)
na_mm <- load_na_mat_list(ids_mm)


# Example of describing a single gene coverage over all experiments

gene_hg <- "ASCL1"

gene_mat_hg <- gene_vec_to_mat(na_hg, gene_hg)

# If the gene itself has NAs equal to the count of cell types for the dataset,
# then that gene is not expressed at all in the data set

gene_mat_hg[gene_hg, ]
n_cts <- filter(sc_meta, ID %in% ids_hg)$N_celltypes



