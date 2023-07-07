## Sample genes from experiments and calculate the topk intersect across 
## experiments to generate a null. Save out the result as a list.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 1000 
n_samps <- 1000

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Saved list RDS of the sampled topk overlaps
out_hg <- "/space/scratch/amorin/R_objects/06-07-2023_sampled_topk_intesect_human.RDS"
out_mm <- "/space/scratch/amorin/R_objects/06-07-2023_sampled_topk_intesect_mouse.RDS"



# NOTE: The following implementation is very slow when calling multiple times,
# as each dataset is loaded each call. This could be made faster by generating
# the sampled genes first and loading the sampled genes for each dataset once. 
# This can quickly into memory issues however as each sampled matrix needs to be
# held in memory.


# For the given ids, sample a single gene that is measured from each experiment.
# Return a dataframe of the unique pairs and their topk intersect.

sample_topk_intersect <- function(ids, genes, msr_mat, k = 1000) {
  
  # Sample one gene that is measured in each data set
  sample_genes <- lapply(ids, function(x) {
    msr_genes <- msr_mat[msr_mat[, x] == 1, x]
    names(sample(msr_genes, 1))
  })
  
  # Load the sampled gene for each experiment and bind into a matrix
  sample_mat <- lapply(1:length(ids), function(x) {
    load_agg_mat_list(ids = ids[x], sub_genes = sample_genes[[x]], genes = genes)[[1]]
  })
  sample_mat <- do.call(cbind, sample_mat)
  
  # Get topk overlap between sampled genes
  sample_topk <- colwise_topk_intersect(sample_mat)
  sample_df <- mat_to_df(sample_topk, symmetric = TRUE, value_name = "Topk")
  
  return(sample_df)
}





# Use all experiments, or choose a subset in which a a given gene is measured
# gene_hg <- "ASCL1"
# gene_mm <- "Ascl1"
# ids_hg <- names(which(msr_hg[gene_hg, ] == 1))
# ids_mm <- names(which(msr_mm[gene_mm, ] == 1))


set.seed(5)


if (!file.exists(out_hg)) {
  
  topk_hg <- lapply(1:n_samps, function(x) {
    
    message(paste("Human sample #", x, Sys.time()))
    
    sample_topk_intersect(ids = ids_hg,
                          genes = pc_hg$Symbol,
                          msr_mat = msr_hg,
                          k = k)
  })
  
  saveRDS(topk_hg, out_hg)
  
} else {
  
  topk_hg <- readRDS(out_hg)
  
}



if (!file.exists(out_mm)) {
  
  topk_mm <- lapply(1:n_samps, function(x) {
    
    message(paste("Mouse sample #", x, Sys.time()))
    
    sample_topk_intersect(ids = ids_mm,
                          genes = pc_mm$Symbol,
                          msr_mat = msr_mm,
                          k = k)
  })
  
  saveRDS(topk_mm, out_mm)
  
} else {
  
  topk_mm <- readRDS(out_mm)
  
}



# Summarize

med_hg <- lapply(topk_hg, function(x) median(x$Topk))
med_mm <- lapply(topk_mm, function(x) median(x$Topk))
hist(unlist(med_hg), breaks = 20)
hist(unlist(med_mm), breaks = 20)


sd_hg <- lapply(topk_hg, function(x) sd(x$Topk))
sd_mm <- lapply(topk_mm, function(x) median(x$Topk))
hist(unlist(sd_hg), breaks = 20)
hist(unlist(sd_mm), breaks = 20)
