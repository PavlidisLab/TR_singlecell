## TODO
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
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Load aggregate matrix into list
agg_hg <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_mm <- load_agg_mat_list(ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

stopifnot(all(unlist(lapply(agg_hg, function(x) identical(rownames(x), pc_hg$Symbol)))))
stopifnot(all(unlist(lapply(agg_mm, function(x) identical(rownames(x), pc_mm$Symbol)))))

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Saved list RDS of the summaries
out_hg <- "/space/scratch/amorin/R_objects/05-07-2023_TF_summary_human.RDS"
out_mm <- "/space/scratch/amorin/R_objects/05-07-2023_TF_summary_mouse.RDS"


# For each gene in genes, generate a summary dataframe that includes:
# 1) The gene ranking of coexpressed partners (average aggregate rank)
# 2) The count of times the gene pair was mutually measured
# 3) The count of NAs (no mutual measurement)
# 4) For the best correlation partner, what was its best rank across datasets?
# 5) And in what dataset did it achieve this best rank?

all_rank_summary <- function(agg_l, msr_mat, genes) {
  
  summ_l <- lapply(genes, function(x) {
    
    message(paste(x))
    
    gene_mat <- gene_vec_to_mat(agg_l, x)
    gene_mat[x, ] <- NA
    gene_mat <- gene_mat[, which(msr_mat[x, ] == 1), drop = FALSE]
    
    if (length(gene_mat) == 0) {
      return(NA)
    }
    
    colrank_gene_mat <- colrank_mat(gene_mat, ties_arg = "min")
    
    avg_rsr <- rowMeans(gene_mat, na.rm = TRUE)
    avg_colrank <- rowMeans(colrank_gene_mat, na.rm = TRUE)
    
    rank_rsr <- rank(-avg_rsr, ties.method = "min")
    rank_colrank <- rank(avg_colrank, ties.method = "min")
    
    best_rank <- apply(colrank_gene_mat, 1, min)
    
    data.frame(
      Symbol = rownames(gene_mat),
      Avg_RSR = avg_rsr,
      Avg_colrank = avg_colrank,
      Rank_RSR = rank_rsr,
      Rank_colrank = rank_colrank,
      Best_rank = best_rank)
  })
  
  names(summ_l) <- genes
  
  return(summ_l)
}




if (!file.exists(out_hg)) {
  summ_hg <- all_rank_summary(agg_l = agg_hg, msr_mat = msr_hg, genes = tfs_hg$Symbol)
  saveRDS(summ_hg, out_hg)
} else {
  summ_hg <- readRDS(out_hg)
}



if (!file.exists(out_mm)) {
  summ_mm <- all_rank_summary(agg_l = agg_mm, msr_mat = msr_mm, genes = tfs_mm$Symbol)
  saveRDS(summ_mm, out_mm)
} else {
  summ_mm <- readRDS(out_mm)
}
