library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

out_dir <- "/space/scratch/amorin/R_objects/cpm_vs_lognorm.RDS"
ids <- read.delim(file.path(sc_dir, "Batch_cpm.tsv"), header = FALSE)$V1
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
msr_mat <- readRDS(msr_mat_hg_path)

dat_log <- load_agg_mat_list(ids, genes = pc_hg$Symbol)
dat_cpm <- dat_log
# dat_cpm <- load_agg_mat_list(ids, genes = pc_hg$Symbol, pattern = "_RSR_allrank_cpm.tsv")
k <- 1000


check_k <- function(vec_sort, k) {
  
  if (vec_sort[k] == vec_sort[k - 1]) {
    tie_start <- head(sort(table(vec_sort), decreasing = TRUE))[1]
    k <- sum(vec_sort > as.numeric(names(tie_start)), na.rm = TRUE)
  }
  
  return(k)
}



topk_sort <- function(vec, k) {
  vec_sort <- sort(vec, decreasing = TRUE)
  k <- check_k(vec_sort, k)
  names(vec_sort[1:k])
}



topk_intersect <- function(vec1, vec2) length(intersect(vec1, vec2))



cor_l <- lapply(1:length(dat_log), function(id) {
  
  lapply(pc_hg$Symbol, function(gene) {
    
    if (msr_mat[gene, id] == 0) {
      return(NA)
    }
    
    WGCNA::cor(dat_log[[id]][, gene], dat_cpm[[id]][, gene], method = "spearman")
  })
})



topk_l <- lapply(1:length(dat_log), function(id) {
  
  lapply(pc_hg$Symbol, function(gene) {
    
    if (msr_mat[gene, id] == 0) {
      return(NA)
    }
    
    vec1 <- topk_sort(dat_log[[id]][, gene], k = 1000)
    vec2 <- topk_sort(dat_cpm[[id]][, gene], k = 1000)
    topk_intersect(vec1, vec2)
  })
})


saveRDS(list(Cor = cor_l, Topk = topk_l), file = out_dir)
