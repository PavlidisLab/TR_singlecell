library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

out_dir <- "/space/scratch/amorin/R_objects/cpm_vs_lognorm.RDS"
ids <- read.delim(file.path(sc_dir, "Batch_cpm.tsv"), header = FALSE)$V1
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
msr_mat <- readRDS(msr_mat_hg_path)

dat_log <- load_agg_mat_list(ids, genes = pc_hg$Symbol)
dat_cpm <- load_agg_mat_list(ids, genes = pc_hg$Symbol, pattern = "_RSR_allrank_cpm.tsv")
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
  
  l <- lapply(pc_hg$Symbol, function(gene) {
    
    if (msr_mat[gene, id] == 0) {
      return(NA)
    }
    
    WGCNA::cor(dat_log[[id]][, gene], dat_cpm[[id]][, gene], method = "spearman")
    
  })
  
  names(l) <- pc_hg$Symbol
  return(l)
})



topk_l <- lapply(1:length(dat_log), function(id) {
  
  l <- lapply(pc_hg$Symbol, function(gene) {
    
    if (msr_mat[gene, id] == 0) {
      return(NA)
    }
    
    vec1 <- topk_sort(dat_log[[id]][, gene], k = 1000)
    vec2 <- topk_sort(dat_cpm[[id]][, gene], k = 1000)
    topk_intersect(vec1, vec2)
    
  })
  
  names(l) <- pc_hg$Symbol
  return(l)
})


saveRDS(list(Cor = cor_l, Topk = topk_l), file = out_dir)


dat_l <- readRDS(out_dir)
cor_l <- dat_l$Cor
topk_l <- dat_l$Topk



plot_hist <- function(vec, xlab) {
  
  data.frame(Value = vec) %>% 
    ggplot(., aes(x = Value)) +
    geom_histogram(bins = 100) +
    xlab(xlab) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20))
}


plot_hist(vec = unlist(cor_l[[1]]), xlab = "Scor")



hist_cor <- lapply(cor_l, function(x) plot_hist(vec = unlist(x), xlab = "Scor"))
cowplot::plot_grid(plotlist = hist_cor)

hist_topk <- lapply(topk_l, function(x) plot_hist(vec = unlist(x), xlab = "Topk"))
cowplot::plot_grid(plotlist = hist_topk)
