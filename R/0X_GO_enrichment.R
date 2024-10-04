# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(parallel)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
source("R/00_config.R")


go_coexpr_hg_path <- "/space/scratch/amorin/R_objects/coexpr_goenrich_hg.RDS"
go_coexpr_mm_path <- "/space/scratch/amorin/R_objects/coexpr_goenrich_mm.RDS"
go_coexpr_ortho_path <- "/space/scratch/amorin/R_objects/coexpr_goenrich_ortho.RDS"



# Assumes dat_l as named list of TF coexpression rankings
# Assumes summary df has Rank_aggr_coexpr and Symbol in column names


go_enrich_list <- function(dat_l, species, topn = 50, ncores = 1) {
  
  tfs <- names(dat_l)
  
  universe <- dat_l[[1]]$Symbol
  
  gene_entrez <- bitr(universe, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = species)
  
  
  go_l <- mclapply(tfs, function(tf) {
    
    top_genes <- slice_min(dat_l[[tf]], Rank_aggr_coexpr, n = topn)[["Symbol"]]
    top_genes <- filter(gene_entrez, SYMBOL %in% top_genes)
    
    enrichGO(gene = top_genes$ENTREZID,
             universe = gene_entrez$ENTREZID,
             OrgDb = species,
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.05,
             readable = TRUE)

    
  }, mc.cores = ncores)
  names(go_l) <- tfs
  
  return(go_l)
}



# Note that ortho (ranking placing equal weight on mouse and human aggregate)
# uses human symbols

if (!file.exists(go_coexpr_hg_path)) {
  coexpr_hg <- readRDS(rank_tf_hg_path)
  go_coexpr_hg <- go_enrich_list(dat_l = coexpr_hg, species = "org.Hs.eg.db", ncores = ncore)
  saveRDS(go_hg, go_coexpr_hg_path)
}
  

if (!file.exists(go_coexpr_mm_path)) {
  coexpr_mm <- readRDS(rank_tf_mm_path)
  go_coexpr_mm <- go_enrich_list(dat_l = coexpr_mm, species = "org.Mm.eg.db", ncores = ncore)
  saveRDS(go_mm, go_coexpr_mm_path)
}


if (!file.exists(go_coexpr_ortho_path)) {
  coexpr_ortho <- readRDS(rank_tf_ortho_path)
  go_coexpr_ortho <- go_enrich_list(dat_l = coexpr_ortho, species = "org.Hs.eg.db", ncores = ncore)
  saveRDS(go_ortho, go_coexpr_ortho_path)
}
