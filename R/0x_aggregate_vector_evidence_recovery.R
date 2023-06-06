## TODO: point to missing genes
## TODO: missing genes excluded from rank df makes a difference?
## TODO: common process for curated and top k evidence for get_rank_df
## TODO: remove tf?
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
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
ids <- union(ids_hg, ids_mm)

# Load aggregate matrix into list
agg_hg <- lapply(ids_hg, function(x) lowertri_to_symm(readRDS(file.path(amat_dir, x, paste0(x, "_RSR_allrank.RDS")))))
names(agg_hg) <- ids_hg
agg_mm <- lapply(ids_mm, function(x) lowertri_to_symm(readRDS(file.path(amat_dir, x, paste0(x, "_RSR_allrank.RDS")))))
names(agg_mm) <- ids_mm
gc()

# Ranked targets from paper
evidence_l <- readRDS(evidence_path)
tfs_mm <- names(evidence_l$Mouse)
tfs_hg <- names(evidence_l$Human)

#
genes_hg <- rownames(agg_hg[[1]])
genes_mm <- rownames(agg_mm[[1]])



# Functions
# ------------------------------------------------------------------------------


# TODO: smarter

get_rank_df <- function(gene_vec, 
                        evidence_df, 
                        evidence_col,
                        k = 1000) {
  
  
  rank_df <- if (evidence_col == "Curated_target") {
    
    data.frame(Score = gene_vec) %>%
      rownames_to_column(var = "Symbol") %>%
      left_join(evidence_df[, c("Symbol", evidence_col)], by = "Symbol") %>%
      # filter(Symbol != tf) %>%
      arrange(!!sym(evidence_col)) %>%
      mutate(Label = Symbol %in% filter(evidence_df, Curated_target)$Symbol) %>%
      arrange(desc(Score))
    
  } else {
    
    data.frame(Score = gene_vec) %>%
      rownames_to_column(var = "Symbol") %>%
      left_join(evidence_df[, c("Symbol", evidence_col)], by = "Symbol") %>%
      # filter(Symbol != tf) %>%
      arrange(!!sym(evidence_col)) %>%
      mutate(Label = c(rep(1, k), rep(0, nrow(.) - k))) %>%
      arrange(desc(Score))
    
  }
  
  return(rank_df)
}


# TODO:

get_tf_performance <- function(agg_l,
                               evidence_df,
                               evidence_col,
                               tf,
                               k = 1000) {
  
  # TODO:
  sc_tf_mat <- gene_vec_to_mat(agg_l, tf)
  sc_tf_mat <- cbind(sc_tf_mat, Aggregate_rank = rowMeans(sc_tf_mat))
  
  # TODO:
  auprc_l <- lapply(colnames(sc_tf_mat), function(x) {
    gene_vec <- sc_tf_mat[, x]
    rank_df <- get_rank_df(gene_vec, evidence, evidence_col, k)
    get_au_perf(rank_df, label_col = "Label", score_col = "Score", measure = "AUPRC")
  })
  names(auprc_l) <- colnames(sc_tf_mat)
  auprc_l <- unlist(auprc_l)
  
  return(auprc_l)
}



# TODO: 

get_all_performance <- function(agg_l,
                                evidence_df,
                                evidence_col,
                                genes,
                                k = 1000,
                                ncores = 1) {
  
  # TODO:
  
  perf_l <- lapply(agg_l, function(agg_mat) {
    
    # TODO:
    
    perf <- mclapply(genes, function(x) {
      
      gene_vec <- agg_mat[, x]
      rank_df <- get_rank_df(gene_vec, evidence_df, evidence_col, k)
      auprc <- get_au_perf(rank_df, label_col = "Label", score_col = "Score", measure = "AUPRC")
      topk <- sum(rank_df$Label[1:k])
      data.frame(AUPRC = auprc, Topk = topk)
      
    }, mc.cores = ncores)
    
    perf_df <- data.frame(Symbol = genes, do.call(rbind, perf)) %>%
      arrange(desc(AUPRC))
    
    return(perf_df)
    
  })
  
  names(perf_l) <- names(agg_l)
  return(perf_l)
}



# Individual versus aggregate recovery
# ------------------------------------------------------------------------------


# Human

tf_performance_hg <- lapply(tfs_hg, function(tf) {
  
  get_tf_performance(agg_l = agg_hg, 
                     evidence_df = evidence_l$Human[[tf]], 
                     evidence_col = "Rank_binding", 
                     tf = tf, 
                     k = 1000)
})
names(tf_performance_hg) <- tfs_hg



# Mouse

tf_performance_mm <- lapply(tfs_mm, function(tf) {
  
  get_tf_performance(agg_l = agg_mm, 
                     evidence_df = evidence_l$Mouse[[tf]], 
                     evidence_col = "Rank_binding", 
                     tf = tf, 
                     k = 1000)
})
names(tf_performance_mm) <- tfs_mm





# Every gene recovery of evidence (slow!)
# ------------------------------------------------------------------------------


outfile_hg1 <- "/space/scratch/amorin/R_objects/05-06-2023_coexpr_agg_recovery_integrated_evidence_human.RDS"
outfile_hg2 <- "/space/scratch/amorin/R_objects/05-06-2023_coexpr_agg_recovery_curated_evidence_human.RDS"
outfile_mm1 <- "/space/scratch/amorin/R_objects/05-06-2023_coexpr_agg_recovery_integrated_evidence_mouse.RDS"
outfile_mm2 <- "/space/scratch/amorin/R_objects/05-06-2023_coexpr_agg_recovery_curated_evidence_mouse.RDS"


# Human

if (!file.exists(outfile_hg1)) {
  
  hg1 <- lapply(tfs_hg, function(tf) {
    
    get_all_performance(agg_l = agg_hg,
                        evidence_df = evidence_l$Human[[tf]],
                        evidence_col = "Rank_integrated",
                        genes = genes_hg,
                        k = 1000,
                        ncores = ncore)
  })
  names(hg1) <- tfs_hg
  
  saveRDS(hg1, outfile_hg1)
}



if (!file.exists(outfile_hg2)) {
  
  hg2 <- lapply(tfs_hg, function(tf) {
    
    get_all_performance(agg_l = agg_hg,
                        evidence_df = evidence_l$Human[[tf]],
                        evidence_col = "Curated_target",
                        genes = genes_hg,
                        k = 1000,
                        ncores = ncore)
  })
  names(hg2) <- tfs_hg
  
  saveRDS(hg2, outfile_hg2)
}



# Mouse


if (!file.exists(outfile_mm1)) {
  
  mm1 <- lapply(tfs_mm, function(tf) {
    
    get_all_performance(agg_l = agg_mm,
                        evidence_df = evidence_l$Mouse[[tf]],
                        evidence_col = "Rank_integrated",
                        genes = genes_mm,
                        k = 1000,
                        ncores = ncore)
  })
  names(mm1) <- tfs_mm
  
  saveRDS(mm1, outfile_mm1)
}



if (!file.exists(outfile_mm2)) {
  
  mm2 <- lapply(tfs_mm, function(tf) {
    
    get_all_performance(agg_l = agg_mm,
                        evidence_df = evidence_l$Mouse[[tf]],
                        evidence_col = "Curated_target",
                        genes = genes_mm,
                        k = 1000,
                        ncores = ncore)
  })
  names(mm2) <- tfs_mm
  
  saveRDS(mm2, outfile_mm2)
}



# which_tf_auprc <- lapply(auprc_l3, function(x) which(x$Symbol == tf))
# 
# 
# i <- 5
# 
# hist(auprc_l3[[i]]$Topk, breaks = 100)
# abline(v = filter(auprc_l3[[i]], Symbol == tf)$Topk, col = "red")
# 
# hist(auprc_l3[[i]]$AUPRC, breaks = 100)
# abline(v = filter(auprc_l3[[i]], Symbol == tf)$AUPRC, col = "red")
# 
# 
# plot(auprc_l3$GSE180928$Topk, auprc_l3$GSE180928$AUPRC)
# cor(auprc_l3$GSE180928$Topk, auprc_l3$GSE180928$AUPRC, method = "spearman")


# Tally the genes whose aggregate expression vector was among the top k at 
# recovering evidence (also top k)

# tally_topk <- lapply(auprc_l3, function(x) x$Symbol[1:k]) %>% 
#   unlist() %>% 
#   table() %>% 
#   sort(decreasing = TRUE)




# Notice that while RUNX1 coexpression itself is not highly ranked, the genes
# that are highly ranked (by their coexpression recovry of RUNX1 targets) are
# also highly ranked RUNX1 targets...

# filter(evidence, Symbol %in% names(tally_topk)[1:100]

       