## This script generates metrics of similarity between a query gene vector from 
## one aggregated correlation matrix and all other gene vectors from another
## network. The main idea idea is to find the rank position of the gene of
## interest with itself in another network.
## TODO: function argument of main call
## TODO: collapse main function (general structure of nested loop to mat)
## TODO: lapply call of main either made into a function or reduced to relevant only
## TODO: process/save/load function
## TODO: double check k+1
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

#
sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
pc_ortho <- read.delim(pc_ortho_path)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))

#
genes_hg <- rownames(agg_hg[[1]])
genes_mm <- rownames(agg_mm[[1]])

# TODO: this step should be done upstream
agg_hg <- lapply(agg_hg, function(x) x[genes_hg, genes_hg])
agg_mm <- lapply(agg_mm, function(x) x[genes_mm, genes_mm])



# Functions
# ------------------------------------------------------------------------------


#

which_ranked_cor <- function(rank_vec, 
                             rank_mat, 
                             gene,
                             cor_method = "spearman",
                             ncores = 1) {
  
  rank_cor <-  WGCNA::cor(x = rank_vec, 
                          y = rank_mat,
                          method = cor_method,
                          nThreads = ncores)
  
  rank_cor <- sort(rank_cor[1, ], decreasing = TRUE)
  rank_ix <- which(names(rank_cor) == gene)
  
  return(list(Rank = rank_ix, List = rank_cor))
}




#

which_ranked_auprc <- function(rank_vec, 
                               rank_mat, 
                               gene, 
                               k = 1000,
                               ncores = 1) {
  
  genes <- rownames(rank_mat)
  x <- sort(rank_vec, decreasing = TRUE)
  x <- x[names(x) != gene]
  
  auprc_l <- mclapply(genes, function(g) {
    
    y <- sort(rank_mat[, g], decreasing = TRUE)[1:(k+1)]
    y <- y[names(y) != g]

    rank_df <- data.frame(
      Symbol = names(x),
      Rank = x,
      Label = names(x) %in% names(y)
    )
    
    if (all(rank_df$Label)) return(1)
    if (all(!rank_df$Label)) return(0)
    
    get_au_perf(rank_df, label_col = "Label", score_col = "Rank", measure = "AUPRC")
    
  }, mc.cores = ncores)
  names(auprc_l) <- genes
  
  rank_auprc <- sort(unlist(auprc_l), decreasing = TRUE)
  rank_ix <- which(names(rank_auprc) == gene)
  
  return(list(Rank = rank_ix, List = rank_auprc))
}



#

which_ranked_intersect <- function(rank_vec,
                                   rank_mat,
                                   gene,
                                   k = 1000,
                                   ncores = 1) {
  
  genes <- rownames(rank_mat)
  x <- sort(rank_vec, decreasing = TRUE)[1:(k+1)]
  x <- x[names(x) != gene]
  
  intersect_l <- mclapply(genes, function(i) {
    
    y <- sort(rank_mat[, i], decreasing = TRUE)[1:(k+1)]
    y <- y[names(y) != i]

    length(intersect(names(x), names(y)))
    
  }, mc.cores = ncores)
  names(intersect_l) <- genes
  
  rank_intersect <- sort(unlist(intersect_l), decreasing = TRUE)
  rank_ix <- which(names(rank_intersect) == gene)
  
  return(list(Rank = rank_ix, List = rank_intersect))
}




# Return an n x n matrix where n is equal to the count of aggregate matrices
# in agg_l. For the given gene

rank_similarity_matrix <- function(agg_l, 
                                   gene, 
                                   msr,  # AUPRC|Intersect|Cor
                                   k = 1000,
                                   ncores = 1) {
  
  ids <- names(agg_l)
  genes <- rownames(agg_l[[1]])
  
  stopifnot(length(ids) > 0, gene %in% genes)
  
  rank_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  rownames(rank_mat) <- colnames(rank_mat) <- ids
  
  for (i in ids) {
    for (j in ids) {
      if (i == j) {
        next
      }
      
      rank_ix <- if (msr == "AUPRC") {
        
        which_ranked_auprc(
          rank_vec = agg_l[[i]][, gene],
          rank_mat = agg_l[[j]],
          gene = gene,
          k = k,
          ncores = ncores)$Rank
        
      } else if (msr == "Intersect") {
        
        which_ranked_intersect(
          rank_vec = agg_l[[i]][, gene],
          rank_mat = agg_l[[j]],
          gene = gene,
          k = k,
          ncores = ncores)$Rank
        
      } else {
        
        which_ranked_cor(
          rank_vec = agg_l[[i]][, gene],
          rank_mat = agg_l[[j]],
          gene = gene,
          ncores = ncores)$Rank
      }

      rank_mat[i, j] <- rank_ix
    }
  }
  
  return(rank_mat)
}




similarity_matrix <- function(agg_l, 
                              gene, 
                              msr,  # AUPRC|Intersect
                              k = 1000,
                              ncores = 1) {
  
  ids <- names(agg_l)
  genes <- rownames(agg_l[[1]])
  
  stopifnot(length(ids) > 0, gene %in% genes)
  
  sim_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  rownames(sim_mat) <- colnames(sim_mat) <- ids
  
  for (i in ids) {
    for (j in ids) {
      if (i == j) {
        next
      }
      
      x <- sort(agg_l[[i]][, gene], decreasing = TRUE)
      x <- x[names(x) != gene]
      
      y <- sort(agg_l[[j]][, gene], decreasing = TRUE)[1:(k+1)]
      y <- y[names(y) != i]
      
      sim <- if (msr == "AUPRC") {
        
        rank_df <- data.frame(
          Symbol = names(x),
          Rank = x,
          Label = names(x) %in% names(y)
        )
        
        if (all(rank_df$Label)) return(1)
        if (all(!rank_df$Label)) return(0)
        
        get_au_perf(rank_df, label_col = "Label", score_col = "Rank", measure = "AUPRC")
        
      } else if (msr == "Intersect") {
       
        x <- x[1:k]
        length(intersect(names(x), names(y)))
        
      }
        
      sim_mat[i, j] <- sim
      
    }
  }
  
  return(sim_mat)
}






# Inspecting single genes
# ------------------------------------------------------------------------------


gene_hg <- "RPL3"  # RPL3  RPL19
gene_mm <- "Rpl3"


# Cor

gene_mat_hg <- gene_vec_to_mat(agg_hg, gene_hg)
gene_cor_hg <- cor(gene_mat_hg, method = "spearman", use = "pairwise.complete.obs")
cor_heatmap(gene_cor_hg)

gene_mat_mm <- gene_vec_to_mat(agg_mm, gene_mm)
gene_cor_mm <- cor(gene_mat_mm, method = "spearman", use = "pairwise.complete.obs")
cor_heatmap(gene_cor_mm)


# AUPRC

gene_auprc_hg <- similarity_matrix(agg_l = agg_hg, gene = gene_hg, k = 1000, msr = "AUPRC")
gene_auprc_mm <- similarity_matrix(agg_l = agg_mm, gene = gene_mm, k = 1000, msr = "AUPRC")

pheatmap(gene_auprc_hg,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_col = "black",
         color = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 20,
         cellwidth = 50,
         cellheight = 50)


# Topk

gene_topk_hg <- similarity_matrix(agg_l = agg_hg, gene = gene_hg, k = 1000, msr = "Intersect")
gene_topk_mm <- similarity_matrix(agg_l = agg_mm, gene = gene_mm, k = 1000, msr = "Intersect")




# Comparing metrics


comp_df <- data.frame(
  Cor = c(
    mat_to_df(gene_cor_hg, symmetric = TRUE)$Value,
    mat_to_df(gene_cor_mm, symmetric = TRUE)$Value
  ),
  Topk = c(
    mat_to_df(gene_topk_hg, symmetric = TRUE)$Value,
    mat_to_df(gene_topk_mm, symmetric = TRUE)$Value
  ),
  AUPRC = c(
    mat_to_df(gene_auprc_hg, symmetric = TRUE)$Value,
    mat_to_df(gene_auprc_mm, symmetric = TRUE)$Value
  )
)



ggplot(comp_df, aes(x = Cor, y = AUPRC)) +
  geom_point(shape = 19, size = 3) +
  ggtitle(gene_hg) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


cor(comp_df)


# Calling which rank over numerous genes (ribosomal + TF) - prohibitively slow
# ------------------------------------------------------------------------------


tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

test_genes_hg <- intersect(c(str_to_upper(tfs), ribo_genes$Symbol_hg), genes_hg)
test_genes_mm <- intersect(c(str_to_title(tfs), ribo_genes$Symbol_mm), genes_mm)

outfile_int_hg <- "/space/scratch/amorin/R_objects/03-06-2023_aggregate_vector_comparison_intersect_human.RDS"
outfile_auprc_hg <- "/space/scratch/amorin/R_objects/03-06-2023_aggregate_vector_comparison_auprc_human.RDS"
outfile_cor_hg <- "/space/scratch/amorin/R_objects/03-06-2023_aggregate_vector_comparison_cor_human.RDS"
outfile_int_mm <- "/space/scratch/amorin/R_objects/03-06-2023_aggregate_vector_comparison_intersect_mouse.RDS"
outfile_auprc_mm <- "/space/scratch/amorin/R_objects/03-06-2023_aggregate_vector_comparison_auprc_mouse.RDS"
outfile_cor_mm <- "/space/scratch/amorin/R_objects/03-06-2023_aggregate_vector_comparison_cor_mouse.RDS"


# Intersects


if (!file.exists(outfile_int_hg)) {
  
  int_l_hg <- lapply(test_genes_hg, function(gene) {
    
    rank_similarity_matrix(agg_l = agg_hg,
                           gene = gene,
                           msr = "Intersect",
                           k = 1000,
                           ncores = ncore)
  })
  names(int_l_hg) <- test_genes_hg
  saveRDS(int_l_hg, outfile_int_hg)
} else {
  int_l_hg <- readRDS(outfile_int_hg)
}




if (!file.exists(outfile_int_mm)) {
  
  int_l_mm <- lapply(test_genes_mm, function(gene) {
    
    rank_similarity_matrix(agg_l = agg_mm,
                           gene = gene,
                           msr = "Intersect",
                           k = 1000,
                           ncores = ncore)
  })
  names(int_l_mm) <- test_genes_mm
  saveRDS(int_l_mm, outfile_int_mm)
} else {
  int_l_mm <- readRDS(outfile_int_mm)
}



# AUPRCs


# if (!file.exists(outfile_auprc_hg)) {
#   
#   auprc_l_hg <- lapply(test_genes_hg, function(gene) {
#     
#     rank_similarity_matrix(agg_l = agg_hg,
#                            gene = gene,
#                            msr = "AUPRC",
#                            k = 1000,
#                            ncores = ncore)
#   })
#   names(auprc_l_hg) <- test_genes_hg
#   saveRDS(auprc_l_hg, outfile_auprc_hg)
# } 
# 
# 
# if (!file.exists(outfile_auprc_mm)) {
#   
#   auprc_l_mm <- lapply(test_genes_mm, function(gene) {
#     
#     rank_similarity_matrix(agg_l = agg_mm,
#                            gene = gene,
#                            msr = "AUPRC",
#                            k = 1000,
#                            ncores = ncore)
#   })
#   names(auprc_l_mm) <- test_genes_mm
#   saveRDS(auprc_l_mm, outfile_auprc_mm)
# }
# 
# 
# 
# # Cors
# 
# 
# if (!file.exists(outfile_cor_hg)) {
#   
#   cor_l_hg <- lapply(test_genes_hg, function(gene) {
#     
#     rank_similarity_matrix(agg_l = agg_hg,
#                            gene = gene,
#                            msr = "Cor",
#                            ncores = ncore)
#   })
#   names(cor_l_hg) <- test_genes_hg
#   saveRDS(cor_l_hg, outfile_cor_hg)
# }
# 
# 
# 
# if (!file.exists(outfile_cor_mm)) {
#   
#   cor_l_mm <- lapply(test_genes_mm, function(gene) {
#     
#     rank_similarity_matrix(agg_l = agg_mm,
#                            gene = gene,
#                            msr = "Cor",
#                            ncores = ncore)
#   })
#   names(cor_l_mm) <- test_genes_mm
#   saveRDS(cor_l_mm, outfile_cor_mm)
# }


##


int_l_hg[[gene_hg]]



tt1 <- which_ranked_intersect(rank_vec = agg_hg[["GSE216019"]][, gene_hg], 
                              rank_mat = agg_hg[["GSE180928"]],
                              gene = gene_hg, 
                              k = 1000, 
                              ncores = ncore)
tt1$Rank

data.frame(Topk = tt1$List) %>% 
  ggplot(aes(x = Topk)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = tt1$List[gene_hg], col = "red") +
  ylab("Count") +
  xlab("Topk intersect") +
  ggtitle("GSE216019 RPL3 in GSE180928") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))



# Collapse the ribosomal gene rankings: First standardize each rank (0 worst, 
# 1 best), and then take the median of [exp1, exp2] ranks across ribo genes.

med_standard_rank <- function(mat_l, max_rank) {
  
  # because current ranks are 1=best to max_rank=worst, reverse for standard
  mat_l <- lapply(mat_l, function(x) (max_rank + 1 - x) / max_rank)
  med_mat <- apply(simplify2array(mat_l), 1:2, median)
  return(med_mat)
}


med_ribo_srank_hg <- med_standard_rank(int_l_hg[ribo_genes$Symbol_hg], max_rank = length(genes_hg))
med_ribo_srank_mm <- med_standard_rank(int_l_mm[ribo_genes$Symbol_mm], max_rank = length(genes_mm))


med_ribo_srank_df <- rbind(mat_to_df(diag_to_na(med_ribo_srank_hg)),
                           mat_to_df(diag_to_na(med_ribo_srank_mm)))


hist(med_ribo_srank_df$Value, breaks = 100)


# Inspecting example of RPL19 being the top hit between datasets


gene_topk_rank_hg <- which_ranked_intersect(
  rank_vec = agg_hg$GSE202352[, "RPL19"],
  rank_mat = agg_hg$GSE180928,
  gene = "RPL19",
  k = 1000,
  ncores = ncore)


gene_topk_rank_hg$List %>% head


# Inspect TFs


med_tf_srank_hg <- med_standard_rank(int_l_hg[tfs_hg], max_rank = length(genes_hg))
med_tf_srank_mm <- med_standard_rank(int_l_mm[tfs_mm], max_rank = length(genes_mm))


med_tf_srank_df <- rbind(mat_to_df(diag_to_na(med_tf_srank_hg)),
                         mat_to_df(diag_to_na(med_tf_srank_mm)))


hist(med_tf_srank_df$Value, breaks = 10)
