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

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Loading genes of interest
sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
pc_ortho <- read.delim(pc_ortho_path)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

# Genes to focus/subset on when loading aggregate coexpression matrices
subset_hg <- c(str_to_upper(tfs), "RPL3", "EGR1")   # NULL
subset_mm <- c(str_to_title(tfs), "Rpl3", "Egr1")   # NULL

# Load aggregate matrix into list
agg_hg <- load_agg_mat_list(ids = ids_hg, sub_genes = subset_hg)
agg_mm <- load_agg_mat_list(ids = ids_mm, sub_genes = subset_mm)

genes_hg <- rownames(agg_hg[[1]])
genes_mm <- rownames(agg_mm[[1]])

# TODO: this needs to be replaced with upstream ordering
agg_hg <- lapply(agg_hg, function(x) x[genes_hg, ])
agg_mm <- lapply(agg_mm, function(x) x[genes_mm, ])


stopifnot(all(unlist(lapply(agg_hg, function(x) identical(rownames(x), genes_hg)))))
stopifnot(all(unlist(lapply(agg_mm, function(x) identical(rownames(x), genes_mm)))))



# Inspecting single genes
# ------------------------------------------------------------------------------


gene_hg <- "ASCL1"  # RPL3  RPL19
gene_mm <- "Ascl1"  # Rpl3  Rpl19


# Cor

gene_mat_hg <- gene_vec_to_mat(agg_hg, gene_hg)
gene_cor_hg <- colwise_cor(gene_mat_hg)
cor_heatmap(gene_cor_hg)

gene_mat_mm <- gene_vec_to_mat(agg_mm, gene_mm)
gene_cor_mm <- colwise_cor(gene_mat_mm)
cor_heatmap(gene_cor_mm)


# AUPRC

gene_auprc_hg <- colwise_topk_auprc(gene_mat_hg)
gene_auprc_mm <- colwise_topk_auprc(gene_mat_mm)



# pheatmap(gene_auprc_hg,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          border_col = "black",
#          color = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
#          display_numbers = TRUE,
#          number_color = "black",
#          fontsize = 20,
#          cellwidth = 50,
#          cellheight = 50)


# Topk

gene_topk_hg <- colwise_topk_intersect(gene_mat_hg)
gene_topk_mm <- colwise_topk_intersect(gene_mat_mm)

# Jaccard of top and bottom k

gene_jacc_hg <- binarize_topk_btmk(gene_mat_hg) %>% colwise_jaccard()
gene_jacc_mm <- binarize_topk_btmk(gene_mat_mm) %>% colwise_jaccard()



# TODO: better
# Comparing metrics


comp_df <- data.frame(
  Scor = c(
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
  ),
  Jaccard = c(
    mat_to_df(gene_jacc_hg, symmetric = TRUE)$Value,
    mat_to_df(gene_jacc_mm, symmetric = TRUE)$Value
  )
)


comp_df <- cbind(comp_df,
                 rbind(
                   mat_to_df(gene_jacc_hg, symmetric = TRUE)[, c("Row", "Col")],
                   mat_to_df(gene_jacc_mm, symmetric = TRUE)[, c("Row", "Col")])
)



# TODO: Replace hacky NA filter with formal binary detection matrix

na_hg <- gene_cor_hg %>% 
  diag_to_na() %>% 
  apply(., 2, function(x) any(x == 1))

na_hg <- names(na_hg[!is.na(na_hg)])


na_mm <- gene_cor_mm %>% 
  diag_to_na() %>% 
  apply(., 2, function(x) any(x == 1))

na_mm <- names(na_mm[!is.na(na_mm)])


comp_df2 <- filter(comp_df, !(Row %in% c(na_hg, na_mm) | Col %in% c(na_hg, na_mm)))


# ggplot(comp_df2, aes(x = AUPRC, y = Topk)) +
ggplot(comp_df2, aes(x = Scor, y = Topk)) +
  geom_point(shape = 19, size = 3) +
  ggtitle(gene_hg) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


cor(select_if(comp_df2, is.numeric))



rank_mat <- readRDS("~/scratch/R_objects/27-06-2023/ASCL1.RDS")


rank_mat <- rank_mat[setdiff(rownames(rank_mat), na_hg), setdiff(colnames(rank_mat), na_hg)]
diag(rank_mat) <- NA


pheatmap(rank_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_col = "black",
         color = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
         display_numbers = TRUE,
         number_format = "%1.0f",
         number_color = "black",
         na_col = "black",
         fontsize = 22,
         cellwidth = 50,
         cellheight = 50,
         height = 18,
         width = 18,
         filename = "tt.png")



F1 <- query_gene_rank_topk(query_vec = agg_hg$NowickiOsuch2023[, gene_hg],
                           subject_mat = load_agg_mat_list("Ravindra2021")[[1]],
                           gene = gene_hg,
                           ncores = ncore)


data.frame(Topk = F1$Topk) %>% 
  ggplot(aes(x = Topk)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = F1$Topk[gene_hg], col = "red") +
  ylab("Count of genes") +
  xlab("Topk intersect") +
  ggtitle("GSE216019 ASCL1 in GSE180928") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


head(sort(table(F1$Topk), decreasing = TRUE))



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
