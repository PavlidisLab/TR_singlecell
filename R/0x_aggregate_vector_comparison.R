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

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Loading genes of interest
sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

# Genes to focus/subset on when loading aggregate coexpression matrices
subset_hg <- c(str_to_upper(tfs), "RPL3", "EGR1", "FOS", "FOSB")   # NULL
subset_mm <- c(str_to_title(tfs), "Rpl3", "Egr1", "Fos", "Fosb")   # NULL

# Load aggregate matrix into list
agg_hg <- load_agg_mat_list(ids = ids_hg, sub_genes = subset_hg, genes = pc_hg$Symbol)
agg_mm <- load_agg_mat_list(ids = ids_mm, sub_genes = subset_mm, genes = pc_mm$Symbol)

stopifnot(all(unlist(lapply(agg_hg, function(x) identical(rownames(x), pc_hg$Symbol)))))
stopifnot(all(unlist(lapply(agg_mm, function(x) identical(rownames(x), pc_mm$Symbol)))))

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)


# Inspecting single genes
# ------------------------------------------------------------------------------


gene_hg <- "ASCL1"
gene_mm <- str_to_title(gene_hg)
stopifnot(gene_mm %in% pc_mm$Symbol)

# Bind individual gene vectors into a matrix, removing datasets with no expression
gene_mat_hg <- gene_vec_to_mat(agg_hg, gene_hg)[, which(msr_hg[gene_hg, ] == 1)]
gene_mat_mm <- gene_vec_to_mat(agg_mm, gene_mm)[, which(msr_mm[gene_mm, ] == 1)]

# Cor
gene_cor_hg <- colwise_cor(gene_mat_hg)
gene_cor_mm <- colwise_cor(gene_mat_mm)

# Topk
gene_topk_hg <- colwise_topk_intersect(gene_mat_hg)
gene_topk_mm <- colwise_topk_intersect(gene_mat_mm)

# Jaccard of top and bottom k
gene_jacc_hg <- binarize_topk_btmk(gene_mat_hg) %>% colwise_jaccard()
gene_jacc_mm <- binarize_topk_btmk(gene_mat_mm) %>% colwise_jaccard()

# AUPRC - skipped: slower and not symmetric. Highly cor with topk
# gene_auprc_hg <- colwise_topk_auprc(gene_mat_hg)
# gene_auprc_mm <- colwise_topk_auprc(gene_mat_mm)


# Collect similarity metrics into a dataframe


gene_sim_df <-
  rbind(
    mat_to_df(gene_cor_hg, symmetric = TRUE, value_name = "Scor"),
    mat_to_df(gene_cor_mm, symmetric = TRUE, value_name = "Scor")
  ) %>%
  cbind(
    Topk = c(
      mat_to_df(gene_topk_hg, symmetric = TRUE)$Value,
      mat_to_df(gene_topk_mm, symmetric = TRUE)$Value
    ),
    Jaccard = c(
      mat_to_df(gene_jacc_hg, symmetric = TRUE)$Value,
      mat_to_df(gene_jacc_mm, symmetric = TRUE)$Value
    )
  )



ggplot(gene_sim_df, aes(x = Scor, y = Topk)) +
  geom_point(shape = 19, size = 3) +
  ggtitle(gene_hg) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


cor(select_if(gene_sim_df, is.numeric))




# Loading the matrix showing a gene's similarity rank (by topk) with itself
# in all pairs of datasets

# rank_mat <- readRDS(paste0("~/scratch/R_objects/27-06-2023/", gene_hg, ".RDS"))
rank_mat <- readRDS(paste0("~/scratch/R_objects/04-07-2023/", gene_hg, ".RDS"))

# TODO: remove when ranks regenerated
keep_exps <- intersect(rownames(rank_mat), names(which(msr_hg[gene_hg, ] == 1)))
rank_mat <- rank_mat[keep_exps, keep_exps]
gene_topk_hg <- gene_topk_hg[keep_exps, keep_exps]
gene_sim_df <- filter(gene_sim_df, Row %in% keep_exps & Col %in% keep_exps)

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
         width = 18)



# Pile up at the left of the histogram (lower/better ranks) suggests a signal,
# as random ranking would be uniform

hist(rank_mat, breaks = 20)
hist(replicate(length(rank_mat), sample(1:length(pc_hg$Symbol), 1)), breaks = 20)



# Inspecting top pairs: by topk magnitude, and by rank. The highest value/most
# similar pair across all experiment pairs may not actually have the highest
# gene similarity rank across all genes within the subject dataset. If there
# are multiple pairs with the best rank, take the highest topk magnitude pair.

best_value_id <- slice_max(gene_sim_df, Topk, n = 1)

best_rank_ix <- which(rank_mat == min(rank_mat, na.rm = TRUE), arr.ind = TRUE)
best_rank_ix <- best_rank_ix[which.max(gene_topk_hg[best_rank_ix, drop = FALSE]), ]
best_rank_id <- data.frame(Row = rownames(rank_mat)[best_rank_ix[1]], Col = colnames(rank_mat)[best_rank_ix[2]])


best_value <- query_gene_rank_topk(
  query_vec = agg_hg[[best_value_id$Row]][, gene_hg],
  subject_mat = load_agg_mat_list(best_value_id$Col, genes = pc_hg$Symbol)[[1]],
  gene = gene_hg,
  ncores = ncore)


best_rank <- query_gene_rank_topk(
  query_vec = agg_hg[[best_rank_id$Row]][, gene_hg],
  subject_mat = load_agg_mat_list(best_rank_id$Col, genes = pc_hg$Symbol)[[1]],
  gene = gene_hg,
  ncores = ncore)



pxa <- data.frame(Topk = best_value$Topk) %>% 
  ggplot(aes(x = Topk)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = best_value$Topk[gene_hg], col = "red") +
  ylab("Count of genes") +
  xlab("Topk intersect") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


pxb <- data.frame(Topk = best_rank$Topk) %>% 
  ggplot(aes(x = Topk)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = best_rank$Topk[gene_hg], col = "red") +
  ylab("Count of genes") +
  xlab("Topk intersect") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))



# Generate a null similarity matrix, pulling random genes from datasets with
# matched expression with the gene of interest


keep_exps <- intersect(rownames(rank_mat), names(which(msr_hg[gene_hg, ] == 1)))

set.seed(77)

sample_genes <- lapply(keep_exps, function(x) {
  msr_genes <- msr_hg[msr_hg[, x] == 1, x]
  names(sample(msr_genes, 1))
})


sample_mat <- lapply(1:length(keep_exps), function(x) {
  load_agg_mat_list(ids = keep_exps[x], sub_genes = sample_genes[[x]], genes = pc_hg$Symbol)[[1]]
})
sample_mat <- do.call(cbind, sample_mat)
colnames(sample_mat) <- paste0(keep_exps, "_", unlist(sample_genes))


sample_topk <- colwise_topk_intersect(sample_mat)
sample_df <- mat_to_df(sample_topk, symmetric = TRUE, value_name = "Topk")


plot_df <- data.frame(
  Topk = c(gene_sim_df$Topk, sample_df$Topk),
  Group = c(rep(gene_hg, length(gene_sim_df$Topk)), rep("Sampled", length(sample_df$Topk)))
)


ggplot(plot_df, aes(x = Topk, fill = Group)) +
  geom_density(alpha = 0.8) +
  # geom_histogram(position = "stack", bins = 100) + 
  theme_classic()




# cor_heatmap(gene_cor_hg)
# cor_heatmap(gene_cor_mm)

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


# Collapse the ribosomal gene rankings: First standardize each rank (0 worst, 
# 1 best), and then take the median of [exp1, exp2] ranks across ribo genes.

# med_standard_rank <- function(mat_l, max_rank) {
#   
#   # because current ranks are 1=best to max_rank=worst, reverse for standard
#   mat_l <- lapply(mat_l, function(x) (max_rank + 1 - x) / max_rank)
#   med_mat <- apply(simplify2array(mat_l), 1:2, median)
#   return(med_mat)
# }
# 
# 
# med_ribo_srank_hg <- med_standard_rank(int_l_hg[ribo_genes$Symbol_hg], max_rank = length(genes_hg))
# med_ribo_srank_mm <- med_standard_rank(int_l_mm[ribo_genes$Symbol_mm], max_rank = length(genes_mm))
# 
# 
# med_ribo_srank_df <- rbind(mat_to_df(diag_to_na(med_ribo_srank_hg)),
#                            mat_to_df(diag_to_na(med_ribo_srank_mm)))
# 
# 
# hist(med_ribo_srank_df$Value, breaks = 100)



# # Inspect TFs
# 
# 
# med_tf_srank_hg <- med_standard_rank(int_l_hg[tfs_hg], max_rank = length(genes_hg))
# med_tf_srank_mm <- med_standard_rank(int_l_mm[tfs_mm], max_rank = length(genes_mm))
# 
# 
# med_tf_srank_df <- rbind(mat_to_df(diag_to_na(med_tf_srank_hg)),
#                          mat_to_df(diag_to_na(med_tf_srank_mm)))
# 
# 
# hist(med_tf_srank_df$Value, breaks = 10)
