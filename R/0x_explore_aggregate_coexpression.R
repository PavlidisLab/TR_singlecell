##
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

# 
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

#
pc_ortho <- read.delim(pc_ortho_path)

# Ranked targets from paper
evidence_l <- readRDS(evidence_path)
tfs_mm <- names(evidence_l$Mouse)
tfs_hg <- names(evidence_l$Human)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID
ids <- union(ids_hg, ids_mm)

# Load aggregate matrix into list
# agg_files <- list.files(amat_dir, pattern = "_RSR_allrank.RDS", recursive = TRUE, full.names = TRUE)
# agg_l <- lapply(agg_files, readRDS)
agg_hg <- lapply(ids_hg, function(x) readRDS(file.path(amat_dir, x, paste0(x, "_RSR_allrank.RDS"))))
names(agg_hg) <- ids_hg
agg_hg <- lapply(agg_hg, lowertri_to_symm)
agg_mm <- lapply(ids_mm, function(x) readRDS(file.path(amat_dir, x, paste0(x, "_RSR_allrank.RDS"))))
names(agg_mm) <- ids_mm
agg_mm <- lapply(agg_mm, lowertri_to_symm)

gc()

# Load mat and meta
# dat_files <- list.files(amat_dir, pattern = "_clean_mat_and_meta.RDS", recursive = TRUE, full.names = TRUE)
# dat_l <- lapply(dat_files, readRDS)
# dat_mm <- lapply(ids_mm, function(x) readRDS(file.path(amat_dir, x, paste0(x, "_clean_mat_and_meta.RDS")))
# names(dat_mm) <- ids_mm
# dat_hg <- lapply(ids_hg, function(x) readRDS(file.path(amat_dir, x, paste0(x, "_clean_mat_and_meta.RDS")))
# names(dat_hg) <- ids_hg


genes_hg <- rownames(agg_hg[[1]])
genes_mm <- rownames(agg_mm[[1]])



gene_vec_to_mat <- function(agg_l, gene) {
  do.call(cbind, lapply(agg_l, function(x) x[gene, ]))
}


gene <- "RPL13"


gene_mat <- gene_vec_to_mat(agg_hg, gene)


# Correlation of gene's rankings across studies
gene_rank_cor <- cor(gene_mat, method = "spearman", use = "pairwise.complete.obs")


# For the given gene, retrieve the rank position of its similarity for the 
# other studies


rank_similarity_matrix <- function(agg_l, gene, ncores = 1) {
  
  ids <- names(agg_l)
  genes <- rownames(agg_l[[1]])
  
  stopifnot(length(ids) > 0, gene %in% genes)
  
  rank_cor_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  rownames(rank_cor_mat) <- colnames(rank_cor_mat) <- ids
  
  for (i in ids) {
    for (j in ids) {
      if (i == j) {
        next
      }
      rank_cor <-  WGCNA::cor(x = agg_l[[i]][, gene], y = agg_l[[j]], nThreads = ncores) 
      rank_cor <- sort(rank_cor[1, ], decreasing = TRUE)
      rank_cor_mat[i, j] <- which(names(rank_cor) == gene)
    }
  }
  
  return(rank_cor_mat)
}



rmat_hg <- rank_similarity_matrix(agg_hg, gene, ncores = ncore)
# rmat_mm <- rank_similarity_matrix(agg_mm, gene, ncores = ncore)



# TODO:

get_au_perf <- function(rank_df, label_col, measure = NULL) {
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("AUROC", "AUPRC"))
  
  # positives/has evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(rank_df[[label_col]]))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  
  if (measure == "AUROC") {
    perf <- ROCR::performance(pred, measure = "auc")@y.values[[1]]
  } else {
    perf <- ROCR::performance(pred, measure = "aucpr")@y.values[[1]]
  }
  
  return(perf)
}




### Isolating inner comparison of two experiments and one gene

k <- 1000  # length(genes_hg)
# x <- sort(agg_hg$GSE180928[, gene], decreasing = TRUE)
x <- sort(agg_hg$GSE180928[, gene], decreasing = TRUE)[1:(k+1)]
x <- x[names(x) != gene]
y <- sort(agg_hg$GSE216019[, gene], decreasing = TRUE)[1:(k+1)]
y <- y[names(y) != gene]


rank_df <- data.frame(
  Symbol = names(x),
  Rank = 1:length(x),
  Label = names(x) %in% names(y)
  # Label = factor(names(x) %in% names(y), levels = c(TRUE, FALSE))
)

perf_df <- get_perf_df(rank_df, label_col = "Label", measure = "PR")
perf <- get_au_perf(rank_df, label_col = "Label", measure = "AUPRC")
##


### Isolating comparison of one gene list versus the top k of all other genes
auprc_l <- lapply(genes_hg, function(i) {
  
  x <- sort(agg_hg$GSE180928[, gene], decreasing = TRUE)
  x <- x[names(x) != gene]
  y <- sort(agg_hg$GSE216019[, i], decreasing = TRUE)[1:(k+1)]
  y <- y[names(y) != gene]
  
  rank_df <- data.frame(
    Symbol = names(x),
    Rank = 1:length(x),
    Label = names(x) %in% names(y)
  )
  
  perf <- get_au_perf(rank_df, label_col = "Label", measure = "AUPRC")
  
})

names(auprc_l) <- genes_hg
auprc_l <- sort(unlist(auprc_l), decreasing = TRUE)
which(names(auprc_l) == gene)
hist(unlist(auprc_l), breaks = 100)
abline(v = auprc_l[gene], col = "red")
##


### Isolating comparison of top k of one gene list versus the top k of all other genes
auprc_l2 <- lapply(genes_hg, function(i) {
  
  x <- sort(agg_hg$GSE180928[, gene], decreasing = TRUE)[1:(k+1)]
  x <- x[names(x) != gene]
  y <- sort(agg_hg$GSE216019[, i], decreasing = TRUE)[1:(k+1)]
  y <- y[names(y) != gene]
  
  rank_df <- data.frame(
    Symbol = names(x),
    Rank = 1:length(x),
    Label = factor(names(x) %in% names(y), levels = c(TRUE, FALSE))
  )
  
  # if (all(!rank_df$Label)) return(0)
  
  perf <- get_au_perf(rank_df, label_col = "Label", measure = "AUPRC")
  perf <- get_perf_df(rank_df, label_col = "Label", measure = "PR")
  
})

names(auprc_l2) <- genes_hg
auprc_l2 <- sort(unlist(auprc_l2), decreasing = TRUE)
which(names(auprc_l2) == gene)
hist(unlist(auprc_l2), breaks = 100)
abline(v = auprc_l2[gene], col = "red")
##


## Compare the two lists
cor(auprc_l, auprc_l2[names(auprc_l)], method = "spearman")



## Into a function
rank_auprc_matrix <- function(agg_l, gene, k) {
  
  ids <- names(agg_l)
  genes <- rownames(agg_l[[1]])
  
  stopifnot(length(ids) > 0, gene %in% genes)
  
  rank_auprc_mat <- matrix(1, nrow = length(agg_l), ncol = length(agg_l))
  rownames(rank_auprc_mat) <- colnames(rank_auprc_mat) <- ids
  
  for (i in ids) {
    for (j in ids) {
      if (i == j) {
        next
      }
      
      x <- sort(agg_hg[[i]][, gene], decreasing = TRUE)[1:(k+1)]
      x <- x[names(x) != gene]
      y <- sort(agg_hg$GSE216019[, i], decreasing = TRUE)[1:(k+1)]
      y <- y[names(y) != gene]
      
      rank_df <- data.frame(
        Symbol = names(x),
        Rank = 1:length(x),
        Label = names(x) %in% names(y)
      )
      
      if (all(!rank_df$Label)) return(0)
      
      perf <- get_au_perf(rank_df, label_col = "Label", measure = "AUPRC")
    }
  }
  
  return(rank_cor_mat)
  
}
##

# Using ribosomal genes as a sanity check
lapply(agg_hg, function(x) head(sort(x[, "RPL3"], decreasing = TRUE)))


# List of max pairs for each data set

max_pair_l <- function(agg_l, ncores = 1) {
  mclapply(agg_l, function(x) {
    diag(x) <- NA
    which(x == max(x, na.rm = TRUE), arr.ind = TRUE)
  }, mc.cores = ncores)
}


max_pair_hg <- max_pair_l(agg_hg, ncore)
max_pair_mm <- max_pair_l(agg_mm, ncore)


dat_files <- list.files(amat_dir, pattern = "_clean_mat_and_meta.RDS", recursive = TRUE, full.names = TRUE)
dat <- readRDS("/space/scratch/amorin/TR_singlecell/Tabula_Sapiens/Tabula_Sapiens_clean_mat_and_meta.RDS")
mat <- dat$Mat
meta <- dat$Meta
gene1 <- "RPL13"
gene2 <- "RPS18"
ct_cors <- all_celltype_cor(mat, meta, gene1, gene2)


plot_scatter <- function(df, cell_type) {
  
  ggplot(df, aes(x = df[, 1], y = df[, 2])) + 
    geom_point(size = 2.5, shape = 21, fill = "#756bb1", alpha = 0.3) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
    xlab(colnames(df)[1]) +
    ylab(colnames(df)[2]) +
    # ggtitle(paste0(x, ": ", round(cor_l[[x]][tf, gene], 3))) +
    ggtitle(cell_type) +
    theme_classic() +
    theme(plot.title = element_text(size = 15),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15))
}


cor_heatmap <- function(cor_vec, cell_size = 10) {
  
  pheatmap(t(cor_vec),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           cellwidth = cell_size,
           cellheight = cell_size)
}



all_celltype_scatter <- function(mat,
                                 meta,
                                 gene1,
                                 gene2,
                                 min_cell = 20,
                                 cor_method = "pearson") {
  
  
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta),
            c(gene1, gene2) %in% rownames(mat))
  
  cts <- unique(meta$Cell_type)
  
  plot_l <- lapply(cts, function(ct) {
    
    ct_mat <- t(mat[c(gene1, gene2), filter(meta, Cell_type == ct)$ID])
    
    # if (sum(ct_mat[, 1] != 0) < min_cell || sum(ct_mat[, 2] != 0) < min_cell) {
    #   return(NA)
    # } else {
    # WGCNA::cor(ct_mat[, gene1], ct_mat[, gene2], method = cor_method)
    # }
    
    p <- plot_scatter(data.frame(ct_mat), cell_type = ct)
    
  })
  names(plot_l) <- cts
  
  plot_l <- plot_l[!is.na(plot_l)]
  
  # plot_l <- cowplot::plot_grid(plotlist = plot_l)
  
  # cor_vec <- sort(unlist(cor_l), decreasing = TRUE)
  
  return(plot_l)
}


# ct_cors <- all_celltype_cor(mat, meta, gene1, gene2)
ct_scatter_l <- all_celltype_scatter(mat, meta, gene1, gene2)
ct_scatter_plot <- plot_grid(plotlist = ct_scatter_l)
ct_heatmap_h <- cor_heatmap(ct_cors)
ct_heatmap_v <- cor_heatmap(t(ct_cors))

px1_h <- plot_grid(ct_scatter_l[[names(ct_cors[1])]], ct_heatmap$gtable, nrow = 2)
px2_h <- plot_grid(ct_scatter_plot, ct_heatmap$gtable, nrow = 2)
px1_v <-  plot_grid(ct_scatter_l[[names(ct_cors[1])]], ct_heatmap_v$gtable, ncol = 2)
px2_v <-  plot_grid(ct_scatter_l[[names(ct_cors[1])]], ct_heatmap_v$gtable, ncol = 2)


ct_scatter_plot_3 <- plot_grid(
  plotlist = c(ct_scatter_l[names(ct_cors)[1]],
               ct_scatter_l[names(ct_cors)[floor(length(ct_cors) / 2)]],
               ct_scatter_l[names(ct_cors)[length(ct_cors)]]),
  ncol = 1)


px3_v <- plot_grid(ct_scatter_plot_3, ct_heatmap_v$gtable, ncol = 2)
px3_h <- plot_grid(ct_scatter_plot_3, ct_heatmap_h$gtable, nrow = 2, rel_heights = c(2, 1))
