## Organize examples like the most correlated gene pairs across datasets, such as
## the ribosomal L/S genes, and provide plots of the underlying correlations
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
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Using ribosomal genes as a sanity check
sribo_hg <- read.csv("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.csv("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))

# List of max pairs for each data set

max_pair_l <- function(agg_l, ncores = 1) {
  
  l <- mclapply(agg_l, function(x) {
    diag(x) <- NA
    ix <- which(x == max(x, na.rm = TRUE), arr.ind = TRUE)
    data.frame(Gene1 = rownames(ix)[1], Gene2 = rownames(ix)[2])
  }, mc.cores = ncores)
  
  # data.frame(Study = names(agg_l), do.call(rbind, l))
  do.call(rbind, l)
}


max_pair_hg <- max_pair_l(agg_hg, ncore)
max_pair_mm <- max_pair_l(agg_mm, ncore)


dat_files <- list.files(amat_dir, pattern = "_clean_mat_and_meta.RDS", recursive = TRUE, full.names = TRUE)
dat <- readRDS("/space/scratch/amorin/TR_singlecell/GSE180928/GSE180928_clean_mat_and_meta.RDS")
mat <- dat$Mat
meta <- dat$Meta
gene1 <- "HBA1"
gene2 <- "HBA2"
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


px1_h <- plot_grid(ct_scatter_l[[names(ct_cors[1])]], ct_heatmap_h$gtable, nrow = 2)
px2_h <- plot_grid(ct_scatter_plot, ct_heatmap_h$gtable, nrow = 2)
px1_v <-  plot_grid(ct_scatter_l[[names(ct_cors[1])]], ct_heatmap_v$gtable, ncol = 2)
px2_v <-  plot_grid(ct_scatter_plot, ct_heatmap_v$gtable, ncol = 2, rel_widths = c(2, 1))


ct_scatter_plot_3 <- plot_grid(
  plotlist = c(ct_scatter_l[names(ct_cors)[1]],
               ct_scatter_l[names(ct_cors)[ceiling(length(ct_cors) / 2)]],
               ct_scatter_l[names(ct_cors)[length(ct_cors)]]),
  ncol = 1)


px3_h <- plot_grid(ct_scatter_plot_3, ct_heatmap_h$gtable, nrow = 2, rel_heights = c(2, 1))
px3_v <- plot_grid(ct_scatter_plot_3, ct_heatmap_v$gtable, ncol = 2, rel_widths = c(2, 1))
