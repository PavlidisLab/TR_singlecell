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



# Get the average aggregate coexpression across all datasets
# ------------------------------------------------------------------------------


get_avg_coexpr <- function(ids, genes) {
  
  avg_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
  colnames(avg_mat) <- rownames(avg_mat) <- genes
  
  for (id in ids) {
    mat <- load_agg_mat_list(x, genes = genes)[[1]]
    avg_mat <- avg_mat + mat
  }
  
  avg_mat <- avg_mat / length(ids)
  
  avg_df <- mat_to_df(avg_mat, symmetric = TRUE, value_name = "Avg_coexpr")
  avg_df <- arrange(avg_df, desc(Avg_coexpr))
  
  return(avg_df)
  
}




if (!file.exists(avg_coexpr_hg_path)) {
  avg_coexpr_hg <- get_avg_coexpr(ids_hg, genes = pc_hg$Symbol)
  saveRDS(avg_coexpr_hg, avg_coexpr_hg_path)
} else {
  avg_coexpr_hg <- readRDS(avg_coexpr_hg_path)
}



if (!file.exists(avg_coexpr_mm_path)) {
  avg_coexpr_mm <- get_avg_coexpr(ids_mm, genes = pc_mm$Symbol)
  saveRDS(avg_coexpr_mm, avg_coexpr_mm_path)
} else {
  avg_coexpr_mm <- readRDS(avg_coexpr_mm_path)
}



# Most common gene in top pairs: ribosomal genes in general top heavy. 
# Human: HSPA1A most common, always paired with HSPA1B or HSP90AA1. RPL13 also
# common, often with RPRS18.
# Mouse: Rps29 most common, paired with Rps28, Rps39... 
# Also note GSE160512: Ttll3-Arpc4 neighbouring genes described by refseq as 
# having run-through transcription 
# https://pubmed.ncbi.nlm.nih.gov/20967262/  
# https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=10093
# ------------------------------------------------------------------------------


head(avg_coexpr_hg)
head(avg_coexpr_mm)


# Visualize top cor pairs 
# ------------------------------------------------------------------------------


# Akin to a forest plot, show the spread of cell-type correlations for a given
# gene pair across all experiments



outfile <- "/space/scratch/amorin/R_objects/top_cor_pair_l.RDS"


if (!file.exists(outfile)) {
  
  cor_l <- list(
    Human = get_all_cor_l(ids_hg, gene1 = "HSPA1A", gene2 = "HSPA1B"),
    Mouse = get_all_cor_l(ids_mm, gene1 = "Rpl13", gene2 = "Rpl32"))
  
  saveRDS(cor_l, outfile)
  
} else {
  
  cor_l <- readRDS(outfile)
  
}



# TODO: finalize plotting
# cor_forest_plot <- function(cor_l) { }


cor_df <- do.call(
  rbind, 
  lapply(names(cor_l), function(x) data.frame(Cor = cor_l[[x]], ID = x))
)


cor_summ <- do.call(rbind, lapply(cor_l, summary)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  arrange(Median) %>% 
  mutate(ID = factor(ID, levels = unique(ID)))



# Cor of every cell type for each dataset as points

ggplot(cor_df, aes(x = Cor, y = reorder(ID, Cor, FUN = median))) +
  geom_point(alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  xlab("Pearson's correlation across cell types") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10))


ggplot(cor_df, aes(x = Cor, y = reorder(ID, Cor, FUN = median))) +
  geom_point() +
  xlab("Pearson's correlation across cell types") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10))



# Cor of every cell type for each dataset as a bar (IRQ or min/max?)

ggplot(cor_summ) +
  geom_errorbarh(aes(xmin = `Min.`, xmax = `Max.`, y = ID)) +
  xlab("Pearson's correlation across cell types") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10))



ggplot(cor_summ) +
  geom_errorbarh(aes(xmin = `1st Qu.`, xmax = `3rd Qu.`, y = ID)) +
  xlab("Pearson's correlation across cell types") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10))



# Cell-type scatter plots + cor heatmap of representative experiment and gene pair



demo_id <- "GSE163252"  # "GSE160512"
gene1 <- "Rpl13"       # "Ttll3"
gene2 <- "Rpl32"       # "Arpc4"
dat <- load_dat_list(demo_id)[[1]]
mat <- dat$Mat
meta <- dat$Meta

ct_cors <- all_celltype_cor(mat, meta, gene1, gene2)




plot_scatter <- function(df, cell_type) {
  
  ggplot(df, aes(x = df[, 1], y = df[, 2])) + 
    geom_point(size = 2.5, shape = 21, fill = "#756bb1", alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
    xlab(colnames(df)[1]) +
    ylab(colnames(df)[2]) +
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
                                 min_cell = 20) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta),
            c(gene1, gene2) %in% rownames(mat))
  
  cts <- unique(meta$Cell_type)
  
  plot_l <- lapply(cts, function(ct) {
    ct_mat <- t(mat[c(gene1, gene2), filter(meta, Cell_type == ct)$ID])
    plot_scatter(data.frame(ct_mat), cell_type = ct)
  })
  
  names(plot_l) <- cts
  plot_l <- plot_l[!is.na(plot_l)]
  return(plot_l)
}



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
px3_v <- plot_grid(ct_scatter_plot_3, ct_heatmap_v$gtable, ncol = 2, rel_widths = c(0.5, 0.25))
