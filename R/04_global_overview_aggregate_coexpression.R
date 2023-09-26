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

force_resave <- TRUE

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)



# Get the average aggregate coexpression across all datasets
# Highest average RSR: 
# Ribosomal genes in general top heavy for both species. Similar for AP-1
# (FOS/B, JUN/B), EGR1 (FOS/B, IER2), DUSP1 (FOS/B, JUN/B), CALR-HSPA5
# ------------------------------------------------------------------------------


get_avg_coexpr <- function(ids, genes) {
  
  avg_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
  colnames(avg_mat) <- rownames(avg_mat) <- genes
  
  for (id in ids) {
    mat <- load_agg_mat_list(id, genes = genes)[[1]]
    avg_mat <- avg_mat + mat
  }
  
  avg_mat <- avg_mat / length(ids)
  
  avg_df <- mat_to_df(avg_mat, symmetric = TRUE, value_name = "Avg_coexpr")
  avg_df <- arrange(avg_df, desc(Avg_coexpr))
  
  return(avg_df)
  
}




if (!file.exists(avg_coexpr_hg_path) || force_resave) {
  avg_coexpr_hg <- get_avg_coexpr(ids_hg, genes = pc_hg$Symbol)
  saveRDS(avg_coexpr_hg, avg_coexpr_hg_path)
} else {
  avg_coexpr_hg <- readRDS(avg_coexpr_hg_path)
}



if (!file.exists(avg_coexpr_mm_path) || force_resave) {
  avg_coexpr_mm <- get_avg_coexpr(ids_mm, genes = pc_mm$Symbol)
  saveRDS(avg_coexpr_mm, avg_coexpr_mm_path)
} else {
  avg_coexpr_mm <- readRDS(avg_coexpr_mm_path)
}




head(avg_coexpr_hg, 30)
head(avg_coexpr_mm, 30)



# Generate a list of gene-gene cor for every cell type in each dataset
# ------------------------------------------------------------------------------



outfile <- "/space/scratch/amorin/R_objects/top_cor_pair_list.RDS"


if (!file.exists(outfile) || force_resave) {
  
  cor_l <- list(
    Human = get_all_cor_l(ids_hg, gene1 = "FOS", gene2 = "FOSB"),
    Mouse = get_all_cor_l(ids_mm, gene1 = "Fos", gene2 = "Fosb"))
  
  saveRDS(cor_l, outfile)
  
} else {
  
  cor_l <- readRDS(outfile)
  
}


stop()

# Plotting
# ------------------------------------------------------------------------------


# Akin to a forest plot, show the spread of cell-type correlations for a given
# gene pair across all experiments


# TODO: finalize plotting
# cor_forest_plot <- function(cor_l) { }


# species <- "Human"
species <- "Mouse"


cor_df <- do.call(
  rbind, 
  lapply(names(cor_l[[species]]), function(x) data.frame(Cor = cor_l[[species]][[x]], ID = x))
)


cor_summ <- do.call(rbind, lapply(cor_l[[species]], summary)) %>%
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



demo_id <- "GSE212606Human"  
gene1 <- "RPL3" 
gene2 <- "RPL13"
dat <- load_dat_list(demo_id)[[1]]
mat <- dat$Mat
meta <- dat$Meta



ct_cors <- all_celltype_cor(mat, meta, gene1, gene2)


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



# Inspecting lowest cor

# dat <- load_dat_list("GSE212606Human")[[1]]
# mat <- dat$Mat
# meta <- dat$Meta
# agg <- load_agg_mat_list("GSE212606Human", genes = pc_hg$Symbol)[[1]]
# agg_df <- mat_to_df(agg, symmetric = TRUE)
# head(arrange(agg_df, desc(Value)))
# ct_cors <- all_celltype_cor(mat, meta, "MAP1A", "CLU")



# dat <- load_dat_list("GSE165003")[[1]]
# mat <- dat$Mat
# meta <- dat$Meta
# agg <- load_agg_mat_list("GSE165003", genes = pc_mm$Symbol)[[1]]
# agg_df <- mat_to_df(agg, symmetric = TRUE)
# head(arrange(agg_df, desc(Value)))
# ct_cors <- all_celltype_cor(mat, meta, "Srsf2", "Mettl23")
