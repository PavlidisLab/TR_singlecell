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

force_resave <- FALSE

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




# head(avg_coexpr_hg, 30)
# head(avg_coexpr_mm, 30)
# 
# tail(avg_coexpr_hg, 30)
# tail(avg_coexpr_mm, 30)



# Generate a list of gene-gene cor for every cell type in each dataset
# ------------------------------------------------------------------------------



# topcor_path <- "/space/scratch/amorin/R_objects/top_cor_pair_list.RDS"
topcor_path <- "/space/scratch/amorin/R_objects/top_cor_pair_list_RPS12-RPS27A_Rps24-Rps20.RDS"
btmcor_path <- "/space/scratch/amorin/R_objects/bottom_cor_pair_list.RDS"


if (!file.exists(topcor_path) || force_resave) {

  topcor_l <- list(
    Human = get_all_cor_l(ids_hg, gene1 = "FOS", gene2 = "EGR1"),
    Mouse = get_all_cor_l(ids_mm, gene1 = "Mt1", gene2 = "Mt2"))
  
  saveRDS(topcor_l, topcor_path)
  
} else {
  
  topcor_l <- readRDS(topcor_path)
  
}



if (!file.exists(btmcor_path) || force_resave) {
  
  btmcor_l <- list(
    Human = get_all_cor_l(ids_hg, gene1 = "GAPDH", gene2 = "ZNF644"),
    Mouse = get_all_cor_l(ids_mm, gene1 = "Fth1", gene2 = "Rbm6"))
  
  saveRDS(btmcor_l, btmcor_path)
  
} else {
  
  btmcor_l <- readRDS(btmcor_path)
  
}



stop()


# Plotting
# ------------------------------------------------------------------------------


# Akin to a forest plot, show the spread of cell-type correlations for a given
# gene pair across all experiments


# TODO: finalize plotting and explicitly save out final examples
# cor_forest_plot <- function(cor_l) { }


species <- "Human"
# species <- "Mouse"


cor_df <- do.call(
  rbind, 
  lapply(names(topcor_l[[species]]), function(x) data.frame(Cor = topcor_l[[species]][[x]], ID = x))
)


cor_summ <- do.call(rbind, lapply(topcor_l[[species]], summary)) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "ID") %>% 
  arrange(Median) %>% 
  mutate(ID = factor(ID, levels = unique(ID)))



# Gene-gene cor of every cell type for each dataset as boxplot + points

px <- ggplot(cor_df, aes(x = Cor, y = reorder(ID, Cor, FUN = median))) +
  geom_point(alpha = 0.4, shape = 21) +
  geom_boxplot(outlier.shape = NA, coef = 0, fill = "slategrey") +
  geom_vline(xintercept = 0, colour = "lightgrey") +
  xlab("Pearson's correlation across cell types") +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 10),
        plot.margin = margin(c(10, 20, 10, 10)))



ggsave(px, height = 14, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, "top_cor_pair_RPS12-RPS27A.png"))


# Cell-type scatter plots + cor heatmap of representative experiment and gene pair



demo_id <- "GSE212606Human"
gene1 <- "MAP1A" 
gene2 <- "MAP1B"
dat <- load_dat_list(demo_id)[[1]]
mat <- dat$Mat
meta <- dat$Meta

# All cell type cors
ct_cors <- all_celltype_cor(mat, meta, gene1, gene2)

# All cell type scatter plots
ct_scatter_l <- all_celltype_scatter(mat, meta, gene1, gene2)
ct_scatter_plot <- plot_grid(plotlist = ct_scatter_l)

# Heatmap of cors as vertical/horizontal
ct_heatmap_h <- cor_heatmap(ct_cors)
ct_heatmap_v <- cor_heatmap(t(ct_cors))


# Demonstrating the top, median, and bottom cor scatterplots
ct_scatter_plot_3 <- plot_grid(
  plotlist = c(ct_scatter_l[names(ct_cors)[1]],
               ct_scatter_l[names(ct_cors)[ceiling(length(ct_cors) / 2)]],
               ct_scatter_l[names(ct_cors)[length(ct_cors)]]),
  ncol = 3)


# Combining top/med/bottom scatter with heatmap
px3_h <- plot_grid(ct_scatter_plot_3, ct_heatmap_h$gtable, nrow = 2, rel_heights = c(2, 1))
px3_v <- plot_grid(ct_scatter_plot_3, ct_heatmap_v$gtable, ncol = 2, rel_widths = c(0.5, 0.25))



# Inspecting dataset that had the lowest cor


dat <- load_dat_list("GSE212606Human")[[1]]
mat <- dat$Mat
meta <- dat$Meta
agg <- load_agg_mat_list("GSE212606Human", genes = pc_hg$Symbol)[[1]]
agg_df <- mat_to_df(agg, symmetric = TRUE)
head(arrange(agg_df, desc(Value)))
ct_cors <- all_celltype_cor(mat, meta, "MAP1A", "MAP1B")



# dat <- load_dat_list("GSE165003")[[1]]
# mat <- dat$Mat
# meta <- dat$Meta
# agg <- load_agg_mat_list("GSE165003", genes = pc_mm$Symbol)[[1]]
# agg_df <- mat_to_df(agg, symmetric = TRUE)
# head(arrange(agg_df, desc(Value)))
# ct_cors <- all_celltype_cor(mat, meta, "Srsf2", "Mettl23")
