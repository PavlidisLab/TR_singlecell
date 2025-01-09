## The all correlation 'forest plots' used in the paper were exported
## using the CLI script 'save_all_corplot.R'. This process also saves out a list
## of gene-gene cors across cell types and datasets, which can be used for
## other plots. Here I am spelling that process out, which also provides a 
## sandbox for inspecting these correlations.
## -----------------------------------------------------------------------------

library(tidyverse)
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


# Using human ASCL1-DLL3 as example of the all cor plot + summarized info
# ------------------------------------------------------------------------------


species <- "Human"
gene1 <- "ASCL1"
gene2 <- "DLL3"
ids <- filter(sc_meta, Species == species)[["ID"]]

corplot_dir <- file.path(plot_dir, "All_corplots")
dir.create(corplot_dir, showWarnings = FALSE)
cor_path <- file.path(corplot_dir, paste0(species, "_", gene1, "_", gene2, ".RDS"))


# Generating list of correlations is slow.

if (!file.exists(cor_path)) {
  cor_l <- get_all_cor_l(ids = ids, gene1 = gene1, gene2 = gene2)
  saveRDS(cor_l, cor_path)
} else {
  cor_l <- readRDS(cor_path)
}


# The all cor plot
p1 <- all_corplot(cor_l)


# Organize a df of these correlations

cor_df <- lapply(names(cor_l), function(x) {
  data.frame(ID = x,
             Cell_type = names(cor_l[[x]]),
             Cor = cor_l[[x]])
}) %>% 
  do.call(rbind, .) %>% 
  arrange(Cor)


# Average cor by dataset

cor_group_df <- data.frame(
  ID = names(cor_l),
  Mean_cor = unlist(lapply(cor_l, mean))) %>% 
  arrange(Mean_cor)


# Heatmaps of these correlations
# ------------------------------------------------------------------------------


pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(100)
col_breaks <- seq(-1, 1, length.out =  100)


pheatmap(t(cor_df[, "Cor", drop = FALSE]),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = pal,
         breaks = col_breaks,
         border_color = NA,
         cellheight = 30,
         height = 3)



pheatmap(t(cor_group_df[, "Mean_cor", drop = FALSE]),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = pal,
         breaks = col_breaks,
         border_color = "black",
         cellheight = 30,
         height = 3)



# Here, generating scatter plots of gene-gene counts for each cell type in a
# dataset. Exporting the max cor ASCL1-DLL3 cell types in mouse and human, as
# used in the paper
# ------------------------------------------------------------------------------


# Human
id_hg <- "GarciaAlonso2022Human"
gene1_hg <- "ASCL1" 
gene2_hg <- "DLL3"
dat_hg <- load_dat_list(id_hg)[[1]]
mat_hg <- dat_hg$Mat
meta_hg <- dat_hg$Meta

# Mouse
id_mm <- "GSE200202"
gene1_mm <- "Ascl1"
gene2_mm <- "Dll3"
dat_mm <- load_dat_list(id_mm)[[1]]
mat_mm <- dat_mm$Mat
meta_mm <- dat_mm$Meta


# All cell type scatter plots in a list
ct_scatter_hg <- all_celltype_scatter(mat_hg, meta_hg, gene1_hg, gene2_hg)
ct_scatter_mm <- all_celltype_scatter(mat_mm, meta_mm, gene1_mm, gene2_mm)



max_ct_hg <- ct_scatter_hg$`neural cell` +
  ggtitle(label = "GarciaAlonso2022Human: neural cell",
          subtitle = paste("r=0.54")) +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 18))


max_ct_mm <- ct_scatter_mm$`GABAergic INs` +
  ggtitle(label = "GSE200202: GABAergic INs",
          subtitle = paste("r=0.63")) +
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 18))



# Exporting max cor cell types
ggsave(max_ct_hg, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0(gene1_hg, "_", gene2_hg, "_", id_hg, "_maxcor.png")))

ggsave(max_ct_mm, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0(gene1_mm, "_", gene2_mm, "_", id_mm, "_maxcor.png")))
