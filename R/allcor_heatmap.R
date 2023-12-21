library(pheatmap)
library(tidyverse)

cor_l <- readRDS("~/Plots/TR_singlecell/All_corplots/Human_ASCL1_DLL3.RDS")
# cor_l <- readRDS("~/Plots/TR_singlecell/All_corplots/Mouse_Ascl1_Dll3.RDS")



# Why does this not work?
# cor_df <- rbind(do.call, lapply(names(cor_l), function(x) {
#   data.frame(ID = x,
#              Cell_type = cor_l[[x]],
#              Cor = cor_l[[x]])
# }))


cor_df <- lapply(names(cor_l), function(x) {
  data.frame(ID = x,
             Cell_type = names(cor_l[[x]]),
             Cor = cor_l[[x]])
}) %>% 
  do.call(rbind, .) %>% 
  arrange(Cor)


cor_group_df <- data.frame(
  ID = names(cor_l),
  Mean_cor = unlist(lapply(cor_l, mean))) %>% 
  arrange(Mean_cor)



pheatmap(t(cor_df[, "Cor", drop = FALSE]),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("#0571b0", "white", "#ca0020"))(100),
         breaks = seq(-1, 1, length.out =  100),
         border_color = NA,
         cellheight = 30,
         # cellwidth = 30,
         height = 3,
         filename = file.path(plot_dir, "demo_allcor_heatmap.png"))



pheatmap(t(cor_group_df[, "Mean_cor", drop = FALSE]),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = colorRampPalette(c("#0571b0", "white", "#ca0020"))(100),
         breaks = seq(-1, 1, length.out =  100),
         border_color = "black",
         cellheight = 30,
         # cellwidth = 30,
         height = 3,
         filename = file.path(plot_dir, "demo_all_mean_cor_heatmap.png"))



# demo_id <- "GSE200202"
# gene1 <- "Ascl1" 
# gene2 <- "Dll3"
# dat <- load_dat_list(demo_id)[[1]]
# mat <- dat$Mat
# meta <- dat$Meta
# 
# # All cell type scatter plots
# ct_scatter_l <- all_celltype_scatter(mat, meta, gene1, gene2)



demo_id <- "GarciaAlonso2022Human"
gene1 <- "ASCL1" 
gene2 <- "DLL3"
dat <- load_dat_list(demo_id)[[1]]
mat <- dat$Mat
meta <- dat$Meta

# All cell type scatter plots
ct_scatter_l <- all_celltype_scatter(mat, meta, gene1, gene2)

