library(tidyverse)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)

gse93374_rsr <- readRDS("~/scratch/R_objects/GSE93374_ranksumrank_celltype_cor.RDS")
gse93374_binthresh <- readRDS("~/scratch/R_objects/GSE93374_binthresh_celltype_cor.RDS")


tfs <- names(rank_l$Mouse)


# Which ranking columns to use


keep_cols <- c(
  "Rank_integrated",
  "Rank_sum_rank",
  "Rank_binthresh")


plot_clrs <- c(
  "Rank_integrated" = "black",
  "Rank_sum_rank" = "royalblue",
  "Rank_binthresh" = "gold1"
)


label_col <- "Group"


perf_l <- lapply(tfs, function(tf) {
  
  
  cor_df <- data.frame(
    Symbol = rownames(gse93374_binthresh),
    Rank_sum_rank = gse93374_rsr[, tf],
    Binthresh = gse93374_binthresh[, tf],
    Rank_binthresh = rank(-gse93374_binthresh[, tf])
  )
  
  
  evidence <- left_join(rank_l$Mouse[[tf]], cor_df, by = "Symbol") %>% 
    filter(Symbol != tf) %>% 
    mutate(Group = Rank_integrated <= 500)
  
  
  stopifnot(all(keep_cols %in% colnames(evidence)))
  
  
  # Precision recall dfs and AURPC values
  pr_df <- all_perf_df(evidence, keep_cols, label_col = label_col, measure = "PR")
  auprc <- all_au_perf(evidence, keep_cols, label_col = label_col, measure = "AUPRC")
  
  # ROC dfs and AUROC values
  roc_df <- all_perf_df(evidence, keep_cols, label_col = label_col, measure = "ROC")
  auroc <- all_au_perf(evidence, keep_cols, label_col = label_col, measure = "AUROC")
  
  
  stopifnot(all(names(cols) %in% keep_cols))

  
  plot_roc <- plot_perf(df = roc_df, auc_l = auroc, measure = "ROC", cols = plot_clrs, title = tf, ncol_legend = 1)
  plot_pr <- plot_perf(df = pr_df, auc_l = auprc, measure = "PR", cols = plot_clrs, title = tf, ncol_legend = 1)
  
  return(list(
    Evidence = evidence,
    Plot_ROC = plot_roc,
    Plot_PR = plot_pr,
    AUPRC = as.numeric(auprc),
    AUROC = as.numeric(auroc)))
  
})
names(perf_l) <- tfs


# AUPRC heatmap

all_auprc_mat <- do.call(cbind, lapply(perf_l, `[[`, "AUPRC"))
rownames(all_auprc_mat) <- keep_cols
all_auprc_mat <- all_auprc_mat[order(rowMeans(all_auprc_mat)), ]


pheatmap(all_auprc_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
         # display_numbers = round(all_auprc_mat, 3),
         number_color = "black",
         border_color = NA,
         fontsize = 20,
         cellwidth = 50,
         cellheight = 50,
         height = 9,
         width = 9)





# AUROC heatmap

all_auroc_mat <- do.call(cbind, lapply(perf_l, `[[`, "AUROC"))
rownames(all_auroc_mat) <- keep_cols
all_auroc_mat <- all_auroc_mat[order(rowMeans(all_auroc_mat)), ]


pheatmap(all_auroc_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(brewer.pal(n = 9, name = "YlGn"))(100),
         # display_numbers = round(all_auproc_mat, 3),
         number_color = "black",
         border_color = NA,
         fontsize = 20,
         cellwidth = 50,
         cellheight = 50,
         height = 9,
         width = 9)
