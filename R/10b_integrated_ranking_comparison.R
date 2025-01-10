## Compare the raw AUC values of the integrated aggregates versus the positive
## coexpression and binding aggregates
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(gplots)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# The minimum count of curated targets for a TF for reporting
min_targets <- 5

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Aggregate coexpr profiles versus null 
coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)

# Integrated ranking versus null 
int_auc_hg <- readRDS(int_auc_hg_path)
int_auc_mm <- readRDS(int_auc_mm_path)

# Aggregate binding profiles versus null
bind_auc_hg <- readRDS(unibind_auc_hg_path)
bind_auc_mm <- readRDS(unibind_auc_mm_path)


# Join the summary dataframes from the list of benchmark results
# ------------------------------------------------------------------------------


join_auc_df <- function(coexpr_l, bind_l, int_l) {
  
  coexpr_df <- do.call(rbind, lapply(coexpr_l, `[[`, "Perf_df"))
  bind_df <- do.call(rbind, lapply(bind_l, `[[`, "Perf_df"))
  int_df <- do.call(rbind, lapply(int_l, `[[`, "Perf_df"))
  
  colnames(coexpr_df)[3:ncol(coexpr_df)] <- paste0(colnames(coexpr_df)[3:ncol(coexpr_df)], "_coexpr")
  colnames(bind_df)[3:ncol(bind_df)] <- paste0(colnames(bind_df)[3:ncol(bind_df)], "_bind")
  colnames(int_df)[3:ncol(int_df)] <- paste0(colnames(int_df)[3:ncol(int_df)], "_int")
  
  join_df <- 
    left_join(coexpr_df, bind_df, by = c("Symbol", "N_targets")) %>% 
    left_join(., int_df, by = c("Symbol", "N_targets"))
  
  return(join_df)
}



agg_df_hg <- join_auc_df(coexpr_l = coexpr_auc_hg, 
                         bind_l = bind_auc_hg, 
                         int_l = int_auc_hg)


agg_df_mm <- join_auc_df(coexpr_l = coexpr_auc_mm, 
                         bind_l = bind_auc_mm,
                         int_auc_mm)


agg_df_hg <- agg_df_hg %>% 
  mutate(N_msr_coexpr = rowSums(msr_hg[agg_df_hg$Symbol, ])) %>% 
  filter(N_targets >= min_targets & !is.na(AUROC_bind))


agg_df_mm <- agg_df_mm %>% 
  mutate(N_msr_coexpr = rowSums(msr_mm[agg_df_mm$Symbol, ])) %>% 
  filter(N_targets >= min_targets & !is.na(AUROC_bind))



# Wilcoxon tests for whether integrated ranking outperforms coexpr/bind

test_integrated <- function(summary_df) {
  
  test_coexpr_df <- data.frame(
    AUROC = c(summary_df$AUROC_int, summary_df$AUROC_coexpr),
    AUPRC = c(summary_df$AUPRC_int, summary_df$AUPRC_coexpr),
    Group = factor(rep(c("Integrated", "Coexpression"), each = nrow(summary_df)),
                   levels = c("Integrated", "Coexpression"))
  )

  test_bind_df <- data.frame(
    AUROC = c(summary_df$AUROC_int, summary_df$AUROC_bind),
    AUPRC = c(summary_df$AUPRC_int, summary_df$AUPRC_bind),
    Group = factor(rep(c("Integrated", "Coexpression"), each = nrow(summary_df)),
                   levels = c("Integrated", "Coexpression"))
  )
  
  c(
    AUROC_coexpr = wilcox.test(AUROC ~ Group, data = test_coexpr_df, alternative = "greater")$p.value,
    AUPRC_coexpr = wilcox.test(AUPRC ~ Group, data = test_coexpr_df, alternative = "greater")$p.value,
    AUPRC_bind = wilcox.test(AUROC ~ Group, data = test_bind_df, alternative = "greater")$p.value,
    AUPRC_bind = wilcox.test(AUPRC ~ Group, data = test_bind_df, alternative = "greater")$p.value
  )
}


test_hg <- test_integrated(agg_df_hg)
test_mm <- test_integrated(agg_df_mm)



# Plots
# ------------------------------------------------------------------------------


# Organize long dfs of AUCs for distribution plots, as well as dfs summarizing
# median AUC to overlay on density

plot_df_hg <- data.frame(
  Value = c(agg_df_hg$AUROC_coexpr, agg_df_hg$AUROC_bind, agg_df_hg$AUROC_int),
  Group = rep(c("Coexpression", "Binding", "Integrated"), each = nrow(agg_df_hg))) %>%
  mutate(Group = factor(Group, levels = c("Coexpression", "Binding", "Integrated")))

med_df_hg <- plot_df_hg %>% 
  group_by(Group) %>% 
  summarize(Med = median(Value))



plot_df_mm <- data.frame(
  Value = c(agg_df_mm$AUROC_coexpr, agg_df_mm$AUROC_bind, agg_df_mm$AUROC_int),
  Group = rep(c("Coexpression", "Binding", "Integrated"), each = nrow(agg_df_mm))) %>%
  mutate(Group = factor(Group, levels = c("Coexpression", "Binding", "Integrated")))

med_df_mm <- plot_df_mm %>% 
  group_by(Group) %>% 
  summarize(Med = median(Value))



# Violin+box plots

vbplot <- function(plot_df, ylab, title) {
  
  ggplot(plot_df, aes(x = Group, y = Value)) +
    geom_violin(fill = "slategrey") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    ylab(ylab) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}


auroc_vbox_hg <- vbplot(plot_df_hg, "AUROC", "Human")
auroc_vbox_mm <- vbplot(plot_df_mm, "AUROC", "Mouse")


ggsave(auroc_vbox_hg, height = 9, width = 9, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "auroc_integrated_vbplot_hg.png")))

ggsave(auroc_vbox_mm, height = 9, width = 9, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "auroc_integrated_vbplot_mm.png")))



# Density plots overlaid with median

dplot <- function(plot_df, med_df, xlab, title) {
  
  ggplot(plot_df, aes(colour = Group, x = Value)) +
    geom_density(linewidth = 1.4) +
    geom_vline(data = med_df, 
               aes(xintercept = Med, colour = Group), 
               linetype = "dashed", size = 1.5) +
    xlab(xlab) +
    ylab("Density") +
    ggtitle(title) +
    scale_colour_manual(values = c('#e41a1c','#377eb8','#4daf4a')) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          legend.position = c(0.2, 0.75),
          legend.text = element_text(size = 20),
          legend.title = element_blank(),
          plot.margin = margin(c(10, 20, 10, 10)))
}


auroc_density_hg <- dplot(plot_df_hg, med_df_hg, "AUROC", "Human")
auroc_density_mm <- dplot(plot_df_mm, med_df_mm, "AUROC", "Mouse")


ggsave(auroc_density_mm, height = 6, width = 9, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "auroc_integrated_density_mm.png")))

ggsave(auroc_density_mm, height = 6, width = 9, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "auroc_integrated_density_mm.png")))



# Heatmaps of AUC values across all three

auc_heatmap <- function(summary_df, title) {
  
  mat <- summary_df[, c("AUROC_coexpr", "AUROC_bind", "AUROC_int")]
  rownames(mat) <- summary_df$Symbol
  colnames(mat) <- c("Coexpression", "Binding", "Integrated")
  mat <- mat[order(rowMeans(mat)), ]
  
  pheatmap(mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = FALSE,
           color = gplots::bluered(11),
           breaks = seq(min(mat), max(mat), length.out = 11),
           border_color = "black",
           gaps_col = c(1:3),
           fontsize_col = 15,
           main = title,
           width = 4)
  
  
}


heatmap_hg <- auc_heatmap(agg_df_hg, "Human")
heatmap_mm <- auc_heatmap(agg_df_mm, "Mouse")
