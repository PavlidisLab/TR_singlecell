## TODO:
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# The minimum count of curated targets for a TF for reporting
min_targets <- 5

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Aggregate coexpr profiles versus null 
coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)

# Integrated ranking versus null 
int_auc_hg_path <- "/space/scratch/amorin/R_objects/integrated_recover_curated_hg.RDS"
int_auc_mm_path <- "/space/scratch/amorin/R_objects/integrated_recover_curated_mm.RDS"
int_auc_hg <- readRDS(int_auc_hg_path)
int_auc_mm <- readRDS(int_auc_mm_path)

# Aggregate binding profiles versus null
bind_auc_hg <- readRDS(unibind_auc_hg_path)
bind_auc_mm <- readRDS(unibind_auc_mm_path)


# Functions
# ------------------------------------------------------------------------------


# Join the summary dataframes from the list of coexpression and bind AUCs

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



#
# ------------------------------------------------------------------------------

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





auroc_box_hg <- data.frame(
  Value = c(agg_df_hg$AUROC_coexpr, agg_df_hg$AUROC_bind, agg_df_hg$AUROC_int),
  Group = rep(c("Coexpr", "Bind", "Int"), each = nrow(agg_df_hg))) %>%
  mutate(Group = factor(Group, levels = c("Coexpr", "Bind", "Int"))) %>% 
  ggplot(., aes(x = Group, y = Value)) +
  geom_boxplot() +
  ylab("AUROC") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))




auroc_box_mm <- data.frame(
  Value = c(agg_df_mm$AUROC_coexpr, agg_df_mm$AUROC_bind, agg_df_mm$AUROC_int),
  Group = rep(c("Coexpr", "Bind", "Int"), each = nrow(agg_df_mm))) %>%
  mutate(Group = factor(Group, levels = c("Coexpr", "Bind", "Int"))) %>% 
  ggplot(., aes(x = Group, y = Value)) +
  geom_boxplot() +
  ylab("AUROC") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))




auroc_mat_hg <- agg_df_hg[, c("AUROC_coexpr", "AUROC_bind", "AUROC_int")]
rownames(auroc_mat_hg) <- agg_df_hg$Symbol
colnames(auroc_mat_hg) <- c("Coexpression", "Binding", "Integrated")
auroc_mat_hg <- auroc_mat_hg[order(rowMeans(auroc_mat_hg)), ]


auroc_mat_mm <- agg_df_mm[, c("AUROC_coexpr", "AUROC_bind", "AUROC_int")]
rownames(auroc_mat_mm) <- agg_df_mm$Symbol
colnames(auroc_mat_mm) <- c("Coexpression", "Binding", "Integrated")
auroc_mat_mm <- auroc_mat_mm[order(rowMeans(auroc_mat_mm)), ]




# color = c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')
# breaks = seq(0, 1, length.out = 9)

# color = viridis::inferno(10)
# color = hcl.colors(10)
# color = rev(brewer.pal(10, "RdYlBu"))
# breaks = seq(min(tt1), max(tt1), length.out = 10)
# breaks = seq(0, 1, length.out = 11)


pheatmap(auroc_mat_hg,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         color = gplots::bluered(11),
         breaks = seq(min(auroc_mat_hg), max(auroc_mat_hg), length.out = 11),
         border_color = "black",
         gaps_col = c(1:3),
         fontsize_col = 15)


# pheatmap(t(auroc_mat_hg),
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          # show_rownames = FALSE,
#          color = gplots::bluered(11),
#          breaks = seq(min(auroc_mat_hg), max(auroc_mat_hg), length.out = 11),
#          border_color = "black",
#          gaps_row = c(1:3),
#          fontsize_row = 15,
#          fontsize_col = 8)



pheatmap(auroc_mat_mm,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         color = gplots::bluered(11),
         breaks = seq(min(auroc_mat_mm), max(auroc_mat_mm), length.out = 11),
         border_color = "black",
         gaps_col = c(1:3),
         fontsize_col = 15)




plot(density(agg_df_hg$AUROC_int), col = "blue")
lines(density(agg_df_hg$AUROC_coexpr), col = "red")
lines(density(agg_df_hg$AUROC_bind), col = "black")




plot(density(agg_df_hg$AUROC_int), col = "blue")
lines(density(agg_df_hg$AUROC_coexpr), col = "red")
lines(density(agg_df_hg$AUROC_bind), col = "black")

plot(density(agg_df_hg$AUROC_quantile_int), col = "blue")
lines(density(agg_df_hg$AUROC_quantile_coexpr), col = "red")
lines(density(agg_df_hg$AUROC_quantile_bind), col = "black")




plot(density(agg_df_hg$AUPRC_int), col = "blue")
lines(density(agg_df_hg$AUPRC_coexpr), col = "red")
lines(density(agg_df_hg$AUPRC_bind), col = "black")
