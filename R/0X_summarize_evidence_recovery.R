## Organize and describe the ability of the coexpression and binding rankings
## to recover literature curated targets
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

n_samps <- 1000

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Saved list RDS of the summaries
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# List of AUROC/AUPRC for TF lists ability to recover curated targets
unibind_auc_hg <- readRDS(unibind_auc_hg_path)
unibind_auc_mm <- readRDS(unibind_auc_mm_path)
coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)
avi_auc_hg <- readRDS(avg_vs_ind_auc_hg_path)
avi_auc_mm <- readRDS(avg_vs_ind_auc_mm_path)
rev_coexpr_auc_hg <- readRDS("/space/scratch/amorin/R_objects/reverse_coexpr_recover_curated_hg.RDS")
rev_coexpr_auc_mm <- readRDS("/space/scratch/amorin/R_objects/reverse_coexpr_recover_curated_mm.RDS")

# The minimum count of curated targets for a TF for reporting
min_targets <- 5


# Functions
# ------------------------------------------------------------------------------


# Give the proportion of data within the supplied numeric columns of df that are
# in the extreme ranges

extreme_proportions <- function(df, stat_cols, upper = 0.9, lower = 0.1) {
  
  res <- lapply(df[, stat_cols], function(col) {
    
    col <- col[!is.na(col)]
    n <- length(col)
    
    stats <- c(
      N = n,
      Eq1 = sum(col == 1) / n,
      Gt_upper = sum(col > upper) / n,
      Lt_lower = sum(col < lower) / n
    )
    
    stats <- round(stats, 3)
  })
  
  df <- do.call(rbind, res)
  row.names(df) <- stat_cols
  colnames(df) <- c("N", "Eq1", paste0("Gt", upper), paste0("Lt", lower))
  
  return(df)
}



# Join the summary dataframes from the list of coexpression and unibind AUCs

join_auc_df <- function(coexpr_l, unibind_l, min_targets) {
  
  coexpr_df <- do.call(rbind, lapply(coexpr_l, `[[`, "Perf_df"))
  unibind_df <- do.call(rbind, lapply(unibind_l, `[[`, "Perf_df"))
  
  left_join(coexpr_df, unibind_df,
            by = c("Symbol", "N_targets"),
            suffix = c("_coexpr", "_unibind")) %>% 
    filter(N_targets >= min_targets)
}



# The AUC percentile of the averaged/final network relative to the distribution
# of the individual TFs
# ------------------------------------------------------------------------------


avi_df_hg <-  lapply(avi_auc_hg, `[[`, "Summary_df") %>% 
  do.call(rbind, .) %>%
  rownames_to_column(var = "Symbol") %>% 
  filter(N_targets >= min_targets)


avi_df_mm <-  lapply(avi_auc_mm, `[[`, "Summary_df") %>% 
  do.call(rbind, .) %>%
  rownames_to_column(var = "Symbol") %>% 
  filter(N_targets >= min_targets)


# Summaries

avi_summ_hg <- summary(Filter(is.numeric, avi_df_hg))
avi_summ_mm <- summary(Filter(is.numeric, avi_df_mm))


avi_cols <- c("AUROC_percentile", "AUPRC_percentile")

ext_avi_hg <- extreme_proportions(avi_df_hg, avi_cols)
ext_avi_mm <- extreme_proportions(avi_df_mm, avi_cols)


# Inspect TFs where the average was outperformed by most individual
# TODO: Why are LYL1 and NRF1 poor in both species? 

# view(filter(avi_df_hg, AUROC_percentile < 0.1 | AUPRC_percentile < 0.1))
# view(filter(avi_df_mm, AUROC_percentile < 0.1 | AUPRC_percentile < 0.1))


# Join and summarize the aggregated coexpression and binding AUCs
# ------------------------------------------------------------------------------


agg_df_hg <- join_auc_df(coexpr_l = coexpr_auc_hg, 
                         unibind_l = unibind_auc_hg, 
                         min_targets = min_targets)


agg_df_mm <- join_auc_df(coexpr_l = coexpr_auc_mm, 
                         unibind_l = unibind_auc_mm,
                         min_targets = min_targets)


# Subset of TFs with data for both binding and coexpr along with min targets

agg_df_common_hg <- filter(agg_df_hg, !is.na(AUPRC_coexpr) & !is.na(AUPRC_unibind))
agg_df_common_mm <- filter(agg_df_mm, !is.na(AUPRC_coexpr) & !is.na(AUPRC_unibind))


# Summaries

agg_summ_hg <- summary(Filter(is.numeric, agg_df_hg))
agg_summ_mm <- summary(Filter(is.numeric, agg_df_mm))

agg_summ_common_hg <- summary(Filter(is.numeric, agg_df_common_hg))
agg_summ_common_mm <- summary(Filter(is.numeric, agg_df_common_mm))


agg_cols <- c(
  "AUPRC_percentile_observed_coexpr",
  "AUPRC_percentile_observed_unibind",
  "AUROC_percentile_observed_coexpr",
  "AUROC_percentile_observed_unibind"
)


ext_agg_hg <- extreme_proportions(agg_df_hg, agg_cols)
ext_agg_mm <- extreme_proportions(agg_df_mm, agg_cols)

ext_agg_common_hg <- extreme_proportions(agg_df_common_hg, agg_cols)
ext_agg_common_mm <- extreme_proportions(agg_df_common_mm, agg_cols)


# Relationship between the count of datasets, targets, and performance
# TODO: this should also look at # of coexpr and binding datasets

cor(select_if(agg_df_common_hg, is.numeric), method = "spearman")
cor(select_if(agg_df_common_mm, is.numeric), method = "spearman")
cor.test(agg_df_common_hg$AUPRC_percentile_observed_coexpr, agg_df_common_hg$N_targets, method = "spearman")



# Isolating TFs that performed well for for both aggregations
# ------------------------------------------------------------------------------


both_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_percentile_observed_coexpr > 0.9 & AUPRC_percentile_observed_unibind > 0.9) |
    (AUROC_percentile_observed_coexpr > 0.9 & AUROC_percentile_observed_unibind > 0.9)
)


both_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_percentile_observed_coexpr > 0.9 & AUPRC_percentile_observed_unibind > 0.9) |
    (AUROC_percentile_observed_coexpr > 0.9 & AUROC_percentile_observed_unibind > 0.9)
)


both_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% both_hg$Symbol & Symbol_mm %in% both_mm$Symbol) %>% 
  pull(Symbol_hg)



both_prop_hg <- length(both_ortho) / nrow(both_hg)
both_prop_mm <- length(both_ortho) / nrow(both_mm)


# TODO: explain when percentile is high for one AUC but not the other
filter(both_hg, sum())
both_hg$AUPRC_percentile_observed_coexpr - both_hg$AUROC_percentile_observed_coexpr



# Looking at TFs that performed well in one aggregation but not the other
# ------------------------------------------------------------------------------


# Coexpression

only_coexpr_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_percentile_observed_coexpr > 0.9 & AUPRC_percentile_observed_unibind < 0.5) |
  (AUROC_percentile_observed_coexpr > 0.9 & AUROC_percentile_observed_unibind < 0.5)
)


only_coexpr_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_percentile_observed_coexpr > 0.9 & AUPRC_percentile_observed_unibind < 0.5) |
    (AUROC_percentile_observed_coexpr > 0.9 & AUROC_percentile_observed_unibind < 0.5)
)


only_coexpr_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% only_coexpr_hg$Symbol & Symbol_mm %in% only_coexpr_mm$Symbol) %>% 
  pull(Symbol_hg)



# Binding

only_binding_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_percentile_observed_unibind > 0.9 & AUPRC_percentile_observed_coexpr < 0.5) |
    (AUROC_percentile_observed_unibind > 0.9 & AUROC_percentile_observed_coexpr < 0.5)
)


only_binding_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_percentile_observed_unibind > 0.9 & AUPRC_percentile_observed_coexpr < 0.5) |
    (AUROC_percentile_observed_unibind > 0.9 & AUROC_percentile_observed_coexpr < 0.5)
)


only_binding_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% only_binding_hg$Symbol & Symbol_mm %in% only_binding_mm$Symbol) %>% 
  pull(Symbol_hg)



# Inspecting behavior of null AUCs
# ------------------------------------------------------------------------------


# Instances where a TFs aggregation was not performant relative to 
# the null, and yet its AUROC is greater than 0.5. This would suggest that its
# aggregate ranking is somewhat "generic" in its ability to recover targets.


check_coexpr_auroc <- agg_df_hg %>% 
  filter(AUROC_percentile_observed_coexpr < 0.5 & AUROC_coexpr > 0.5) %>% 
  slice_max(AUROC_coexpr)


check_binding_auroc <- agg_df_hg %>% 
  filter(AUROC_percentile_observed_unibind < 0.5 & AUROC_unibind > 0.5) %>% 
  slice_max(AUROC_unibind)


# Summarize and order null AUROCs


null_coexpr_auroc <- 
  lapply(coexpr_auc_hg, function(x) summary(unlist(lapply(x$Null, `[[`, "AUROC")))) %>% 
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol %in% agg_df_hg$Symbol) %>% 
  arrange(desc(Median))



null_binding_auroc <- 
  lapply(unibind_auc_hg, function(x) summary(unlist(lapply(x$Null, `[[`, "AUROC")))) %>% 
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol %in% agg_df_hg$Symbol) %>% 
  arrange(desc(Median))


# Note that the TFs with the highest null AUROCs still exceeded the null
# expectation in their respective aggregations.

top_nulls <- filter(agg_df_hg, Symbol %in% c("GATA6", "FOXC1"))



# Are there TFs where reversing ranks (neg coexpr) increased performance?
# ------------------------------------------------------------------------------



rev_coexpr_hg <- lapply(rev_coexpr_auc_hg, `[[`, "Perf_df") %>%
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rename_with(~paste0(., "_reverse_coexpr"), -c("Symbol", "N_targets")) %>% 
  left_join(agg_df_hg, by = c("Symbol", "N_targets")) %>% 
  filter(N_targets >= min_targets)



rev_coexpr_mm <- lapply(rev_coexpr_auc_mm, `[[`, "Perf_df") %>%
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rename_with(~paste0(., "_reverse_coexpr"), -c("Symbol", "N_targets")) %>% 
  left_join(agg_df_mm, by = c("Symbol", "N_targets")) %>% 
  filter(N_targets >= min_targets)



only_rev_hg <- filter(
  rev_coexpr_hg, 
  (AUPRC_percentile_observed_reverse_coexpr > 0.9 & AUPRC_percentile_observed_coexpr < 0.5) |
    (AUROC_percentile_observed_reverse_coexpr > 0.9 & AUROC_percentile_observed_coexpr < 0.5)
)


only_rev_mm <- filter(
  rev_coexpr_mm, 
  (AUPRC_percentile_observed_reverse_coexpr > 0.9 & AUPRC_percentile_observed_coexpr < 0.5) |
    (AUROC_percentile_observed_reverse_coexpr > 0.9 & AUROC_percentile_observed_coexpr < 0.5)
)


only_rev_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% only_rev_hg$Symbol & Symbol_mm %in% only_rev_mm$Symbol) %>% 
  pull(Symbol_hg)



# Elevated performance in reverse coexpression and binding: repression?

rev_binding_hg <- filter(only_rev_hg, Symbol %in% only_binding_hg$Symbol)
rev_binding_mm <- filter(only_rev_mm, Symbol %in% only_binding_mm$Symbol)


# Performant in coexpression and reverse coexpression... suggests targets are
# at extremes of rankings


both_coexpr_hg <- filter(
  rev_coexpr_hg, 
  (AUPRC_percentile_observed_reverse_coexpr > 0.9 & AUPRC_percentile_observed_coexpr > 0.9) |
    (AUROC_percentile_observed_reverse_coexpr > 0.9 & AUROC_percentile_observed_coexpr > 0.9)
)


both_coexpr_mm <- filter(
  rev_coexpr_mm, 
  (AUPRC_percentile_observed_reverse_coexpr > 0.9 & AUPRC_percentile_observed_coexpr > 0.9) |
    (AUROC_percentile_observed_reverse_coexpr > 0.9 & AUROC_percentile_observed_coexpr > 0.9)
)



check_tf <- "PAX7"
targets <- filter(curated, TF_Symbol == check_tf)$Target_Symbol
targets_rank <- rank_tf_hg[[check_tf]] %>% filter(Symbol %in% targets) %>% arrange(Rank_RSR)




# Plotting
# ------------------------------------------------------------------------------


# Histograms of Average versus individual AUCs

plot_grid(
  plot_hist(avi_df_hg, stat_col = paste0("AUPRC", "_percentile")),
  plot_hist(avi_df_hg, stat_col = paste0("AUROC", "_percentile")),
  plot_hist(avi_df_mm, stat_col = paste0("AUPRC", "_percentile")),
  plot_hist(avi_df_mm, stat_col = paste0("AUROC", "_percentile")),
  nrow = 2)


# Density plot of a TF's individual AUCs overlaid with the average AUC

tf <- "ASCL1"
auc_df <- arrange(avi_auc_hg[[tf]]$AUC_df, AUROC)
auc_df_no_avg <- filter(auc_df, ID != "Average")
auc_avg <- filter(auc_df, ID == "Average")


ggplot(auc_df_no_avg, aes(x = AUROC)) +
  geom_density(fill = "lightgrey") +
  geom_vline(xintercept = auc_avg$AUROC, col = "firebrick") +
  ylab("Density") +
  ggtitle(tf) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))



# Histogram of percentiles of observered AUCs versus null


stat <- "AUPRC"

# Using all TFs available for coexpression
# plot_df_hg <- agg_df_hg 
# plot_df_mm <- agg_df_mm

# Using only TFs common to coexpression and unibind
plot_df_hg <- agg_df_common_hg
plot_df_mm <- agg_df_common_mm


p_auprc_hist_hg <- list(
  plot_hist(plot_df_hg, stat_col = paste0(stat, "_percentile_observed_coexpr")),
  plot_hist(plot_df_hg, stat_col = paste0(stat, "_percentile_observed_unibind")))


p_auprc_hist_mm <- list(
  plot_hist(plot_df_mm, stat_col = paste0(stat, "_percentile_observed_coexpr")),
  plot_hist(plot_df_mm, stat_col = paste0(stat, "_percentile_observed_unibind")))


plot_grid(
  plot_grid(plotlist = p_auprc_hist_hg), 
  plot_grid(plotlist = p_auprc_hist_mm),
  nrow = 2)


plot_grid(
  plot_hist(plot_df_hg, stat_col = paste0("AUPRC", "_percentile_observed_coexpr")),
  plot_hist(plot_df_hg, stat_col = paste0("AUROC", "_percentile_observed_coexpr")),
  plot_hist(plot_df_mm, stat_col = paste0("AUPRC", "_percentile_observed_coexpr")),
  plot_hist(plot_df_mm, stat_col = paste0("AUROC", "_percentile_observed_coexpr")),
  nrow = 4)




# Demo a single TF's auc performance relative to null distn
# TODO: function and clean

tf_hg <- "ASCL1"
tf_mm <- "Ascl1"

null <- unlist(lapply(unibind_auc_hg[[tf_hg]]$Null, `[[`, stat))
obs <- unibind_auc_hg[[tf_hg]]$Perf_df[[stat]]

xmin <- min(obs, min(null)) * 0.95
xmax <- max(obs, max(null)) * 1.05
hist(null, xlim = c(xmin, xmax), main = tf_hg, breaks = 100)
abline(v = obs, col = "red")

null <- unlist(lapply(unibind_auc_mm[[tf_mm]]$Null, `[[`, stat))
obs <- unibind_auc_mm[[tf_mm]]$Perf_df[[stat]]

xmin <- min(obs, min(null)) * 0.95
xmax <- max(obs, max(null)) * 1.05
hist(null, xlim = c(xmin, xmax), main = tf_mm, breaks = 100)
abline(v = obs, col = "red")



# Scatterplot of coexpression versus unibind performance for each TF


qplot(agg_df_common_hg, xvar = "AUPRC_coexpr", yvar = "AUPRC_unibind")
qplot(agg_df_common_mm, xvar = "AUPRC_coexpr", yvar = "AUPRC_unibind")
qplot(agg_df_common_hg, xvar = "AUROC_coexpr", yvar = "AUROC_unibind")
qplot(agg_df_common_mm, xvar = "AUROC_coexpr", yvar = "AUROC_unibind")


ggplot(agg_df_common_hg, 
       aes(x = AUPRC_percentile_observed_coexpr, 
           y = AUPRC_percentile_observed_unibind)) +
  
  # geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = 0.9, ymax = 1.05), fill = NA, colour = "forestgreen") +
  # geom_rect(aes(xmin = -0.05, xmax = 0.1, ymin = 0.9, ymax = 1.05), fill = NA, colour = "grey") +
  # geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = -0.05, ymax = 0.1), fill = NA, colour = "grey") +
  # geom_rect(aes(xmin = -0.05, xmax = 0.5, ymin = -0.05, ymax = 0.5), fill = NA, colour = "firebrick") +
  
  geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = 0.9, ymax = 1.05), fill = "forestgreen", colour = NA, alpha = 0.006) +
  geom_rect(aes(xmin = -0.05, xmax = 0.1, ymin = 0.9, ymax = 1.05), fill = "grey", colour = NA, alpha = 0.006) +
  geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = -0.05, ymax = 0.1), fill = "grey", colour = NA, alpha = 0.006) +
  geom_rect(aes(xmin = -0.05, xmax = 0.5, ymin = -0.05, ymax = 0.5), fill = "firebrick", colour = NA, alpha = 0.006) +
  
  geom_point(shape = 21, size = 3.4) +
  xlab("AUPRC percentile coexpression") +
  ylab("AUPRC percentile binding") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))


# Scatterplot of count of targets and performance

qplot(agg_df_common_hg, yvar = "AUPRC_coexpr", xvar = "N_targets")
qplot(agg_df_common_hg, yvar = "AUPRC_percentile_observed_coexpr", xvar = "N_targets")

qplot(agg_df_common_hg, yvar = "AUROC_coexpr", xvar = "N_targets")
qplot(agg_df_common_hg, yvar = "AUROC_percentile_observed_coexpr", xvar = "N_targets")



# Looking at the distribution of observed percentiles for coexpr/binding
# TODO: function


plot_perc_df <- agg_df_common_hg
# plot_perc_df <- agg_df_hg

plot_perc_df <- data.frame(
  Percentile = c(
    plot_perc_df$AUPRC_percentile_observed_coexpr,
    plot_perc_df$AUPRC_percentile_observed_unibind,
    plot_perc_df$AUROC_percentile_observed_coexpr,
    plot_perc_df$AUROC_percentile_observed_unibind
  ),
  Group = c(
    rep("Coexpr_AUPRC", length(plot_perc_df$AUPRC_percentile_observed_coexpr)),
    rep("Binding_AUPRC", length(plot_perc_df$AUPRC_percentile_observed_unibind)),
    rep("Coexpr_AUROC", length(plot_perc_df$AUROC_percentile_observed_coexpr)),
    rep("Binding_AUROC", length(plot_perc_df$AUROC_percentile_observed_unibind))
  )
)

plot_perc_df$Group <- factor(plot_perc_df$Group, levels = unique(plot_perc_df$Group))




pb <- ggplot(plot_perc_df, aes(y = Group, x = Percentile)) +
  geom_violin(fill = "slategrey") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  xlab("Percentile observed") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 20))


# ggsave(pb, height = 4, width = 8, device = "png", dpi = 600,
#        filename = file.path(paste0(plot_dir, "binding_vs_agg_AUC_all_human.png")))

# ggsave(pb, height = 4, width = 8, device = "png", dpi = 600,
#        filename = file.path(paste0(plot_dir, "binding_vs_agg_AUC_common_human.png")))



# Plotting null versus observed AUROCs


plot(density(agg_df_hg$AUROC_coexpr), ylim = c(0, 20))
lines(density(null_coexpr_auroc$Median), col = "red")


left_join(agg_df_hg, null_coexpr_auroc, by = "Symbol") %>% 
  qplot(., xvar = "Median", yvar = "AUROC_coexpr") +
  geom_hline(yintercept = 0.5, col = "red") +
  geom_vline(xintercept = 0.5, col = "red") +
  ylab("Observed coexpression AUROC") +
  xlab("Median null coexpression AUROC")



left_join(agg_df_hg, null_binding_auroc, by = "Symbol") %>% 
  qplot(., xvar = "Median", yvar = "AUROC_unibind") +
  geom_hline(yintercept = 0.5, col = "red") +
  geom_vline(xintercept = 0.5, col = "red") +
  ylab("Observed binding AUROC") +
  xlab("Median null binding AUROC")



# Scatter of coexpr versus reverse coexpr
plot(rev_coexpr_hg$AUPRC_percentile_observed_coexpr, rev_coexpr_hg$AUPRC_percentile_observed_reverse_coexpr)
plot(rev_coexpr_mm$AUPRC_percentile_observed_coexpr, rev_coexpr_mm$AUPRC_percentile_observed_reverse_coexpr)

plot(rev_coexpr_hg$AUPRC_coexpr, rev_coexpr_hg$AUPRC_reverse_coexpr)
plot(rev_coexpr_mm$AUPRC_coexpr, rev_coexpr_mm$AUPRC_reverse_coexpr)

plot(rev_coexpr_hg$AUROC_percentile_observed_coexpr, rev_coexpr_hg$AUROC_percentile_observed_reverse_coexpr)
plot(rev_coexpr_mm$AUROC_percentile_observed_coexpr, rev_coexpr_mm$AUROC_percentile_observed_reverse_coexpr)

plot(rev_coexpr_hg$AUROC_coexpr, rev_coexpr_hg$AUROC_reverse_coexpr)
plot(rev_coexpr_mm$AUROC_coexpr, rev_coexpr_mm$AUROC_reverse_coexpr)


# Scatter of reverse coexpr versus binding
plot(rev_coexpr_hg$AUPRC_percentile_observed_reverse_coexpr, rev_coexpr_hg$AUPRC_percentile_observed_unibind)
plot(rev_coexpr_mm$AUPRC_percentile_observed_reverse_coexpr, rev_coexpr_mm$AUPRC_percentile_observed_unibind)

plot(rev_coexpr_hg$AUPRC_reverse_coexpr, rev_coexpr_hg$AUPRC_unibind)
plot(rev_coexpr_mm$AUPRC_reverse_coexpr, rev_coexpr_mm$AUPRC_unibind)

plot(rev_coexpr_hg$AUROC_percentile_observed_reverse_coexpr, rev_coexpr_hg$AUROC_percentile_observed_unibind)
plot(rev_coexpr_mm$AUROC_percentile_observed_reverse_coexpr, rev_coexpr_mm$AUROC_percentile_observed_unibind)

plot(rev_coexpr_hg$AUROC_reverse_coexpr, rev_coexpr_hg$AUROC_unibind)
plot(rev_coexpr_mm$AUROC_reverse_coexpr, rev_coexpr_mm$AUROC_unibind)


# Distribution of performances for coexpr, reverse coexpr, and binding

plot_df <- data.frame(
  AUROC = c(
    rev_coexpr_hg$AUROC_reverse_coexpr,
    rev_coexpr_hg$AUROC_coexpr,
    rev_coexpr_hg$AUROC_unibind
  ),
  Group = c(
    rep("Coexpr_reverse", length(rev_coexpr_hg$AUROC_reverse_coexpr)),
    rep("Coexpr", length(rev_coexpr_hg$AUROC_coexpr)),
    rep("Binding", length(rev_coexpr_hg$AUROC_unibind))
  )
)

plot_df <- filter(plot_df, !is.na(AUROC))

boxplot(plot_df$AUROC ~ plot_df$Group)


# plot_df <- data.frame(
#   AUPRC = c(
#     rev_coexpr_hg$AUPRC_reverse_coexpr,
#     rev_coexpr_hg$AUPRC_coexpr,
#     rev_coexpr_hg$AUPRC_unibind
#   ),
#   Group = c(
#     rep("Coexpr_reverse", length(rev_coexpr_hg$AUPRC_reverse_coexpr)),
#     rep("Coexpr", length(rev_coexpr_hg$AUPRC_coexpr)),
#     rep("Binding", length(rev_coexpr_hg$AUPRC_unibind))
#   )
# )
# 
# plot_df <- filter(plot_df, !is.na(AUPRC))
# 
# boxplot(plot_df$AUPRC ~ plot_df$Group)



plot_df <- data.frame(
  AUPRC = c(
    rev_coexpr_hg$AUPRC_percentile_observed_reverse_coexpr,
    rev_coexpr_hg$AUPRC_percentile_observed_coexpr,
    rev_coexpr_hg$AUPRC_percentile_observed_unibind
  ),
  Group = c(
    rep("Coexpr_reverse", length(rev_coexpr_hg$AUPRC_percentile_observed_reverse_coexpr)),
    rep("Coexpr", length(rev_coexpr_hg$AUPRC_percentile_observed_coexpr)),
    rep("Binding", length(rev_coexpr_hg$AUPRC_percentile_observed_unibind))
  )
)

plot_df <- filter(plot_df, !is.na(AUPRC))

boxplot(plot_df$AUPRC ~ plot_df$Group)



## TODO: Candidate to move to test?


# tf <- "Pax6"
# measure <- "both"
# 
# 
# rank_df <- rank_tf_mm[[tf]] %>% 
#   filter(Symbol != tf) %>% 
#   arrange(desc(Avg_RSR))
# 
# 
# score_vec <- setNames(rank_df$Avg_RSR, rank_df$Symbol)
# 
# 
# curated_tf <- curated %>%
#   filter(str_to_upper(TF_Symbol) == str_to_upper(tf) &
#            !(str_to_upper(Target_Symbol) == str_to_upper(tf)))
# 
# 
# labels <- unique(str_to_title(curated_tf$Target_Symbol))
# labels <- labels[labels %in% pc_mm$Symbol]
# label_vec <- names(score_vec) %in% labels
# 
# 
# auc <- get_auc(score_vec = score_vec, label_vec = label_vec, measure = measure)
# 
# 
# null <- get_null_aucormance(score_vec = score_vec,
#                              label_all = targets_curated_mm,
#                              measure = measure,
#                              n_target = length(labels),
#                              n_samps = 1000,
#                              ncores = 8)


# unlist(lapply(null, `[[`, "AUROC"))
# unlist(lapply(null, `[[`, "AUPRC"))
# hist(unlist(lapply(null, `[[`, "AUPRC")), xlim = c(0, 0.05))
# abline(v = auc$AUPRC, col = "red")



# tt <- summarize_obs_and_null_auc(tf = tf,
#                                  score_vec = score_vec,
#                                  label_vec = label_vec,
#                                  label_all = targets_curated_mm,
#                                  n_samps = 1000,
#                                  ncores = 8)
# 
# 
# 
# 
# tt <- curated_obs_and_null_auc(tf = tf,
#                                rank_df = rank_tf_mm[[tf]],
#                                score_col = "Avg_RSR",
#                                curated_df = curated,
#                                label_all = targets_curated_mm,
#                                pc_df = pc_mm,
#                                species = "Mouse",
#                                n_samps = 1000,
#                                ncores = 8)
# 
# 
# tt <- curated_obs_and_null_auc_list(tfs = tfs_curated_mm[6:8],
#                                     rank_l = rank_tf_mm,
#                                     score_col = "Avg_RSR",
#                                     curated_df = curated,
#                                     label_all = targets_curated_mm,
#                                     pc_df = pc_mm,
#                                     species = "Mouse",
#                                     n_samps = 1000,
#                                     ncores = 8,
#                                     verbose = TRUE)
# 
# 
# 
# save_curated_auc_list(path = "tt.RDS",
#                       tfs = tfs_curated_mm[6:8],
#                       rank_l = rank_tf_mm,
#                       score_col = "Avg_RSR",
#                       curated_df = curated,
#                       label_all = targets_curated_mm,
#                       pc_df = pc_mm,
#                       species = "Mouse",
#                       n_samps = 1000,
#                       ncores = 8,
#                       verbose = TRUE,
#                       force_resave = TRUE)
