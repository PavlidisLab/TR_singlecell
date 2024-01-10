## Organize and describe the ability of the coexpression and binding rankings
## to recover literature curated targets
# TODO: smarter organization of aggregate tables
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

# The minimum count of curated targets for a TF for reporting
min_targets <- 5

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
# TODO: finalize pathing
collection <- "Permissive"
unibind_auc_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_recover_curated_hg.RDS"))
unibind_auc_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_recover_curated_mm.RDS"))

coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)

avi_auc_hg <- readRDS(avg_vs_ind_auc_hg_path)
avi_auc_mm <- readRDS(avg_vs_ind_auc_mm_path)

rev_coexpr_auc_hg <- readRDS("/space/scratch/amorin/R_objects/reverse_coexpr_recover_curated_hg.RDS")
rev_coexpr_auc_mm <- readRDS("/space/scratch/amorin/R_objects/reverse_coexpr_recover_curated_mm.RDS")

bind_dat <- readRDS(bind_dat_path)
bind_summary <- readRDS(bind_summary_path)


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

avi_underperf_hg <- filter(avi_df_hg, AUROC_percentile < 0.1 | AUPRC_percentile < 0.1)
avi_underperf_mm <- filter(avi_df_mm, AUROC_percentile < 0.1 | AUPRC_percentile < 0.1)

avi_underperf_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% avi_underperf_hg$Symbol & Symbol_mm %in% avi_underperf_mm$Symbol) %>% 
  pull(Symbol_hg)



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
  "AUPRC_percentile_coexpr",
  "AUPRC_percentile_unibind",
  "AUROC_percentile_coexpr",
  "AUROC_percentile_unibind"
)


ext_agg_hg <- extreme_proportions(agg_df_hg, agg_cols)
ext_agg_mm <- extreme_proportions(agg_df_mm, agg_cols)

ext_agg_common_hg <- extreme_proportions(agg_df_common_hg, agg_cols)
ext_agg_common_mm <- extreme_proportions(agg_df_common_mm, agg_cols)


# Difference of coexpression aggregate AUCs to binding AUCs


# # More positive values mean that coexpression AUC was more performant
agg_auc_delta <- list(
  AUPRC_human = agg_df_common_hg$AUPRC_coexpr - agg_df_common_hg$AUPRC_unibind,
  AUPRC_mouse = agg_df_common_mm$AUPRC_coexpr - agg_df_common_mm$AUPRC_unibind,
  AUROC_human = agg_df_common_hg$AUROC_coexpr - agg_df_common_hg$AUROC_unibind,
  AUROC_mouse = agg_df_common_mm$AUROC_coexpr - agg_df_common_mm$AUROC_unibind
)



# More positive values mean that coexpression better separated from its null than binding
# agg_auc_delta <- list(
#   AUPRC_human = agg_df_common_hg$AUPRC_diff_coexpr - agg_df_common_hg$AUPRC_diff_unibind,
#   AUPRC_mouse = agg_df_common_mm$AUPRC_diff_coexpr - agg_df_common_mm$AUPRC_diff_unibind,
#   AUROC_human = agg_df_common_hg$AUROC_diff_coexpr - agg_df_common_hg$AUROC_diff_unibind,
#   AUROC_mouse = agg_df_common_mm$AUROC_diff_coexpr - agg_df_common_mm$AUROC_diff_unibind
# )


delta_auc_summ <- lapply(agg_auc_delta, summary)



# Relationship between the count of datasets, targets, and performance
# TODO: this should also look at # of coexpr and binding datasets
# ------------------------------------------------------------------------------


agg_df_hg$N_coexpr_msr <- rowSums(msr_hg[agg_df_hg$Symbol, ])
agg_df_mm$N_coexpr_msr <- rowSums(msr_mm[agg_df_mm$Symbol, ])

# Spearman cor test AUPRC (raw and percentile) versus count of data
cor.test(agg_df_hg$AUPRC_coexpr, agg_df_hg$N_targets, method = "spearman")
cor.test(agg_df_hg$AUPRC_coexpr, agg_df_hg$N_coexpr_msr, method = "spearman")
cor.test(agg_df_hg$AUPRC_percentile_coexpr, agg_df_hg$N_targets, method = "spearman")
cor.test(agg_df_hg$AUPRC_percentile_coexpr, agg_df_hg$N_coexpr_msr, method = "spearman")

# Spearman cor test AUROC (raw and percentile) versus count of data
cor.test(agg_df_hg$AUROC_coexpr, agg_df_hg$N_targets, method = "spearman")
cor.test(agg_df_hg$AUROC_coexpr, agg_df_hg$N_coexpr_msr, method = "spearman")
cor.test(agg_df_hg$AUROC_percentile_coexpr, agg_df_hg$N_targets, method = "spearman")
cor.test(agg_df_hg$AUROC_percentile_coexpr, agg_df_hg$N_coexpr_msr, method = "spearman")

# Linear model of performance ~ n_targets + n_datasets_measuring_TF
n_model_hg1 <- lm(AUPRC_coexpr ~ N_targets + N_coexpr_msr, data = agg_df_hg)
n_model_hg2 <- lm(AUPRC_percentile_coexpr ~ N_targets + N_coexpr_msr, data = agg_df_hg)
n_model_hg3 <- lm(AUROC_coexpr ~ N_targets + N_coexpr_msr, data = agg_df_hg)
n_model_hg4 <- lm(AUROC_percentile_coexpr ~ N_targets + N_coexpr_msr, data = agg_df_hg)

# summary(n_model_hg1)
# summary(n_model_hg2)
# summary(n_model_hg3)
# summary(n_model_hg4)


# Isolating TFs that performed well for for both aggregations
# ------------------------------------------------------------------------------


both_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_percentile_coexpr > 0.9 & AUPRC_percentile_unibind > 0.9) |
    (AUROC_percentile_coexpr > 0.9 & AUROC_percentile_unibind > 0.9)
)


both_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_percentile_coexpr > 0.9 & AUPRC_percentile_unibind > 0.9) |
    (AUROC_percentile_coexpr > 0.9 & AUROC_percentile_unibind > 0.9)
)


both_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% both_hg$Symbol & Symbol_mm %in% both_mm$Symbol) %>% 
  pull(Symbol_hg)



both_prop_hg <- length(both_ortho) / nrow(both_hg)
both_prop_mm <- length(both_ortho) / nrow(both_mm)


# TODO: explain when percentile is high for one AUC but not the other
# filter(both_hg, sum())
# both_hg$AUPRC_percentile_coexpr - both_hg$AUROC_percentile_coexpr



# Looking at TFs that performed well in one aggregation but not the other
# ------------------------------------------------------------------------------


# Coexpression

only_coexpr_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_percentile_coexpr > 0.9 & AUPRC_percentile_unibind < 0.5) |
  (AUROC_percentile_coexpr > 0.9 & AUROC_percentile_unibind < 0.5)
)


only_coexpr_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_percentile_coexpr > 0.9 & AUPRC_percentile_unibind < 0.5) |
    (AUROC_percentile_coexpr > 0.9 & AUROC_percentile_unibind < 0.5)
)


only_coexpr_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% only_coexpr_hg$Symbol & Symbol_mm %in% only_coexpr_mm$Symbol) %>% 
  pull(Symbol_hg)



# Binding

only_binding_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_percentile_unibind > 0.9 & AUPRC_percentile_coexpr < 0.5) |
    (AUROC_percentile_unibind > 0.9 & AUROC_percentile_coexpr < 0.5)
)


only_binding_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_percentile_unibind > 0.9 & AUPRC_percentile_coexpr < 0.5) |
    (AUROC_percentile_unibind > 0.9 & AUROC_percentile_coexpr < 0.5)
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
  filter(AUROC_percentile_coexpr < 0.5 & AUROC_coexpr > 0.5) %>% 
  slice_max(AUROC_coexpr)


check_binding_auroc <- agg_df_hg %>% 
  filter(AUROC_percentile_unibind < 0.5 & AUROC_unibind > 0.5) %>% 
  slice_max(AUROC_unibind)


# Summarize and order null AUROCs
# TODO: function


null_coexpr_auroc_hg <- 
  lapply(coexpr_auc_hg, function(x) summary(unlist(lapply(x$Null, `[[`, "AUROC")))) %>% 
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol %in% agg_df_hg$Symbol) %>% 
  arrange(desc(Median))


null_coexpr_auroc_mm <- 
  lapply(coexpr_auc_mm, function(x) summary(unlist(lapply(x$Null, `[[`, "AUROC")))) %>% 
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol %in% agg_df_mm$Symbol) %>% 
  arrange(desc(Median))


null_binding_auroc_hg <- 
  lapply(unibind_auc_hg, function(x) summary(unlist(lapply(x$Null, `[[`, "AUROC")))) %>% 
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol %in% agg_df_hg$Symbol) %>% 
  arrange(desc(Median))


null_binding_auroc_mm <- 
  lapply(unibind_auc_mm, function(x) summary(unlist(lapply(x$Null, `[[`, "AUROC")))) %>% 
  do.call(rbind, .) %>%
  as.data.frame() %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol %in% agg_df_mm$Symbol) %>% 
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
  (AUPRC_percentile_reverse_coexpr > 0.9 & AUPRC_percentile_coexpr < 0.5) |
    (AUROC_percentile_reverse_coexpr > 0.9 & AUROC_percentile_coexpr < 0.5)
)


only_rev_mm <- filter(
  rev_coexpr_mm, 
  (AUPRC_percentile_reverse_coexpr > 0.9 & AUPRC_percentile_coexpr < 0.5) |
    (AUROC_percentile_reverse_coexpr > 0.9 & AUROC_percentile_coexpr < 0.5)
)


only_rev_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% only_rev_hg$Symbol & Symbol_mm %in% only_rev_mm$Symbol) %>% 
  pull(Symbol_hg)



# Elevated performance in reverse coexpression and binding: repression?

rev_binding_hg <- filter(only_rev_hg, Symbol %in% only_binding_hg$Symbol)

rev_binding_mm <- filter(only_rev_mm, Symbol %in% only_binding_mm$Symbol)

rev_binding_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% rev_binding_hg$Symbol & Symbol_mm %in% rev_binding_mm$Symbol) %>% 
  pull(Symbol_hg)



# Performant in coexpression and reverse coexpression... suggests targets are
# at extremes of rankings


both_coexpr_hg <- filter(
  rev_coexpr_hg, 
  (AUPRC_percentile_reverse_coexpr > 0.9 & AUPRC_percentile_coexpr > 0.9) |
    (AUROC_percentile_reverse_coexpr > 0.9 & AUROC_percentile_coexpr > 0.9)
)


both_coexpr_mm <- filter(
  rev_coexpr_mm, 
  (AUPRC_percentile_reverse_coexpr > 0.9 & AUPRC_percentile_coexpr > 0.9) |
    (AUROC_percentile_reverse_coexpr > 0.9 & AUROC_percentile_coexpr > 0.9)
)



# Curated target ranks
# ------------------------------------------------------------------------------

check_tf <- "LEF1"

# filter(bind_dat$Permissive_hg$Meta, Symbol == check_tf)

targets <- get_curated_labels(tf = check_tf,
                              curated_df = curated,
                              ortho_df = pc_ortho,
                              pc_df = pc_hg,
                              species = "Human",
                              remove_self = TRUE)


filter(agg_df_common_hg, Symbol == check_tf)
filter(agg_df_common_mm, Symbol == str_to_title(check_tf))


targets_rank_coexpr <- rank_tf_hg[[check_tf]] %>% 
  filter(Symbol %in% targets) %>% 
  arrange(Rank_RSR)


targets_rank_binding <- 
  data.frame(Bind = bind_summary$Human_TF[, check_tf]) %>% 
  rownames_to_column(var = "Symbol") %>% 
  mutate(Rank_bind = rank(-Bind, ties.method = "min")) %>% 
  filter(Symbol %in% targets) %>% 
  arrange(Rank_bind)


sum(targets_rank_coexpr$Rank_RSR <= 500)
sum(targets_rank_binding$Rank_bind <= 500)



# Plotting
# ------------------------------------------------------------------------------


# Histograms of Average versus individual AUCs

p1a <- plot_hist(avi_df_hg, stat_col = paste0("AUPRC", "_percentile"), xlab = "Percentile AUPRC", title = "Human")
p1b <- plot_hist(avi_df_hg, stat_col = paste0("AUROC", "_percentile"), xlab = "Percentile AUROC", title = "Human")
p1c <- plot_hist(avi_df_mm, stat_col = paste0("AUPRC", "_percentile"), xlab = "Percentile AUPRC", title = "Mouse")
p1d <- plot_hist(avi_df_mm, stat_col = paste0("AUROC", "_percentile"), xlab = "Percentile AUROC", title = "Mouse")

p1 <- plot_grid(p1a, p1b, p1c, p1d, nrow = 2)

ggsave(p1, height = 9, width = 12, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_average_vs_individual_hists.png")))



# Density plot of a TF's individual AUCs overlaid with the average AUC


plot_auc_density <- function(plot_df, vline, stat, tf) {
  
  xlim <- c(min(plot_df[[stat]], vline), max(plot_df[[stat]], vline))
  
  ggplot(plot_df, aes(x = !!sym(stat))) +
    geom_density(fill = "lightgrey") +
    geom_vline(xintercept = vline, col = "firebrick", linewidth = 1.4) +
    ylab("Density") +
    ggtitle(tf) +
    xlim(xlim) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          plot.margin = margin(c(10, 20, 10, 10)))
}



plot_tf2 <- "ASCL1"
plot_stat2 <- "AUROC"
plot_df2 <- avi_auc_hg[[plot_tf2]]$AUC_df
vline2 <- filter(plot_df2, ID == "Average")[[plot_stat2]]
plot_df2 <- filter(plot_df2, ID != "Average")

p2 <- plot_auc_density(plot_df2, vline = vline2, stat = plot_stat2, tf = plot_tf2)

ggsave(p2, height = 4, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, plot_tf2, "_coexpr_average_vs_individual_density.png")))



# Histogram of percentiles of aggregate AUCs versus null


# Using all TFs available for coexpression
plot_df3_hg <- agg_df_hg
plot_df3_mm <- agg_df_mm

# Using only TFs common to coexpression and unibind
# plot_df3_hg <- agg_df_common_hg
# plot_df3_mm <- agg_df_common_mm


p3a <- plot_hist(plot_df3_hg, stat_col = paste0("AUPRC_percentile_coexpr"), xlab = "AUPRC percentile coexpression", title = "Human") 
p3b <- plot_hist(plot_df3_mm, stat_col = paste0("AUPRC_percentile_coexpr"), xlab = "AUPRC percentile coexpression", title = "Mouse")

p3c <- plot_hist(plot_df3_hg, stat_col = paste0("AUROC_percentile_coexpr"), xlab = "AUROC percentile coexpression", title = "Human")
p3d <- plot_hist(plot_df3_mm, stat_col = paste0("AUROC_percentile_coexpr"), xlab = "AUROC percentile coexpression", title = "Mouse")

p3e <- plot_hist(plot_df3_hg, stat_col = paste0("AUPRC_percentile_unibind"), xlab = "AUPRC percentile binding", title = "Human")
p3f <- plot_hist(plot_df3_mm, stat_col = paste0("AUPRC_percentile_unibind"), xlab = "AUPRC percentile binding", title = "Mouse")

p3g <- plot_hist(plot_df3_hg, stat_col = paste0("AUROC_percentile_unibind"), xlab = "AUROC percentile binding", title = "Human")
p3h <- plot_hist(plot_df3_mm, stat_col = paste0("AUROC_percentile_unibind"), xlab = "AUROC percentile binding", title = "Mouse")


# TODO: fussing with fit in figure. when finalized clean up save 
# # Focus on AUPRC, showing both species and both aggregations
# p3 <- plot_grid(p3a, p3b, p3e, p3f, ncol = 1)


ggsave(p3a, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_perc_auprc_human.png")))

ggsave(p3b, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_perc_auprc_mouse.png")))

ggsave(p3c, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_perc_auroc_human.png")))

ggsave(p3d, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_perc_auroc_mouse.png")))

ggsave(p3e, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "binding_perc_auprc_human.png")))

ggsave(p3f, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "binding_perc_auprc_mouse.png")))

ggsave(p3g, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "binding_perc_auroc_human.png")))

ggsave(p3h, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "binding_perc_auroc_mouse.png")))




# Density plot of a TF's null AUCs overlaid with the aggregate AUC

plot_tf4 <- "ASCL1"
plot_stat4 <- "AUROC"
plot_l4 <- coexpr_auc_hg
plot_df4 <- data.frame(AUC = unlist(lapply(plot_l4[[plot_tf4]]$Null, `[[`, plot_stat4)))
colnames(plot_df4) <- plot_stat4
vline4 <- plot_l4[[plot_tf4]]$Perf_df[[plot_stat4]]

p4 <- plot_auc_density(plot_df4, vline = vline4, stat = plot_stat4, tf = plot_tf4)

ggsave(p4, height = 6, width = 9, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, plot_tf4, "_coexpr_aggregate_vs_null_density.png")))


# Scatterplot of coexpression versus binding percentile


p5_stat <- "AUROC"


p5 <- 
  ggplot(
  agg_df_common_hg,
  aes(x = !!sym(paste0(p5_stat, "_percentile_coexpr")),
      y = !!sym(paste0(p5_stat, "_percentile_unibind")))
) +
  
  # geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = 0.9, ymax = 1.05), fill = NA, colour = "forestgreen") +
  # geom_rect(aes(xmin = -0.05, xmax = 0.1, ymin = 0.9, ymax = 1.05), fill = NA, colour = "grey") +
  # geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = -0.05, ymax = 0.1), fill = NA, colour = "grey") +
  # geom_rect(aes(xmin = -0.05, xmax = 0.5, ymin = -0.05, ymax = 0.5), fill = NA, colour = "firebrick") +

  geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = 0.9, ymax = 1.05), fill = "forestgreen", colour = NA, alpha = 0.006) +
  geom_rect(aes(xmin = -0.05, xmax = 0.1, ymin = 0.9, ymax = 1.05), fill = "grey", colour = NA, alpha = 0.006) +
  geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = -0.05, ymax = 0.1), fill = "grey", colour = NA, alpha = 0.006) +
  geom_rect(aes(xmin = -0.05, xmax = 0.5, ymin = -0.05, ymax = 0.5), fill = "firebrick", colour = NA, alpha = 0.006) +

  geom_point(shape = 21, size = 3.4) +
  xlab(paste0(p5_stat, " percentile coexpression")) +
  ylab(paste0(p5_stat, " percentile binding")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))


ggsave(p5, height = 9, width = 9, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_vs_binding_", p5_stat, "_perc_scatter_human.png")))



# Relationship between AUC and percentile AUC
qplot(agg_df_hg, xvar = "AUPRC_coexpr", yvar = "AUPRC_percentile_coexpr")
qplot(agg_df_hg, xvar = "AUROC_coexpr", yvar = "AUROC_percentile_coexpr")

# Relationship between coexpression and binding AUC
qplot(agg_df_common_hg, xvar = "AUPRC_coexpr", yvar = "AUPRC_unibind")
qplot(agg_df_common_mm, xvar = "AUPRC_coexpr", yvar = "AUPRC_unibind")
qplot(agg_df_common_hg, xvar = "AUROC_coexpr", yvar = "AUROC_unibind")
qplot(agg_df_common_mm, xvar = "AUROC_coexpr", yvar = "AUROC_unibind")

# Scatterplot of count of targets/datasets and performance

pl_auprc_n <- list(
  qplot(agg_df_hg, yvar = "AUPRC_coexpr", xvar = "N_targets"),
  qplot(agg_df_hg, yvar = "AUPRC_coexpr", xvar = "N_coexpr_msr"),
  qplot(agg_df_hg, yvar = "AUPRC_percentile_coexpr", xvar = "N_targets"),
  qplot(agg_df_hg, yvar = "AUPRC_percentile_coexpr", xvar = "N_coexpr_msr")
)


pl_auroc_n <- list(
  qplot(agg_df_hg, yvar = "AUROC_coexpr", xvar = "N_targets"),
  qplot(agg_df_hg, yvar = "AUROC_coexpr", xvar = "N_coexpr_msr"),
  qplot(agg_df_hg, yvar = "AUROC_percentile_coexpr", xvar = "N_targets"),
  qplot(agg_df_hg, yvar = "AUROC_percentile_coexpr", xvar = "N_coexpr_msr")
)

p_auroc_n <- plot_grid(plotlist = pl_auroc_n, nrow = 2)


# Vbplot of AUC percentiles for coexpression and binding aggregations 
# TODO: better joining, function and save

# plot_df6 <- filter(rev_coexpr_hg, Symbol %in% agg_df_common_hg$Symbol)
plot_df6 <- filter(rev_coexpr_mm, Symbol %in% agg_df_common_mm$Symbol)


plot_df6 <- data.frame(
  Percentile = c(
    plot_df6$AUPRC_percentile_coexpr,
    plot_df6$AUPRC_percentile_unibind,
    plot_df6$AUPRC_percentile_reverse_coexpr,
    plot_df6$AUROC_percentile_coexpr,
    plot_df6$AUROC_percentile_unibind,
    plot_df6$AUROC_percentile_reverse_coexpr
  ),
  Group1 = c(
    rep("(+) Coexpression", length(plot_df6$AUPRC_percentile_coexpr)),
    rep("Binding", length(plot_df6$AUPRC_percentile_unibind)),
    rep("(-) Coexpression", length(plot_df6$AUPRC_percentile_unibind)),
    rep("(+) Coexpression", length(plot_df6$AUROC_percentile_coexpr)),
    rep("Binding", length(plot_df6$AUROC_percentile_unibind)),
    rep("(-) Coexpression", length(plot_df6$AUROC_ratio_reverse_coexpr))
  )
)

plot_df6$Group2 <- c(rep("AUPRC", nrow(plot_df6) / 2), rep("AUROC", nrow(plot_df6) / 2))
plot_df6$Group1 <- factor(plot_df6$Group1, levels = unique(plot_df6$Group1))


p6 <- ggplot(plot_df6, aes(x = Group1, y = Percentile)) +
  facet_wrap(~Group2, scales = "free") +
  geom_violin(fill = "slategrey") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ylab("AUC Percentile") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title.x = element_blank(),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 25),
        plot.title = element_text(size = 25))


ggsave(p6, height = 9, width = 18, device = "png", dpi = 300,
       # filename = file.path(paste0(plot_dir, "coexpr_vs_binding_percentile_vbplot_human.png")))
       filename = file.path(paste0(plot_dir, "coexpr_vs_binding_percentile_vbplot_mouse.png")))



# Plotting distribution of the median null versus observed AUROCs


p7a <- 
  data.frame(
  AUROC = c(agg_df_hg$AUROC_coexpr, null_coexpr_auroc_hg$Median),
  Group = c(rep("Aggregate", nrow(agg_df_hg)), 
            rep("Median null", nrow(agg_df_hg)))) %>% 
  ggplot(., aes(x = AUROC, fill = Group)) +
  geom_density() +
  ylab("Density") +
  ggtitle("Human coexpression") +
  scale_fill_manual(values = c("royalblue", "lightgrey")) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        plot.margin = margin(c(10, 10, 10, 10)))


p7b <- 
  data.frame(
  AUROC = c(agg_df_mm$AUROC_coexpr, null_coexpr_auroc_mm$Median),
  Group = c(rep("Aggregate", nrow(agg_df_mm)), 
            rep("Median null", nrow(agg_df_mm)))) %>% 
  ggplot(., aes(x = AUROC, fill = Group)) +
  geom_density() +
  ylab("Density") +
  ggtitle("Mouse coexpression") +
  scale_fill_manual(values = c("goldenrod", "lightgrey")) +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25),
        legend.position = c(0.8, 0.8),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        plot.margin = margin(c(10, 10, 10, 10)))


p7 <- plot_grid(p7a, p7b, ncol = 1)



ggsave(p7, height = 12, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, "coexpr_obs_vs_median_null_density.png"))



p_null_coexpr_auroc_scatter <-
  left_join(agg_df_hg, null_coexpr_auroc, by = "Symbol") %>%
  qplot(., xvar = "Median", yvar = "AUROC_coexpr") +
  geom_hline(yintercept = 0.5, col = "red") +
  geom_vline(xintercept = 0.5, col = "red") +
  ylab("Aggregate coexpression AUROC") +
  xlab("Median null coexpression AUROC")


p_null_coexpr_perc_auroc_scatter <-
  left_join(agg_df_hg, null_coexpr_auroc, by = "Symbol") %>%
  qplot(., xvar = "Median", yvar = "AUROC_percentile_coexpr") +
  geom_vline(xintercept = 0.5, col = "red") +
  ylab("Aggregate coexpression percentile AUROC") +
  xlab("Median null coexpression AUROC")



p_null_binding_auroc_scatter <-
  left_join(agg_df_hg, null_binding_auroc, by = "Symbol") %>%
  qplot(., xvar = "Median", yvar = "AUROC_unibind") +
  geom_hline(yintercept = 0.5, col = "red") +
  geom_vline(xintercept = 0.5, col = "red") +
  ylab("Aggregate binding AUROC") +
  xlab("Median null binding AUROC")



pl_delta_auc_hist <- list(
  plot_hist(data.frame(AUPRC_human = agg_auc_delta$AUPRC_human), stat_col = "AUPRC_human", xlab = "Coexpression - Binding AUPRC", title = "Human"),
  plot_hist(data.frame(AUROC_human = agg_auc_delta$AUROC_human), stat_col = "AUROC_human", xlab = "Coexpression - Binding AUROC", title = "Human"),
  plot_hist(data.frame(AUPRC_mouse = agg_auc_delta$AUPRC_mouse), stat_col = "AUPRC_mouse", xlab = "Coexpression - Binding AUPRC", title = "Mouse"),
  plot_hist(data.frame(AUROC_mouse = agg_auc_delta$AUROC_mouse), stat_col = "AUROC_mouse", xlab = "Coexpression - Binding AUROC", title = "Mouse")
)


p_delta_auc_hist <- plot_grid(plotlist = pl_delta_auc_hist, nrow = 2)



ggsave(p_delta_auc_hist, height = 9, width = 12, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "delta_coexpr_bind_histograms.png")))


# diff from null
pl_diff_auc_hist <- list(
  plot_hist(plot_df3_hg, stat_col = paste0("AUROC_diff_coexpr"), xlab = "AUROC diff coexpression", title = "Human") + geom_vline(xintercept = 0),
  plot_hist(plot_df3_hg, stat_col = paste0("AUROC_diff_unibind"), xlab = "AUROC diff binding", title = "Human") + geom_vline(xintercept = 0),
  plot_hist(plot_df3_mm, stat_col = paste0("AUROC_diff_coexpr"), xlab = "AUROC diff coexpression", title = "Mouse") + geom_vline(xintercept = 0),
  plot_hist(plot_df3_mm, stat_col = paste0("AUROC_diff_unibind"), xlab = "AUROC diff binding", title = "Mouse") + geom_vline(xintercept = 0)
)


p_diff_auc_hist <- plot_grid(plotlist = pl_diff_auc_hist, nrow = 2)



# Scatter of coexpr versus reverse coexpr
qplot(rev_coexpr_hg, xvar = "AUPRC_percentile_coexpr", yvar = "AUPRC_percentile_reverse_coexpr")
qplot(rev_coexpr_mm, xvar = "AUPRC_percentile_coexpr", yvar = "AUPRC_percentile_reverse_coexpr")
qplot(rev_coexpr_hg, xvar = "AUPRC_coexpr", yvar = "AUPRC_reverse_coexpr")
qplot(rev_coexpr_mm, xvar = "AUPRC_coexpr", yvar = "AUPRC_reverse_coexpr")

qplot(rev_coexpr_hg, xvar = "AUROC_percentile_coexpr", yvar = "AUROC_percentile_reverse_coexpr")
qplot(rev_coexpr_mm, xvar = "AUROC_percentile_coexpr", yvar = "AUROC_percentile_reverse_coexpr")
qplot(rev_coexpr_hg, xvar = "AUROC_coexpr", yvar = "AUROC_reverse_coexpr")
qplot(rev_coexpr_mm, xvar = "AUROC_coexpr", yvar = "AUROC_reverse_coexpr")


# Scatter of reverse coexpr versus binding
qplot(rev_coexpr_hg, xvar = "AUPRC_percentile_unibind", yvar = "AUPRC_percentile_reverse_coexpr")
qplot(rev_coexpr_mm, xvar = "AUPRC_percentile_unibind", yvar = "AUPRC_percentile_reverse_coexpr")
qplot(rev_coexpr_hg, xvar = "AUPRC_unibind", yvar = "AUPRC_reverse_coexpr")
qplot(rev_coexpr_mm, xvar = "AUPRC_unibind", yvar = "AUPRC_reverse_coexpr")

qplot(rev_coexpr_hg, xvar = "AUROC_percentile_unibind", yvar = "AUROC_percentile_reverse_coexpr")
qplot(rev_coexpr_mm, xvar = "AUROC_percentile_unibind", yvar = "AUROC_percentile_reverse_coexpr")
qplot(rev_coexpr_hg, xvar = "AUROC_unibind", yvar = "AUROC_reverse_coexpr")
qplot(rev_coexpr_mm, xvar = "AUROC_unibind", yvar = "AUROC_reverse_coexpr")


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


plot_df <- data.frame(
  AUPRC = c(
    rev_coexpr_hg$AUPRC_reverse_coexpr,
    rev_coexpr_hg$AUPRC_coexpr,
    rev_coexpr_hg$AUPRC_unibind
  ),
  Group = c(
    rep("Coexpr_reverse", length(rev_coexpr_hg$AUPRC_reverse_coexpr)),
    rep("Coexpr", length(rev_coexpr_hg$AUPRC_coexpr)),
    rep("Binding", length(rev_coexpr_hg$AUPRC_unibind))
  )
)

plot_df <- filter(plot_df, !is.na(AUPRC))

boxplot(plot_df$AUPRC ~ plot_df$Group)



plot_df <- data.frame(
  AUPRC = c(
    rev_coexpr_hg$AUPRC_percentile_reverse_coexpr,
    rev_coexpr_hg$AUPRC_percentile_coexpr,
    rev_coexpr_hg$AUPRC_percentile_unibind
  ),
  Group = c(
    rep("Coexpr_reverse", length(rev_coexpr_hg$AUPRC_percentile_reverse_coexpr)),
    rep("Coexpr", length(rev_coexpr_hg$AUPRC_percentile_coexpr)),
    rep("Binding", length(rev_coexpr_hg$AUPRC_percentile_unibind))
  )
)

plot_df <- filter(plot_df, !is.na(AUPRC))

boxplot(plot_df$AUPRC ~ plot_df$Group)




# Sketching heatmap of percentiles

a1 <- left_join(agg_df_common_hg, pc_ortho, by = c("Symbol" = "Symbol_hg"))
a2 <- left_join(agg_df_common_mm, pc_ortho, by = c("Symbol" = "Symbol_mm"))
a3 <- left_join(a1, a2, by = "ID", suffix = c("_human", "_mouse"))

keep_cols <- c("AUROC_percentile_coexpr_human",
               "AUROC_percentile_unibind_human",
               "AUROC_percentile_coexpr_mouse",
               "AUROC_percentile_unibind_mouse")


# keep_cols <- c("AUROC_diff_coexpr_human", 
#                "AUROC_diff_unibind_human", 
#                "AUROC_diff_coexpr_mouse", 
#                "AUROC_diff_unibind_mouse")

  
a4 <- a3[,keep_cols]
rownames(a4) <- a3$Symbol_human
a5 <- a4[complete.cases(a4), ]
a6 <- a5[order(rowMeans(a5), decreasing = TRUE), ]
colnames(a6) <- c("Coexpression_human", "Binding_human", "Coexpression_mouse", "Binding_mouse")


# pheatmap(t(a6), 
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          gaps_row = c(2),
#          border_col = "black",
#          color = colorRampPalette(c("#0571b0", "white", "#ca0020"))(100),
#          cellwidth = 8,
#          cellheight = 8,
#          filename = paste0(plot_dir, "ortho_performance_heatmap.png")
# )


pheatmap(a6, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(2),
         border_col = "black",
         color = colorRampPalette(c("#0571b0", "white", "#ca0020"))(100),
         cellwidth = 10,
         cellheight = 10,
         filename = paste0(plot_dir, "ortho_performance_heatmap.png")
)



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
