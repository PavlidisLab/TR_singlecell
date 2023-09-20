## TODO:
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

# TODO: formalize unibind loading
# Processed list of meta and matrices
# bind_dat_path <- "/space/scratch/amorin/R_objects/processed_unibind_data.RDS"
# bind_dat <- readRDS(bind_dat_path)

# Average bind scores and output of binding specificity model
bind_summary_path <- "/space/scratch/amorin/R_objects/unibind_bindscore_summary.RDS"
bind_model_path <- "/space/scratch/amorin/R_objects/unibind_bindscore_modelfit.RDS"
bind_summary <- readRDS(bind_summary_path)
bind_model <- readRDS(bind_model_path)

# Curated TFs with ChIP-seq and all targets for null
# tfs_curated_hg <- intersect(colnames(bind_summary$Human_TF), str_to_upper(curated$TF_Symbol))
# tfs_curated_mm <- intersect(colnames(bind_summary$Mouse_TF), str_to_upper(curated$TF_Symbol))
# targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
# targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))

# For each dataset, load the subset gene x TF aggregation matrix 
# agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
# agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

# List of AUROC/AUPRC for TF lists ability to recover curated targets
unibind_perf_hg <- readRDS(unibind_recover_curated_hg_path)
unibind_perf_mm <- readRDS(unibind_recover_curated_mm_path)
coexpr_perf_hg <- readRDS(coexpr_recover_curated_hg_path)
coexpr_perf_mm <- readRDS(coexpr_recover_curated_mm_path)




#
# ------------------------------------------------------------------------------



# TODO: can keep be done upstream?

join_perf_df <- function(coexpr_l, unibind_l) {
  
  keep_coexpr <- names(which(unlist(lapply(coexpr_l, function(x) "Perf_df" %in% names(x)))))
  coexpr_df <- do.call(rbind, lapply(coexpr_l[keep_coexpr], `[[`, "Perf_df"))
  
  keep_unibind <- names(which(unlist(lapply(unibind_l, function(x) "Perf_df" %in% names(x)))))
  unibind_df <- do.call(rbind, lapply(unibind_l[keep_unibind], `[[`, "Perf_df"))
  
  left_join(coexpr_df, unibind_df,
            by = c("Symbol", "N_targets"),
            suffix = c("_coexpr", "_unibind"))
}




perf_hg <- join_perf_df(coexpr_l = coexpr_perf_hg, unibind_l = unibind_perf_hg)
perf_mm <- join_perf_df(coexpr_l = coexpr_perf_mm, unibind_l = unibind_perf_mm)


# Keep only TFs with a minimum count of curated targets

min_targets <- 5
perf_sub_hg <- filter(perf_hg, N_targets >= min_targets)
perf_sub_mm <- filter(perf_mm, N_targets >= min_targets)


summary(Filter(is.numeric, perf_sub_hg))
summary(Filter(is.numeric, perf_sub_mm))


# Minimum count of curated targets common to unibind and coexpr

perf_common_hg <- filter(perf_sub_hg, !is.na(AUPRC_coexpr) & !is.na(AUPRC_unibind))
perf_common_mm <- filter(perf_sub_mm, !is.na(AUPRC_coexpr) & !is.na(AUPRC_unibind))


summary(Filter(is.numeric, perf_common_hg))
summary(Filter(is.numeric, perf_common_mm))


# Proportion with notably good or bad performance

stat_cols <- c(
  "AUPRC_percentile_observed_coexpr",
  "AUPRC_percentile_observed_unibind",
  "AUROC_percentile_observed_coexpr",
  "AUROC_percentile_observed_unibind"
)


# TODO: function that returns proportion of data in percentile bins



top_and_bottom_prop <- function(df, stat_cols) {
  
  lapply(df[, stat_cols], function(col) {
    
    col <- col[!is.na(col)]
    n <- length(col)
    
    round(c(
      N = n,
      Eq1 = sum(col == 1) / n,
      Gt09 = sum(col > 0.9) / n,
      Lt01 = sum(col < 0.1) / n
      ), 3)
  })
}


top_and_bottom_prop(perf_sub_hg, stat_cols)
top_and_bottom_prop(perf_sub_mm, stat_cols)

top_and_bottom_prop(perf_common_hg, stat_cols)
top_and_bottom_prop(perf_common_mm, stat_cols)






plot_hist <- function(df, stat_col, title = NULL, xlab = NULL) {
  
  if (is.null(xlab)) xlab <- stat_col
  
  ggplot(df, aes(x = !!sym(stat_col))) +
    geom_histogram(bins = 100) +
    ggtitle(title) +
    xlab(xlab) +
    ylab("Frequency") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20))
  
}


# This uses all TFs for coexpr and unibind

stat <- "AUPRC"
plot_df_hg <- perf_sub_hg  # perf_common_hg 
plot_df_mm <- perf_sub_mm  # perf_common_mm

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






# Demo a single TF's performance relative to null distn


tf_hg <- "ASCL1"
tf_mm <- "Ascl1"

xmin <- 0.01
xmax <- 0.06

null <- unlist(lapply(unibind_perf_hg[[tf_hg]]$Null, `[[`, stat))
obs <- unibind_perf_hg[[tf_hg]]$Perf_df[[stat]]
hist(null, xlim = c(xmin, xmax), main = tf_hg, breaks = 100)
abline(v = obs, col = "red")


null <- unlist(lapply(unibind_perf_mm[[tf_mm]]$Null, `[[`, stat))
obs <- unibind_perf_mm[[tf_mm]]$Perf_df[[stat]]
hist(null, xlim = c(xmin, xmax), main = tf_mm, breaks = 100)
abline(v = obs, col = "red")



# 

ggplot(perf_common_hg, aes(x = AUPRC_percentile_observed_coexpr, y = AUPRC_percentile_observed_unibind)) +
  geom_point(shape = 21, size = 2.4) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))



# plot(perf_common_hg$AUPRC_percentile_observed_coexpr, perf_common_hg$AUPRC_percentile_observed_unibind)
# cor(perf_common_hg$AUPRC_percentile_observed_coexpr, perf_common_hg$AUPRC_percentile_observed_unibind, use = "pairwise.complete.obs", method = "spearman")
# 
# plot(perf_common_hg$AUROC_percentile_observed_coexpr, perf_common_hg$AUROC_percentile_observed_unibind)
# cor(perf_common_hg$AUROC_percentile_observed_coexpr, perf_common_hg$AUROC_percentile_observed_unibind, use = "pairwise.complete.obs", method = "spearman")
# 
# plot(perf_common_mm$AUPRC_percentile_observed_coexpr, perf_common_mm$AUPRC_percentile_observed_unibind)
# cor(perf_common_mm$AUPRC_percentile_observed_coexpr, perf_common_mm$AUPRC_percentile_observed_unibind, use = "pairwise.complete.obs", method = "spearman")
# 
# plot(perf_common_mm$AUPRC_percentile_observed_coexpr, perf_common_mm$AUPRC_percentile_observed_unibind)
# cor(perf_common_mm$AUPRC_percentile_observed_coexpr, perf_common_mm$AUPRC_percentile_observed_unibind, use = "pairwise.complete.obs", method = "spearman")



top_both_hg <- perf_common_hg %>% 
  filter(AUPRC_percentile_observed_coexpr > 0.9 &
           AUPRC_percentile_observed_unibind > 0.9)


top_both_mm <- perf_common_mm %>% 
  filter(AUPRC_percentile_observed_coexpr > 0.9 &
           AUPRC_percentile_observed_unibind > 0.9)


top_ortho <- intersect(top_both_hg$Symbol, str_to_upper(top_both_mm$Symbol))




diff_hg <- perf_common_hg %>%
  mutate(Diff_AUPRC = abs(
    AUPRC_percentile_observed_coexpr - AUPRC_percentile_observed_unibind
  )) 
# %>% filter(Diff_AUPRC > 0.8)


diff_mm <- perf_common_mm %>%
  mutate(Diff_AUPRC = abs(
    AUPRC_percentile_observed_coexpr - AUPRC_percentile_observed_unibind
  )) 
# %>% filter(Diff_AUPRC > 0.8)



# Individual vs average
# ------------------------------------------------------------------------------


df1 <- do.call(rbind, lapply(avg_vs_ind_recover_curated_hg, `[[`, "Summary_df"))
df2 <- do.call(rbind, lapply(avg_vs_ind_recover_curated_mm, `[[`, "Summary_df"))

df1_sub <- filter(df1, N_targets >= 5)
df2_sub <- filter(df2, N_targets >= 5)

summary(df1_sub$AUROC_percentile)
summary(df1_sub$AUPRC_percentile)
hist(df1_sub$AUROC_percentile, breaks = 50)
hist(df1_sub$AUPRC_percentile, breaks = 50)

summary(df2_sub$AUROC_percentile)
summary(df2_sub$AUPRC_percentile)
hist(df2_sub$AUROC_percentile, breaks = 50)
hist(df2_sub$AUPRC_percentile, breaks = 50)


tf <- "ASCL1"
tf <- "ASCL1"
auc_df <- avg_vs_ind_recover_curated_hg[[tf]]$AUC_df
auprc_all <- auc_df$AUPRC


auprc_no_agg <- sort(auprc_all[names(auprc_all) != "Average"])
auprc_agg <- auprc_all[["Average"]]
hist(auprc_no_agg, breaks = 100, xlim = c(0, max(auprc_agg, max(auprc_no_agg)) * 1.2))
abline(v = auprc_agg, col = "red")


auroc_all <- unlist(lapply(avg_vs_ind_recover_curated_hg[[tf]], `[[`, "AUROC"))
auroc_no_agg <- sort(auroc_all[names(auroc_all) != "Average"])
auroc_agg <- auroc_all[["Average"]]
hist(auroc_no_agg, xlim = c(0, max(auroc_agg, max(auroc_no_agg)) * 1.2))
abline(v = auroc_agg, col = "red")


# Inspecting retrieved ranks for a given TF


score_mat <- gene_vec_to_mat(agg_l, tf)
score_mat <- subset_to_measured(score_mat, msr_mat = msr_mat, gene = tf)
score_mat <- cbind(score_mat, Average = rowMeans(score_mat))


rank_l <- lapply(1:ncol(score_mat), function(x) {
  score_rank <- rank(-score_mat[, x], ties.method = "min")
  # score_rank[labels_curated]
  # median(score_rank[labels_curated])
  # mean(score_rank[labels_curated])
  # sort(score_rank[labels_curated])[1:10]
  median(sort(score_rank[labels_curated])[1:10])
})
names(rank_l) <- colnames(score_mat)

rank_l[c(names(auprc_no_agg), "Average")]
rank_l[c(names(auroc_no_agg), "Average")]

plot(unlist(rank_l[c(names(auprc_no_agg), "Average")]))
plot(unlist(rank_l[c(names(auroc_no_agg), "Average")]))

# Inspecting performance of a specific dataset

id <- "GSE222956"
score_vec <- sort(score_mat[, id], decreasing = TRUE)
score_rank <- rank(-score_mat[, id], ties.method = "min")
label_vec <- names(score_vec) %in% labels_curated
auprc_all[id]
auroc_all[id]
get_auc(score_vec, label_vec, "both")
perf_df <- get_performance_df(score_vec, label_vec, measure = "both")


score_vec[labels_curated]
score_rank[labels_curated]






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
# null <- get_null_performance(score_vec = score_vec,
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
