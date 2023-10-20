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


evidence_l <- readRDS(evidence_path)

# Curated TFs with ChIP-seq and all targets for null
# tfs_curated_hg <- intersect(colnames(bind_summary$Human_TF), str_to_upper(curated$TF_Symbol))
# tfs_curated_mm <- intersect(colnames(bind_summary$Mouse_TF), str_to_upper(curated$TF_Symbol))
# targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
# targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))

# For each dataset, load the subset gene x TF aggregation matrix 
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
# agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

# List of AUROC/AUPRC for TF lists ability to recover curated targets
unibind_auc_hg <- readRDS(unibind_auc_hg_path)
unibind_auc_mm <- readRDS(unibind_auc_mm_path)
coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)
avg_vs_ind_auc_hg <- readRDS(avg_vs_ind_auc_hg_path)
avg_vs_ind_auc_mm <- readRDS(avg_vs_ind_auc_mm_path)



# Functions
# ------------------------------------------------------------------------------

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




# Individual TF vector AUC vs averaged/final AUC
# ------------------------------------------------------------------------------


min_targets <- 5


avg_vs_ind_df_hg <-  lapply(avg_vs_ind_auc_hg, `[[`, "Summary_df") %>% 
  do.call(rbind, .) %>%
  rownames_to_column(var = "Symbol") %>% 
  filter(N_targets >= min_targets)


avg_vs_ind_df_mm <-  lapply(avg_vs_ind_auc_mm, `[[`, "Summary_df") %>% 
  do.call(rbind, .) %>%
  rownames_to_column(var = "Symbol") %>% 
  filter(N_targets >= min_targets)



summary(Filter(is.numeric, avg_vs_ind_df_hg))
summary(Filter(is.numeric, avg_vs_ind_df_mm))

top_and_bottom_prop(avg_vs_ind_df_hg, c("AUROC_percentile", "AUPRC_percentile"))
top_and_bottom_prop(avg_vs_ind_df_mm, c("AUROC_percentile", "AUPRC_percentile"))


plot_grid(
  plot_hist(avg_vs_ind_df_hg, stat_col = paste0("AUPRC", "_percentile")),
  plot_hist(avg_vs_ind_df_hg, stat_col = paste0("AUROC", "_percentile")),
  plot_hist(avg_vs_ind_df_mm, stat_col = paste0("AUPRC", "_percentile")),
  plot_hist(avg_vs_ind_df_mm, stat_col = paste0("AUROC", "_percentile")),
  nrow = 2)



tf <- "ASCL1"
auc_df <- arrange(avg_vs_ind_auc_hg[[tf]]$AUC_df, AUROC)
auc_df_no_avg <- filter(auc_df, ID != "Average")
auc_avg <- filter(auc_df, ID == "Average")


# hist(auc_df_no_avg$AUPRC, breaks = 100, xlim = c(0, max(auprc_agg, max(auprc_no_agg)) * 1.2))
hist(auc_df_no_avg$AUPRC, breaks = 100, xlim = c(0, 0.07))
abline(v = auc_avg$AUPRC, col = "red")

hist(auc_df_no_avg$AUROC, breaks = 20, xlim = c(0.3, 0.8))
abline(v = auc_avg$AUROC, col = "red")


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


plot_hist(auc_df_no_avg, stat_col = "AUROC") +
  geom_vline(xintercept = auc_avg$AUROC, col = "firebrick")


plot(density(auc_df_no_avg$AUROC))
abline(v = auc_avg$AUROC, col = "red")




#
# ------------------------------------------------------------------------------



# TODO: can keep be done upstream?

join_auc_df <- function(coexpr_l, unibind_l) {
  
  coexpr_df <- do.call(rbind, lapply(coexpr_l, `[[`, "Perf_df"))
  unibind_df <- do.call(rbind, lapply(unibind_l, `[[`, "Perf_df"))
  
  left_join(coexpr_df, unibind_df,
            by = c("Symbol", "N_targets"),
            suffix = c("_coexpr", "_unibind"))
}




auc_hg <- join_auc_df(coexpr_l = coexpr_auc_hg, unibind_l = unibind_auc_hg)
auc_mm <- join_auc_df(coexpr_l = coexpr_auc_mm, unibind_l = unibind_auc_mm)


# Keep only TFs with a minimum count of curated targets

min_targets <- 5
auc_sub_hg <- filter(auc_hg, N_targets >= min_targets)
auc_sub_mm <- filter(auc_mm, N_targets >= min_targets)

# Minimum count of curated targets common to unibind and coexpr

auc_common_hg <- filter(auc_sub_hg, !is.na(AUPRC_coexpr) & !is.na(AUPRC_unibind))
auc_common_mm <- filter(auc_sub_mm, !is.na(AUPRC_coexpr) & !is.na(AUPRC_unibind))

summary(Filter(is.numeric, auc_sub_hg))
summary(Filter(is.numeric, auc_sub_mm))

summary(Filter(is.numeric, auc_common_hg))
summary(Filter(is.numeric, auc_common_mm))


# Proportion with notably good or bad auc performance

stat_cols <- c(
  "AUPRC_percentile_observed_coexpr",
  "AUPRC_percentile_observed_unibind",
  "AUROC_percentile_observed_coexpr",
  "AUROC_percentile_observed_unibind"
)





top_and_bottom_prop(auc_sub_hg, stat_cols)
top_and_bottom_prop(auc_sub_mm, stat_cols)

top_and_bottom_prop(auc_common_hg, stat_cols)
top_and_bottom_prop(auc_common_mm, stat_cols)




# Histogram of percentiles of observered AUCs versus null


stat <- "AUPRC"

# Using all TFs available for coexpression
plot_df_hg <- auc_sub_hg 
plot_df_mm <- auc_sub_mm

# Using only TFs common to coexpression and unibind
# plot_df_hg <- auc_common_hg 
# plot_df_mm <- auc_common_mm


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

ggplot(auc_common_hg, aes(x = AUPRC_percentile_observed_coexpr, 
                          y = AUPRC_percentile_observed_unibind)) +
  geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = 0.9, ymax = 1.05), fill = NA, colour = "forestgreen") +
  geom_rect(aes(xmin = -0.05, xmax = 0.1, ymin = 0.9, ymax = 1.05), fill = NA, colour = "grey") +
  geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = -0.05, ymax = 0.1), fill = NA, colour = "grey") +
  geom_rect(aes(xmin = -0.05, xmax = 0.5, ymin = -0.05, ymax = 0.5), fill = NA, colour = "firebrick") +
  geom_point(shape = 21, size = 3.4) +
  xlab("AUPRC percentile coexpression") +
  ylab("AUPRC percentile binding") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))



ggplot(auc_common_hg, aes(x = AUPRC_percentile_observed_coexpr, 
                          y = AUPRC_percentile_observed_unibind)) +
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



# Looking at the distribution of observed percentiles for coexpr/binding
# TODO: function

plot(density(auc_common_hg$AUPRC_percentile_observed_coexpr))
lines(density(auc_common_hg$AUPRC_percentile_observed_unibind), col = "red")


plot(density(auc_sub_hg$AUPRC_coexpr))
lines(density(auc_sub_hg$AUPRC_unibind, na.rm = TRUE), col = "red")



# plot_perc_df <- auc_common_hg
plot_perc_df <- auc_sub_hg

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


ggsave(pb, height = 4, width = 8, device = "png", dpi = 600,
       filename = file.path(paste0(plot_dir, "binding_vs_agg_AUC_all_human.png")))

# ggsave(pb, height = 4, width = 8, device = "png", dpi = 600,
#        filename = file.path(paste0(plot_dir, "binding_vs_agg_AUC_common_human.png")))


# plot(auc_common_hg$AUPRC_percentile_observed_coexpr, auc_common_hg$AUPRC_percentile_observed_unibind)
# cor(auc_common_hg$AUPRC_percentile_observed_coexpr, auc_common_hg$AUPRC_percentile_observed_unibind, use = "pairwise.complete.obs", method = "spearman")
# 
# plot(auc_common_hg$AUROC_percentile_observed_coexpr, auc_common_hg$AUROC_percentile_observed_unibind)
# cor(auc_common_hg$AUROC_percentile_observed_coexpr, auc_common_hg$AUROC_percentile_observed_unibind, use = "pairwise.complete.obs", method = "spearman")
# 
# plot(auc_common_mm$AUPRC_percentile_observed_coexpr, auc_common_mm$AUPRC_percentile_observed_unibind)
# cor(auc_common_mm$AUPRC_percentile_observed_coexpr, auc_common_mm$AUPRC_percentile_observed_unibind, use = "pairwise.complete.obs", method = "spearman")
# 
# plot(auc_common_mm$AUPRC_percentile_observed_coexpr, auc_common_mm$AUPRC_percentile_observed_unibind)
# cor(auc_common_mm$AUPRC_percentile_observed_coexpr, auc_common_mm$AUPRC_percentile_observed_unibind, use = "pairwise.complete.obs", method = "spearman")


# Isolating TFs that performed well for coexpression and binding


top_both_hg <- auc_common_hg %>% 
  filter(AUPRC_percentile_observed_coexpr > 0.9 &
           AUPRC_percentile_observed_unibind > 0.9)


top_both_mm <- auc_common_mm %>% 
  filter(AUPRC_percentile_observed_coexpr > 0.9 &
           AUPRC_percentile_observed_unibind > 0.9)


top_ortho <- intersect(top_both_hg$Symbol, str_to_upper(top_both_mm$Symbol))



# Looking at TFs that performed well in one but not both 

diff_hg <- auc_common_hg %>%
  mutate(Diff_AUPRC = abs(
    AUPRC_percentile_observed_coexpr - AUPRC_percentile_observed_unibind
  )) 
# %>% filter(Diff_AUPRC > 0.8)


diff_mm <- auc_common_mm %>%
  mutate(Diff_AUPRC = abs(
    AUPRC_percentile_observed_coexpr - AUPRC_percentile_observed_unibind
  )) 
# %>% filter(Diff_AUPRC > 0.8)





# Inspecting retrieved ranks for a given TF


agg_l <- agg_tf_hg
msr_mat <- msr_hg
score_mat <- gene_vec_to_mat(agg_l, tf)
score_mat <- subset_to_measured(score_mat, msr_mat = msr_mat, gene = tf)
score_mat <- cbind(score_mat, Average = rowMeans(score_mat))
labels_curated <- get_curated_labels(tf = tf, curated_df = curated, pc_df = pc_hg, species = "Human", remove_self = TRUE)


rank_l <- lapply(1:ncol(score_mat), function(x) {
  score_rank <- rank(-score_mat[, x], ties.method = "min")
  # score_rank[labels_curated]
  # median(score_rank[labels_curated])
  # mean(score_rank[labels_curated])
  # sort(score_rank[labels_curated])[1:10]
  median(sort(score_rank[labels_curated])[1:10])
})
names(rank_l) <- colnames(score_mat)

rank_l[auc_df$ID]



# Inspecting aucormance of a specific dataset

id <- "GSE222956"
score_vec <- sort(score_mat[, id], decreasing = TRUE)
score_rank <- rank(-score_mat[, id], ties.method = "min")
label_vec <- names(score_vec) %in% labels_curated
auprc_all[id]
auroc_all[id]
get_auc(score_vec, label_vec, "both")
auc_df <- get_aucormance_df(score_vec, label_vec, measure = "both")


score_vec[labels_curated]
score_rank[labels_curated]


# Demonstrating retrieval of genomic evidence at topk as labels
# ------------------------------------------------------------------------------


topk_labels <- evidence_l$Human[[tf]] %>%
  filter(Symbol != tf) %>%
  slice_min(Rank_binding, n = 500) %>%
  pull(Symbol)


genomic_auc <- get_colwise_auc(score_mat,
                               topk_labels,
                               ncores = 8)


auc_df <- arrange(genomic_auc, AUROC)
auc_df_no_avg <- filter(auc_df, ID != "Average")
auc_avg <- filter(auc_df, ID == "Average")

hist(auc_df_no_avg$AUPRC, breaks = 100, xlim = c(0, 0.07))
abline(v = auc_avg$AUPRC, col = "red")

hist(auc_df_no_avg$AUROC, breaks = 100, xlim = c(0.3, 0.8))
abline(v = auc_avg$AUROC, col = "red")


# Demo a single TF
# ------------------------------------------------------------------------------


tf <- "ASCL1"
agg_l <- agg_tf_hg
msr_mat <- msr_hg
labels_curated <- get_curated_labels(tf = tf, curated_df = curated, pc_df = pc_hg, species = "Human", remove_self = TRUE)

score_mat <- gene_vec_to_mat(agg_l, tf)
score_mat <- subset_to_measured(score_mat, msr_mat = msr_mat, gene = tf)
score_mat <- score_mat[setdiff(rownames(score_mat), tf), ]
score_mat <- cbind(score_mat, Average = rowMeans(score_mat))

# Scoring where ties have been broken
score_mat_noties <- colrank_mat(-score_mat, ties_arg = "random")

# Single dataset
id <- "Average"
score_vec <- sort(score_mat[, id], decreasing = TRUE)
score_vec_noties <- sort(score_mat_noties[, id], decreasing = TRUE)
stopifnot(identical(names(score_vec)[1:100], names(score_vec_noties)[1:100]))

label_vec <- names(score_vec) %in% labels_curated
label_vec_noties <- names(score_vec_noties) %in% labels_curated

auc <- get_auc(score_vec, label_vec, measure = "both")
auc_noties <- get_auc(score_vec_noties, label_vec_noties, measure = "both")

auc_df <- get_performance_df(score_vec, label_vec, measure = "both")
auc_df_noties <- get_performance_df(score_vec_noties, label_vec_noties, measure = "both")


# Directly accessing ROCR prediction/plots
pred <- ROCR::prediction(predictions = score_vec, labels = label_vec)
pred_noties <- ROCR::prediction(predictions = score_vec_noties, labels = label_vec_noties)
roc <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
roc_noties <- ROCR::performance(pred_noties, measure = "tpr", x.measure = "fpr")
plot(roc)
plot(roc_noties)


# All datasets
auc_all <- get_colwise_auc(score_mat, labels = labels_curated, ncores = 8)
auc_all_noties <- get_colwise_auc(score_mat_noties, labels = labels_curated, ncores = 8)

auc_df_all <- get_colwise_performance_df(score_mat, labels = labels_curated, ncores = 8)
auc_df_all_noties <- get_colwise_performance_df(score_mat_noties, labels = labels_curated, ncores = 8)

auc_l <- split(auc_df_all, auc_df_all$ID)
auc_noties_l <- split(auc_df_all_noties, auc_df_all_noties$ID)



pdf("auc_l.pdf")


for (i in names(auc_l)) {
  id <- i
  score_vec <- sort(score_mat[, id], decreasing = TRUE)
  label_vec <- names(score_vec) %in% labels_curated
  pred <- ROCR::prediction(predictions = score_vec, labels = label_vec)
  roc <- ROCR::aucormance(pred, measure = "tpr", x.measure = "fpr")
  plot(roc, main = i)
  # plot(y = tt[[i]]$TPR, x = tt[[i]]$FPR, main = i)
}
graphics.off()


pdf("auc_noties_l.pdf")

for (i in names(auc_noties_l)) {
  id <- i
  score_vec_noties <- sort(score_mat_noties[, id], decreasing = TRUE)
  label_vec_noties <- names(score_vec_noties) %in% labels_curated
  pred_noties <- ROCR::prediction(predictions = score_vec_noties, labels = label_vec_noties)
  roc_noties <- ROCR::aucormance(pred_noties, measure = "tpr", x.measure = "fpr")
  plot(roc_noties, main = i)
  # plot(y = tt2[[i]]$TPR, x = tt2[[i]]$FPR, main = i)
}
graphics.off()




cols <- c("Average" = "black",
          # "Binding" = "darkblue",
          # "Perturbation" = "firebrick",
          # "Integrated" = "darkgreen",
          "Single_coexpression" = "lightgrey")


# cols <- names()

plot_df <- auc_df_all
plot_df$ID <- ifelse(plot_df$ID == "Average", plot_df$ID, "Single_coexpression")


plot_df_noties <- auc_df_all_noties
plot_df_noties$ID <- ifelse(plot_df_noties$ID == "Average", plot_df_noties$ID, "Single_coexpression")


plot_df$ID <- factor(plot_df$ID, levels = rev(unique(names(cols))))
plot_df_noties$ID <- factor(plot_df_noties$ID, levels = rev(unique(names(cols))))




p_all <- plot_perf(df = plot_df, measure = "ROC", cols = cols, title = tf)
p_all_noties <- plot_perf(df = plot_df_noties, measure = "ROC", cols = cols, title = tf)


pdf("myplot.pdf")

p_l <- lapply(names(auc_l), function(x) {
  
  p <- plot_auc(df = filter(auc_df_all, ID == x), measure = "ROC", cols = cols, title = x)
  print(p)
  
})

graphics.off()


pdf("myplot_noties.pdf")

p_l <- lapply(names(auc_noties_l), function(x) {
  
  p <- plot_auc(df = filter(auc_df_all_noties, ID == x), measure = "ROC", cols = cols, title = x)
  print(p)
  
})

graphics.off()






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
