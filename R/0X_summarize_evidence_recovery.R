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
unibind_perf_mm <- readRDS(unibind_recover_curated_hg_path)
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


ggplot(perf_sub_hg, aes(x = Percentile_observed_PR)) +
  geom_histogram(bins = 100) +
  ggtitle("Human") +
  xlab("AUPRC percentile observed greater than null") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




ggplot(perf_sub_hg, aes(x = Percentile_observed_ROC)) +
  geom_histogram(bins = 100) +
  ggtitle("Human") +
  xlab("AUROC percentile observed greater than null") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




hist(perf_sub_hg$Percentile_observed_PR, breaks = 100)
hist(perf_sub_hg$Percentile_observed_ROC, breaks = 100)

hist(perf_sub_mm$Percentile_observed_PR, breaks = 100)
hist(perf_sub_mm$Percentile_observed_ROC, breaks = 100)


sum(perf_sub_hg$Percentile_observed_PR == 1) / nrow(perf_sub_hg)
sum(perf_sub_mm$Percentile_observed_PR == 1) / nrow(perf_sub_mm)

sum(perf_sub_hg$Percentile_observed_ROC == 1) / nrow(perf_sub_hg)
sum(perf_sub_mm$Percentile_observed_ROC == 1) / nrow(perf_sub_mm)

sum(perf_sub_hg$Percentile_observed_PR > 0.9) / nrow(perf_sub_hg)
sum(perf_sub_mm$Percentile_observed_PR > 0.9) / nrow(perf_sub_mm)

sum(perf_sub_hg$Percentile_observed_ROC > 0.9) / nrow(perf_sub_hg)
sum(perf_sub_mm$Percentile_observed_ROC > 0.9) / nrow(perf_sub_mm)

sum(perf_sub_hg$Percentile_observed_PR < 0.1) / nrow(perf_sub_hg)
sum(perf_sub_mm$Percentile_observed_PR < 0.1) / nrow(perf_sub_mm)

sum(perf_sub_hg$Percentile_observed_ROC < 0.1) / nrow(perf_sub_hg)
sum(perf_sub_mm$Percentile_observed_ROC < 0.1) / nrow(perf_sub_mm)



# Demo a single TF


tf_hg <- "SREBF2"
tf_mm <- "Sp1"

# hist(curated_auprc_hg[[tf_hg]]$Null, xlim = c(0, curated_auprc_hg[[tf_hg]]$Perf_df$AUC * 1.5), main = tf_hg)
hist(curated_auprc_hg[[tf_hg]]$Null, main = tf_hg)
abline(v = curated_auprc_hg[[tf_hg]]$Perf_df$AUC, col = "red")


# hist(curated_auprc_mm[[tf_mm]]$Null, xlim = c(0, curated_auprc_mm[[tf_mm]]$Perf_df$AUC * 1.5), main = tf_mm)
hist(curated_auprc_mm[[tf_mm]]$Null, main = tf_mm)
abline(v = curated_auprc_mm[[tf_mm]]$Perf_df$AUC, col = "red")


perf_sub_hg %>% 
  mutate(Diff = abs(Percentile_observed_PR - Percentile_observed_ROC)) %>% 
  arrange(desc(Diff)) %>% 
  head()


ggplot(perf_sub_hg, aes(x = Percentile_observed_PR, y = Percentile_observed_ROC)) +
  geom_point() +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




##


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



tt <- summarize_obs_and_null_auc(tf = tf,
                                 score_vec = score_vec,
                                 label_vec = label_vec,
                                 label_all = targets_curated_mm,
                                 n_samps = 1000,
                                 ncores = 8)




tt <- curated_obs_and_null_auc(tf = tf,
                               rank_df = rank_tf_mm[[tf]],
                               score_col = "Avg_RSR",
                               curated_df = curated,
                               label_all = targets_curated_mm,
                               pc_df = pc_mm,
                               species = "Mouse",
                               n_samps = 1000,
                               ncores = 8)


tt <- curated_obs_and_null_auc_list(tfs = tfs_curated_mm[6:8],
                                    rank_l = rank_tf_mm,
                                    score_col = "Avg_RSR",
                                    curated_df = curated,
                                    label_all = targets_curated_mm,
                                    pc_df = pc_mm,
                                    species = "Mouse",
                                    n_samps = 1000,
                                    ncores = 8,
                                    verbose = TRUE)



save_curated_auc_list(path = "tt.RDS",
                      tfs = tfs_curated_mm[6:8],
                      rank_l = rank_tf_mm,
                      score_col = "Avg_RSR",
                      curated_df = curated,
                      label_all = targets_curated_mm,
                      pc_df = pc_mm,
                      species = "Mouse",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)


##




unibind_auprc_hg <- readRDS(unibind_auprc_hg_path)
unibind_auroc_hg <- readRDS(unibind_auroc_hg_path)
unibind_auprc_mm <- readRDS(unibind_auprc_mm_path)
unibind_auroc_mm <- readRDS(unibind_auroc_mm_path)


join_perf_df <- function(auprc_l, auroc_l) {
  
  auprc_df <- do.call(rbind, lapply(auprc_l, `[[`, "Perf_df"))
  auroc_df <- do.call(rbind, lapply(auroc_l, `[[`, "Perf_df"))
  
  left_join(auprc_df, auroc_df,
            by = c("Symbol", "N_targets"),
            suffix = c("_PR", "_ROC"))
}


unibind_df_hg <- join_perf_df(unibind_auprc_hg, unibind_auroc_hg)
unibind_df_mm <- join_perf_df(unibind_auprc_mm, unibind_auroc_mm)

unibind_df_sub_hg <- filter(unibind_df_hg, N_targets >= 5)
unibind_df_sub_mm <- filter(unibind_df_mm, N_targets >= 5)

summary(Filter(is.numeric, unibind_df_sub_hg))
summary(Filter(is.numeric, unibind_df_sub_mm))


ggplot(unibind_df_sub_hg, aes(x = Percentile_observed_PR)) +
  geom_histogram(bins = 100) +
  ggtitle("Human") +
  xlab("AUPRC percentile observed greater than null") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




ggplot(unibind_df_sub_hg, aes(x = Percentile_observed_ROC)) +
  geom_histogram(bins = 100) +
  ggtitle("Human") +
  xlab("AUROC percentile observed greater than null") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




### from unibind script



curated_auprc_hg_path <- "/space/scratch/amorin/R_objects/curated_auprc_hg.RDS"
curated_auprc_mm_path <- "/space/scratch/amorin/R_objects/curated_auprc_mm.RDS"
curated_auroc_hg_path <- "/space/scratch/amorin/R_objects/curated_auroc_hg.RDS"
curated_auroc_mm_path <- "/space/scratch/amorin/R_objects/curated_auroc_mm.RDS"


curated_auprc_hg <- readRDS(curated_auprc_hg_path)
curated_auroc_hg <- readRDS(curated_auroc_hg_path)
curated_auprc_mm <- readRDS(curated_auprc_mm_path)
curated_auroc_mm <- readRDS(curated_auroc_mm_path)


curated_df_hg <- join_perf_df(curated_auprc_hg, curated_auroc_hg)
curated_df_mm <- join_perf_df(curated_auprc_mm, curated_auroc_mm)


df_hg <- left_join(curated_df_hg, 
                   unibind_df_hg,
                   by = c("Symbol", "N_targets"),
                   suffix = c("_coexpr", "_binding"))



df_mm <- left_join(curated_df_mm, 
                   unibind_df_mm,
                   by = c("Symbol", "N_targets"),
                   suffix = c("_coexpr", "_binding"))



plot(df_hg$Percentile_observed_PR_coexpr, df_hg$Percentile_observed_PR_binding)
cor(df_hg$Percentile_observed_PR_coexpr, df_hg$Percentile_observed_PR_binding, use = "pairwise.complete.obs", method = "spearman")


plot(df_mm$Percentile_observed_PR_coexpr, df_mm$Percentile_observed_PR_binding)
cor(df_mm$Percentile_observed_PR_coexpr, df_mm$Percentile_observed_PR_binding, use = "pairwise.complete.obs", method = "spearman")


df_hg %>% 
  filter(N_targets >= 5 &
           Percentile_observed_PR_binding > 0.9 &
           Percentile_observed_PR_coexpr > 0.9)


df_mm %>% 
  filter(N_targets >= 5 &
           Percentile_observed_PR_binding > 0.9 &
           Percentile_observed_PR_coexpr > 0.9)



diff_hg <- df_hg %>%
  filter(N_targets >= 5) %>%
  mutate(Diff_PR = abs(
    Percentile_observed_PR_binding - Percentile_observed_PR_coexpr
  )) %>%
  filter(Diff_PR > 0.8)




diff_mm <- df_mm %>%
  filter(N_targets >= 5) %>%
  mutate(Diff_PR = abs(
    Percentile_observed_PR_binding - Percentile_observed_PR_coexpr
  )) %>%
  filter(Diff_PR > 0.8)



diff_ortho <- intersect(diff_hg$Symbol, str_to_upper(diff_mm$Symbol))



filter(df_hg, Symbol %in% diff_ortho)
filter(df_mm, Symbol %in% str_to_title(diff_ortho))
