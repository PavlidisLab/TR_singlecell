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

# For each dataset, load the subset gene x TF aggregation matrix 
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

# Saved list RDS of the summaries
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Curated TFs with ChIP-seq and all targets for null
tfs_curated_hg <- intersect(tfs_hg$Symbol, str_to_upper(curated$TF_Symbol))
tfs_curated_mm <- intersect(tfs_mm$Symbol, str_to_title(curated$TF_Symbol))
targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))

# Rankings from Morin 2023
evidence_l <- readRDS(evidence_path)

# List of AUROC/AUPRC for TF lists ability to recover curated targets
unibind_auc_hg <- readRDS(unibind_auc_hg_path)
unibind_auc_mm <- readRDS(unibind_auc_mm_path)
coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)
avi_auc_hg <- readRDS(avg_vs_ind_auc_hg_path)
avi_auc_mm <- readRDS(avg_vs_ind_auc_mm_path)
rev_coexpr_auc_hg <- readRDS("/space/scratch/amorin/R_objects/reverse_coexpr_recover_curated_hg.RDS")
rev_coexpr_auc_mm <- readRDS("/space/scratch/amorin/R_objects/reverse_coexpr_recover_curated_mm.RDS")



# Demo a single TF
# ------------------------------------------------------------------------------


tf <- "ASCL1"

agg_l <- agg_tf_hg
msr_mat <- msr_hg
pc_df <- pc_hg
species <- "Human"

labels_curated <- get_curated_labels(tf = tf, 
                                     curated_df = curated, 
                                     pc_df = pc_df, 
                                     species = species, 
                                     remove_self = TRUE)


score_mat <- gene_vec_to_mat(agg_l, tf)
score_mat <- subset_to_measured(score_mat, msr_mat = msr_mat, gene = tf)
score_mat <- cbind(score_mat, Average = rowMeans(score_mat))

# No ties version: used for smoother plotting
score_mat_noties <- colrank_mat(-score_mat, ties_arg = "random")



## TODO: Looking at effect of imputation

score_mat_raw <- score_mat
rm(score_mat)

score_mat <- gene_vec_to_mat(agg_l, tf)
score_mat[tf, ] <- NA
score_mat <- subset_to_measured(score_mat, msr_mat = msr_mat, gene = tf)


# Logical matrix of whether TF-genes were co-measured in an experiment
comsr <- lapply(rownames(score_mat), function(y) {
  msr_mat[tf, colnames(score_mat)] & msr_mat[y, colnames(score_mat)]
})
comsr <- do.call(rbind, comsr)
rownames(comsr) <- rownames(score_mat)


# When a TF-gene was not co-measured, impute to the median NA
med <- median(score_mat, na.rm = TRUE)
score_mat_med <- score_mat
score_mat_med[comsr == FALSE] <- med

# NA to 0.5
score_mat05 <- score_mat
score_mat05[comsr == FALSE] <- 0.5


score_mat_med <- cbind(score_mat_med, Average = rowMeans(score_mat_med))
score_mat05 <- cbind(score_mat05, Average = rowMeans(score_mat05))

##

## TODO: binary heatmap

a1 <- rank_tf_hg$ASCL1 %>% 
  mutate(Curated = factor(Symbol %in% labels_curated)) %>% 
  arrange(desc(Rank_RSR)) %>% 
  mutate(Symbol = factor(Symbol, levels = unique(Symbol)))


pheatmap(a2, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         color = c("white", "forestgreen"), 
         cellwidth = 70,
         legend = FALSE,
         filename = "ascl1_binary_curated_heatmap.png",
         height = 30)


px <- 
  ggplot(a1, aes(y = Symbol, x = factor(1), fill = Curated)) +
  geom_tile(width = 2) +
  scale_fill_manual(values = c("white", "forestgreen")) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(c(10, 20, 10, 10)))
  

ggsave(px, height = 9, width = 3, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "ascl1_binary_curated_heatmap.png")))



##




# Inspecting retrieved coexpression ranks for the given TF for all networks


rank_l <- lapply(1:ncol(score_mat), function(x) {
  score_rank <- rank(-score_mat[, x], ties.method = "min")
  score_rank[labels_curated]
})

rank_mat <- do.call(cbind, rank_l)
colnames(rank_mat) <- colnames(score_mat)
rank_mat <- rank_mat[names(sort(rowSums(rank_mat))), ]


## TODO: different imputations


rank_mat_l <- 
  lapply(list(score_mat_raw, score_mat_med, score_mat05), function(mat) {
  
  rank_l <- lapply(1:ncol(mat), function(x) {
    score_rank <- rank(-mat[, x], ties.method = "min")
    score_rank[labels_curated]
  })
  
  rank_mat <- do.call(cbind, rank_l)
  colnames(rank_mat) <- colnames(mat)
  rank_mat <- rank_mat[names(sort(rowSums(rank_mat))), ]
  
  return(rank_mat)
})
names(rank_mat_l) <- c("Raw", "Imp_median", "Imp_0.5")



# view(do.call(cbind, lapply(rank_mat_l, function(x) x[, "Salcher2022"])))



##




# Inspect performance of a single dataset
# ------------------------------------------------------------------------------


id <- "Average"

score_vec <- sort(score_mat[, id], decreasing = TRUE)
score_rank <- rank(-score_mat[, id], ties.method = "min")
score_vec_noties <- sort(score_mat_noties[, id], decreasing = TRUE)
# stopifnot(identical(names(score_vec)[1:100], names(score_vec_noties)[1:100]))

label_vec <- names(score_vec) %in% labels_curated
label_vec_noties <- names(score_vec_noties) %in% labels_curated

# AUPRC and AUROC
auc <- get_auc(score_vec, label_vec, measure = "both")
auc_noties <- get_auc(score_vec_noties, label_vec_noties, measure = "both")

# Performance at every step
perf_df <- get_performance_df(score_vec, label_vec, measure = "both")
perf_df_noties <- get_performance_df(score_vec_noties, label_vec_noties, measure = "both")

# RSR and rank RSR of curated targets
rank_df <- data.frame(
  Symbol = labels_curated,
  RSR = score_vec[labels_curated],
  Rank_RSR = score_rank[labels_curated])


# Directly accessing ROCR prediction/plots
pred <- ROCR::prediction(predictions = score_vec, labels = label_vec)
pred_noties <- ROCR::prediction(predictions = score_vec_noties, labels = label_vec_noties)
roc <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
roc_noties <- ROCR::performance(pred_noties, measure = "tpr", x.measure = "fpr")
plot(roc)
plot(roc_noties)


# Inspect performance of all datasets
# ------------------------------------------------------------------------------


auc_all <- get_colwise_auc(score_mat, labels = labels_curated, ncores = 8)
auc_all_noties <- get_colwise_auc(score_mat_noties, labels = labels_curated, ncores = 8)

# Saved copy
auc_saved <- avi_auc_hg[[tf]]$AUC_df

auc_df_all <- get_colwise_performance_df(score_mat, labels = labels_curated, ncores = 8)
auc_df_all_noties <- get_colwise_performance_df(score_mat_noties, labels = labels_curated, ncores = 8)


cols <- c("Average" = "black",
          "Single_coexpression" = "lightgrey")


plot_df <- auc_df_all
plot_df$ID <- ifelse(plot_df$ID == "Average", plot_df$ID, "Single_coexpression")
plot_df$ID <- factor(plot_df$ID, levels = rev(unique(names(cols))))


plot_df_noties <- auc_df_all_noties
plot_df_noties$ID <- ifelse(plot_df_noties$ID == "Average", plot_df_noties$ID, "Single_coexpression")
plot_df_noties$ID <- factor(plot_df_noties$ID, levels = rev(unique(names(cols))))



p_avi_roc <- plot_perf(df = plot_df, measure = "ROC", cols = cols, title = tf)
p_avi_noties_roc <- plot_perf(df = plot_df_noties, measure = "ROC", cols = cols, title = tf)

p_avi_roc <- plot_perf(df = plot_df, measure = "PR", cols = cols, title = tf)
p_avi_noties_roc <- plot_perf(df = plot_df_noties, measure = "PR", cols = cols, title = tf)



# Demonstrate coexpr agg + perturb + binding + integrated recovery of curated
# ------------------------------------------------------------------------------


rank_cols <- c("Rank_perturbation", "Rank_binding", "Rank_integrated")


genomic_perf <- lapply(rank_cols, function(x) {
  df <- arrange(evidence_l$Human[[tf]], !!sym(x))
  score_vec <- 1 / (1:nrow(df))
  label_vec <- df$Symbol %in% labels_curated
  df <- get_performance_df(score_vec, label_vec, measure = "both")
  df$ID <- x
  return(df)
})


genomic_perf_df <- do.call(rbind, genomic_perf)

plot_df_genomic <- rbind(plot_df_noties, genomic_perf_df)


cols <- c("Average" = "black",
          "Rank_binding" = "darkblue",
          "Rank_perturbation" = "firebrick",
          "Rank_integrated" = "darkgreen",
          "Single_coexpression" = "lightgrey")


p_genomic_roc <- plot_perf(df = plot_df_genomic, measure = "ROC", cols = cols, title = tf)
p_genomic_pr <- plot_perf(df = plot_df_genomic, measure = "PR", cols = cols, title = tf)




# Demonstrating retrieval of genomic evidence at topk as labels
# ------------------------------------------------------------------------------


topk_labels <- evidence_l$Human[[tf]] %>%
  filter(Symbol != tf) %>%
  slice_min(Rank_integrated, n = 500) %>%
  pull(Symbol)


genomic_auc <- get_colwise_auc(score_mat,
                               topk_labels,
                               ncores = 8)


genomic_auc <- arrange(genomic_auc, AUROC)
genomic_auc_no_avg <- filter(genomic_auc, ID != "Average")
genomic_auc_avg <- filter(genomic_auc, ID == "Average")

hist(genomic_auc_no_avg$AUPRC, breaks = 100, xlim = c(0, 0.07))
abline(v = genomic_auc_avg$AUPRC, col = "red")

hist(genomic_auc_no_avg$AUROC, breaks = 100, xlim = c(0.3, 0.8))
abline(v = genomic_auc_avg$AUROC, col = "red")
