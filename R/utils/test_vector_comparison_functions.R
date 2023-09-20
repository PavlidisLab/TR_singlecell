# Example of ability of individual aggregate coexpression vectors + 
# averaged vectors to recovery curated targets


tf <- tfs_curated_mm[15]
agg_l <- agg_tf_mm
msr_mat <- msr_mm


score_mat <- gene_vec_to_mat(agg_l, tf)
score_mat <- subset_to_measured(score_mat, msr_mat = msr_mat, gene = tf)
score_mat <- cbind(score_mat, Average = rowMeans(score_mat))


# curated <- data.frame(TF_Symbol = tf, Target_Symbol = topk_sort(score_mat[, 1], 40))
labels_curated <- get_curated_labels(tf = tf, curated_df = curated, pc_df = pc_mm, species = "Mouse", remove_self = TRUE)

all_agg_auc <- get_colwise_auc(score_mat, labels_curated, ncores = 8)


all_agg_auc2 <- get_colwise_curated_auc_list(tfs = tf,
                                             agg_l = agg_tf_mm,
                                             msr_mat = msr_mm,
                                             curated_df = curated,
                                             pc_df = pc_mm,
                                             species = "Mouse",
                                             ncores = 8,
                                             verbose = TRUE)


identical(all_agg_auc, all_agg_auc2[[tf]]$AUC_df)
