bind_to_df <- function(auprc_l, auroc_l) {
  
  auprc_df <- do.call(rbind, lapply(auprc_l, `[[`, "Perf_df"))
  auroc_df <- do.call(rbind, lapply(auroc_l, `[[`, "Perf_df"))
  
  left_join(auprc_df, auroc_df,
            by = c("Symbol", "N_targets"),
            suffix = c("_PR", "_ROC"))
}


curated_df_hg <- bind_to_df(curated_auprc_hg, curated_auroc_hg)
curated_df_mm <- bind_to_df(curated_auprc_mm, curated_auroc_mm)

curated_df_sub_hg <- filter(curated_df_hg, N_targets >= 5)
curated_df_sub_mm <- filter(curated_df_mm, N_targets >= 5)

summary(Filter(is.numeric, curated_df_sub_hg))
summary(Filter(is.numeric, curated_df_sub_mm))


ggplot(curated_df_sub_hg, aes(x = Percentile_observed_PR)) +
  geom_histogram(bins = 100) +
  ggtitle("Human") +
  xlab("AUPRC percentile observed greater than null") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




ggplot(curated_df_sub_hg, aes(x = Percentile_observed_ROC)) +
  geom_histogram(bins = 100) +
  ggtitle("Human") +
  xlab("AUROC percentile observed greater than null") +
  ylab("Frequency") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




hist(curated_df_sub_hg$Percentile_observed_PR, breaks = 100)
hist(curated_df_sub_hg$Percentile_observed_ROC, breaks = 100)

hist(curated_df_sub_mm$Percentile_observed_PR, breaks = 100)
hist(curated_df_sub_mm$Percentile_observed_ROC, breaks = 100)


sum(curated_df_sub_hg$Percentile_observed_PR == 1) / nrow(curated_df_sub_hg)
sum(curated_df_sub_mm$Percentile_observed_PR == 1) / nrow(curated_df_sub_mm)

sum(curated_df_sub_hg$Percentile_observed_ROC == 1) / nrow(curated_df_sub_hg)
sum(curated_df_sub_mm$Percentile_observed_ROC == 1) / nrow(curated_df_sub_mm)

sum(curated_df_sub_hg$Percentile_observed_PR > 0.9) / nrow(curated_df_sub_hg)
sum(curated_df_sub_mm$Percentile_observed_PR > 0.9) / nrow(curated_df_sub_mm)

sum(curated_df_sub_hg$Percentile_observed_ROC > 0.9) / nrow(curated_df_sub_hg)
sum(curated_df_sub_mm$Percentile_observed_ROC > 0.9) / nrow(curated_df_sub_mm)

sum(curated_df_sub_hg$Percentile_observed_PR < 0.1) / nrow(curated_df_sub_hg)
sum(curated_df_sub_mm$Percentile_observed_PR < 0.1) / nrow(curated_df_sub_mm)

sum(curated_df_sub_hg$Percentile_observed_ROC < 0.1) / nrow(curated_df_sub_hg)
sum(curated_df_sub_mm$Percentile_observed_ROC < 0.1) / nrow(curated_df_sub_mm)



# Demo a single TF


tf_hg <- "SREBF2"
tf_mm <- "Sp1"

# hist(curated_auprc_hg[[tf_hg]]$Null, xlim = c(0, curated_auprc_hg[[tf_hg]]$Perf_df$AUC * 1.5), main = tf_hg)
hist(curated_auprc_hg[[tf_hg]]$Null, main = tf_hg)
abline(v = curated_auprc_hg[[tf_hg]]$Perf_df$AUC, col = "red")


# hist(curated_auprc_mm[[tf_mm]]$Null, xlim = c(0, curated_auprc_mm[[tf_mm]]$Perf_df$AUC * 1.5), main = tf_mm)
hist(curated_auprc_mm[[tf_mm]]$Null, main = tf_mm)
abline(v = curated_auprc_mm[[tf_mm]]$Perf_df$AUC, col = "red")


curated_df_sub_hg %>% 
  mutate(Diff = abs(Percentile_observed_PR - Percentile_observed_ROC)) %>% 
  arrange(desc(Diff)) %>% 
  head()


ggplot(curated_df_sub_hg, aes(x = Percentile_observed_PR, y = Percentile_observed_ROC)) +
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


##




unibind_auprc_hg <- readRDS(unibind_auprc_hg_path)
unibind_auroc_hg <- readRDS(unibind_auroc_hg_path)
unibind_auprc_mm <- readRDS(unibind_auprc_mm_path)
unibind_auroc_mm <- readRDS(unibind_auroc_mm_path)


bind_to_df <- function(auprc_l, auroc_l) {
  
  auprc_df <- do.call(rbind, lapply(auprc_l, `[[`, "Perf_df"))
  auroc_df <- do.call(rbind, lapply(auroc_l, `[[`, "Perf_df"))
  
  left_join(auprc_df, auroc_df,
            by = c("Symbol", "N_targets"),
            suffix = c("_PR", "_ROC"))
}


unibind_df_hg <- bind_to_df(unibind_auprc_hg, unibind_auroc_hg)
unibind_df_mm <- bind_to_df(unibind_auprc_mm, unibind_auroc_mm)

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


curated_df_hg <- bind_to_df(curated_auprc_hg, curated_auroc_hg)
curated_df_mm <- bind_to_df(curated_auprc_mm, curated_auroc_mm)


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
