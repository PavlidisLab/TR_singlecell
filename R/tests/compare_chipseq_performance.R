## 
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

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Saved list RDS of the summaries
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Binding summaries
bind_dat <- readRDS("/space/scratch/amorin/R_objects/processed_unibind_data.RDS")
# bind_summary_perm <- readRDS("/space/scratch/amorin/R_objects/unibind_Permissive_bindscore_summary.RDS")
# bind_summary_rob <- readRDS("/space/scratch/amorin/R_objects/unibind_Robust_bindscore_summary.RDS")

# List of AUROC/AUPRC for TF lists ability to recover curated targets
unibind_auc_perm_hg <- readRDS("/space/scratch/amorin/R_objects/unibind_Permissive_recover_curated_hg.RDS")
unibind_auc_rob_hg <- readRDS("/space/scratch/amorin/R_objects/unibind_Robust_recover_curated_hg.RDS")
unibind_auc_perm_mm <- readRDS("/space/scratch/amorin/R_objects/unibind_Permissive_recover_curated_mm.RDS")
unibind_auc_rob_mm <- readRDS("/space/scratch/amorin/R_objects/unibind_Robust_recover_curated_mm.RDS")

coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)

# Rankings from Morin 2023
evidence_l <- readRDS(evidence_path)

# Curated TFs with ChIP-seq and all targets for null
tfs_curated_hg <- intersect(tfs_hg$Symbol, str_to_upper(curated$TF_Symbol))
tfs_curated_mm <- intersect(tfs_mm$Symbol, str_to_title(curated$TF_Symbol))
targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))



# 1: Count of TFs/experiments exclusive to permissive
# ------------------------------------------------------------------------------


perm_only_hg <- setdiff(bind_dat$Permissive_hg$Meta$ID, bind_dat$Robust_hg$Meta$ID)
perm_only_mm <- setdiff(bind_dat$Permissive_mm$Meta$ID, bind_dat$Robust_mm$Meta$ID)

setdiff(filter(bind_dat$Permissive_hg$Meta, Symbol == "RUNX1")$ID, 
        filter(bind_dat$Robust_hg$Meta, Symbol == "RUNX1")$ID)


intersect(filter(bind_dat$Permissive_hg$Meta, Symbol == "RUNX1")$File, 
          filter(bind_dat$Robust_hg$Meta, Symbol == "RUNX1")$File)


identical(bind_dat$Permissive_hg$Mat_raw[, "GSE64862.BONE.RUNX1.MA0002.2.damo.bed"],
          bind_dat$Robust_hg$Mat_raw[, "GSE64862.BONE.RUNX1.MA0002.2.damo.bed"])


count_hg <- left_join(
  count(bind_dat$Permissive_hg$Meta, Symbol, name = "Permissive"),
  count(bind_dat$Robust_hg$Meta, Symbol, name = "Robust"),
  by = "Symbol") %>% 
  mutate(Diff = Permissive - Robust)

sum(!is.na(count_hg$Robust))

count_mm <- left_join(
  count(bind_dat$Permissive_mm$Meta, Symbol, name = "Permissive"),
  count(bind_dat$Robust_mm$Meta, Symbol, name = "Robust"),
  by = "Symbol") %>% 
  mutate(Diff = Permissive - Robust)



# 2: Compare the aggregate rankings: does the inclusion of the permissive set
# of experiment result in an appreciable drop in group performance?
# Conclusion: TFs with data in both sets don't see a hit in performance, despite
# the permissive set allowing more experiments. TFs exclusive to the permissive 
# set show appreciably decreased performance as group relative to the robust TFs. 
# ------------------------------------------------------------------------------




auc_hg <- left_join(
  do.call(rbind, lapply(unibind_auc_perm_hg, `[[`, "Perf_df")),
  do.call(rbind, lapply(unibind_auc_rob_hg, `[[`, "Perf_df")),
  by = c("Symbol", "N_targets"),
  suffix = c("_permissive", "_robust")
) %>%
  filter(N_targets >= min_targets)



auc_hg <- lapply(coexpr_auc_hg, `[[`, "Perf_df") %>%
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rename_with(~paste0(., "_coexpr"), -c("Symbol", "N_targets")) %>% 
  left_join(auc_hg, by = c("Symbol", "N_targets")) %>% 
  filter(N_targets >= min_targets)


sum(!is.na(auc_hg$AUPRC_robust))



auc_mm <- left_join(
  do.call(rbind, lapply(unibind_auc_perm_mm, `[[`, "Perf_df")),
  do.call(rbind, lapply(unibind_auc_rob_mm, `[[`, "Perf_df")),
  by = c("Symbol", "N_targets"),
  suffix = c("_permissive", "_robust")
) %>%
  filter(N_targets >= min_targets)



auc_mm <- lapply(coexpr_auc_mm, `[[`, "Perf_df") %>%
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  rename_with(~paste0(., "_coexpr"), -c("Symbol", "N_targets")) %>% 
  left_join(auc_mm, by = c("Symbol", "N_targets")) %>% 
  filter(N_targets >= min_targets)



auc_hg$Group <- ifelse(is.na(auc_hg$AUPRC_robust), "Permissive_only", "Robust_only")
auc_mm$Group <- ifelse(is.na(auc_mm$AUPRC_robust), "Permissive_only", "Robust_only")


summary(Filter(is.numeric, filter(auc_hg, Group == "Permissive_only")))
summary(Filter(is.numeric, filter(auc_hg, Group == "Robust_only")))
summary(Filter(is.numeric, auc_hg))
summary(Filter(is.numeric, auc_mm))


qplot(auc_hg, xvar = "AUPRC_permissive", yvar = "AUPRC_robust")
qplot(auc_mm, xvar = "AUPRC_permissive", yvar = "AUPRC_robust")

qplot(auc_hg, xvar = "AUROC_permissive", yvar = "AUROC_robust")
qplot(auc_mm, xvar = "AUROC_permissive", yvar = "AUROC_robust")

qplot(auc_hg, xvar = "AUPRC_diff_permissive", yvar = "AUPRC_diff_robust")
qplot(auc_mm, xvar = "AUPRC_diff_permissive", yvar = "AUPRC_diff_robust")

qplot(auc_hg, xvar = "AUROC_diff_permissive", yvar = "AUROC_diff_robust")
qplot(auc_mm, xvar = "AUROC_diff_permissive", yvar = "AUROC_diff_robust")



a1 <- filter(auc_hg, Group == "Permissive_only")$AUROC_diff_permissive
a2 <- filter(auc_hg, Group == "Robust_only")$AUROC_diff_permissive
a3 <- filter(auc_hg, Group == "Robust_only")$AUROC_diff_robust



plot_df <- data.frame(
  Group = c(rep("Permissive_only", length(a1)),
            rep("Robust_1", length(a2)),
            rep("Robust_2", length(a3))),
  Delta_AUROC = c(a1, a2, a3)
)



ggplot(plot_df, aes(x = Group, y = Delta_AUROC)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



# Does considering only the robust set change the claim of coexpression being
# "on par" with the Unibind aggregates?



auc_delta_perm <- list(
  AUPRC_human = auc_hg$AUPRC_diff_coexpr - auc_hg$AUPRC_diff_permissive,
  AUPRC_mouse = auc_mm$AUPRC_diff_coexpr - auc_mm$AUPRC_diff_permissive,
  AUROC_human = auc_hg$AUROC_diff_coexpr - auc_hg$AUROC_diff_permissive,
  AUROC_mouse = auc_mm$AUROC_diff_coexpr - auc_mm$AUROC_diff_permissive
)


auc_delta_rob <- list(
  AUPRC_human = auc_hg$AUPRC_diff_coexpr - auc_hg$AUPRC_diff_robust,
  AUPRC_mouse = auc_mm$AUPRC_diff_coexpr - auc_mm$AUPRC_diff_robust,
  AUROC_human = auc_hg$AUROC_diff_coexpr - auc_hg$AUROC_diff_robust,
  AUROC_mouse = auc_mm$AUROC_diff_coexpr - auc_mm$AUROC_diff_robust
)


summ_delta_perm <- lapply(auc_delta_perm, summary)
summ_delta_rob <- lapply(auc_delta_rob, summary)


a4 <- filter(auc_hg, Group == "Permissive_only")$AUROC_diff_coexpr
a5 <- filter(auc_hg, Group == "Robust_only")$AUROC_diff_coexpr


plot_df2 <- data.frame(
  Group = c(rep("Coexpr_1", length(a4)),
            rep("Coexpr_2", length(a5))),
  Delta_AUROC = c(a4, a5)
)


plot_df2 <- rbind(plot_df, plot_df2)
plot_df2$Group <- factor(plot_df2$Group, 
                         levels = c("Coexpr_1", "Coexpr_2", "Permissive_only", "Robust_1", "Robust_2"))


ggplot(plot_df2, aes(x = Group, y = Delta_AUROC)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


qplot(auc_hg, xvar = "AUROC_diff_coexpr", yvar = "AUROC_diff_robust")
qplot(auc_hg, xvar = "AUROC_diff_coexpr", yvar = "AUROC_diff_permissive")





# 3: Performance of Morin 2023 ChIP-seq aggregates versus Unibind. Validation
# of ENCODE pipeline?
# ------------------------------------------------------------------------------


tfs_hg <- names(evidence_l$Human)
tfs_mm <- names(evidence_l$Mouse)


set.seed(5)


gr2023_hg <- lapply(tfs_hg, function(x) {
  
  rank_df <- mutate(evidence_l$Human[[x]], 
                    Rank_integrated = -Rank_integrated,
                    Rank_binding = -Rank_binding,
                    Rank_perturbation = -Rank_perturbation)

  int <- curated_obs_and_null_auc(tf = x,
                                  rank_df = rank_df,
                                  score_col = "Rank_integrated",
                                  curated_df = curated,
                                  label_all = targets_curated_hg,
                                  pc_df = pc_hg,
                                  species = "Human",
                                  n_samps = 1000,
                                  ncores = 8)
  
  bind <- curated_obs_and_null_auc(tf = x,
                                   rank_df = rank_df,
                                   score_col = "Rank_binding",
                                   curated_df = curated,
                                   label_all = targets_curated_hg,
                                   pc_df = pc_hg,
                                   species = "Human",
                                   n_samps = 1000,
                                   ncores = 8)
  
  pert <- curated_obs_and_null_auc(tf = x,
                                   rank_df = rank_df,
                                   score_col = "Rank_perturbation",
                                   curated_df = curated,
                                   label_all = targets_curated_hg,
                                   pc_df = pc_hg,
                                   species = "Human",
                                   n_samps = 1000,
                                   ncores = 8)
  
  data.frame(Integrated = int$Perf_df$AUROC_diff,
             Binding = bind$Perf_df$AUROC_diff,
             Perturbation = pert$Perf_df$AUROC_diff)
  
  
})
names(gr2023_hg) <- tfs_hg



gr2023_mm <- lapply(tfs_mm, function(x) {
  
  rank_df <- mutate(evidence_l$Mouse[[x]], 
                    Rank_integrated = -Rank_integrated,
                    Rank_binding = -Rank_binding,
                    Rank_perturbation = -Rank_perturbation)
  
  int <- curated_obs_and_null_auc(tf = x,
                                  rank_df = rank_df,
                                  score_col = "Rank_integrated",
                                  curated_df = curated,
                                  label_all = targets_curated_mm,
                                  pc_df = pc_mm,
                                  species = "Mouse",
                                  n_samps = 1000,
                                  ncores = 8)
  
  bind <- curated_obs_and_null_auc(tf = x,
                                   rank_df = rank_df,
                                   score_col = "Rank_binding",
                                   curated_df = curated,
                                   label_all = targets_curated_mm,
                                   pc_df = pc_mm,
                                   species = "Mouse",
                                   n_samps = 1000,
                                   ncores = 8)
  
  pert <- curated_obs_and_null_auc(tf = x,
                                   rank_df = rank_df,
                                   score_col = "Rank_perturbation",
                                   curated_df = curated,
                                   label_all = targets_curated_mm,
                                   pc_df = pc_mm,
                                   species = "Mouse",
                                   n_samps = 1000,
                                   ncores = 8)
  
  data.frame(Integrated = int$Perf_df$AUROC_diff,
             Binding = bind$Perf_df$AUROC_diff,
             Perturbation = pert$Perf_df$AUROC_diff)
  
  
})
names(gr2023_mm) <- tfs_mm




gr2023_df_hg <- do.call(rbind, gr2023_hg) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Symbol") %>% 
  # rename_with(~paste0(., "_GR2023"), -c("Symbol")) %>%
  left_join(auc_hg, by = "Symbol")



gr2023_df_mm <- do.call(rbind, gr2023_mm) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Symbol") %>% 
  # rename_with(~paste0(., "_GR2023"), -c("Symbol")) %>%
  left_join(auc_mm, by = "Symbol")





keep_cols <- c(
  "Integrated",
  "Binding",
  "Perturbation",
  "AUROC_diff_coexpr",
  "AUROC_diff_permissive",
  "AUROC_diff_robust"
)



rename_cols <- c(
  "Integrated",
  "Binding",
  "Perturbation",
  "Coexpression",
  "Unibind_permissive",
  "Unibind_robust"
)


mat_hg <- gr2023_df_hg[, keep_cols]
rownames(mat_hg) <- gr2023_df_hg$Symbol
colnames(mat_hg) <- rename_cols


pheatmap(mat_hg, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         na_col = "black",
         border_color = NA,
         gaps_col = c(3, 4),
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 15)



mat_mm <- gr2023_df_mm[, keep_cols]
rownames(mat_mm) <- gr2023_df_mm$Symbol
colnames(mat_mm) <- rename_cols

pheatmap(mat_mm, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         na_col = "black",
         border_color = NA,
         gaps_col = c(3, 4),
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 15)





# gr2023_df_hg <- lapply(gr2023_hg, `[[`, "Perf_df") %>%
#   do.call(rbind, .) %>% 
#   as.data.frame() %>% 
#   rename_with(~paste0(., "_GR2023"), -c("Symbol")) %>% 
#   left_join(auc_hg, by = "Symbol")
# 
# 
# 
# gr2023_df_mm <- lapply(gr2023_mm, `[[`, "Perf_df") %>%
#   do.call(rbind, .) %>% 
#   as.data.frame() %>% 
#   rename_with(~paste0(., "_GR2023"), -c("Symbol")) %>% 
#   left_join(auc_mm, by = "Symbol")

# qplot(gr2023_df_hg, xvar = "AUROC_GR2023", yvar = "AUROC_robust")
# qplot(gr2023_df_hg, xvar = "AUROC_diff_GR2023", yvar = "AUROC_diff_robust")
# 
# 
# auc_delta_gr2023 <- list(
#   AUPRC_human = gr2023_df_hg$AUPRC_diff_GR2023 - gr2023_df_hg$AUPRC_diff_robust,
#   AUPRC_mouse = gr2023_df_mm$AUPRC_diff_GR2023 - gr2023_df_mm$AUPRC_diff_robust,
#   AUROC_human = gr2023_df_hg$AUROC_diff_GR2023 - gr2023_df_hg$AUROC_diff_robust,
#   AUROC_mouse = gr2023_df_mm$AUROC_diff_GR2023 - gr2023_df_mm$AUROC_diff_robust
# )
