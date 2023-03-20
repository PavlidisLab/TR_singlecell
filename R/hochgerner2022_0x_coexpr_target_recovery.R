library(tidyverse)
library(Seurat)
library(cowplot)
library(parallel)
library(WGCNA)
library(ggrepel)
library(pheatmap)
library(ROCR)
source("R/utils/functions.R")
source("R/00_config.R")

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)

# Load correlation matrices generated per cell type and across all cells
cor_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_celltype.RDS")
# cor_all <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_all.RDS")
cor_all <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_amygdala_biorvx_cormat.RDS")

# Ribosomal genes as positive control for coexpr.
# Something wonky with table shifts columns over - just want gene symbols
ribo_genes <- read.delim("/space/grp/amorin/Metadata/MGI_GO_term_ribosomal_genes_15-03-2023.txt", row.names = NULL, stringsAsFactors = FALSE)
ribo_genes <- intersect(ribo_genes$MGI.Gene.Marker.ID, rownames(sdat))


# 
# ------------------------------------------------------------------------------


# # Given a list of cor matrices per cell type and a specified TF, get that TFs
# # gene cors in a gene x cell type matrix
# 
# tf_by_ct_cmat <- function(cor_list, tf, rm_tf = TRUE) {
#   
#   stopifnot(identical(rownames(cor_list[[1]]), rownames(cor_list[[2]])))
#   
#   cor_tf <- do.call(cbind, lapply(cor_ct, function(x) x[tf, ]))
#   
#   if (rm_tf) {
#     cor_tf <- cor_tf[rownames(cor_tf) != tf, ]
#   }
#   
#   return(cor_tf)
# }
# 
# 
# # Convert cor matrix to column-wise (cell-type) ranks such that 1 is best
# 
# rank_cormat <- function(cmat) {
#   
#   rank_mat <- apply(-cmat, 2, rank, ties.method = "min")
#   
#   return(rank_mat)
# }
# 
# 
# # 
# 
# join_evidence <- function(rank_list, cmat, tf) {
#   
#   stopifnot(tf %in% rownames(cmat), tf %in% names(rank_list))
#   
#   cor_df <- data.frame(Cor = cmat[tf, ]) %>% rownames_to_column(var = "Symbol")
#   
#   target_df <- left_join(target_df, cor_df, by = "Symbol")
#   
#   return(target_df)
# }



get_perf_df <- function(rank_df, label_col, measure = NULL) {
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("ROC", "PR"))
  
  # positives/has evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(rank_df[[label_col]]))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  
  if (measure == "ROC") {
    perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
    perf_df <- data.frame(TPR = unlist(perf@y.values),
                          FPR = unlist(perf@x.values))
  } else {
    perf <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
    perf_df <- data.frame(Precision = unlist(perf@y.values),
                          Recall = unlist(perf@x.values))
  }
  
  return(perf_df)
}



get_au_perf <- function(rank_df, label_col, measure = NULL) {
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("AUC", "AUPRC"))
  
  # positives/has evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(rank_df[[label_col]]))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  
  if (measure == "AUC") {
    perf <- ROCR::performance(pred, measure = "auc")@y.values[[1]]
  } else {
    perf <- ROCR::performance(pred, measure = "aucpr")@y.values[[1]]
  }
  
  return(perf)
}


plot_auc <- function(auc_df_list) {
  
  # auc_labels <- c(
  #   Integrated = paste0("Integrated (AUC=", signif(auc_df[tf, "Integrated"], 3), ")"),
  #   Perturbation = paste0("Perturbation (AUC=", signif(auc_df[tf, "Perturbation"], 3), ")"),
  #   Binding = paste0("Binding (AUC=", signif(auc_df[tf, "Binding"], 3), ")"))
  
  p <- 
    ggplot() +
    geom_path(data = auc_df_list[[1]], aes(x = FPR, y = TPR), linewidth = 1) +
    # geom_path(data = pr_list$Perturbation, aes(x = Recall, y = Precision, col = "Perturbation"), size = 1) +
    # geom_path(data = pr_list$Binding, aes(x = Recall, y = Precision, col = "Binding"), size = 1) +
    # ggtitle(tf) +
    # scale_color_manual(labels = auc_labels, values = colours) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 30),
          legend.title = element_blank(),
          legend.text = element_text(size = 25),
          legend.position = c(0.55, 0.85))
  
  if (length(auc_df_list) > 1) {
    for (i in 2:length(auc_df_list)) {
      p <- p + geom_path(data = auc_df_list[[i]], 
                         aes(x = FPR, y = TPR), 
                         linewidth = 0.5,
                         colour = "lightgrey")
    }
  }
  
  return(p)
}


plot_auprc <- function(auprc_df) {
  
  # auc_labels <- c(
  #   Integrated = paste0("Integrated (AUC=", signif(auc_df[tf, "Integrated"], 3), ")"),
  #   Perturbation = paste0("Perturbation (AUC=", signif(auc_df[tf, "Perturbation"], 3), ")"),
  #   Binding = paste0("Binding (AUC=", signif(auc_df[tf, "Binding"], 3), ")"))
  
  ggplot() +
    geom_path(data = auprc_df, aes(x = Recall, y = Precision), linewidth = 1) +
    # geom_path(data = pr_list$Perturbation, aes(x = Recall, y = Precision, col = "Perturbation"), size = 1) +
    # geom_path(data = pr_list$Binding, aes(x = Recall, y = Precision, col = "Binding"), size = 1) +
    # ggtitle(tf) +
    # scale_color_manual(labels = auc_labels, values = colours) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 30),
          legend.title = element_blank(),
          legend.text = element_text(size = 25),
          legend.position = c(0.55, 0.85))
  
}




tf <- "Mecp2"


# All cor

cor_all_tf <- data.frame(Cor_all = cor_all[, tf]) %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol != tf) %>%
  mutate(
    Cor_all_abs = abs(Cor_all),
    Rank_cor_all = rank(-Cor_all, ties.method = "min"),
    Rank_cor_all_abs = rank(-Cor_all_abs, ties.method = "min"))
  

evidence <- left_join(cor_all_tf, rank_l$Mouse[[tf]], by = "Symbol")

evidence <- filter(evidence, !is.na(Rank_integrated))



# Precision recall dfs

rank_int_pr <- get_perf_df(
  rank_df = arrange(evidence, Rank_integrated),
  label_col = "Curated_target",
  measure = "PR")

cor_all_pr <- get_perf_df(
  rank_df = arrange(evidence, Rank_cor_all),
  label_col = "Curated_target",
  measure = "PR")

cor_all_abs_pr <- get_perf_df(
  rank_df = arrange(evidence, Rank_cor_all_abs),
  label_col = "Curated_target",
  measure = "PR")


# ROC dfs

rank_int_roc <- get_perf_df(
  rank_df = arrange(evidence, Rank_integrated),
  label_col = "Curated_target",
  measure = "ROC")

cor_all_roc <- get_perf_df(
  rank_df = arrange(evidence, Rank_cor_all),
  label_col = "Curated_target",
  measure = "ROC")

cor_all_abs_roc <- get_perf_df(
  rank_df = arrange(evidence, Rank_cor_all_abs),
  label_col = "Curated_target",
  measure = "ROC")


# Area under the PR curve

rank_int_auprc <- get_au_perf(
  rank_df = arrange(evidence, Rank_integrated),
  label_col = "Curated_target",
  measure = "AUPRC")

cor_all_auprc <- get_au_perf(
  rank_df = arrange(evidence, Rank_cor_all),
  label_col = "Curated_target",
  measure = "AUPRC")

cor_all_abs_auprc <- get_au_perf(
  rank_df = arrange(evidence, Rank_cor_all_abs),
  label_col = "Curated_target",
  measure = "AUPRC")


# Area under the ROC

rank_int_auc <- get_au_perf(
  rank_df = arrange(evidence, Rank_integrated),
  label_col = "Curated_target",
  measure = "AUC")

cor_all_auc <- get_au_perf(
  rank_df = arrange(evidence, Rank_cor_all),
  label_col = "Curated_target",
  measure = "AUC")

cor_all_abs_auc <- get_au_perf(
  rank_df = arrange(evidence, Rank_cor_all_abs),
  label_col = "Curated_target",
  measure = "AUC")



plot_auprc(rank_int_pr)
plot_auprc(cor_all_pr)
plot_auprc(cor_all_abs_pr)


# plot_auc(rank_int_roc)
# plot_auc(cor_all_roc)
# plot_auc(cor_all_abs_roc)


plot_auc(list(rank_int_roc, cor_all_roc, cor_all_abs_roc))


# Inspecting NAs
sum(is.na(evidence$Rank_integrated))
filter(evidence, is.na(Rank_integrated)) %>% arrange(desc(Cor_all_abs)) %>% head(20)
sum(is.na(evidence$Cor_all))

# Example of high cor that is suspicious
FeaturePlot(sdat, features = c("Runx1", "Mnd1"))

plot(sdat@assays$RNA@data["Runx1", sdat$Cell_type == "GABA-11-Adora2a-Id4"],
     sdat@assays$RNA@data["Mnd1", sdat$Cell_type == "GABA-11-Adora2a-Id4"])

sum(sdat$Cell_type == "GABA-11-Adora2a-Id4")


tf <- "Runx1"
# Cell type cor
tt <- tf_by_ct_cmat(cor_ct, tf)
tt2 <- rank_cormat(tt)
cor_tf <- do.call(cbind, lapply(cor_ct, function(x) x[tf, ]))
# cor_tf <- abs(do.call(cbind, lapply(cor_ct, function(x) x[tf, ])))



rank_tf <- apply(-cor_tf, 2, rank, ties.method = "min")
rank_agg <- rowSums(rank_tf)

rank_order <- data.frame(Cor_rank = sort(rank(rank_agg))) %>% 
  rownames_to_column(var = "Symbol")


target_rank <- left_join(rank_l$Mouse[[tf]], rank_order, by = "Symbol")

boxplot(target_rank$Cor_rank ~ target_rank$Curated_target)


# Some NAs refer to genes not expressed in that cell type
gene_nas <- apply(cor_tf, 1, function(x) sum(is.na(x)))
ct_nas <- apply(cor_tf, 2, function(x) sum(is.na(x)))

# All NAs correspond to cell types where TF is not expressed
all_nas <- names(ct_nas[ct_nas == nrow(sdat@assays$RNA)])
all(sdat@assays$RNA@data[tfs[1], sdat$Cell_type %in% all_nas] == 0)




cor_ribo_in <- cor_all[ribo_genes, ribo_genes]
cor_ribo_out <- cor_all[ribo_genes, setdiff(rownames(cor_all), ribo_genes)]

cor_ribo_df <- rbind(
  data.frame(mat_to_df(cor_ribo_in), Group = "In"),
  data.frame(mat_to_df(cor_ribo_out), Group = "Out")
)

boxplot(cor_ribo_df$Value ~ cor_ribo_df$Group)


px1 <- 
  ggplot(cor_ribo_df, aes(x = Value, fill = Group)) +
  geom_density(alpha = 0.6) +
  theme_classic() +
  ylab("Density") +
  xlab("Pearson's correlation") +
  ggtitle("Ribosomal gene cor") +
  scale_fill_manual(values = c("deepskyblue4", "lightgrey")) +
  theme(
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 30),
    plot.title = element_text(hjust = 0.5, size = 30),
    legend.position = c(0.75, 0.90),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    plot.margin = margin(10, 20, 10, 10))


px2 <- 
  ggplot(cor_ribo_df, aes(x = Group, y = Value)) +
  geom_violin(width = 0.4, fill = "lightslategrey") +
  geom_boxplot(width = 0.1, fill = "white") +
  ylab("Pearson's correlation") +
  ggtitle("Ribosomal gene cor") +
  theme_classic() +
  theme(axis.title = element_text(size = 30),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 27),
        plot.title = element_text(size = 30, hjust = 0.5))


pheatmap(cor_ribo_in,
         clustering_distance_cols = as.dist(1 - cor_ribo_in),
         clustering_distance_rows = as.dist(1 - cor_ribo_in))



cor_all["Ascl1", c("Dll1", "Dll3", "Dll4")]


gene <- "Pax6"
nongene <- setdiff(names(rank_l$Mouse), gene)
# head(sort(cor_all[gene, ], decreasing = TRUE), 20)


cor_df <- data.frame(SC_cor = cor_all[gene, ]) %>% 
  rownames_to_column(var = "Symbol") %>% 
  left_join(rank_l$Mouse[[gene]], by = "Symbol") %>% 
  filter(Symbol != gene) %>% 
  filter(!is.na(Rank_integrated) & !is.na(SC_cor)) %>% 
  arrange(Rank_integrated) %>% 
  mutate(Group = Rank_integrated <= 500)


mismatch_target_cor <- lapply(nongene, function(x) {
  mismatch_genes <- intersect(rank_l$Mouse[[x]]$Symbol[1:500], cor_df$Symbol)
  cor_all[gene, mismatch_genes]
})
names(mismatch_target_cor) <- nongene



cor_plot_df <- data.frame(
  Group = c(
    rep(paste0("Top500_", gene), nrow(filter(cor_df, Group))),
    rep(paste0("Out_", gene), nrow(filter(cor_df,!Group))),
    unlist(lapply(nongene, function(x) {
      rep(paste0("Top500_", x), length(mismatch_target_cor[[x]]))
    }))
  ),
  Pcor = c(
    filter(cor_df, Group)$SC_cor,
    filter(cor_df, !Group)$SC_cor,
    unlist(mismatch_target_cor, use.names = FALSE)
  )
)


px3 <- ggplot() +
  geom_density(data = filter(cor_plot_df, Group == paste0("Top500_", gene)),
               aes(x = Pcor),
               fill = "red", colour = NA, alpha = 0.4) +
  geom_density(data = filter(cor_plot_df, Group == paste0("Out_", gene)),
               aes(x = Pcor),
               colour = "black", linewidth = 1.5) +
  xlab("Pearson's cor") +
  ylab("Density") +
  ggtitle(gene) +
  theme_classic() +
  theme(axis.title = element_text(size = 30),
        axis.text = element_text(size = 25),
        plot.title = element_text(size = 30, hjust = 0.5))


for (x in nongene) {
  px3 <- px3 +
    geom_density(data = filter(cor_plot_df, Group == paste0("Top500_", x)),
                 aes(x = Pcor),
                 colour = "grey", linewidth = 1.5)
}



boxplot(cor_df$SC_cor ~ cor_df$Group)
boxplot(abs(cor_df$SC_cor) ~ cor_df$Group)


boxplot(cor_plot_df$Pcor ~ cor_plot_df$Group)
boxplot(abs(cor_plot_df$Pcor) ~ cor_plot_df$Group)





bin_cor_df <- data.frame(SC_cor = cor_all[gene, ]) %>% 
  rownames_to_column(var = "Symbol") %>% 
  left_join(rank_l$Mouse[[gene]], by = "Symbol") %>% 
  filter(Symbol != gene) %>% 
  filter(!is.na(Rank_integrated) & !is.na(SC_cor)) %>% 
  arrange(Rank_integrated) %>% 
  # mutate(Group_integrated = cut(Rank_integrated, breaks = seq(1, length(Rank_integrated), by = 500), include.lowest = TRUE)) %>% 
  mutate(Group_integrated = cut(Rank_integrated, breaks = seq(1, 20000, by = 500), include.lowest = TRUE)) %>% 
  arrange(Rank_perturbation) %>% 
  # mutate(Group_perturbation = cut(Rank_perturbation, breaks = seq(1, length(Rank_perturbation), by = 500), include.lowest = TRUE)) %>% 
  mutate(Group_perturbation = cut(Rank_perturbation, breaks = seq(1, 20000, by = 500), include.lowest = TRUE)) %>% 
  arrange(Rank_binding) %>% 
  # mutate(Group_binding = cut(Rank_binding, breaks = seq(1, length(Rank_binding), by = 500), include.lowest = TRUE)) %>% 
  mutate(Group_binding = cut(Rank_binding, breaks = seq(1, 20000, by = 500), include.lowest = TRUE)) %>% 
  arrange(Rank_integrated) 


ggplot(bin_cor_df) +
  geom_boxplot(aes(x = Group_integrated, y = abs(SC_cor))) +
  # geom_boxplot(aes(x = Group_perturbation, y = abs(SC_cor))) +
  # geom_boxplot(aes(x = Group_binding, y = abs(SC_cor))) +
  ylab("Absolute Pearson's correlation") +
  xlab("Binned rankings") +
  ggtitle(gene) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))
