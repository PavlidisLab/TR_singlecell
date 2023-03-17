library(tidyverse)
library(Seurat)
library(cowplot)
library(parallel)
library(future)
library(WGCNA)
library(ggrepel)
library(pheatmap)
source("R/utils/functions.R")
source("R/00_config.R")


seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
de_top_qntl <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_amygdala_biorvx_TR_DEA_quantile.RDS")
de_top_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_amygdala_biorvx_TR_DEA_celltype.RDS")
rank_path <- "/space/scratch/amorin/R_objects/Apr2022_ranked_target_list.RDS"

# Visualize batch, samples, and provided labels on PC space
DimPlot(sdat, reduction = "pca", group.by = "Batch")
DimPlot(sdat, reduction = "pca", group.by = "Sample")
DimPlot(sdat, reduction = "pca", group.by = "FC_time")
DimPlot(sdat, reduction = "pca", group.by = "Cell_type") + NoLegend()


# Visualize clusters, batch, samples, and provided labels on UMAP space
DimPlot(sdat, reduction = "umap", group.by = "seurat_clusters")
DimPlot(sdat, reduction = "umap", group.by = "Batch")
DimPlot(sdat, reduction = "umap", group.by = "Sample")
DimPlot(sdat, reduction = "umap", group.by = "Cell_type", label = TRUE, label.size = 6, pt.size = 0.5) + NoLegend()


# Plot expression of genes on UMAP
plot_gene <- c("Pax6", "Neurod1")
FeaturePlot(sdat, features = plot_gene)



# Get the average gene expression per cell type
# ------------------------------------------------------------------------------


ct_avg <- get_ct_avg(sdat)
ct_avg_scale <- get_ct_avg(sdat, scale = TRUE)

head(sort(ct_avg[, "microglia"], decreasing = TRUE), 15)
head(sort(ct_avg_scale[, "microglia"], decreasing = TRUE), 15)




# Look at expression of top TR-targets in cells expressing the TR
# ------------------------------------------------------------------------------


gene <- "Pax6"
top_targets <- rank_l$Mouse[[gene]]$Symbol[1:500]
top_targets <- unique(c(gene, top_targets[top_targets %in% genes]))

# Create metadata column of status in top expressed celltype/quantile
sdat <- top_expr_quantile(sdat, gene = gene)
sdat <- top_expr_celltype(sdat, avg_mat = ct_avg, gene = gene)
# sdat <- top_expr_celltype(sdat, avg_mat = ct_avg_scale, gene = gene)

# Plot expression of genes on UMAP
pxa <- FeaturePlot(sdat, features = gene)
pxb <- FeaturePlot(sdat, features = "Top_expr_quantile", cols = c("lightgrey", "red")) + ggtitle("Top 10% expressed cells") + NoLegend()
pxc <- FeaturePlot(sdat, features = "Top_expr_celltype", cols = c("lightgrey", "red")) + ggtitle("Top expressed cell type") + NoLegend()
px <- plot_grid(pxa, pxb, pxc, nrow = 1)



# DE status by group status
# ------------------------------------------------------------------------------


# wilcox.test(sdat@assays$RNA@data["Ctsz", ] ~ sdat$Top_expr_quantile)$p.value
# fit <- glm(sdat@assays$RNA@data["Ctsz", ] ~ sdat$Top_expr_quantile, family = "poisson")
# summary(fit)

# Default much faster but only uses subset of genes
# de_top_qntl <- FindMarkers(sdat, group.by = "Top_expr_quantile", ident.1 = "TRUE")
# de_top_ct <- FindMarkers(sdat, group.by = "Top_expr_celltype", ident.1 = "TRUE")



# Using all genes (so very slow)
# de_top_qntl <- FindMarkers(sdat,
#                            group.by = "Top_expr_quantile",
#                            ident.1 = "TRUE",
#                            min.cells.group = 1,
#                            min.cells.feature = 1,
#                            min.pct = 0,
#                            logfc.threshold = 0,
#                            only.pos = FALSE)



# de_top_qntl <- FindMarkers(sdat,
#                            group.by = "Top_expr_celltype",
#                            ident.1 = "TRUE",
#                            min.cells.group = 1,
#                            min.cells.feature = 1,
#                            min.pct = 0,
#                            logfc.threshold = 0,
#                            only.pos = FALSE)


de_top_qntl <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_amygdala_biorvx_TR_DEA_quantile.RDS")
de_top_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_amygdala_biorvx_TR_DEA_celltype.RDS")


bind_mat <- function(list, col) {
  
  stopifnot(col %in% colnames(list[[1]]))
  list <- lapply(list, function(x) x[rownames(list[[1]]), ])  # common order
  
  mat <- do.call(cbind, lapply(list, `[[`, col))
  rownames(mat) <- rownames(list[[1]])
  colnames(mat) <- names(list)
  
  return(mat)
}


score_pval_mat <- function(pval_mat) {
  
  pval_mat <- apply(pval_mat, 2, function(x) {
    nonzero_min <- min(x[x != 0])
    ifelse(x == 0, nonzero_min, x)
  })
  pval_mat <- -log10(pval_mat)
  
  return(pval_mat)  
}


group_evidence <- function(df, topn = 500) {
  
  df <- mutate(df, Group = case_when(
    Curated_target & Rank_integrated <= topn ~ "Top",
    Rank_integrated <= 500 ~ "Genomic",
    Curated_target ~ "Low-throughput",
    TRUE ~ "Out"
  ))
  
  df$Group <- 
    factor(df$Group, levels = c("Out", "Low-throughput", "Genomic", "Top"))
  
  return(df)
}


fc_mat_qntl <- bind_mat(de_top_qntl, "avg_log2FC")
fc_mat_ct <- bind_mat(de_top_ct, "avg_log2FC")

score_mat_qntl <- score_pval_mat(bind_mat(de_top_qntl, "p_val"))
score_mat_ct <- score_pval_mat(bind_mat(de_top_ct, "p_val"))


gene <- "Pax6"


fc_df <- data.frame(SC_FC = fc_mat_qntl[, gene],
                    Pval_score = score_mat_qntl[, gene]) %>%
  rownames_to_column(var = "Symbol") %>%
  left_join(rank_l$Mouse[[gene]], by = "Symbol") %>%
  # filter(Symbol != gene) %>%
  filter(!is.na(Rank_integrated) & !is.na(SC_FC)) %>%
  arrange(Rank_integrated) %>% 
  group_evidence()




plot_volcano <- function(df) {
  
  ggplot() +
    geom_jitter(data = df[df$Group == "Out", ],
                aes(x = SC_FC, y = Pval_score, fill = Group),
                alpha = 1, shape = 19, colour = "grey", size = 2, height = 1, width = 0.01) +
    geom_jitter(data = df[df$Group == "Low-throughput", ],
                aes(x = SC_FC, y = Pval_score, fill = Group),
                alpha = 1, shape = 21, size = 3, col = "black", height = 1, width = 0.01) +
    geom_jitter(data = df[df$Group == "Genomic", ],
                aes(x = SC_FC, y = Pval_score, fill = Group),
                alpha = 1, shape = 21, size = 3, col = "black", height = 1, width = 0.01) +
    geom_jitter(data = df[df$Group == "Top", ],
                aes(x = SC_FC, y = Pval_score, fill = Group),
                alpha = 1, shape = 21, size = 3, col = "black", height = 1, width = 0.01) +
    xlab("Single cell average log2FC") +
    ylab("Single cell -log10(p-value)") +
    scale_fill_manual(name = "Evidence",
                      values = c("Out" = "grey",
                                 "Low-throughput" = "royalblue",
                                 "Genomic" = "firebrick",
                                 "Top" = "goldenrod"),
                      labels = c("None",
                                 "Curated",
                                 "Top 500 integrated",
                                 "Curated and top 500 integrated"),
                      breaks = levels(df$Group)) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 30),
      plot.title = element_text(size = 30),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      legend.position = c(0.85, 0.2),
      plot.margin = margin(10, 10, 10, 10)
    )
}


plot_volcano(fc_df)


cor(fc_df$SC_FC, 
    fc_df[, c("Rank_integrated", "Rank_perturbation", "Rank_binding", "Count_DE")],
    method = "spearman", 
    use = "pairwise.complete.obs")


cor(abs(fc_df$SC_FC), 
    fc_df[, c("Rank_integrated", "Rank_perturbation", "Rank_binding", "Count_DE")],
    method = "spearman", 
    use = "pairwise.complete.obs")



plot(density(abs(fc_df$SC_FC[1:500])), col = "red")
lines(density(abs(fc_df$SC_FC[10000:10500])), col = "black")

boxplot(fc_df$SC_FC ~ fc_df$Group)
boxplot(abs(fc_df$SC_FC) ~ fc_df$Group)



plot(density(fc_mat_qntl[, "Ascl1"]))
fc_mat_qntl[c("Dll1", "Dll3", "Dll4"), "Ascl1"]
head(sort(fc_mat_qntl[, "Ascl1"], decreasing = TRUE), 20)


# Cor of TR-targets
# ------------------------------------------------------------------------------



# cor_all <- WGCNA::cor(x = t(mat),
#                       nThreads = 8,
#                       use = "pairwise.complete.obs")
# 
# 
# cor_qntl <- WGCNA::cor(x = t(mat[, sdat$Top_expr_quantile, drop = FALSE]),
#                        nThreads = 8,
#                        use = "pairwise.complete.obs")
# 
# 
# cor_gene_all <- cor_all[gene, ]
# cor_gene_qntl <- cor_qntl[gene, ]



cor_all <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_amygdala_biorvx_cormat.RDS")
ribo_genes <- read.delim("/space/grp/amorin/Metadata/MGI_GO_term_ribosomal_genes_15-03-2023.txt", row.names = NULL, stringsAsFactors = FALSE)

# something wonky with table shifts columns over
ribo_genes <- intersect(ribo_genes$MGI.Gene.Marker.ID, rownames(sdat))


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




# Histogram of counts for TF and another gene (like a top target)
# ------------------------------------------------------------------------------


expr_hist <- function(sdat, gene1, gene2) {
  
  plot_df <- data.frame(gene1 = sdat@assays$RNA@data[gene1, ],
                        gene2 = sdat@assays$RNA@data[gene2, ])
  
  ggplot() +
    geom_histogram(aes(x = gene1), data = plot_df,
                   bins = 30, colour = "black", fill = "red", alpha = 0.5) +
    geom_histogram(aes(x = gene2), data = plot_df,
                   bins = 30, colour = "black", fill = "black", alpha = 0.5) +
    xlab("Expression counts") +
    ylab("Cell counts") +
    ggtitle(paste(gene1, gene2)) +
    theme_classic() +
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          plot.title = element_text(size = 25))
}


gene1 <- "Ascl1"
gene2 <- "Dll1"
sdat <- top_expr_quantile(sdat, gene = gene1)
sdat <- top_expr_celltype(sdat, avg_mat = ct_avg, gene = gene1)

# All cells
expr_hist(sdat, gene1 = gene1, gene2 = gene2)

# Top Quantile
expr_hist(subset(sdat, subset = Top_expr_quantile), gene1 = gene1, gene2 = gene2)

# Top cell type
expr_hist(subset(sdat, subset = Top_expr_celltype), gene1 = gene1, gene2 = gene2)



# Boxplot of gene expression across cell types
# ------------------------------------------------------------------------------


ct_expr_boxplot <- function(sdat, gene) {
  
  plot_df <- mutate(sdat@meta.data, Expr = sdat@assays$RNA@data[gene,])
  
  ggplot(plot_df, aes(x = reorder(Cell_type, Expr, median), y = Expr)) +
    geom_boxplot() +
    xlab("Cell type") +
    ylab("Expression counts") +
    ggtitle(gene) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
          axis.text.y = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25))
}



# Looking at expression of paper TRs


tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

px6 <- lapply(tfs, function(x) ct_expr_boxplot(sdat, x))
names(px6) <- tfs


# Select gene examples
head(VariableFeatures(sdat), 10)
ct_expr_boxplot(sdat, "Cst3") # Cst3  top microglia expr
ct_expr_boxplot(sdat, "Lyl1") # Lyl1 RUNX1 target
ct_expr_boxplot(sdat, "Csf1r") # Csf1r microglia marker
ct_expr_boxplot(sdat, "C1qc") # C1qc top Runx1 cor
ct_expr_boxplot(sdat, "Acta2") # Acta2 top var

