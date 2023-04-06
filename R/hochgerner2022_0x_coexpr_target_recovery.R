## Exploring TR-gene correlation and ability to recover annotated targets
## TODO: standardize cmat/cormat naming
## TODO: decide if prep df handles ranking
## TODO: better collection of ranked dfs when finalize which to use
## -----------------------------------------------------------------------------

library(plyr)
library(tidyverse)
library(Seurat)
library(cowplot)
library(parallel)
library(ROCR)
library(pheatmap)
library(RColorBrewer)
library(ggridges)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Load Seurat object
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
sdat <- readRDS(seurat_path)

# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)

# Load correlation matrices generated per cell type and across all cells
cor_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_celltype.RDS")
cor_all <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_cormat_all.RDS")

#
de_top_qntl <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_TR_DEA_quantile.RDS")
de_top_ct <- readRDS("/space/scratch/amorin/R_objects/Hochgerner2022_TR_DEA_celltype.RDS")

# Ribosomal genes as positive control for coexpr.
# Something wonky with table shifts columns over - just want gene symbols
ribo_genes <- read.delim("/space/grp/amorin/Metadata/MGI_GO_term_ribosomal_genes_15-03-2023.txt", row.names = NULL, stringsAsFactors = FALSE)
ribo_genes <- intersect(ribo_genes$MGI.Gene.Marker.ID, rownames(sdat))

tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")


# 
# ------------------------------------------------------------------------------


# TODO:
# TODO: Keep rank?

get_cor_all_df <- function(cmat, tf, rm_tf = TRUE) {
  
  if (rm_tf) {
    cmat <- cmat[setdiff(rownames(cmat), tf), ]
  }
  
  df <- 
    data.frame(Cor_all = cmat[, tf]) %>%
    rownames_to_column(var = "Symbol") %>%
    mutate(
      Cor_all_abs = abs(Cor_all),
      Rank_cor_all = rank(-Cor_all, ties.method = "min"),
      Rank_cor_all_abs = rank(-Cor_all_abs, ties.method = "min")
    )
  
  return(df)
}


# Given a list of cor matrices per cell type and a specified TF, get that TFs
# gene cors in a gene x cell type matrix

tf_by_ct_cmat <- function(cmat_list, tf, rm_tf = TRUE) {

  stopifnot(identical(rownames(cmat_list[[1]]), rownames(cmat_list[[2]])))

  cor_tf <- do.call(cbind, lapply(cmat_list, function(x) x[, tf]))

  if (rm_tf) {
    cor_tf <- cor_tf[rownames(cor_tf) != tf, ]
  }

  return(cor_tf)
}



# TODO

threshold_cmat <- function(mat, top_qtl = 0.95, btm_qtl = NULL) {
  
  mat <- apply(mat, 2, function(x) {
    
    if (!is.null(btm_qtl)) {
      
      top_qtl <- quantile(x, top_qtl, na.rm = TRUE)
      btm_qtl <- quantile(x, btm_qtl, na.rm = TRUE)
      x[x >= top_qtl] <- 1
      x[x <= btm_qtl] <- 1
      x[x < top_qtl & x > btm_qtl] <- 0

    } else {
     
      top_qtl <- quantile(x, top_qtl, na.rm = TRUE)
      x[x >= top_qtl] <- 1
      x[x < top_qtl] <- 0
    }
    
    return(x)
  })

  return(mat)
}


# TODO:

prep_df <- function(vec, name, add_rank = TRUE) {
  
  df <- data.frame(vec)
  colnames(df) <- name
  df <- rownames_to_column(df, var = "Symbol")
  
  if (add_rank) {
    df[[paste0("Rank_", name)]] <- rank(-vec, ties.method = "min")
  }
  
  return(df)
}


# General stats 
# ------------------------------------------------------------------------------


# All genes and cell types
genes <- rownames(sdat)
cts <- unique(sdat$Cell_type)

# Matrix of genes x TFs where elements are the count of NA cors across cell types
n_na <- count_nas(cor_ct)


# Cor across all cells
# ------------------------------------------------------------------------------


tf <- "Mef2c"
rank_df <- rank_l$Mouse[[tf]]

# Dataframe of TF-gene correlation (+/- abs) and their ranking
cor_all_tf <- get_cor_all_df(cor_all, tf, rm_tf = FALSE)

# Inspect TF-gene pairs with appreciable cor but are not in the rankings
top_cor_missing <- filter(cor_all_tf, Symbol %in% setdiff(rownames(sdat), rank_df$Symbol))

# Inspect top ranked genes missing from single cell dataset
top_ranked_missing <- filter(rank_df, Symbol %in% setdiff(rank_df$Symbol, rownames(sdat)))

# Only keep common genes
rank_df <- filter(rank_df, Symbol %in% intersect(rank_df$Symbol, rownames(sdat)))


# Cor by cell type
# ------------------------------------------------------------------------------


# Matrix of TF-gene cor per cell type for the given TF
cor_tf <- tf_by_ct_cmat(cor_ct, tf, rm_tf = TRUE)

# Examine TF-gene NA counts across cell types. NAs result when not enough cells
# passed non-zero expression filter for the given cell type and TF-gene pair.

gene_nas <- n_na[, tf]
ct_nas <- apply(cor_tf, 2, function(x) sum(is.na(x)))
all_gene_nas <- names(which(gene_nas == length(cts)))
all_ct_nas <- names(which(ct_nas == nrow(cor_tf)))

# Remove all NAs
cor_tf <- cor_tf[setdiff(rownames(cor_tf), all_gene_nas), setdiff(colnames(cor_tf), all_ct_nas)]

# Harris 2021 single cell coexpression aggregation (rank sum rank)
cor_ct_rm_tf <- lapply(cor_ct, function(x) x[setdiff(rownames(x), tf), ])
cor_agg_rsr <- aggregate_cor(cor_ct_rm_tf)[rownames(cor_tf), tf]

# Threshold aggregate: Binarize top top 0.5% of cor, sum, rank
cor_agg_thresh <- threshold_cmat(cor_tf, top_qtl = 0.995)[rownames(cor_tf), ]
cor_agg_thresh2 <- rowSums(cor_agg_thresh, na.rm = TRUE)
cor_agg_thresh3 <- rank(-cor_agg_thresh2, ties.method = "min")

# Min rank scheme, where a TF-gene pair is assigned its best rank across all
# cell types. Break ties using average cor across across cell types.
cor_agg_mean <- rowMeans(cor_tf, na.rm = TRUE)
cor_agg_min <- colrank_mat(cor_tf)
cor_agg_min2  <- apply(cor_agg_min, 1, min, na.rm = TRUE)
cor_agg_min3 <- data.table::frank(list(cor_agg_min2, -cor_agg_mean), ties.method = "min")
names(cor_agg_min4) <- names(cor_agg_min3) <- names(cor_agg_min2) <- rownames(cor_tf)


# Collect into a data frame

stopifnot(identical(rownames(cor_tf), names(cor_agg_rsr)))
stopifnot(identical(names(cor_agg_rsr), names(cor_agg_thresh3)))
stopifnot(identical(names(cor_agg_rsr), names(cor_agg_min4)))


cor_ct_df <- data.frame(
  Symbol = rownames(cor_tf),
  Count_NA = gene_nas[rownames(cor_tf)],
  Mean_cor = cor_agg_mean,
  Rank_sum_rank = cor_agg_rsr,
  Thresh_count = cor_agg_thresh2,
  Thresh_rank = cor_agg_thresh3,
  Min_rank = cor_agg_min2,
  Tiebreak_rank = cor_agg_min3
)



# DE genes
# ------------------------------------------------------------------------------

# TODO: to function

de_ct <- de_top_ct[[tf]] %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol %in% rank_df$Symbol) %>% 
  arrange(p_val, -avg_log2FC) %>% 
  dplyr::select(Symbol, p_val, avg_log2FC) %>% 
  dplyr::rename(DE_CT_pval = p_val, DE_CT_log2FC = avg_log2FC) %>% 
  mutate(Rank_DE_CT = 1:nrow(.))
  

de_qntl <- de_top_qntl[[tf]] %>% 
  rownames_to_column(var = "Symbol") %>% 
  filter(Symbol %in% rank_df$Symbol) %>% 
  arrange(p_val, -avg_log2FC) %>% 
  dplyr::select(Symbol, p_val, avg_log2FC) %>% 
  dplyr::rename(DE_qntl_pval = p_val, DE_qntl_log2FC = avg_log2FC) %>% 
  mutate(Rank_DE_qntl = 1:nrow(.))


# Null sets
# ------------------------------------------------------------------------------


# Rank genes by their average expression across all cells

avg_all <- rowMeans(sdat@assays$RNA@data)
avg_all <- avg_all[rank_df$Symbol]
rank_avg_all <- rank(-avg_all, ties.method = "min")


# Sampled targets matched to expression level of observed targets

set.seed(14)

sampled_targets <- sample_expression_level(
  sdat = subset(sdat, features = rank_df$Symbol), 
  targets = filter(rank_df, Curated_target)$Symbol, 
  rank_window = 10)

# Targets from mismatched TRs

mismatched_targets <- lapply(rank_l$Mouse[setdiff(tfs, tf)], function(x) {
  intersect(filter(x, Curated_target)$Symbol, rank_df$Symbol)
})


mismatched_targets[["Matched_expression"]] <- sampled_targets


# Collect data for joint comparison/plotting
# ------------------------------------------------------------------------------


# prep_df(agg_tf1, "Agg1") %>% head
# plyr::join_all(list(de_table, bind_summary), by = "Symbol")
# cor_all_tf
# agg_tf1
# agg_tf2
# agg_tf3
# rank_min_tf
# rank_mean_tf
# rank_tiebreak_tf
# rank_avg_all


evidence <- 
  left_join(rank_df, cor_all_tf, by = "Symbol") %>% 
  left_join(prep_df(agg_tf1, "Agg1"), by = "Symbol") %>% 
  left_join(prep_df(agg_tf2, "Agg2"), by = "Symbol") %>% 
  left_join(prep_df(agg_tf3, "Agg3"), by = "Symbol") %>% 
  left_join(prep_df(rank_min_tf, "Rank_min", add_rank = FALSE), by = "Symbol") %>% 
  left_join(prep_df(rank_mean_tf, "Rank_mean", add_rank = FALSE), by = "Symbol") %>% 
  left_join(prep_df(rank_tiebreak_tf, "Rank_tiebreak", add_rank = FALSE), by = "Symbol") %>% 
  left_join(de_ct, by = "Symbol") %>% 
  left_join(de_qntl, by = "Symbol") %>% 
  left_join(prep_df(rank_avg_all, "Rank_expression", add_rank = FALSE), by = "Symbol")


# Which ranking columns to use

keep_cols <- c(
  "Rank_integrated",
  "Rank_cor_all",
  "Rank_cor_all_abs",
  "Rank_Agg1",
  "Rank_Agg2",
  "Rank_Agg3",
  "Rank_min",
  "Rank_mean",
  "Rank_tiebreak",
  "Rank_DE_CT",
  "Rank_DE_qntl",
  "Rank_expression"
)


stopifnot(all(keep_cols %in% colnames(evidence)))


# List of ROC dfs
roc_l <- lapply(keep_cols, function(x) {
  
  get_perf_df(
    rank_df = dplyr::arrange(evidence, !!sym(x)),
    label_col = "Curated_target",
    measure = "ROC") %>% 
    mutate(Group = x)
})
names(roc_l) <- keep_cols


# List of AUROCs
auc_l <- lapply(keep_cols, function(x) {
  
  get_au_perf(
    rank_df = dplyr::arrange(evidence, !!sym(x)),
    label_col = "Curated_target",
    measure = "AUROC")
})
names(auc_l) <- keep_cols



# List of PR dfs
pr_l <- lapply(keep_cols, function(x) {
  
  get_perf_df(
    rank_df = dplyr::arrange(evidence, !!sym(x)),
    label_col = "Curated_target",
    measure = "PR") %>% 
    mutate(Group = x)
})
names(pr_l) <- keep_cols


# List of AUPRCs
auprc_l <- lapply(keep_cols, function(x) {
  
  get_au_perf(
    rank_df = dplyr::arrange(evidence, !!sym(x)),
    label_col = "Curated_target",
    measure = "AUPRC")
})
names(auprc_l) <- keep_cols




# Prepare for plotting
# ------------------------------------------------------------------------------


# Long df for plotting

roc_df <- do.call(rbind, roc_l) %>% 
  mutate(Group = factor(Group, levels = keep_cols))

# Labels for plotting

auc_labels <- paste0(keep_cols, " AUC=", round(unlist(auc_l), 3))

# Colours

cols <- c("Rank_integrated" = "black",
          "Rank_cor_all" = "dodgerblue1",
          "Rank_cor_all_abs" = "dodgerblue2",
          "Rank_Agg1" = "orangered1",
          "Rank_Agg2" = "orangered2",
          "Rank_Agg3" = "orangered3",
          "Rank_min" = "seagreen2",
          "Rank_mean" = "seagreen3",
          "Rank_tiebreak" = "seagreen4",
          "Rank_DE_qntl" = "purple",
          "Rank_DE_CT" = "green",
          "Rank_expression" = "lightgrey")



plot_auc(roc_df, auc_labels, cols, tf)



# Example of comparing observed ranking to null targets
# Here using integrated rank

roc_null_l <- lapply(names(mismatched_targets), function(x) {
  
  null_df <- evidence %>% 
    mutate(Curated_target = Symbol %in% mismatched_targets[[x]]) %>% 
    arrange(Rank_integrated)
  
  get_perf_df(
    rank_df = null_df,
    label_col = "Curated_target",
    measure = "ROC") %>% 
    mutate(Group = x)
})
names(roc_null_l) <- names(mismatched_targets)

# Add observed integrated rank
roc_null_l[[tf]] <- mutate(roc_l$Rank_integrated, Group = tf)

# Plot df

roc_null_df <- do.call(rbind, roc_null_l) %>%
  mutate(Group = factor(Group, levels = names(roc_null_l)))


# Null AUCs 
auc_null_l <- lapply(names(mismatched_targets), function(x) {
  
  null_df <- evidence %>% 
    mutate(Curated_target = Symbol %in% mismatched_targets[[x]]) %>% 
    arrange(Rank_integrated)
  
  get_au_perf(
    rank_df = null_df,
    label_col = "Curated_target",
    measure = "AUROC")
})
names(auc_null_l) <- names(mismatched_targets)


auc_null_l[[tf]] <- auc_l$Rank_integrated


auc_null_labels <- paste0(names(auc_null_l), " AUC=", round(unlist(auc_null_l), 3))


cols_null <- c(rep("lightgrey", 7), "red", "black")
names(cols_null) <- names(auc_null_l)

plot_auc(roc_null_df, auc_null_labels, cols_null, tf)




# across all TFs
# ------------------------------------------------------------------------------




all_auprc_l <- mclapply(tfs, function(tf) {
  
  # rank_df <- rank_l$Mouse[[tf]] %>% 
  #   filter(Symbol %in% intersect(Symbol, rownames(sdat))) %>% 
  #   filter(Symbol != tf)
  
  
  rank_df <- rank_l$Mouse[[tf]] %>% 
    filter(Symbol %in% intersect(Symbol, rownames(sdat)))
  
  
  cor_all_tf <- get_cor_all_df(cor_all, tf, rm_tf = FALSE)
  
  cor_tf <- tf_by_ct_cmat(cor_ct, tf, rm_tf = FALSE)
  all_ct_nas <- apply(cor_tf, 2, function(x) all(is.na(x)))
  all_ct_nas <- names(all_ct_nas[all_ct_nas])
  stopifnot(all(sdat@assays$RNA@data[tf, sdat$Cell_type %in% all_ct_nas] == 0))
  cor_tf <- cor_tf[, setdiff(colnames(cor_tf), all_ct_nas)]
  
  
  thresh_tf1 <- threshold_cormat(cor_tf, top_qtl = 0.995, btm_qtl = 0.005)
  agg_tf1 <- rowSums(thresh_tf1, na.rm = TRUE)
  thresh_tf2 <- threshold_cormat(cor_tf, top_qtl = 0.995)
  agg_tf2 <- rowSums(thresh_tf2, na.rm = TRUE)
  thresh_tf3 <- threshold_cormat(abs(cor_tf), top_qtl = 0.995)
  agg_tf3 <- rowSums(thresh_tf3, na.rm = TRUE)
  
  
  rank_tf <- rank_cormat(cor_tf)
  rank_min_tf <- apply(rank_tf, 1, min)
  rank_mean_tf <- rowMeans(rank_tf)
  rank_tiebreak_tf <- data.table::frank(list(rank_min_tf, rank_mean_tf), ties.method = "min")
  names(rank_min_tf) <- names(rank_mean_tf) <- names(rank_tiebreak_tf) <- rownames(rank_tf)
  
  
  de_ct <- de_top_ct[[tf]] %>% 
    rownames_to_column(var = "Symbol") %>% 
    filter(Symbol %in% rank_df$Symbol) %>% 
    arrange(p_val, -avg_log2FC) %>% 
    dplyr::select(Symbol, p_val, avg_log2FC) %>% 
    dplyr::rename(DE_CT_pval = p_val, DE_CT_log2FC = avg_log2FC) %>% 
    mutate(Rank_DE_CT = 1:nrow(.))
  
  
  de_qntl <- de_top_qntl[[tf]] %>% 
    rownames_to_column(var = "Symbol") %>% 
    filter(Symbol %in% rank_df$Symbol) %>% 
    arrange(p_val, -avg_log2FC) %>% 
    dplyr::select(Symbol, p_val, avg_log2FC) %>% 
    dplyr::rename(DE_qntl_pval = p_val, DE_qntl_log2FC = avg_log2FC) %>% 
    mutate(Rank_DE_qntl = 1:nrow(.))
  
  
  evidence <- 
    left_join(rank_df, cor_all_tf, by = "Symbol") %>% 
    left_join(prep_df(agg_tf1, "Agg1"), by = "Symbol") %>% 
    left_join(prep_df(agg_tf2, "Agg2"), by = "Symbol") %>% 
    left_join(prep_df(agg_tf3, "Agg3"), by = "Symbol") %>% 
    left_join(prep_df(rank_min_tf, "Rank_min", add_rank = FALSE), by = "Symbol") %>% 
    left_join(prep_df(rank_mean_tf, "Rank_mean", add_rank = FALSE), by = "Symbol") %>% 
    left_join(prep_df(rank_tiebreak_tf, "Rank_tiebreak", add_rank = FALSE), by = "Symbol") %>% 
    left_join(de_ct, by = "Symbol") %>% 
    left_join(de_qntl, by = "Symbol") %>% 
    left_join(prep_df(rank_avg_all, "Rank_expression", add_rank = FALSE), by = "Symbol")
  
  
  auprc_l <- lapply(keep_cols, function(x) {
    
    get_au_perf(
      rank_df = dplyr::arrange(evidence, !!sym(x)),
      label_col = "Curated_target",
      measure = "AUPRC")
  })
  names(auprc_l) <- keep_cols
  
  return(unlist(auprc_l))
  
}, mc.cores = 8)



all_auprc_df <- round(do.call(rbind, all_auprc_l), 3)
rownames(all_auprc_df) <- tfs


pheatmap(all_auprc_df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(1, 3, 6, 9, 11),
         color = colorRampPalette(brewer.pal(n = 7, name = "OrRd"))(100),
         display_numbers = all_auprc_df,
         number_color = "black",
         fontsize = 20,
         cellwidth = 50,
         cellheight = 50)




all_auroc_l <- mclapply(tfs, function(tf) {
  
  rank_df <- rank_l$Mouse[[tf]] %>%
    filter(Symbol %in% intersect(Symbol, rownames(sdat))) %>%
    filter(Symbol != tf)
  
  # rank_df <- rank_l$Mouse[[tf]] %>% 
  #   filter(Symbol %in% intersect(Symbol, rownames(sdat)))
  
  
  cor_all_tf <- get_cor_all_df(cor_all, tf, rm_tf = FALSE)
  
  cor_tf <- tf_by_ct_cmat(cor_ct, tf, rm_tf = FALSE)
  all_ct_nas <- apply(cor_tf, 2, function(x) all(is.na(x)))
  all_ct_nas <- names(all_ct_nas[all_ct_nas])
  stopifnot(all(sdat@assays$RNA@data[tf, sdat$Cell_type %in% all_ct_nas] == 0))
  cor_tf <- cor_tf[, setdiff(colnames(cor_tf), all_ct_nas)]
  
  
  thresh_tf1 <- threshold_cormat(cor_tf, top_qtl = 0.995, btm_qtl = 0.005)
  agg_tf1 <- rowSums(thresh_tf1, na.rm = TRUE)
  thresh_tf2 <- threshold_cormat(cor_tf, top_qtl = 0.995)
  agg_tf2 <- rowSums(thresh_tf2, na.rm = TRUE)
  thresh_tf3 <- threshold_cormat(abs(cor_tf), top_qtl = 0.995)
  agg_tf3 <- rowSums(thresh_tf3, na.rm = TRUE)
  
  
  rank_tf <- rank_cormat(cor_tf)
  rank_min_tf <- apply(rank_tf, 1, min)
  rank_mean_tf <- rowMeans(rank_tf)
  rank_tiebreak_tf <- data.table::frank(list(rank_min_tf, rank_mean_tf), ties.method = "min")
  names(rank_min_tf) <- names(rank_mean_tf) <- names(rank_tiebreak_tf) <- rownames(rank_tf)
  
  
  de_ct <- de_top_ct[[tf]] %>% 
    rownames_to_column(var = "Symbol") %>% 
    filter(Symbol %in% rank_df$Symbol) %>% 
    arrange(p_val, -avg_log2FC) %>% 
    dplyr::select(Symbol, p_val, avg_log2FC) %>% 
    dplyr::rename(DE_CT_pval = p_val, DE_CT_log2FC = avg_log2FC) %>% 
    mutate(Rank_DE_CT = 1:nrow(.))
  
  
  de_qntl <- de_top_qntl[[tf]] %>% 
    rownames_to_column(var = "Symbol") %>% 
    filter(Symbol %in% rank_df$Symbol) %>% 
    arrange(p_val, -avg_log2FC) %>% 
    dplyr::select(Symbol, p_val, avg_log2FC) %>% 
    dplyr::rename(DE_qntl_pval = p_val, DE_qntl_log2FC = avg_log2FC) %>% 
    mutate(Rank_DE_qntl = 1:nrow(.))
  
  
  evidence <- 
    left_join(rank_df, cor_all_tf, by = "Symbol") %>% 
    left_join(prep_df(agg_tf1, "Agg1"), by = "Symbol") %>% 
    left_join(prep_df(agg_tf2, "Agg2"), by = "Symbol") %>% 
    left_join(prep_df(agg_tf3, "Agg3"), by = "Symbol") %>% 
    left_join(prep_df(rank_min_tf, "Rank_min", add_rank = FALSE), by = "Symbol") %>% 
    left_join(prep_df(rank_mean_tf, "Rank_mean", add_rank = FALSE), by = "Symbol") %>% 
    left_join(prep_df(rank_tiebreak_tf, "Rank_tiebreak", add_rank = FALSE), by = "Symbol") %>% 
    left_join(de_ct, by = "Symbol") %>% 
    left_join(de_qntl, by = "Symbol") %>% 
    left_join(prep_df(rank_avg_all, "Rank_expression", add_rank = FALSE), by = "Symbol")
  
  
  auroc_l <- lapply(keep_cols, function(x) {
    
    get_au_perf(
      rank_df = dplyr::arrange(evidence, !!sym(x)),
      label_col = "Curated_target",
      measure = "AUROC")
  })
  names(auroc_l) <- keep_cols
  
  return(unlist(auroc_l))
  
}, mc.cores = 8)



all_auroc_df <- round(do.call(rbind, all_auroc_l), 3)
rownames(all_auroc_df) <- tfs


pheatmap(all_auroc_df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         gaps_col = c(1, 3, 6, 9, 11),
         color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100),
         display_numbers = all_auroc_df,
         number_color = "black",
         fontsize = 20,
         cellwidth = 50,
         cellheight = 50)




tt <- lapply(tfs, function(x) {
  
  cor_tf <- tf_by_ct_cmat(cor_ct, x, rm_tf = FALSE)
  
  thresh_tf <- threshold_cormat(cor_tf)
  
  agg_tf <- rowSums(thresh_tf, na.rm = TRUE)
  
})

tt <- do.call(cbind, tt)
colnames(tt) <- tfs


for (tf in tfs) {
  tt[tf, tf] <- 0
}

all_0 <- which(rowSums(tt) == 0)

tt <- tt[-all_0, ]


tt_z <- t(scale(t(tt)))


view(data.frame(tt))
view(data.frame(tt_z))



# Misc
# ------------------------------------------------------------------------------



# Example of high cor that is suspicious
FeaturePlot(sdat, features = c("Runx1", "Mnd1"))
plot_scatter(subset(sdat, subset = Cell_type == "GABA-11-Adora2a-Id4"), gene1 = "Runx1", gene2 = "Mnd1")
plot_scatter(subset(sdat, subset = Cell_type == "microglia"), gene1 = "Ascl1", gene2 = "Gm1992")




# Inspect ribosomal gene cor
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



# Ridge plot
# ------------------------------------------------------------------------------


evidence <- left_join(rank_df, cor_all_tf, by = "Symbol") %>% 
  filter(Symbol != tf) %>% 
  mutate(
    Group = cut(Rank_integrated, breaks = seq(1, nrow(.), by = 100), include.lowest = TRUE),
    Group = ifelse(Rank_integrated > 500, "Bottom", Group),
    Group = factor(Group, levels = rev(unique(Group))))


ggplot(evidence, aes(y = Group, x = Cor_all)) +
  geom_density_ridges() +
  geom_vline(xintercept = 0) +
  theme_classic()


ggplot(evidence, aes(y = Group, x = Cor_all)) +
  geom_boxplot() +
  theme_classic()


ggplot(evidence, aes(y = Group, x = Cor_all_abs)) +
  geom_density_ridges() +
  geom_vline(xintercept = 0) +
  theme_classic()


ggplot(evidence, aes(y = Group, x = Cor_all_abs)) +
  geom_boxplot() +
  theme_classic()



ggplot(evidence, aes(y = Curated_target, x = Cor_all)) +
  geom_density_ridges() +
  geom_vline(xintercept = 0) +
  theme_classic()
