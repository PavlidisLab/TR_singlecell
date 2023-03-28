##
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(pheatmap)
source("R/00_config.R")

#
expr_l <- readRDS(expr_mat_l_path)
cor_l <- readRDS(cor_mat_l_path)


# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)


#
# ------------------------------------------------------------------------------


# Get a data frame of gene cor values for requested gene across all data sets 

cor_df <- function(cor_list, gene) {
  
  df_list <- lapply(cor_list, function(x) {
    
    if (!(gene %in% rownames(x))) {
      return(data.frame(Cor = NA, Symbol = rownames(x)))
    }
    
    vec <- x[gene, ]
    vec <- vec[names(vec) != gene]
    data.frame(Cor = vec) %>% rownames_to_column(var = "Symbol")
  })
  
  df <- plyr::join_all(df_list, by = "Symbol")
  colnames(df)[2:ncol(df)] <- names(cor_list)
  
  return(df)
}



# Convert df of TF-gene cors to ranks

rank_cor <- function(cor_df) {
  
  rank_l <- lapply(2:ncol(cor_df), function(x) rank(-cor_df[, x]))  # TODO this doesn't need to be lapply, no?
  mat <- do.call(cbind, rank_l)
  rownames(mat) <- cor_df$Symbol
  colnames(mat) <- colnames(cor_df)[2:ncol(cor_df)]
  mat <- cbind(mat, Average_rank = rowMeans(mat, na.rm = TRUE))
  mat <- mat[order(mat[, "Average_rank"]), ]
  
  return(mat)
}



#
# ------------------------------------------------------------------------------



tf <- "PAX6"
tf_cor <- cor_df(cor_l, tf)
tf_cor_rank <- rank_cor(tf_cor)


# Get the "cor of cors" - how consistent TR-gene coexpr is across data sets

all_cor <- cor(dplyr::select_if(tf_cor, is.numeric), 
               use = "pairwise.complete.obs")



#
# ------------------------------------------------------------------------------



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


# TODO:

get_au_perf <- function(rank_df, label_col, measure = NULL) {
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("AUROC", "AUPRC"))
  
  # positives/has evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(rank_df[[label_col]]))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  
  if (measure == "AUROC") {
    perf <- ROCR::performance(pred, measure = "auc")@y.values[[1]]
  } else {
    perf <- ROCR::performance(pred, measure = "aucpr")@y.values[[1]]
  }
  
  return(perf)
}


# TODO
# TODO: labels should just be generated inside

plot_auc <- function(roc_df, auc_labels, cols, tf) {
  
  p <- 
    # ggplot(roc_df, aes(x = Recall, y = Precision, col = Group)) +
    ggplot(roc_df, aes(x = FPR, y = TPR, col = Group)) +
    geom_path() +
    ggtitle(tf) +
    scale_color_manual(labels = auc_labels, values = cols) +
    guides(colour = guide_legend(ncol = 2)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 30),
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.position = c(0.7, 0.2))
  
  return(p)
}


rank_df <- left_join(
  rank_l$Human[[tf]][, c("Symbol", "Rank_integrated", "Curated_target")],
  rownames_to_column(data.frame(tf_cor_rank), var = "Symbol"),
  by = "Symbol")
  

keep_cols <- c("Rank_integrated", "Average_rank", names(expr_l))


# List of ROC dfs
roc_l <- lapply(keep_cols, function(x) {
  
  get_perf_df(
    rank_df = dplyr::arrange(rank_df, !!sym(x)),
    label_col = "Curated_target",
    measure = "ROC") %>%
    # measure = "PR") %>% 
    mutate(Group = x)
})
names(roc_l) <- keep_cols


# List of AUROCs
auc_l <- lapply(keep_cols, function(x) {
  
  get_au_perf(
    rank_df = dplyr::arrange(rank_df, !!sym(x)),
    label_col = "Curated_target",
    measure = "AUROC")
    # measure = "AUPRC")
})
names(auc_l) <- keep_cols


roc_df <- do.call(rbind, roc_l) %>% 
  mutate(Group = factor(Group, levels = keep_cols))

auc_labels <- paste0(names(auc_l), " AUC=", round(unlist(auc_l), 3))


cols <- c(rep("lightgrey", length(auc_labels)))
names(cols) <- names(auc_l)
cols["Rank_integrated"] <- "black"
cols["Average_rank"] <- "red"


plot_auc(roc_df, auc_labels, cols, tf)




# Plotting
# ------------------------------------------------------------------------------


# TODO: 

gene <- "CSF1"

p_list <- lapply(names(expr_l), function(x) {
  
  mat <- expr_l[[x]]
  
  if(!all(c(tf, gene) %in% rownames(mat))) {
    return(NA)
  }
  
  df <- data.frame(tf = mat[tf,], gene = mat[gene, ])
  ggplot(df, aes(x = tf, y = gene)) + 
    geom_point(size = 2, shape = 21, fill = "#756bb1") +
    geom_smooth(method = "lm", colour = "black") +
    xlab(tf) +
    ylab(gene) +
    ggtitle(paste0(x, ": ", round(cor_l[[x]][tf, gene], 3))) +
    theme_classic() +
    theme(plot.title = element_text(size = 20))
})


p_list <- p_list[!is.na(p_list)]

cowplot::plot_grid(plotlist = p_list)

head(sort(expr_l$FANTOM[gene, ], decreasing = TRUE), 10)
head(sort(expr_l$FANTOM[tf, ], decreasing = TRUE), 10)


# Heatmap of cor of cor. Pre-cluster to get order for triangular mat

hc <-  hclust(as.dist(1 - all_cor))
all_cor <- all_cor[hc$order, hc$order]
all_cor[lower.tri(all_cor)] <-  NA
diag(all_cor) <- NA
all_cor <- all_cor[1:nrow(all_cor)-1, 2:ncol(all_cor)]


pal_length <- 100
bluered_pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(pal_length)
color_breaks <- seq(min(all_cor, na.rm = TRUE), max(all_cor, na.rm = TRUE), length.out = pal_length)


pheatmap(all_cor, 
         cluster_col = FALSE, 
         cluster_row = FALSE,
         color = bluered_pal,
         breaks = color_breaks,
         na_col = "white",
         border_color = NA,
         fontsize = 18)
