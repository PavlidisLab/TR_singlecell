##
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(pheatmap)
source("R/00_config.R")

expr_l <- readRDS(expr_mat_l_path)
cor_l <- readRDS(cor_mat_l_path)



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






gene <- "CRB1"



# Plotting
# ------------------------------------------------------------------------------


# TODO: 

p_list <- lapply(names(mat_list), function(x) {
  
  mat <- mat_list[[x]]
  
  if(!all(c(tf, gene) %in% rownames(mat))) {
    return(NA)
  }
  
  df <- data.frame(tf = mat[tf,], gene = mat[gene, ])
  ggplot(df, aes(x = tf, y = gene)) + 
    geom_point(size = 2, shape = 21, fill = "#756bb1") +
    geom_smooth(method = "lm", colour = "black") +
    xlab(tf) +
    ylab(gene) +
    ggtitle(paste0(x, ": ", round(cor_list[[x]][tf, gene], 3))) +
    theme_classic() +
    theme(plot.title = element_text(size = 20))
})


p_list <- p_list[!is.na(p_list)]

cowplot::plot_grid(plotlist = p_list)


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
