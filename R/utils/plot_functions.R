## Plotting functions
## -----------------------------------------------------------------------------

library(tidyverse)


# TODO

plot_perf <- function(df, 
                      auc_l, 
                      measure = NULL, 
                      cols, 
                      title,
                      ncol_legend = 1) {
  
  stopifnot(measure %in% c("ROC", "PR"), "Group" %in% colnames(df))
  
  if (measure == "ROC") {
    p <- ggplot(df, aes(x = FPR, y = TPR, col = Group))
  } else {
    p <- ggplot(df, aes(x = Recall, y = Precision, col = Group))
  }
  
  labels <- paste0(names(auc_l), " AUC=", round(unlist(auc_l), 3))
  
  p <- p +
    geom_path() +
    ggtitle(title) +
    scale_color_manual(labels = labels, values = cols) +
    guides(colour = guide_legend(ncol = ncol_legend)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 30),
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.position = c(0.7, 0.2))
  
  return(p)
}


# Plot functions interacting with Seurat object
# ------------------------------------------------------------------------------


# TODO:

plot_scatter <- function(sdat, 
                         gene1, 
                         gene2, 
                         assay = "RNA", 
                         slot = "data",
                         jitter = TRUE) {
  
  stopifnot(assay %in% Assays(sdat), slot %in% slotNames(sdat@assays[[assay]]))
  
  counts <- GetAssayData(object = sdat@assays[[assay]], slot = slot)
  
  plot_df <- data.frame(t(as.matrix(counts[c(gene1, gene2), ])))
  
  p <- ggplot(plot_df, aes(x = !!sym(gene1), y = !!sym(gene2)))
  
  if (jitter) {
    p <- p + geom_jitter(shape = 21, size = 2.4)
  } else {
    p <- p + geom_point(shape = 21, size = 2.4)
  }
  
  p <- p + 
    xlab(gene1) +
    ylab(gene2) +
    theme_classic() +
    theme(axis.title = element_text(size = 25),
          axis.text = element_text(size = 20))
  
  return(p)
}
