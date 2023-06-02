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
          legend.position = c(0.8, 0.25))
  
  return(p)
}



# For QC plots for count matrices
# ------------------------------------------------------------------------------


qc_hist <- function(meta, xvar, xline, xlab, ylab, log10_xvar = TRUE) {
  
  if (log10_xvar) {
    meta[, xvar] <- log10(meta[, xvar])
    xline <- log10(xline)
  }
  
  ggplot(meta, aes(x = !!sym(xvar))) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = xline, colour = "red") +
    xlab(xlab) +
    ylab(ylab) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.margin = margin(10, 20, 10, 10))
}



all_hist <- function(meta) {
  
  pa <- qc_hist(meta, xvar = "UMI_counts", xline = 500, xlab = "Column (cell)-wise counts", ylab = "log10 UMI counts")
  pb <- qc_hist(meta, xvar = "Gene_counts", xline = 250, xlab = "Row (gene)-wise counts", ylab = "log10 Gene counts")
  pc <- qc_hist(meta, xvar = "RNA_novelty", xline = 0.8, xlab = "Column (cell)-wise counts", ylab = "log10 Gene counts / log10 UMI counts", log10_xvar = FALSE)
  
  p <- cowplot::plot_grid(pa, pb, pc, nrow = 2)
  
  if ("MT_ratio" %in% colnames(meta)) {
    pd <- qc_hist(meta, xvar = "MT_ratio", xline = 0.2, xlab = "Column (cell)-wise counts", ylab = "Ratio of mitochondrial counts", log10_xvar = FALSE)
    p <- p + pd
    p <- cowplot::plot_grid(pa, pb, pc, pd, nrow = 2)
  }
  
  return(p)
}



qc_scatter <- function(meta) {
  
  ggplot(meta, aes(x = log10(Gene_counts), y = log10(UMI_counts))) +
    geom_point() +
    geom_vline(xintercept = log10(250), colour = "red") +
    geom_hline(yintercept = log10(500), colour = "red") +
    xlab("log10 Gene counts") +
    ylab("log10 UMI counts") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20))
  
}



# Correlation plots
# ------------------------------------------------------------------------------








# Plot functions interacting with Seurat object
# ------------------------------------------------------------------------------


# Old, typically not working directly with Seurat object

# plot_scatter <- function(sdat, 
#                          gene1, 
#                          gene2, 
#                          assay = "RNA", 
#                          slot = "data",
#                          jitter = TRUE) {
#   
#   stopifnot(assay %in% Assays(sdat), slot %in% slotNames(sdat@assays[[assay]]))
#   
#   counts <- GetAssayData(object = sdat@assays[[assay]], slot = slot)
#   
#   plot_df <- data.frame(t(as.matrix(counts[c(gene1, gene2), ])))
#   
#   p <- ggplot(plot_df, aes(x = !!sym(gene1), y = !!sym(gene2)))
#   
#   if (jitter) {
#     p <- p + geom_jitter(shape = 21, size = 2.4)
#   } else {
#     p <- p + geom_point(shape = 21, size = 2.4)
#   }
#   
#   p <- p + 
#     xlab(gene1) +
#     ylab(gene2) +
#     theme_classic() +
#     theme(axis.title = element_text(size = 25),
#           axis.text = element_text(size = 20))
#   
#   return(p)
# }
