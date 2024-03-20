## Plotting functions
## -----------------------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(scales)



# Generic plots
# ------------------------------------------------------------------------------


# Return ggobject with default config
# TODO: how to connect to ggplot `+` call

pdef <- function(p) {
  p +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}



# TODO:

qplot <- function(df, xvar, yvar, title = NULL) {
  
  ggplot(df, aes(x = !!sym(xvar), y = !!sym(yvar))) +
    geom_point(shape = 21, size = 2.1) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}



# TODO:

qbox <- function(df, xvar, yvar, title = NULL) {
  
  ggplot(df, aes(x = !!sym(xvar), y = !!sym(yvar))) +
    geom_boxplot(width = 0.3, fill = "slategrey") +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}





# TODO:

plot_hist <- function(df, 
                      stat_col, 
                      nbins = 100,
                      title = NULL, 
                      xlab = NULL) {
  
  if (is.null(xlab)) xlab <- stat_col
  
  ggplot(df, aes(x = !!sym(stat_col))) +
    geom_histogram(bins = nbins, colour = "slategrey", fill = "slategrey") +
    ggtitle(title) +
    xlab(xlab) +
    ylab("Frequency") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}




# For AUC
# ------------------------------------------------------------------------------


# TODO

plot_perf <- function(df,
                      measure = NULL,
                      cols,
                      title,
                      ncol_legend = 1) {
  
  
  stopifnot(measure %in% c("ROC", "PR"), "ID" %in% colnames(df))
  
  if (measure == "ROC") {
    p <- ggplot(df, aes(x = FPR, y = TPR, col = ID))
  } else {
    p <- ggplot(df, aes(x = Recall, y = Precision, col = ID))
  }
  
  p <- p +
    geom_path() +
    ggtitle(title) +
    scale_color_manual(values = cols) +
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



# Akin to a forest plot, show the spread of cell-type correlations for a given
# gene pair across all experiments
# cor_l: a named list of named numeric vectors, where list names are the dataset
# IDs, vector names are the cell types, and the values are cell type correlations

all_corplot <- function(cor_l) {
  
  cor_df <- do.call(
    rbind, 
    lapply(names(cor_l), function(x) data.frame(Cor = cor_l[[x]], ID = x))
  )
  
  ggplot(cor_df, aes(x = Cor, y = reorder(ID, Cor, FUN = median))) +
    geom_point(alpha = 0.4, shape = 21, size = 2.4) +
    geom_boxplot(outlier.shape = NA, coef = 0, fill = "slategrey") +
    geom_vline(xintercept = 0, colour = "black", linetype = "dashed") +
    xlab("Pearson's correlation") +
    ylab("Datasets") +
    theme_classic() +
    theme(axis.title.x = element_text(size = 30),
          axis.text.x = element_text(size = 20),
          axis.title.y = element_text(size = 30),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(c(10, 20, 10, 10)))
  
}



# TODO:

plot_scatter <- function(df, cell_type) {
  
  ggplot(df, aes(x = df[, 1], y = df[, 2])) + 
    geom_point(size = 2.5, shape = 21, fill = "#756bb1", alpha = 0.8) +
    geom_smooth(method = "lm", formula = y ~ x, colour = "black") +
    xlab(colnames(df)[1]) +
    ylab(colnames(df)[2]) +
    ggtitle(cell_type) +
    theme_classic() +
    theme(plot.title = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15))
}



# TODO: Filter min cell to NA?

all_celltype_scatter <- function(mat,
                                 meta,
                                 gene1,
                                 gene2,
                                 min_cell = 20) {
  
  stopifnot(c("Cell_type", "ID") %in% colnames(meta),
            c(gene1, gene2) %in% rownames(mat))
  
  cts <- unique(meta$Cell_type)
  
  plot_l <- lapply(cts, function(ct) {
    ids <- filter(meta, Cell_type == ct)$ID
    ct_mat <- as.matrix(t(mat[c(gene1, gene2), ids]))
    plot_scatter(data.frame(ct_mat), cell_type = ct)
  })
  
  names(plot_l) <- cts
  plot_l <- plot_l[!is.na(plot_l)]
  return(plot_l)
}
