## Project functions
## TODO: sample_expression_level() impl requires subset (costly) - try pre-format input mat
## -----------------------------------------------------------------------------


library(parallel)


# TODO:

top_expr_quantile <- function(sdat, gene, qtl = 0.9) {
  
  expr_cutoff <- quantile(sdat@assays$RNA@data[gene, ], qtl)
  expr_cutoff <- ifelse(expr_cutoff == 0L, 1, expr_cutoff)
  sdat$Top_expr_quantile <- sdat@assays$RNA@data[gene, ] >= expr_cutoff
  
  return(sdat)
}



# TODO:

top_expr_celltype <- function(sdat, avg_mat, gene) {
  
  stopifnot("Cell_type" %in% colnames(sdat@meta.data))
  
  top_ct <- names(which.max(avg_mat[gene , ]))
  sdat$Top_expr_celltype <- sdat$Cell_type == top_ct
  
  return(sdat)
}



# TODO:

get_ct_avg <- function(sdat, assay = "RNA", scale = FALSE, ncores = 8) {
  
  stopifnot("Cell_type" %in% colnames(sdat@meta.data))
  
  cts <- unique(sdat$Cell_type)
  
  if (scale) {
    mat <- sdat@assays[[assay]]@scale.data
  } else {
    mat <- sdat@assays[[assay]]@data
  }
  
  ct_avg <- mclapply(cts, function(x) {
    rowMeans(mat[, sdat$Cell_type == x])
  }, mc.cores = ncores)
  
  ct_avg <- do.call(cbind, ct_avg)
  rownames(ct_avg) <- rownames(mat)
  colnames(ct_avg) <- cts
  
  return(ct_avg)
}



# TODO:

mat_to_df <- function(mat, symmetric = TRUE) {
  
  if (symmetric) {
    df <- data.frame(
      Row = rownames(mat)[row(mat)[lower.tri(mat)]],
      Col = colnames(mat)[col(mat)[lower.tri(mat)]],
      Value = mat[lower.tri(mat)],
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Row = rownames(mat)[row(mat)],
      Col = colnames(mat)[col(mat)],
      Value = c(mat),
      stringsAsFactors = FALSE
    )
  }
  return(df)
}



# TODO


sample_expression_level <- function(sdat, targets, rank_window = 10) {
  
  stopifnot(all(targets %in% rownames(sdat)))
  
  # get the ordered average expression of genes across all cells
  avg_all <- sort(rowMeans(sdat@assays$RNA@data))
  
  # for each target sample a gene whose average expression is within rank window
  
  sample_ix <- vapply(targets, function(x) {
    
    ix_orig <- which(names(avg_all) == x)
    ix_new <- ix_orig
    
    while (ix_new == ix_orig) {
      ix_new <- ix_orig + sample(-rank_window:rank_window, 1) 
    }
    
    return(ix_new)
    
  }, FUN.VALUE = integer(1), USE.NAMES = FALSE)
  
  return(names(avg_all[sample_ix]))
}



# Performance of rankings in a data frame. 
# ------------------------------------------------------------------------------


# TODO:

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



# TODO:

all_perf_df <- function(rank_df,
                        keep_cols,
                        label_col,
                        measure = NULL) {
  
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("ROC", "PR"))
 
  perf_l <- lapply(keep_cols, function(x) {
    
    get_perf_df(
      rank_df = dplyr::arrange(rank_df, !!sym(x)),
      label_col,
      measure) %>%
      mutate(Group = x)
  })
  
  perf_df <- do.call(rbind, perf_l) %>% 
    mutate(Group = factor(Group, levels = keep_cols))
   
  return(perf_df)
}


# TODO:

all_au_perf <- function(rank_df,
                        keep_cols,
                        label_col,
                        measure = NULL) {
  
  
  stopifnot(label_col %in% colnames(rank_df), measure %in% c("AUROC", "AUPRC"))
  
  perf_l <- lapply(keep_cols, function(x) {
    
    get_au_perf(
      rank_df = dplyr::arrange(rank_df, !!sym(x)),
      label_col,
      measure)
  })
  
  names(perf_l) <- keep_cols
  
  return(perf_l)
}



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
