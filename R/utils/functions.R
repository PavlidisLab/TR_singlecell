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





# Temp plot functions


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
