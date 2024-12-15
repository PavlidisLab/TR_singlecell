## 
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
library(aggtools)
library(pheatmap)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

set.seed(74)
k <- 200
min_cells <- 100
max_steps <- 10
n_iters <- 100

id <- "GSE180928"
ct <- "Microglia"


tfs <- read.delim(tfs_hg_path, stringsAsFactors = FALSE)


dat <- load_dat_list(id)[[1]]


# 
# ------------------------------------------------------------------------------


# Generate steps of the number of cells to subsample from the total number of cells.
# Steps are exponentially decided by: (prev_step / ratio) where ratio is assumed
# to be less than 1.0.

# TODO: doc max not included; max_steps means end point could not be reached
# max when min_cells=100 is 100 / (0.5 ^ 10) == 102400... maybe work backwards?
# TODO: add to function to remove steps that are less than previous step

dynamic_steps <- function(max_cells,
                          min_cells = 100,
                          max_steps = 10, 
                          ratio = 0.5) {
  
  steps <- c(min_cells)
  
  while (tail(steps, 1) < max_cells) {
    
    next_step <- ceiling(tail(steps, 1) / ratio)

    if (next_step > max_cells || length(steps) >= max_steps) {
      break
    }
    
    steps <- c(steps, next_step)
  }
  
  return(steps)
}



# TODO: doc and reduce

generate_summ_df <- function(topk_l) {
  
  topk_l <- bind_rows(topk_l)
  summ <- sapply(topk_l, summary)
  summ <- t(summ)
  summ <- as.data.frame(summ) %>% rownames_to_column(var = "Symbol")
  summ$Symbol <- factor(summ$Symbol, levels = unique(summ$Symbol))
  
  return(summ)
}



# TODO: Doc
# TODO: inner loop to function

subsample_and_topk_cor <- function(ct_mat, 
                                   steps, 
                                   n_iters,
                                   k,
                                   keep_tfs,
                                   ncores = 1) {
  
  ct_ids <- rownames(ct_mat)
  
  # Generate correlation matrix from all cells as reference
  cor_mat <- sparse_pcor(ct_mat)
  diag(cor_mat) <- 0  # remove self-cor to prevent inflated topk
  
  # Loop through each sample step size, generating cor. matrix n_iter times
  topk_at_step <- lapply(steps, function(step) {
    
    message(paste("Step =", step, Sys.time()))
    
    topk_at_iter <- lapply(1:n_iters, function(i) {
      
      sample_ids <- sample(ct_ids, step, replace = FALSE)
      ct_mat_sub <- ct_mat[sample_ids, ]
      cor_mat_sub <- sparse_pcor(ct_mat_sub)
      
      # NA cors to 0 to allow overlap and remove self cor
      cor_mat_sub[is.na(cor_mat_sub)] <- 0  
      diag(cor_mat_sub) <- 0 
      
      # Only consider TFs for overlap
      topk <- pair_colwise_topk(mat1 = cor_mat[, keep_tfs], 
                                mat2 = cor_mat_sub[, keep_tfs], 
                                k = k, 
                                ncores = ncores)
      
      gc(verbose = FALSE)
      
      return(topk)
      
    })
    
    # Summary of topk for each TF across iterations
    generate_summ_df(topk_at_iter)
    
  })
  names(topk_at_step) <- paste0("step_", steps)
  
  return(topk_at_step)
}



#
# ------------------------------------------------------------------------------


ct_ids <- filter(dat$Meta, Cell_type == ct)$ID
ct_mat <- dat$Mat[, ct_ids]

keep_genes <- names(which(rowSums(ct_mat > 0) >= 20))
keep_tfs <- intersect(tfs$Symbol, keep_genes)
ct_mat <- t(ct_mat[keep_genes, ])

max_cells <- nrow(ct_mat)


steps <- dynamic_steps(max_cells = max_cells, min_cells = min_cells)



file <- file.path("/space/scratch/amorin/R_objects/TRsc", 
                  paste("subsample_cor", id, ct, n_iters, "iters.RDS", sep = "_"))



if (!file.exists(file)) {
  
  subsample_l <- subsample_and_topk_cor(
    ct_mat = ct_mat,
    steps = steps,
    n_iters = n_iters,
    k = k,
    keep_tfs = keep_tfs,
    ncores = ncore)
  
  saveRDS(subsample_l, file)
  
} else {
  
  subsample_l <- readRDS(file)
  
}




topk_barplot <- function(summ_df) {
  
  ggplot(summ_df, aes(x = Symbol, y = Mean)) +
    geom_crossbar(aes(x = Symbol, ymin = `1st Qu.`, ymax = `3rd Qu.`)) +
    geom_point(aes(x = Symbol, y = Mean), shape = 21, colour = "firebrick") +
    ylab("Top200") +
    theme_classic() +
    theme(text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  
}




pl1 <- lapply(subsample_l, topk_barplot)


# TODO: mean vs median

mean_df <- lapply(1:length(steps), function(x) {
  
  data.frame(Group = rep(paste0("Step=", steps[x]), length(keep_tfs)),
             Topk = unlist(subsample_l[[x]]$Mean)) %>% 
    mutate(Group = factor(Group, levels = unique(Group)))
  
}) %>% 
  bind_rows()



ggplot(mean_df, aes(x = Topk, colour = Group)) +
  geom_density()



mean_df <- lapply(subsample_l, function(x) {
  
  x %>% 
    arrange(match(Symbol, keep_tfs)) %>% 
    pull(Mean)
  
}) %>% 
  bind_cols() %>% 
  as.matrix()

rownames(mean_df) <- keep_tfs
colnames(mean_df) <- paste0("Step=", steps)
mean_df <- mean_df[order(mean_df[, length(steps)]), ]


# cols <- c('#f7f4f9','#e7e1ef','#d4b9da','#c994c7','#df65b0','#e7298a','#ce1256','#980043','#67001f')
cols <- viridis::inferno(10)

pheatmap(mean_df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         color = cols)


avg_expr <- colMeans(log2(ct_mat[, keep_tfs] + 1))


summary_df <- data.frame(Symbol = keep_tfs, 
                         Avg_expr = avg_expr,
                         mean_df[keep_tfs, ])



cor(select_if(summary_df, is.numeric), method = "spearman")
# plot(summary_df$Avg_expr, summary_df$Step.25600)


##

# sample_ids <- sample(ct_ids, steps[1], replace = FALSE)
# ct_mat_sub <- ct_mat[sample_ids, ]
# cor_mat_sub <- sparse_pcor(ct_mat_sub)
# 
# # NA cors to 0 to allow overlap and remove self cor
# cor_mat_sub[is.na(cor_mat_sub)] <- 0  
# diag(cor_mat) <- diag(cor_mat_sub) <- 0 
# 
# # Only consider TFs for overlap
# topk <- pair_colwise_topk(mat1 = cor_mat[, keep_tfs], 
#                           mat2 = cor_mat_sub[, keep_tfs], 
#                           k = k, 
#                           ncores = ncore)



# topk_at_step <- lapply(steps, function(step) {
#   
#   message(paste("Step =", step, Sys.time()))
#   
#   topk_l <- lapply(1:10, function(x) {
#     
#     sample_ids <- sample(ct_ids, step, replace = FALSE)
#     ct_mat_sub <- ct_mat[sample_ids, ]
#     cor_mat_sub <- sparse_pcor(ct_mat_sub)
#     
#     # NA cors to 0 to allow overlap and remove self cor
#     cor_mat_sub[is.na(cor_mat_sub)] <- 0  
#     diag(cor_mat) <- diag(cor_mat_sub) <- 0 
#     
#     # Only consider TFs for overlap
#     topk <- pair_colwise_topk(mat1 = cor_mat[, keep_tfs], 
#                               mat2 = cor_mat_sub[, keep_tfs], 
#                               k = k, 
#                               ncores = ncore)
#     
#     gc(verbose = FALSE)
#     
#     return(topk)
#     
#   })
#   
#   generate_summ_df(topk_l)
#   
# })
# names(topk_at_step) <- paste0("step_", steps)