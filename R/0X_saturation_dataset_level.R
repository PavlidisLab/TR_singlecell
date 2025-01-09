## 
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
library(aggtools)
library(ggrepel)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 200
n_iter <- 100

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
# ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes and TFs
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
tfs_hg <- read.delim(tfs_hg_path, stringsAsFactors = FALSE)
# pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
# tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
# msr_mm <- readRDS(msr_mat_mm_path)

comsr_hg <- readRDS(comsr_mat_hg_path)
# comsr_mm <- readRDS(comsr_mat_mm_path)


#
# ------------------------------------------------------------------------------


# When a TF-gene was not co-measured, impute to the median value across profiles

impute_na_to_med <- function(tf_mat, msr_mat) {
  
  ids <- colnames(tf_mat)
  msr_mat <- msr_hg[, ids]
  med <- median(tf_mat, na.rm = TRUE)
  tf_mat[msr_mat == 0] <- med
  
  return(tf_mat)
}



# Prepare a TF's coexpression profile across all datasets in which it is measured

prepare_tf_mat <- function(agg_l, tf, msr_mat) {
  
  # Bind TF profiles into one matrix, set self to NA to prevent inflated overlap
  tf_mat <- gene_vec_to_mat(agg_l, tf, msr_mat)
  tf_mat[tf, ] <- NA
  
  # Impute non-measured TF-gene pairs (tied values) to global median
  tf_ids <- colnames(tf_mat)
  msr_mat <- msr_mat[, tf_ids]
  tf_mat <- impute_na_to_med(tf_mat, msr_mat)
  
  return(tf_mat)
}



# Get the Top K overlap between individual dataset profiles and the global profile

expwise_topk <- function(tf_mat, tf_avg, k) {
  
  tf_mat <- cbind(tf_mat, Avg = tf_avg)
  topk_mat <- colwise_topk_intersect(tf_mat, k = k)
  
  topk_df <- data.frame(
    Topk = topk_mat[setdiff(rownames(topk_mat), "Avg"), "Avg"]
  ) %>% 
    rownames_to_column(var = "ID")
  
  return(topk_df)
}


# TODO:

format_summary_df <- function(topk_l, steps) {
  
  summ_df <- bind_rows(topk_l) %>% 
    as.matrix() %>%   # Needed for table -> integer data type
    as.data.frame() %>% 
    mutate(N_step = steps)
  
  return(summ_df)
}


# TODO:
# Sample procedure
# NOTE: Sampling from the full imputed matrix. I tested imputing NAs from only
# the sampled experiments, which did not make a difference, so using full for
# simplicity!


topk_at_subsample_steps <- function(tf_mat, 
                                    tf_avg,
                                    tf_ids,
                                    steps, 
                                    n_iter, 
                                    ncore) {
  
  topk_l <- mclapply(steps, function(step) {
    
    topk_at_step <- lapply(1:n_iter, function(iter) {
      
      sample_ids <- sample(tf_ids, size = step, replace = FALSE)
      sample_avg <- rowMeans(tf_mat[, sample_ids], na.rm = TRUE)
      
      topk_intersect(
        topk_sort(sample_avg, k = k),
        topk_sort(tf_avg, k = k)
      )
      
    })
    
    summary(unlist(topk_at_step))
    
  }, mc.cores = ncore)
  
  summ_df <- format_summary_df(topk_l, steps)
  
  return(summ_df)
}




# TODO:
# TODO: remember hard coded parallel

subsample_tf_topk <- function(tf, agg_l, msr_mat, k, ncore) {
  
  # Prepare global collection of TF profiles and the global average
  tf_mat <- prepare_tf_mat(agg_l = agg_tf_hg, tf = tf, msr_mat = msr_hg)
  tf_avg <- rowMeans(tf_mat, na.rm = TRUE)
  
  # Get the overlap of individual profiles with global aggregate
  exp_topk <- expwise_topk(tf_mat = tf_mat, tf_avg = tf_avg, k = k)
  
  # Count of samples from 2 to all experiments but one
  tf_ids <- colnames(tf_mat)
  steps <- 2:(length(tf_ids) - 1)
  
  # Calculate the top k overlap of each steps samples compared to global
  topk_summ <- topk_at_subsample_steps(tf_mat = tf_mat,
                                       tf_avg = tf_avg,
                                       tf_ids = tf_ids,
                                       steps = steps,
                                       n_iter = n_iter, 
                                       ncore = 4)
  
  return(list(Experiments = exp_topk, 
              Steps = topk_summ))
  
}



# TODO:
# TODO: where to put parallel

subsample_all_tfs <- function(agg_l, tfs, msr_mat, k, ncore) {
  
  tf_l <- mclapply(tfs, function(tf) {
    
    message(tf, paste(Sys.time()))
    
    tryCatch({
      tf_topk <- subsample_tf_topk(tf = tf, 
                                   agg_l = agg_tf_hg,
                                   msr_mat = msr_hg, 
                                   k = k, 
                                   ncore = ncore)
      }, 
      error = function(e) {
        NA
      }
    )
  }, mc.cores = ncore)
  
  names(tf_l) <- tfs
  
  return(tf_l)
  
}



# 
# ------------------------------------------------------------------------------


file_hg <- "/space/scratch/amorin/TRsc_output/dataset_subsample_cor_hg.RDS"


if (!file.exists(file_hg)) {
  
  agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, 
                                    ids = ids_hg, 
                                    genes = pc_hg$Symbol, 
                                    sub_genes = tfs_hg$Symbol)
  

  tf_l <- subsample_all_tfs(agg_l = agg_tf_hg,
                            tfs = tfs_hg$Symbol,
                            msr_mat = msr_hg, 
                            k = k,
                            ncore = ncore)
  
  saveRDS(tf_l, file_hg)
  
} else {
  
  tf_l <- readRDS(file_hg)
  
}



tf_l <- tf_l[!is.na(tf_l)]




# How many sampled datasets are required to reach 80% overlap on average?
# ------------------------------------------------------------------------------


ndat_for_recovery <- function(tf_l, recover = 0.8) {
  
  tfs <- names(tf_l)
  
  rec_l <- lapply(tfs, function(tf) {
    
    summ_df <- tf_l[[tf]]$Steps
    
    min_step <- summ_df %>% 
      filter(Mean > (k * recover)) %>% 
      slice_min(Mean, n = 1) %>% 
      pull(N_step)
    
    data.frame(Min_step = min_step,
               Prop_total = min_step / (max(summ_df$N_step) + 1))
    
  })
  
  names(rec_l) <- tfs
  
  df <- do.call(rbind, rec_l) %>% rownames_to_column(var = "Symbol")
  
  return(df)
}



rec_df <- ndat_for_recovery(tf_l)


summ_rec <- summary(rec_df)


# Check against average pairwise similarity of experiments

pair_sim_l <- readRDS(paste0("/space/scratch/amorin/TRsc_output/human_mouse_topk=", k, "_similarity_df.RDS"))


rec_df <- left_join(rec_df, pair_sim_l$Human, by = "Symbol")


# TODO: Note Pax4 -- higher mean and high proportion... n=16, not sure what this means

ggplot(rec_df, aes(x = Mean, y = Prop_total)) +
  geom_point(shape = 21) +
  ylab("Proportion recovery") +
  xlab("Average pairwise Top200") +
  theme_classic() +
  theme(text = element_text(size = 20))


cor(rec_df$Prop_total, rec_df$Mean, method = "spearman", use = "pairwise.complete.obs")


ggplot(rec_df, aes(x = Prop_total)) +
  geom_histogram(bins = 100) +
  ylab("Count TRs") +
  xlab("Proportion of recovery") +
  theme_classic() +
  theme(text = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



# Tally datasets that had the highest overlap with global
# ------------------------------------------------------------------------------


rank_expwise_topk <- function(tf_l) {
  
  tfs <- names(tf_l)
  
  rank_l <- lapply(tfs, function(tf) {
    
    summ_df <- tf_l[[tf]]$Experiments
    
    summ_df %>% 
      mutate(Rank = rank(-Topk, ties.method = "random"),
             # RS = rank(Topk) / nrow(summ_df))
             RS = rank(-Rank) / nrow(summ_df))
    
  })
  names(rank_l) <- tfs
  
  return(rank_l)
}




rankstandard_to_mat <- function(rank_l, ids) {
  
  tfs <- names(rank_l)
  exp_vec <- setNames(rep(0, length(ids)), ids)
  
  rs_l <- lapply(tfs, function(tf) {
    tf_rank <- rank_l[[tf]]
    exp_vec[match(tf_rank$ID, ids)] <- tf_rank$RS
    exp_vec
  })
  
  rs_mat <- do.call(rbind, rs_l)
  colnames(rs_mat) <- ids
  rownames(rs_mat) <- tfs
  
  return(rs_mat)
}




#### Tallying #1
# rank_expwise_topk <- function(tf_l) {
#   
#   tfs <- names(tf_l)
#   
#   rank_l <- lapply(tfs, function(tf) {
#     
#     summ_df <- tf_l[[tf]]$Experiments
#     
#     summ_df %>% 
#       slice_max(Topk) %>% 
#       pull(ID)
#     
#   })
#   names(rank_l) <- tfs
#   
#   return(rank_l)
# }
# 
# rank_l <- rank_expwise_topk(tf_l)
# table(unlist(rank_l)) %>% view
####





rank_l <- rank_expwise_topk(tf_l)
rs_mat <- rankstandard_to_mat(rank_l, ids_hg)


# Consider just raw TopK


topk_to_mat <- function(tf_l, ids) {
  
  tfs <- names(tf_l)
  exp_vec <- setNames(rep(0, length(ids)), ids)
  
  topk_l <- lapply(tfs, function(tf) {
    tf_topk <- tf_l[[tf]]$Experiment
    exp_vec[match(tf_topk$ID, ids)] <- tf_topk$Topk
    exp_vec
  })
  
  topk_mat <- do.call(rbind, topk_l)
  colnames(topk_mat) <- ids
  rownames(topk_mat) <- tfs
  
  return(topk_mat)
}



topk_mat <- topk_to_mat(tf_l, ids_hg)


# Average rank standardized exp-wise top K overlap with global profile to see
# which experiments most commonly aligned with global

expwise_mean <- sort(colMeans(rs_mat), decreasing = TRUE)



#
# ------------------------------------------------------------------------------


plot_tf <- "PAX6"


pdf1 <- tf_l[[plot_tf]]$Experiments %>% 
  mutate(Label = ID %in% slice_max(., Topk, n = 1)$ID,
         Group = factor(1, levels = 1))
  


p1 <- 
  ggplot(pdf1, aes(y = Topk, x =  Group)) +
  geom_boxplot(width = 0.1) +
  geom_jitter(aes(x = Group, y = Topk), 
              shape = 21, colour = "darkblue", width = 0.02, height = 0) +
  geom_label_repel(data = filter(pdf1, Label),
                  aes(x = Group, y = Topk, label = ID),
                  nudge_x = 0.3, nudge_y = 5, size = 4) +
  ylim(c(0, k)) +
  ylab("Top200") +
  xlab("Individual experiments") +
  ggtitle(plot_tf) +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text.x = element_blank())
  


pdf2 <- tf_l[[plot_tf]]$Steps
min_step <- filter(rec_df, Symbol == plot_tf)$Min_step


p2 <- 
  ggplot(pdf2, aes(x = N_step, y = Mean)) +
  geom_crossbar(aes(x = N_step, ymin = `1st Qu.`, ymax = `3rd Qu.`)) +
  geom_point(shape = 19, colour = "firebrick") +
  geom_vline(xintercept = min_step, linetype = "dashed", colour = "black") +
  geom_hline(yintercept = (k * 0.8), linetype = "dashed", colour = "black") +
  # geom_smooth() +
  ylim(c(0, k)) +
  ylab("Top200") +
  xlab("Count of sampled experiments") +
  ggtitle(plot_tf) +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))


p2_reduced <- p2 + 
  ylab(NULL) +
  ggtitle(NULL)


egg::ggarrange(p1, p2_reduced, nrow = 1, widths = c(0.3, 1))


heat_rs <- pheatmap::pheatmap(rs_mat,
                              cluster_rows = FALSE,
                              show_rownames = FALSE)


cut_rs <- cutree(heat_rs$tree_col, k = 2)


cut_df <- data.frame(Group = cut_rs) %>% 
  rownames_to_column(var = "ID") %>% 
  left_join(sc_meta, by = "ID") %>% 
  mutate(Platform2 = str_replace_all(Platform, " ", ""),
         Platform2 = ifelse(is.na(Platform2), "Mixed", Platform2),
         Platform2 = ifelse(sapply(str_split(Platform2, ","), length) != 1, "Mixed", Platform2),
         Is_10X = str_detect(Platform, "^10x.*"))



boxplot(cut_df$N_genes ~ cut_df$Group)
boxplot(log10(cut_df$N_cells) ~ cut_df$Group)
boxplot(cut_df$N_celltypes ~ cut_df$Group)


table(cut_df$Platform2, cut_df$Group)
table(cut_df$Is_10X, cut_df$Group)
fisher.test(table(cut_df$Is_10X, cut_df$Group))



# hclust_rs <- hclust(d = rs_mat)
# sil_avg_hg <- cluster::silhouette(cutree(hclust_avg_hg, k = 2), as.dist(1 - cor_avg_hg))
# mean(sil_avg_hg[, "sil_width"])
# plot(sil_avg_hg)


pheatmap::pheatmap(t(data.frame(expwise_mean)), 
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   show_rownames = FALSE,
                   cellheight = 25,
                   fontsize = 15)
