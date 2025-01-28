## Perform an iterative subsampling procedure to see how many datasets on average
## are required to recover the global/aggregate profiles by Topk overlap
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
library(aggtools)
library(ggrepel)
library(egg)
library(pheatmap)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 200
n_iter <- 100  # how many iterations per subsample step
set.seed(5)

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes and TFs
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
tfs_hg <- read.delim(tfs_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# List of dataset pair-wise topk similarity metrics
pair_sim_l <- readRDS(paste0("/space/scratch/amorin/TRsc_output/human_mouse_topk=", k, "_similarity_df.RDS"))

# Output lists
subsample_path_hg <- paste0("/space/scratch/amorin/TRsc_output/dataset_subsample_topk=", k, "_hg.RDS")
subsample_path_mm <- paste0("/space/scratch/amorin/TRsc_output/dataset_subsample_topk=", k, "_mm.RDS")



# Functions
# ------------------------------------------------------------------------------


# Impute gene pairs that are not co-measured (and thus are NA -> 0 -> ties during
# ranking) to the median value of the gene profile matrix, ensuring they have a 
# middling value

impute_na_to_med <- function(gene_mat, msr_mat) {
  
  med <- median(gene_mat, na.rm = TRUE)
  gene_mat[msr_mat == 0] <- med
  
  return(gene_mat)
}



# Collect a gene's coexpression profile across all datasets in which it is 
# measured into a single gene x dataset experiment with unmeasured/NA/tied
# values imputed to the median of the whole matrix (gene profile matrix). Finally,
# append a column of the average of the matrix

prepare_gene_mat <- function(agg_l, gene, msr_mat) {
  
  # Bind gene profiles into one matrix, set self to NA to prevent inflated overlap
  gene_mat <- gene_vec_to_mat(agg_l, gene, msr_mat)
  gene_mat[gene, ] <- NA
  
  # Impute non-measured TF-gene pairs (tied values) to global median
  exp_ids <- colnames(gene_mat)
  msr_mat <- msr_mat[, exp_ids]
  gene_mat <- impute_na_to_med(gene_mat, msr_mat)
  
  # Add the average of the matrix as a column for downstream comparisons
  gene_mat <- cbind(gene_mat, Avg = rowMeans(gene_mat, na.rm = TRUE))
  
  return(gene_mat)
}



# Get the Top K overlap between every individual dataset profiles and the global 
# average, returned as a df

gen_expwise_topk_df <- function(gene_mat, k) {
  
  topk_mat <- colwise_topk_intersect(gene_mat, k = k)
  
  topk_df <- data.frame(
    Topk = topk_mat[setdiff(rownames(topk_mat), "Avg"), "Avg"]
  ) %>% 
    rownames_to_column(var = "ID")
  
  return(topk_df)
}



# Organize the summary of topk overlap across iterations of each subsample step 
# into a df

format_summ_df_across_steps <- function(topk_step_l, steps) {
  
  summ_df <- bind_rows(topk_step_l) %>% 
    as.matrix() %>%   # Needed for table -> integer data type
    as.data.frame() %>% 
    mutate(N_step = steps)
  
  return(summ_df)
}


# Subsampling procedure: for each step in steps, perform n_iter times:
# Sample n step experiments that measure the TF, then compare the top k overlap
# of the profile generated from these sampled experiment to that of the global
# profile. Summarize topK across all iterations per step, returning the result
# as a list, where each list element is the topk summary at the given step

# NOTE: Sampling from the full NA imputed TR profile matrix. I tested imputing 
# NAs using the median from only the sampled experiments, which did not make a
# difference, so using full for simplicity!



# Helper to calculate overlap of one sample iteration

calc_single_sample_overlap <- function(gene_mat, sample_ids, k) {
  
  sample_avg <- rowMeans(gene_mat[, sample_ids], na.rm = TRUE)
  gene_avg <- gene_mat[, "Avg"]
  
  topk_intersect(
    topk_sort(sample_avg, k = k),
    topk_sort(gene_avg, k = k)
  )
}


# Helper to calculate and summarize overlap across all iterations of a step

calc_single_step_overlap <- function(gene_mat, exp_ids, k, step, n_iter) {
  
  step_l <- lapply(1:n_iter, function(iter) {
    sample_ids <- sample(exp_ids, size = step, replace = FALSE)
    calc_single_sample_overlap(gene_mat, sample_ids, k)
  })
  summary(unlist(step_l))
}



# Helper to calculate and summarize overlap across all steps

summ_topk_overlap_across_steps <- function(gene_mat, exp_ids, k, steps, n_iter) {
  
  topk_summ_by_step <- lapply(steps, function(step) {
    calc_single_step_overlap(gene_mat, exp_ids, k, step, n_iter)
  })
  
  format_summ_df_across_steps(topk_summ_by_step, steps)
}



# For a given gene, prepare its gene profile matrix and then perform the topk
# overlap procedure. Produces a list of two results: one is a dataframe of the
# individual dataset overlap to the global average; the second is a list where
# each element summarizes the topk overlap at increasing subsampling steps

analyze_gene_topk_overlap <- function(agg_l, gene, msr_mat, k, n_iter) {
  
  # Collect gene profiles from aggregate list into one matrix
  gene_mat <- prepare_gene_mat(agg_l, gene, msr_mat)

  # Dataset measuring the gene, excluding the global average
  exp_ids <- setdiff(colnames(gene_mat), "Avg")
  
  # Steps as the count of datasets to sample, from 2 to all experiments but one
  steps <- 2:(length(exp_ids) - 1)
  
  # Get the overlap of individual profiles with global aggregate
  expwise_topk <- gen_expwise_topk_df(gene_mat, k)

  # Calculate the top k overlap across subsample steps
  steps_topk <- summ_topk_overlap_across_steps(gene_mat, exp_ids, k, steps, n_iter)
  
  return(list(Experiments = expwise_topk, 
              Steps = steps_topk))
}



# Perform analyze_gene_topk_overlap for every gene in gene_vec, returned as list

analyze_all_gene_topk_overlap <- function(agg_l, 
                                          gene_vec, 
                                          msr_mat, 
                                          k, 
                                          n_iter, 
                                          ncore) {
  
  gene_l <- mclapply(gene_vec, function(gene) {
    
    message(paste(gene, Sys.time()))
    
    tryCatch({
      gene_topk <- analyze_gene_topk_overlap(agg_l, gene, msr_mat, k, n_iter)
      }, error = function(e) NA
    )
  }, mc.cores = ncore)
  
  names(gene_l) <- gene_vec
  
  return(gene_l)
}



# The following functions are for analyzing the output of the topk process


# Create Gene x Dataset matrices where elements represent the standardized rank
# of the top K overlap for the given TR/dataset pair relative to the global TR profile


# Helper that ranks the experiment-wise top k overlap

rank_expwise_topk <- function(gene_l) {
  
  genes <- names(gene_l)
  
  rank_l <- lapply(genes, function(gene) {
    
    summ_df <- gene_l[[gene]]$Experiments %>% 
      mutate(Rank = rank(-Topk, ties.method = "random"),
             RS = rank(-Rank) / nrow(.))
    
  })
  names(rank_l) <- genes
  
  return(rank_l)
}



# Takes the rank list and creates a Gene x Dataset matrix of the std. rank
# values. Note that the ranking is done gene-wise; every gene has an experiment
# vector where a value of 1 means that experiment had the best match to the 
# gene's global profile, 0 the lowest/no match. Set all unmeasured pairs as 0.
# This means that a column/dataset can have multiple values of 1 (it was the 
# best match for multiple genes), but a gene/row will only have a single 1

rank_l_to_mat <- function(rank_l, ids) {
  
  genes <- names(rank_l)
  exp_vec <- setNames(rep(0, length(ids)), ids)

  rs_l <- lapply(genes, function(gene) {
    
    gene_rank <- rank_l[[gene]]
    exp_vec[match(gene_rank$ID, ids)] <- gene_rank$RS
    exp_vec
    
  })
  
  rs_mat <- do.call(rbind, rs_l)
  colnames(rs_mat) <- ids
  rownames(rs_mat) <- genes
  
  return(rs_mat)
}



# For each gene, calculate how many sampled datasets (steps) on average were 
# required to recover the global profile. Here recovery is having 80% of the 
# global profile's Top200 (i.e, overlap of at least 160/200). Also represent
# this count as the proportion of available datasets measuring the gene
# ("proportion of recovery")

ndat_for_recovery <- function(gene_l, recover = 0.8) {
  
  genes <- names(gene_l)
  
  rec_l <- lapply(genes, function(gene) {
    
    summ_df <- gene_l[[gene]]$Steps
    
    min_step <- summ_df %>% 
      filter(Mean >= (k * recover)) %>% 
      slice_min(Mean, n = 1, with_ties = FALSE) %>% 
      pull(N_step)
    
    data.frame(Min_step = min_step,
               Prop_total = min_step / (max(summ_df$N_step) + 1))
    
  })
  
  names(rec_l) <- genes
  
  df <- do.call(rbind, rec_l) %>% rownames_to_column(var = "Symbol")
  
  return(df)
}



# Fit an exponential function between Topk overlap and n_steps, of the form:
# Mean_topk_across_iter ~ TopK * (1 = exp(-b * N_step))
# To extract the b coefficient. High values of b idicate a plateauing effect - 
# the samples are converging on the global Top K with few steps. A small b means
# a linear trend -- saturation is not being reached

fit_exponential <- function(step_df, k = 200, b_init = 0.1) {
  
  tryCatch({
    nls(Mean ~ a * (1 - exp(-b * N_step)), 
        data = step_df, 
        start = list(a = k, b = b_init))  # Initial guesses for parameters
  }, error = function(e) NULL)
  
}



# Fit the exponential function to all gene step dfs, extract the b coeff, and
# organize into a summary df

summ_all_fits <- function(gene_l) {
  
  coef_l <- lapply(gene_l, function(x) {
    fit <- fit_exponential(step_df = x$Steps)
    coef(fit)["b"]
  })
  
  do.call(rbind, coef_l) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Symbol")
  
}



# Run/save/load
# ------------------------------------------------------------------------------


# Human
if (!file.exists(subsample_path_hg)) {
  
  agg_tf_hg <- load_or_generate_agg(
    path = agg_tf_hg_path,
    ids = ids_hg,
    genes = pc_hg$Symbol,
    sub_genes = tfs_hg$Symbol
  )
  
  hg_l <- analyze_all_gene_topk_overlap(
    agg_l = agg_tf_hg,
    gene_vec = tfs_hg$Symbol,
    msr_mat = msr_hg,
    k = k,
    n_iter = n_iter,
    ncore = ncore
  )
  
  saveRDS(hg_l, subsample_path_hg)
} 



# Mouse
if (!file.exists(subsample_path_mm)) {
  
  agg_tf_mm <- load_or_generate_agg(
    path = agg_tf_mm_path,
    ids = ids_mm,
    genes = pc_mm$Symbol,
    sub_genes = tfs_mm$Symbol
  )
  
  mm_l <- analyze_all_gene_topk_overlap(
    agg_l = agg_tf_mm,
    gene_vec = tfs_mm$Symbol,
    msr_mat = msr_mm,
    k = k,
    n_iter = n_iter,
    ncore = ncore
  )
  
  saveRDS(mm_l, subsample_path_mm)
} 



hg_l <- readRDS(subsample_path_hg)
hg_l <- hg_l[!is.na(hg_l)]

mm_l <- readRDS(subsample_path_mm)
mm_l <- mm_l[!is.na(mm_l)]



# Generate summary objects
# ------------------------------------------------------------------------------


# Get the Gene by Dataset ranked topk matrices

topk_rank_hg <- rank_expwise_topk(hg_l) %>% rank_l_to_mat(ids = ids_hg)
topk_rank_mm <- rank_expwise_topk(mm_l) %>% rank_l_to_mat(ids = ids_mm)


# Take the column/dataset-wise average: this metric ("Global_agreement") is
# a means of showing how closely aligned a given dataset is to the global signal,
# across all its measured TRs

agree_hg <- data.frame(
  ID = colnames(topk_rank_hg),
  Global_agreement = colMeans(topk_rank_hg, na.rm = TRUE)) %>% 
  left_join(sc_meta, by = "ID") %>% 
  mutate(Has_10X = str_detect(Platform, "^10x.*")) %>% 
  arrange(Global_agreement)


agree_mm <- data.frame(
  ID = colnames(topk_rank_mm),
  Global_agreement = colMeans(topk_rank_mm, na.rm = TRUE)) %>% 
  left_join(sc_meta, by = "ID") %>% 
  mutate(Has_10X = str_detect(Platform, "^10x.*")) %>% 
  arrange(Global_agreement)



# Summarise and join tables for the mean topk across steps, as well as the 
# mean pairwise similarity

summ_topk_hg <- ndat_for_recovery(hg_l) %>% 
  left_join(pair_sim_l$Human, by = "Symbol") %>% 
  left_join(summ_all_fits(hg_l), by = "Symbol")


summ_topk_mm <- ndat_for_recovery(mm_l) %>% 
  left_join(pair_sim_l$Mouse, by = "Symbol") %>% 
  left_join(summ_all_fits(mm_l), by = "Symbol")



cor_hg <- cor(select_if(summ_topk_hg, is.numeric), 
              method = "spearman", 
              use = "pairwise.complete.obs")


cor_mm <- cor(select_if(summ_topk_mm, is.numeric), 
              method = "spearman", 
              use = "pairwise.complete.obs")



# Plotting
# ------------------------------------------------------------------------------


# Histogram demonstrating most TRs have not reached saturation

pz <- ggplot(summ_topk_hg, aes(x = Prop_total)) +
  geom_histogram(bins = 100, col = "slategrey", fill = "slategrey") +
  ggtitle("Human") +
  ylab("Count TRs") +
  xlab("Proportion of datasets needed for recovery") +
  theme_classic() +
  theme(text = element_text(size = 25),
        plot.margin = margin(c(10, 20, 10, 10)))


ggsave(pz, height = 9, width = 9, device = "png", dpi = 300,
    filename = file.path(plot_dir, "dataset_proportion_recovery_hist_hg.png"))




# Showing TopK overlap of individual experiments to global

plot_topk_ind <- function(gene, gene_l, k) {
  
  # Create a group label for the single experiment that had the best TopK
  plot_df <- gene_l[[gene]]$Experiments %>% 
  mutate(Label = ID %in% slice_max(., Topk, n = 1, with_ties = FALSE)$ID,
         Group = factor(1, levels = 1))
  
  ggplot(plot_df, aes(y = Topk, x = Group)) +
    geom_boxplot(width = 0.1) +
    geom_jitter(aes(x = Group, y = Topk), 
                shape = 21, colour = "slategrey", width = 0.02, height = 0) +
    geom_label_repel(data = filter(plot_df, Label),
                     aes(x = Group, y = Topk, label = ID),
                     nudge_x = -0.3, nudge_y = 5, size = 6) +
    ylim(c(0, k)) +
    ylab(expr("Top"[!!k])) +
    xlab("Individual experiments") +
    ggtitle(gene) +
    theme_classic() +
    theme(text = element_text(size = 30),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}



# Showing spread of TopK across iterations per step

plot_topk_steps <- function(gene, gene_l, summ_df, k) {
  
  # Isolating steps as well as minimum steps to achieve recovery from summary df
  plot_df <- gene_l[[gene]]$Steps
  min_step <- filter(summ_df, Symbol == gene)$Min_step
  
  ggplot(plot_df, aes(x = N_step, y = Mean)) +
    geom_crossbar(aes(x = N_step, ymin = `1st Qu.`, ymax = `3rd Qu.`)) +
    geom_point(shape = 19, colour = "firebrick") +
    geom_vline(xintercept = min_step, linetype = "dashed", colour = "black") +
    geom_hline(yintercept = (k * 0.8), linetype = "dashed", colour = "black") +
    ylim(c(0, k)) +
    ylab(expr("Top"[!!k])) +
    xlab("Count of sampled experiments") +
    ggtitle(gene) +
    theme_classic() +
    theme(text = element_text(size = 30),
          axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
          plot.margin = margin(10, 20, 10, 10))
}




# Paper highlights E2F8 and PAX6 - plot together (reduce common labels)


px1a <- plot_topk_ind("E2F8", hg_l, k)
px1b <- plot_topk_steps("E2F8", hg_l, summ_topk_hg, k) + ylab(NULL) + ggtitle(NULL)
px1 <- egg::ggarrange(px1a, px1b, nrow = 1, widths = c(0.4, 1))

px2a <- plot_topk_ind("PAX6", hg_l, k)
px2b <- plot_topk_steps("PAX6", hg_l, summ_topk_hg, k) + ylab(NULL) + ggtitle(NULL)
px2 <- egg::ggarrange(px2a, px2b, nrow = 1, widths = c(0.4, 1))


ggsave(px1, height = 8, width = 16, device = "png", dpi = 300,
    filename = file.path(plot_dir, "E2F8_dataset_saturation.png"))


ggsave(px2, height = 8, width = 16, device = "png", dpi = 300,
    filename = file.path(plot_dir, paste0("PAX6_dataset_saturation.png")))



# Plotting the relationship between global agreement and experiment meta


plot_agree_scatter <- function(df, xvar, xlab) {
  
  ggplot(df, aes(x = !!sym(xvar), y = Global_agreement)) +
  geom_point(shape = 21, size = 3) +
  geom_smooth(method = "lm") +
  ylim(c(0, 0.8)) +
  xlab(xlab) +
  ylab("Global agreement") +
  theme_classic() +
  theme(text = element_text(size = 20))
}
 


plot_agree_boxplot <- function(df) {
  
  ggplot(df, aes(x = Has_10X, y = Global_agreement)) +
  geom_boxplot(width = 0.2, fill = "slategrey", alpha = 0.4) +
  geom_jitter(width = 0.1, shape = 21, size = 2.1, alpha = 0.6) +
  ylim(c(0, 0.8)) +
  xlab("Contains cells from 10X platform") +
  ylab("Aggreement with global") +
  scale_color_grey() +
  theme_classic() +
  theme(text = element_text(size = 20))
  
}


pya <- plot_agree_scatter(agree_hg, xvar = "N_genes", xlab = "Count of genes")
pyb <- plot_agree_scatter(mutate(agree_hg, N_cells = log10(N_cells)), xvar = "N_cells", xlab = expr("Log"[10] ~ "count of cells"))
pyc <- plot_agree_scatter(mutate(agree_hg, N_celltypes = log10(N_celltypes)), xvar = "N_celltypes", xlab = expr("Log"[10] ~ "count of cell types"))
pyd <- plot_agree_boxplot(agree_hg)
  

# Plot all 4 together in one row with common y lab on leftmost plot

py <- egg::ggarrange(pya, 
                     pyb + theme(axis.title.y = element_blank()), 
                     pyc + theme(axis.title.y = element_blank()), 
                     pyd + theme(axis.title.y = element_blank()),
                     nrow = 1)


ggsave(py, height = 5, width = 20, device = "png", dpi = 300,
    filename = file.path(plot_dir, "global_agreement_vs_exp_meta_hg.png"))




# Heatmap of global agreement vector for demonstration

plot_vec <- setNames(agree_hg$Global_agreement, agree_hg$ID)

pheatmap(t(plot_vec),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         cellheight = 20,
         cellwidth = 12,
         fontsize = 15,
         color = colorRampPalette(c("white", "#ca0020"))(100),
         border_color = "black",
         filename = file.path(plot_dir, "global_agreement_hg.png"))
