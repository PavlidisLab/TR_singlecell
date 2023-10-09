## Summarize the similarity (by topk intersect) of TF genes with themselves
## across datasets, and contrast to null intersects and L/S ribosomal genes
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 200

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# List of paired experiment similarities for TFs
# sim_tf_hg <- readRDS(sim_tf_hg_path)
sim_tf_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_TF_hg_k=", k, ".RDS"))
# sim_tf_mm <- readRDS(sim_tf_mm_path)
sim_tf_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_TF_mm_k=", k, ".RDS"))

# L/S ribo topk overlap
# sim_ribo_hg <- readRDS(sim_ribo_hg_path)
sim_ribo_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_ribo_hg_k=", k, ".RDS"))
# sim_ribo_mm <- readRDS(sim_ribo_mm_path)
sim_ribo_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_ribo_mm_k=", k, ".RDS"))

# Null topk overlap
# sim_null_hg <- readRDS(sim_null_hg_path)
sim_null_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_null_hg_k=", k, ".RDS"))
# sim_null_mm <- readRDS(sim_null_mm_path)
sim_null_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/similarity_null_mm_k=", k, ".RDS"))



# Functions
# ------------------------------------------------------------------------------


# Returns a list of dataframes of summary stats for each TF's similarity

get_summary_df <- function(sim_l, msr_mat = NULL, add_nexp = TRUE) {
  
  stats <- c("Scor", "Topk", "Bottomk", "Jaccard")
  
  df_l <- lapply(stats, function(stat) {
    
    df <- 
      lapply(sim_l, function(x) summary(x[[stat]])) %>% 
      do.call(rbind, .) %>%
      as.data.frame() %>% 
      rownames_to_column(var = "Symbol") %>% 
      arrange(desc(Median)) %>% 
      mutate(Symbol = factor(Symbol, levels = unique(Symbol)))
    
    if (add_nexp) df[["N_exp"]] <-  rowSums(msr_mat[as.character(df$Symbol), ])
    
    return(df)
  })
  
  names(df_l) <- stats
  return(df_l)
}


# Organizing summary dataframes of TF similarity across experiment pairs
# ------------------------------------------------------------------------------


# Summarize each TF's similarity and organize into a df
summ_tf_hg <- get_summary_df(sim_tf_hg, msr_hg)
summ_tf_mm <- get_summary_df(sim_tf_mm, msr_mm)

# Join TF summary with TF family information
summ_tf_hg <- lapply(summ_tf_hg, left_join, tfs_hg[, c("Symbol", "Family")], by = "Symbol")
summ_tf_mm <- lapply(summ_tf_mm, left_join, tfs_mm[, c("Symbol", "Family")], by = "Symbol")


# Require a minimum count of measured experiments for global trends
summ_sub_tf_hg <- lapply(summ_tf_hg, filter, N_exp >= 5)
summ_sub_tf_mm <- lapply(summ_tf_mm, filter, N_exp >= 5)


# Summarize ribo similarity and organize into a df
summ_ribo_hg <- get_summary_df(sim_ribo_hg, msr_hg)
summ_ribo_mm <- get_summary_df(sim_ribo_mm, msr_mm)


# Summary of null topk overlap
summ_null_hg <- get_summary_df(sim_null_hg, add_nexp = FALSE)
summ_null_mm <- get_summary_df(sim_null_mm, add_nexp = FALSE)


# summary(summ_ribo_hg$Topk[, "Mean"])
# summary(summ_ribo_mm$Topk[, "Mean"])

# summary(summ_null_hg$Topk[, "Mean"])
# summary(summ_null_mm$Topk[, "Mean"])

# summary(summ_sub_tf_hg$Topk[, "Mean"])
# summary(summ_sub_tf_mm$Topk[, "Mean"])



# Count of TFs whose expected similarity is greater than the expected null
# ------------------------------------------------------------------------------


calc_prop_gt_null <- function(sim_l, null_l) {
  
  stats <- intersect(names(sim_l), names(null_l))
  n <- nrow(sim_l[[1]])
  
  prop_l <- lapply(stats, function(x) {
    null_mean <- mean(null_l[[x]]$Mean)
    sum(sim_l[[x]]$Mean > null_mean) / n
  })
  
  names(prop_l) <- stats
  return(prop_l)
}


prop_gt_null_hg <- calc_prop_gt_null(sim_l = summ_sub_tf_hg, null_l = summ_null_hg)
prop_gt_null_mm <- calc_prop_gt_null(sim_l = summ_sub_tf_mm, null_l = summ_null_mm)



# Are TFs that are below null in Topk above null in bottom k?
# ------------------------------------------------------------------------------


calc_prop_divergent <- function(topk_df, bottomk_df, null_topk, null_bottomk) {
  
  join_df <- left_join(topk_df, bottomk_df,
                       by = "Symbol",
                       suffix = c("_Topk", "_Bottomk"))
  
  null_topk_mean <- mean(null_topk[, "Mean"])
  null_bottomk_mean <- mean(null_bottomk[, "Mean"])
  
  topk_lt_null <- filter(join_df, Mean_Topk < null_topk_mean)
  bottomk_gt_null <- filter(join_df, Mean_Bottomk > null_bottomk_mean)
  
  div <- intersect(topk_lt_null$Symbol, bottomk_gt_null$Symbol)
  div_df <- filter(join_df, Symbol %in% div)
  
  summary_df <- data.frame(N_topk_lt_null = nrow(topk_lt_null),
                           N_bottomk_gt_null = nrow(bottomk_gt_null),
                           Prop_divergent = nrow(div_df) / nrow(topk_lt_null))
  
  return(list(Summary = summary_df, 
              Below_null = topk_lt_null,
              Divergent = div_df))
}




divergent_hg <- calc_prop_divergent(topk_df = summ_sub_tf_hg$Topk,
                                    bottomk_df = summ_sub_tf_hg$Bottomk,
                                    null_topk = summ_null_hg$Topk,
                                    null_bottomk = summ_null_hg$Bottomk)


divergent_mm <- calc_prop_divergent(topk_df = summ_sub_tf_mm$Topk,
                                    bottomk_df = summ_sub_tf_mm$Bottomk,
                                    null_topk = summ_null_mm$Topk,
                                    null_bottomk = summ_null_mm$Bottomk)



# Examples of TFs with good measurement coverage but low similarity. 
# summ_sub_tf_hg$Topk %>% filter(N_exp > 80) %>% slice_min(Mean, n = 10)
# summ_sub_tf_mm$Topk %>% filter(N_exp > 80) %>% slice_min(Mean, n = 10)
# filter(summ_sub_tf_hg$Bottomk, Symbol == "ZSCAN22")
# filter(summ_sub_tf_mm$Bottomk, Symbol == "Zfp661")



# Examples of most similar individual experiments pairs across TFs
# ------------------------------------------------------------------------------


get_best_pairs <- function(summ_df_l, sim_l, n_tfs = 3, n_pairs = 3) {
  
  stats <- names(summ_df_l)
  
  stat_l <- lapply(stats, function(stat) {
    tfs <- as.character(slice_max(summ_df_l[[stat]], Max., n = n_tfs)$Symbol)
    lapply(sim_l[tfs], slice_max, !!sym(stat), n = n_pairs)
  })
  
  names(stat_l) <- stats
  return(stat_l)
}


best_pairs_tf_hg <- get_best_pairs(summ_df_l = summ_sub_tf_hg, sim_l = sim_tf_hg)
best_pairs_tf_mm <- get_best_pairs(summ_df_l = summ_sub_tf_mm, sim_l = sim_tf_mm)

best_pairs_ribo_hg <- get_best_pairs(summ_df_l = summ_ribo_hg, sim_l = sim_ribo_hg)
best_pairs_ribo_mm <- get_best_pairs(summ_df_l = summ_ribo_mm, sim_l = sim_ribo_mm)



# Rank order of similarity between species
# ------------------------------------------------------------------------------


ortho_rank_sim <- lapply(names(summ_sub_tf_hg), function(x) {
  
  df_hg <- summ_sub_tf_hg[[x]] %>% 
    filter(Symbol %in% pc_ortho$Symbol_hg) %>% 
    left_join(., pc_ortho, by = c("Symbol" = "Symbol_hg"))
  
  df_mm <- summ_sub_tf_mm[[x]] %>% 
    filter(Symbol %in% pc_ortho$Symbol_mm) %>% 
    left_join(., pc_ortho, by = c("Symbol" = "Symbol_mm"))
  
  left_join(df_hg, df_mm, by = "ID", suffix = c("_human", "_mouse"))
})
names(ortho_rank_sim) <- names(summ_sub_tf_hg)



top_ortho_sim <- lapply(names(ortho_rank_sim), function(x) {
  
  ortho_rank_sim[[x]] %>% 
  mutate(Avg_sim = rowMeans(.[, c("Mean_human", "Mean_mouse")])) %>% 
  slice_max(Avg_sim, n = 10)

})
names(top_ortho_sim) <- names(ortho_rank_sim)



ortho_cor <- vapply(ortho_rank_sim, function(x) {
  cor(x$Mean_human, x$Mean_mouse, use = "pairwise.complete.obs", method = "spearman")
}, numeric(1))



# Look at similarity by TF family
# ------------------------------------------------------------------------------


summarize_family_similarity <- function(sim_l, tf_df) {
  
  stat_l <- lapply(sim_l, function(x) {
    
    x %>%
      group_by(Family) %>%
      summarise(
        Mean_stat = mean(Mean),
        N_TFs = n(),
        Mean_exp_msrd = mean(N_exp)
      ) %>%
      arrange(desc(Mean_stat)) %>% 
      mutate(Family = factor(Family, levels = unique(Family)))
  })
  
  names(stat_l) <- names(sim_l)
  return(stat_l)
}


family_sim_hg <- summarize_family_similarity(summ_sub_tf_hg, tfs_hg)
family_sim_mm <- summarize_family_similarity(summ_sub_tf_mm, tfs_mm)




# Plotting
# ------------------------------------------------------------------------------


# Relationship between mean top and bottom k across all TFs

left_join(summ_sub_tf_hg$Topk, 
          summ_sub_tf_hg$Bottomk, 
          by = "Symbol",
          suffix = c("_Topk", "_Bottomk")
          ) %>% 
  qplot(., xvar = "Mean_Topk", yvar = "Mean_Bottomk")
  


left_join(summ_sub_tf_mm$Topk, 
          summ_sub_tf_mm$Bottomk, 
          by = "Symbol",
          suffix = c("_Topk", "_Bottomk")
          ) %>% 
  qplot(., xvar = "Mean_Topk", yvar = "Mean_Bottomk")


# Relationship between similarity for a given TF


gene_hg <- "ASCL1"
gene_mm <- "Ascl1"

qplot(sim_tf_hg[[gene_hg]], xvar = "Topk", yvar = "Bottomk")
qplot(sim_tf_mm[[gene_mm]], xvar = "Topk", yvar = "Bottomk")

qplot(sim_tf_hg[[gene_hg]], xvar = "Scor", yvar = "Jaccard")
qplot(sim_tf_mm[[gene_mm]], xvar = "Scor", yvar = "Jaccard")



# Point plot of TF similarity with text overlay for top n TFs


# Labels to the left
point_plot_similarity1 <- function(summary_df,
                                   null_df,
                                   topn_label = 30,
                                   plot_title,
                                   ylabel) {
  
  null_expected <- mean(null_df[, "Mean"])
  
  topn_genes <- slice_max(summary_df, Mean, n = topn_label)$Symbol
  
  summary_df <- summary_df %>%
    arrange(Mean) %>% 
    mutate(
      Group = Symbol %in% topn_genes,
      Symbol = factor(Symbol, levels = unique(Symbol)))
    
  ggplot(summary_df, aes(y = Mean, x = Symbol)) +
    geom_point() +
    geom_text_repel(data = filter(summary_df, Group),
                    aes(x = Symbol, y = Mean, label = Symbol, fontface = "italic"),
                    max.overlaps = topn_label,
                    force = 0.5,
                    nudge_x = -500,
                    nudge_y = 150,
                    hjust = 1,
                    direction = "y",
                    size = 5,
                    segment.size = 0.2) +
    geom_hline(yintercept = null_expected, colour = "firebrick") +
    annotate("text", x = 60, y = null_expected + 1, label = "Null", colour = "firebrick", size = 5) +
    ggtitle(plot_title) +
    ylab(ylabel) +
    expand_limits(x = nrow(summary_df) + 50) +  # prevent point cut off
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 10, 10, 10)))
}



# Labels to the right
point_plot_similarity2 <- function(summary_df,
                                   null_df,
                                   topn_label = 30,
                                   plot_title,
                                   ylabel) {
  
  null_expected <- mean(null_df[, "Mean"])
  
  topn_genes <- slice_max(summary_df, Mean, n = topn_label)$Symbol
  
  summary_df <- summary_df %>%
    arrange(Mean) %>% 
    mutate(
      Group = Symbol %in% topn_genes,
      Symbol = factor(Symbol, levels = unique(Symbol)))
  
  ggplot(summary_df, aes(y = Mean, x = Symbol)) +
    geom_point() +
    geom_text_repel(data = filter(summary_df, Group),
                    aes(x = Symbol, y = Mean, label = Symbol, fontface = "italic"),
                    max.overlaps = topn_label, 
                    force = 0.5,
                    nudge_x = 300,
                    hjust = 0,
                    direction = "y",
                    size = 5,
                    segment.size = 0.1,
                    segment.color = "grey50") +
    geom_hline(yintercept = null_expected, colour = "firebrick") +
    annotate("text", x = 60, y = null_expected + 1, label = "Null", colour = "firebrick", size = 5) +
    ggtitle(plot_title) +
    ylab(ylabel) +
    expand_limits(x = nrow(summary_df) + 300) +  # prevent point cut off
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 10, 10, 10)))
}



# Using text labels to the left
px_hg1 <- point_plot_similarity1(summ_tf_hg$Topk,
                                 null_df = summ_null_hg$Topk,
                                 plot_title = "Human",
                                 ylabel = paste0("Mean Top k (k=", k, ")"))


px_mm1 <- point_plot_similarity1(summ_tf_mm$Topk,
                                 null_df = summ_null_mm$Topk,
                                 plot_title = "Mouse",
                                 ylabel = paste0("Mean Top k (k=", k, ")"))


# Using text labels to the right
px_hg2 <- point_plot_similarity2(summ_tf_hg$Topk,
                                 null_df = summ_null_hg$Topk,
                                 plot_title = "Human",
                                 ylabel = paste0("Mean Top k (k=", k, ")"))


px_mm2 <- point_plot_similarity2(summ_tf_mm$Topk,
                                 null_df = summ_null_mm$Topk,
                                 plot_title = "Mouse",
                                 ylabel = paste0("Mean Top k (k=", k, ")"))



# Histogram of expected values for TF, null, and ribosomal


hist_expected <- function(summary_df, null_df, ribo_df, xlabel) {
  
  expected_df <- 
    data.frame(Stat = c(summary_df[, "Mean"],
                        null_df[, "Mean"],
                        ribo_df[, "Mean"]),
               Group = c(rep("TF", nrow(summary_df)),
                         rep("Null", nrow(null_df)),
                         rep("Ribosomal", nrow(ribo_df))))
  
  ggplot(expected_df, aes(x = Stat)) + 
    facet_wrap(~Group, nrow = 3, scales = "free") + 
    geom_histogram(bins = 30) +
    xlab(xlabel) +
    ylab("Count") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          strip.background = element_rect(color = "black", fill = c("lightgrey")),
          plot.margin = margin(c(10, 20, 10, 10)))
}



py_hg <- hist_expected(summary_df = summ_tf_hg$Topk,
                       null_df = summ_null_hg$Topk,
                       ribo_df = summ_ribo_hg$Topk,
                       xlabel = paste0("Mean Top k (k=", k, ")"))


py_mm <- hist_expected(summary_df = summ_tf_mm$Topk,
                       null_df = summ_null_mm$Topk,
                       ribo_df = summ_ribo_mm$Topk,
                       xlabel = paste0("Mean Top k (k=", k, ")"))



# Combining point plot and hist of expected (labels left)
pxy_hg1 <- plot_grid(py_hg, px_hg1, rel_widths = c(0.5, 1))
pxy_mm1 <- plot_grid(py_mm, px_mm1, rel_widths = c(0.5, 1))

# Combining point plot and hist of expected (labels right)
pxy_hg2 <- plot_grid(py_hg, px_hg2, rel_widths = c(0.5, 1))
pxy_mm2 <- plot_grid(py_mm, px_mm2, rel_widths = c(0.5, 1))

# Flip mouse histogram for combining with human panels
pxy_mm3 <- plot_grid(px_mm1, py_mm, rel_widths = c(1, 0.5))
pxy_both <- plot_grid(pxy_hg1, pxy_mm3, nrow = 1)


# ggsave(pxy_both, height = 8, width = 20, device = "png", dpi = 300,
#        filename = file.path(plot_dir, "tt.png"))



# Boxplot (show 1st and 3rd IQR) of TF median topk intersect


boxplot_topk_median <- function(summary_df, null_df, plot_title, ylabel) {
  
  null_expected <- mean(null_df[, "Mean"])
  
  summary_df <- summary_df %>%
    arrange(Mean) %>% 
    mutate(Symbol = factor(Symbol, levels = unique(Symbol)))
  
  ggplot(summary_df, aes(x = Symbol)) +
    geom_boxplot(
      # aes(ymin = Min., lower = `1st Qu.`, middle = Median, upper = `3rd Qu.`, ymax = Max.),
      aes(ymin = `1st Qu.`, lower = `1st Qu.`, middle = Median, upper = `3rd Qu.`, ymax = `3rd Qu.`),
      stat = "identity") +
    geom_hline(yintercept = null_expected , colour = "firebrick", linewidth = 1.6) +
    ylab(ylabel) +
    ggtitle(plot_title) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}


pza <- boxplot_topk_median(summary_df = summ_tf_hg$Topk, 
                           null_df = summ_null_hg$Topk,
                           plot_title = "Human",
                           ylabel = paste0("Mean Top k (k=", k, ")"))

pzb <- boxplot_topk_median(summary_df = summ_tf_mm$Topk, 
                           null_df = summ_null_mm$Topk,
                           plot_title = "Mouse",
                           ylabel = paste0("Mean Top k (k=", k, ")"))



# Plot divergent: mean bottom k versus null of TFs with expected topk below null


plot_divergent <- function(divergent_l, null_df, k, plot_title, nlabel = 30) {
  
  plot_df <- data.frame(
    
    Mean_bottomk = c(divergent_l$Below_null$Mean_Bottomk, null_df$Mean),
    
    Symbol = c(divergent_l$Below_null$Symbol, 
               rep("Null", length(null_df$Mean))),
    
    Group = c(rep("TF", length(divergent_l$Below_null$Mean_Bottomk)),
              rep("Null", length(null_df$Mean)))
  )
  
  
  div_symbol <- slice_max(divergent_l$Divergent, Mean_Bottomk, n = nlabel)$Symbol
  
  plot_df$Label <- plot_df$Symbol %in% div_symbol
  
  
  ggplot(plot_df, aes(x = Group, y = Mean_bottomk)) +
    
    geom_violin(width = 0.3, fill = "slategrey") +
    
    geom_boxplot(
      width = 0.1,
      fill = "white",
      # outlier.shape = NA,
      outlier.fill = "white"
    ) +
    geom_point(
      data = filter(plot_df, Label),
      shape = 21,
      alpha = 0.6,
      size = 3.1
    ) +
    geom_text_repel(
      data = filter(plot_df, Label),
      aes(x = Group, y = Mean_bottomk, label = Symbol, fontface = "italic"),
      max.overlaps = 30,
      force = 0.5,
      nudge_x = 500,
      hjust = 0,
      direction = "y",
      size = 5,
      segment.size = 0.1,
      segment.color = "grey50") +
    ggtitle(plot_title) +
    ylab(paste0("Mean Bottom k=", k)) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}


pma <- plot_divergent(divergent_l = divergent_hg,
                      null_df = summ_null_hg$Bottomk,
                      k = k,
                      plot_title = "Human")


pmb <- plot_divergent(divergent_l = divergent_mm,
                      null_df = summ_null_mm$Bottomk,
                      k = k,
                      plot_title = "Mouse")


pm <- plot_grid(pma, pmb, nrow = 1)



# Density plot of topk intersect for select genes + null


density_topk <- function(plot_df, k) {
  
  ggplot(plot_df, aes(x = Topk, fill = Group)) +
    geom_density(alpha = 0.6) +
    theme_classic() +
    ylab("Density") +
    xlab(paste0("Top k=", k)) +
    scale_fill_manual(values = c("lightgrey", "#7570b3", "#d95f02", "#1b9e77")) +
    theme(
      axis.text = element_text(size = 30),
      axis.title = element_text(size = 30),
      plot.title = element_text(hjust = 0.5, size = 30),
      legend.position = c(0.75, 0.90),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(10, 30, 10, 10))
}


# Isolating the max TF by topk similarity, a representative TF, a ribosomal gene,
# and a null

example_tf_hg <- "ASCL1"
example_tf_mm <- "Ascl1"

max_tf_hg <- "E2F8"
max_tf_mm <- "E2f8"
# max_tf_hg <- as.character(slice_max(summ_sub_tf_hg$Topk, Median)$Symbol)
# max_tf_mm <- as.character(slice_max(summ_sub_tf_mm$Topk, Median)$Symbol)

example_ribo_hg <- "RPL32"
example_ribo_mm <- "Rpl32"

max_tf_df_hg <- sim_tf_hg[[max_tf_hg]]
max_tf_df_mm <- sim_tf_mm[[max_tf_mm]]

tf_df_hg <- sim_tf_hg[[example_tf_hg]]
tf_df_mm <- sim_tf_mm[[example_tf_mm]]

ribo_df_hg <- sim_ribo_hg[[example_ribo_hg]]
ribo_df_mm <- sim_ribo_mm[[example_ribo_mm]]

# Pull a random sample
set.seed(154)
rep_null_hg <- sim_null_hg[[sample(1:length(sim_null_hg), 1)]]
rep_null_mm <- sim_null_mm[[sample(1:length(sim_null_mm), 1)]]


plot_df_hg <- data.frame(
  Group = c(
    rep("Null", nrow(rep_null_hg)),
    rep(example_ribo_hg, nrow(ribo_df_hg)),
    rep(max_tf_hg, nrow(max_tf_df_hg)),
    rep(example_tf_hg, nrow(tf_df_hg))
  ), 
  Topk = c(rep_null_hg$Topk, ribo_df_hg$Topk, max_tf_df_hg$Topk, tf_df_hg$Topk))

plot_df_hg$Group <- factor(plot_df_hg$Group, levels = unique(plot_df_hg$Group))


plot_df_mm <- data.frame(
  Group = c(
    rep("Null", nrow(rep_null_mm)),
    rep(example_ribo_mm, nrow(ribo_df_mm)),
    rep(max_tf_mm, nrow(max_tf_df_mm)),
    rep(example_tf_mm, nrow(tf_df_mm))
  ), 
  Topk = c(rep_null_mm$Topk, ribo_df_mm$Topk, max_tf_df_mm$Topk, tf_df_mm$Topk))

plot_df_mm$Group <- factor(plot_df_mm$Group, levels = unique(plot_df_mm$Group))

p4a <- density_topk(plot_df_hg, k = k)
p4b <- density_topk(plot_df_mm, k = k)
p4 <- plot_grid(p4a, p4b)

# boxplot(log10(plot_df_hg$Topk+1) ~ plot_df_hg$Group)
# boxplot(log10(plot_df_mm$Topk+1) ~ plot_df_mm$Group)

# ggsave(p4, height = 9, width = 18, device = "png", dpi = 300,
#        filename = file.path(plot_dir, "density_topk_example.png"))



# Scatterplot of expected similarity of ortho TFs


qplot(ortho_rank_sim$Topk, xvar = "Mean_human", yvar = "Mean_mouse") +
  geom_smooth() +
  xlab(paste0("Human Mean Top k=", k)) +
  ylab(paste0("Mouse Mean Top k=", k)) +
  ggtitle("Top")

qplot(ortho_rank_sim$Bottomk, xvar = "Mean_human", yvar = "Mean_mouse") +
  geom_smooth() +
  xlab(paste0("Human Mean Bottom k=", k)) +
  ylab(paste0("Mouse Mean Bottom k=", k)) +
  ggtitle("Bottom")


qplot(ortho_rank_sim$Jaccard, xvar = "Mean_human", yvar = "Mean_mouse")


# Relationship between topk and count of experiments

topk_na_cor_hg <- suppressWarnings(cor.test(summ_tf_hg$Topk$Mean, summ_tf_hg$Topk$N_exp, method = "spearman"))
topk_na_cor_mm <- suppressWarnings(cor.test(summ_tf_mm$Topk$Mean, summ_tf_mm$Topk$N_exp, method = "spearman"))
qplot(summ_tf_hg$Topk, xvar = "N_exp", yvar = "Mean")
qplot(summ_tf_mm$Topk, xvar = "N_exp", yvar = "Mean")


# Similarity by TF family


plot_family_similarity <- function(summary_df) {
  
  ggplot(summary_df, 
         aes(x = fct_reorder(Family, desc(Mean), .fun = mean), 
             y = Mean)) +
    geom_boxplot(fill = "slategrey", outlier.shape = NA) +
    geom_jitter(width = 0.25, height = 0.1, shape = 21) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = margin(c(10, 20, 10, 10)))
  
}


plot_family_similarity(summ_sub_tf_hg$Topk)
plot_family_similarity(summ_sub_tf_mm$Topk)
