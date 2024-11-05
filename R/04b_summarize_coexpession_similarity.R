## Summarize the similarity (focus on topk intersect) of TF profiles across 
## datasets, and contrast to null intersects and L/S ribosomal genes
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
min_exp <- 5

# Ortho genes and TFs
pc_ortho <- read.delim(pc_ortho_path, stringsAsFactors = FALSE)
tfs_hg <- read.delim(tfs_hg_path, stringsAsFactors = FALSE)
tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Lists of paired experiment similarities for TFs, ribo, and shuffled null
sim_tf_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/TRsc/similarity_TF_hg_k=", k, ".RDS"))
sim_tf_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/TRsc/similarity_TF_mm_k=", k, ".RDS"))
sim_ribo_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/TRsc/similarity_ribo_hg_k=", k, ".RDS"))
sim_ribo_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/TRsc/similarity_ribo_mm_k=", k, ".RDS"))
sim_null_hg <- readRDS(paste0("/space/scratch/amorin/R_objects/TRsc/similarity_null_hg_k=", k, ".RDS"))
sim_null_mm <- readRDS(paste0("/space/scratch/amorin/R_objects/TRsc/similarity_null_mm_k=", k, ".RDS"))

# Export similarity summaries for downstream ortho comparison
out_path <- paste0("/space/scratch/amorin/R_objects/TRsc/human_mouse_topk=", k, "_similarity_df.RDS")



# Functions
# ------------------------------------------------------------------------------


# Assumes that each element of a sim_l is a dataframe with Row/Col (IDs) and 
# similarity stat columns. Generates a new list where each element, one for each
# stat, is a data.frame summarizing the given similarity stat for every element
# in sim_l (eg, get a dataframe of TopK summaries for all TFs)

gen_summary_df_list <- function(sim_l, msr_mat = NULL, add_nexp = TRUE) {
  
  stats <- c("Scor", "Topk", "Bottomk", "Jaccard")
  
  df_l <- lapply(stats, function(stat) {
    
    df <- 
      lapply(sim_l, function(x) summary(x[[stat]])) %>% 
      do.call(rbind, .) %>%
      as.data.frame() %>% 
      rownames_to_column(var = "Symbol") %>% 
      arrange(desc(Mean)) %>% 
      mutate(Symbol = factor(Symbol, levels = unique(Symbol)))
    
    if (add_nexp) df[["N_exp"]] <- rowSums(msr_mat[as.character(df$Symbol), ])
    
    return(df)
  })
  
  names(df_l) <- stats
  return(df_l)
}


# Organizing summary dataframes of TF similarity across experiment pairs
# ------------------------------------------------------------------------------


# Summarize each TF's similarity and organize into a df
summ_tf_hg <- gen_summary_df_list(sim_tf_hg, msr_hg)
summ_tf_mm <- gen_summary_df_list(sim_tf_mm, msr_mm)


# Join TF summary with TF family information
summ_tf_hg <- lapply(summ_tf_hg, left_join, tfs_hg[, c("Symbol", "Family")], by = "Symbol")
summ_tf_mm <- lapply(summ_tf_mm, left_join, tfs_mm[, c("Symbol", "Family")], by = "Symbol")


# Require a minimum count of measured experiments for global trends
summ_sub_tf_hg <- lapply(summ_tf_hg, filter, N_exp >= min_exp)
summ_sub_tf_mm <- lapply(summ_tf_mm, filter, N_exp >= min_exp)


# Summarize ribo similarity and organize into a df
summ_ribo_hg <- gen_summary_df_list(sim_ribo_hg, msr_hg)
summ_ribo_mm <- gen_summary_df_list(sim_ribo_mm, msr_mm)


# Summary of null topk overlap
summ_null_hg <- gen_summary_df_list(sim_null_hg, add_nexp = FALSE)
summ_null_mm <- gen_summary_df_list(sim_null_mm, add_nexp = FALSE)



# Count of TFs whose mean similarity is greater than the maximum mean TopK 
# achieved across the nulls
# ------------------------------------------------------------------------------


calc_prop_gt_null <- function(sim_l, null_l) {
  
  stats <- intersect(names(sim_l), names(null_l))
  n <- nrow(sim_l[[1]])
  
  prop_l <- lapply(stats, function(x) {
    null_max <- max(null_l[[x]]$Mean)
    sum(sim_l[[x]]$Mean > null_max) / n
  })
  
  names(prop_l) <- stats
  return(prop_l)
}


prop_gt_null_hg <- calc_prop_gt_null(sim_l = summ_sub_tf_hg, null_l = summ_null_hg)
prop_gt_null_mm <- calc_prop_gt_null(sim_l = summ_sub_tf_mm, null_l = summ_null_mm)



# Are there TFs that are below typical null in Topk but above null in bottom k?
# ------------------------------------------------------------------------------


# Here, directly select TF who have an average Topk less than the average of the
# typical null Topk. Then intersect with TFs whose average BottomK is greater
# than the typical null BottomK. Note that 'typical' here is taking the mean of
# each null iteration's overlap, then taking the mean of these means. Find that
# the means across each null iteration don't vary much

calc_prop_divergent <- function(tf_l, null_l) {
  
  join_df <- left_join(tf_l$Topk, tf_l$Bottomk,
                       by = "Symbol",
                       suffix = c("_Topk", "_Bottomk"))
  
  null_topk_mean <- mean(null_l$Topk[, "Mean"])
  null_bottomk_mean <- mean(null_l$Bottomk[, "Mean"])
  
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



divergent_hg <- calc_prop_divergent(tf_l = summ_sub_tf_hg, 
                                    null_l = summ_null_hg)


divergent_mm <- calc_prop_divergent(tf_l = summ_sub_tf_mm, 
                                    null_l = summ_null_mm)




# Ratio of mean overlap relative to typical null overlap for top and bottomk

calc_ratio_overlap <- function(tf_l, null_l) {
  
  ratio_df <- left_join(tf_l$Topk[, c("Symbol", "Mean")],
                        tf_l$Bottomk[, c("Symbol", "Mean")],
                        by = "Symbol",
                        suffix = c("_Topk", "_Bottomk")) %>%
    mutate(
      Topk_ratio = Mean_Topk / mean(null_l$Topk$Mean),
      Btmk_ratio = Mean_Bottomk / mean(null_l$Bottomk$Mean),
      Odds_ratio = Topk_ratio / Btmk_ratio,
      LOR = log(Odds_ratio)
    )
  
  return(ratio_df)
}



ratio_hg <- calc_ratio_overlap(tf_l = summ_sub_tf_hg, null_l = summ_null_hg)

ratio_mm <- calc_ratio_overlap(tf_l = summ_sub_tf_mm, null_l = summ_null_mm)




# Examples of most similar individual experiments pairs across TFs and ribo
# ------------------------------------------------------------------------------


# Slicing the globally most similar dataset pairs across all TFs/ribo

most_sim_pairs <- function(summ_df_l, sim_l, n_tfs = 3, n_pairs = 3) {
  
  stats <- names(summ_df_l)
  
  stat_l <- lapply(stats, function(stat) {
    tfs <- as.character(slice_max(summ_df_l[[stat]], Max., n = n_tfs)$Symbol)
    lapply(sim_l[tfs], slice_max, !!sym(stat), n = n_pairs)
  })
  
  names(stat_l) <- stats
  return(stat_l)
}


best_pairs_tf_hg <- most_sim_pairs(summ_df_l = summ_sub_tf_hg, sim_l = sim_tf_hg)
best_pairs_tf_mm <- most_sim_pairs(summ_df_l = summ_sub_tf_mm, sim_l = sim_tf_mm)

best_pairs_ribo_hg <- most_sim_pairs(summ_df_l = summ_ribo_hg, sim_l = sim_ribo_hg)
best_pairs_ribo_mm <- most_sim_pairs(summ_df_l = summ_ribo_mm, sim_l = sim_ribo_mm)



# Rank order of similarity between species
# ------------------------------------------------------------------------------


# Join the stat summaries between species for each stat

ortho_rank_sim <- lapply(names(summ_sub_tf_hg), function(x) {
  
  df_hg <- summ_sub_tf_hg[[x]] %>% 
    filter(Symbol %in% pc_ortho$Symbol_hg) %>% 
    left_join(., pc_ortho, by = c("Symbol" = "Symbol_hg"))
  
  df_mm <- summ_sub_tf_mm[[x]] %>% 
    filter(Symbol %in% pc_ortho$Symbol_mm) %>% 
    left_join(., pc_ortho, by = c("Symbol" = "Symbol_mm"))
  
  df_ortho <- left_join(df_hg, df_mm, 
                        by = "ID", 
                        suffix = c("_human", "_mouse")) %>% 
    filter(!is.na(Mean_human) & !is.na(Mean_mouse) & 
           N_exp_human >= min_exp & N_exp_mouse >= min_exp)
  
  return(df_ortho)
  
})

names(ortho_rank_sim) <- names(summ_sub_tf_hg)



# Slicing the TFs that have the most consistent profiles in mouse and human

top_ortho_sim <- lapply(names(ortho_rank_sim), function(x) {
  
  ortho_rank_sim[[x]] %>% 
    mutate(Avg_sim = rowMeans(.[, c("Mean_human", "Mean_mouse")])) %>% 
    slice_max(Avg_sim, n = 10)

})
names(top_ortho_sim) <- names(ortho_rank_sim)



# Spearman cor of TF's summarized similarities between species

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



# Saving out similarity (within species) for joining with cross-species comparison

out_l <- list(Human = summ_sub_tf_hg$Topk, Mouse = summ_sub_tf_mm$Topk)
saveRDS(out_l, out_path)



# Plotting
# ------------------------------------------------------------------------------


# Relationship between mean top and bottom k across all TFs

overlap_scatter <- function(tf_l) {
  
  pdf <- left_join(tf_l$Topk, tf_l$Bottomk,
                   by = "Symbol",
                   suffix = c("_Topk", "_Bottomk"))
  
  qplot(pdf, xvar = "Mean_Topk", yvar = "Mean_Bottomk")
}


p1a <- overlap_scatter(summ_tf_hg)
p1b <- overlap_scatter(summ_tf_mm)



# Point plot of TF similarity with text overlay for top n TFs

point_plot_similarity <- function(summary_df,
                                  null_df,
                                  topn_label = 35,
                                  plot_title,
                                  ylabel) {
  
  null_line <- max(null_df[, "Mean"])
  
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
                    hjust = 1,
                    direction = "y",
                    size = 5,
                    segment.size = 0.15,
                    segment.color = "black") +
    geom_hline(yintercept = null_line, colour = "firebrick") +
    annotate("text", x = 60, y = null_line + 1, label = "Null", colour = "firebrick", size = 5) +
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



p2a <- point_plot_similarity(summ_tf_hg$Topk,
                             null_df = summ_null_hg$Topk,
                             plot_title = "Human",
                             ylabel = expr("Mean Top"[!!k]))


p2b <- point_plot_similarity(summ_tf_mm$Topk,
                             null_df = summ_null_mm$Topk,
                             plot_title = "Mouse",
                             ylabel = expr("Mean Top"[!!k]))


p2 <- plot_grid(p2a, p2b, nrow = 1)


ggsave(p2, height = 9, width = 20, device = "png", dpi = 600,
       filename = file.path(paste0(plot_dir, "mean_topk=", k, "_human_and_mouse.png")))



# Histogram of mean values for TF, null, and ribosomal

hist_mean <- function(summary_df, null_df, ribo_df, xlabel) {

  mean_df <- data.frame(
    Stat = c(summary_df[, "Mean"],
             null_df[, "Mean"],
             ribo_df[, "Mean"]),
    Group = c(rep("TR", nrow(summary_df)),
              rep("Null", nrow(null_df)),
              rep("Ribosomal", nrow(ribo_df)))
  )

  ggplot(mean_df, aes(x = Stat)) +
    facet_wrap(~Group, nrow = 3, scales = "free_y") +
    geom_histogram(bins = 100, fill = "slategrey", colour = "slategrey") +
    xlab(xlabel) +
    ylab("Count") +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 20),
          strip.text = element_text(size = 20),
          strip.background = element_rect(color = "black", fill = c("lightgrey")),
          plot.margin = margin(c(10, 20, 10, 10)))
}



p3a <- hist_mean(summary_df = summ_tf_hg$Topk,
                 null_df = summ_null_hg$Topk,
                 ribo_df = summ_ribo_hg$Topk,
                 xlabel = expr("Mean Top"[!!k]))


p3b <- hist_mean(summary_df = summ_tf_mm$Topk,
                 null_df = summ_null_mm$Topk,
                 ribo_df = summ_ribo_mm$Topk,
                 xlabel = expr("Mean Top"[!!k]))


# Mean values as a strip plot. We ended up using in the paper, insetting this
# within the point similarity plot in illustrator

strip_mean <- function(summary_df, null_df, ribo_df, ylabel) {
  
  mean_df <- data.frame(
    Stat = c(summary_df[, "Mean"],
             null_df[, "Mean"],
             ribo_df[, "Mean"]),
    Group = c(rep("TR", nrow(summary_df)),
              rep("Null", nrow(null_df)),
              rep("Ribo", nrow(ribo_df)))
  )
  
  ggplot(mean_df, aes(y = Stat, x = Group)) +
    geom_jitter(shape = 21, alpha = 0.4, width = 0.3, fill = "slategrey") + 
    ylab(expr("Mean Top"[!!k])) +
    theme_classic() +
    theme(text = element_text(size = 25),
          axis.title.x = element_blank())
}
  
  
p4a <- strip_mean(summary_df = summ_tf_hg$Topk,
                  null_df = summ_null_hg$Topk,
                  ribo_df = summ_ribo_hg$Topk,
                  ylabel = expr("Mean Top"[!!k])) 


p4b <- strip_mean(summary_df = summ_tf_mm$Topk,
                  null_df = summ_null_mm$Topk,
                  ribo_df = summ_ribo_mm$Topk,
                  ylabel = expr("Mean Top"[!!k]))


ggsave(p4a, height = 6, width = 4, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "strip_mean_topk=", k, "_hg.png")))

ggsave(p4b, height = 6, width = 4, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "strip_mean_topk=", k, "_mm.png")))



# Combining point plot and hist of mean
# Flip mouse to mirror human

p5a <- plot_grid(p3a, p2a, rel_widths = c(0.5, 1))
p5b <- plot_grid(p2b, p3b, rel_widths = c(1, 0.5))
p5 <- plot_grid(p5a, p5b, nrow = 1)

ggsave(p5, height = 10, width = 22, device = "png", dpi = 600,
       filename = file.path(paste0(plot_dir, "mean_topk=", k, "_human_and_mouse_hist.png")))


ggsave(p5, height = 9, width = 22, device = "png", dpi = 600,
       filename = file.path(paste0(plot_dir, "mean_topk=", k, "_human_and_mouse_hist.png")))


# Again for bottomk

p6a <- point_plot_similarity(summ_tf_hg$Bottomk,
                             null_df = summ_null_hg$Bottomk,
                             plot_title = "Human",
                             ylabel = expr("Mean Bottom"[!!k]))


p6b <- point_plot_similarity(summ_tf_mm$Bottomk,
                             null_df = summ_null_mm$Bottomk,
                             plot_title = "Mouse",
                             ylabel = expr("Mean Bottom"[!!k]))


p7a <- hist_mean(summary_df = summ_tf_hg$Bottomk,
                 null_df = summ_null_hg$Bottomk,
                 ribo_df = summ_ribo_hg$Bottomk,
                 xlabel = expr("Mean Bottom"[!!k]))


p7b <- hist_mean(summary_df = summ_tf_mm$Bottomk,
                 null_df = summ_null_mm$Bottomk,
                 ribo_df = summ_ribo_mm$Bottomk,
                 xlabel = expr("Mean Bottom"[!!k]))


# Combining point plot and hist of mean bottomk
p8a <- plot_grid(p7a, p6a, rel_widths = c(0.5, 1))
p8b <- plot_grid(p6b, p7b, rel_widths = c(1, 0.5))
p8 <- plot_grid(p8a, p8b, nrow = 1)

ggsave(p8, height = 9, width = 22, device = "png", dpi = 600,
       filename = file.path(paste0(plot_dir, "mean_bottomk=", k, "_human_and_mouse_hist.png")))


# Density plot of topk intersect for select genes + null

density_topk <- function(plot_df, k, plot_title) {
  
  cols <- c("black", "#1b9e77", "#d95f02", "#7570b3", "lightgrey")
  
  ggplot(plot_df, aes(x = Topk, fill = Group)) +
    geom_density(alpha = 0.6) +
    theme_classic() +
    ylab("Density") +
    xlab(expr("Top"[!!k])) +
    ggtitle(plot_title) +
    scale_fill_manual(values = cols) +
    theme(
      axis.text = element_text(size = 30),
      axis.title = element_text(size = 30),
      plot.title = element_text(size = 30),
      legend.position = c(0.75, 0.90),
      legend.text = element_text(size = 30),
      legend.title = element_blank(),
      plot.margin = margin(5, 20, 5, 5))
}


# In paper compare NEUROD6 and PAX6, using E2F8 and a ribo gene as reference

tf1_hg <- "NEUROD6"
tf2_hg <- "PAX6"
tf3_hg <- "E2F8"
ribo1_hg <- "RPL32"

tf_l_hg <- sim_tf_hg[c(tf1_hg, tf2_hg, tf3_hg)]
ribo_hg <- sim_ribo_hg[[ribo1_hg]]

# Pull a random null iteration
set.seed(154)
rep_null_hg <- sim_null_hg[[sample(1:length(sim_null_hg), 1)]]


pdf9 <- data.frame(
  Group = c(
    rep(ribo1_hg, nrow(ribo_hg)),
    rep(tf1_hg, nrow(tf_l_hg[[1]])),
    rep(tf2_hg, nrow(tf_l_hg[[2]])),
    rep(tf3_hg, nrow(tf_l_hg[[3]])),
    rep("Null", nrow(rep_null_hg))
  ),
  Topk = c(
    ribo_hg$Topk,
    tf_l_hg[[1]]$Topk,
    tf_l_hg[[2]]$Topk,
    tf_l_hg[[3]]$Topk,
    rep_null_hg$Topk
  )
)


pdf9$Group <- factor(pdf9$Group, levels = unique(pdf9$Group))

p9 <- density_topk(pdf9, k = k, plot_title = "Human")

ggsave(p9, height = 8, width = 11, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("density_topk=", k, "_human.png")))



# Violin plot with TF labels for divergent

plot_divergent <- function(divergent_l, null_df, k, plot_title, nlabel = 10) {
  
  plot_df <- data.frame(
    Mean_bottomk = c(divergent_l$Below_null$Mean_Bottomk, 
                     null_df$Mean),
    Symbol = c(divergent_l$Below_null$Symbol, 
               rep("Null", length(null_df$Mean))),
    Group = c(rep("TR", length(divergent_l$Below_null$Mean_Bottomk)),
              rep("Null", length(null_df$Mean)))
  )
  
  # label TFs whose mean Bottomk, but not Topk, exceeds the respectie nulls
  div_symbol <- slice_max(divergent_l$Divergent, Mean_Bottomk, n = nlabel)$Symbol
  plot_df$Label <- plot_df$Symbol %in% div_symbol
  
  
  ggplot(plot_df, aes(x = Group, y = Mean_bottomk)) +
    
    geom_violin(width = 0.3, fill = "slategrey") +
    geom_boxplot(width = 0.1, fill = "white", outlier.fill = "white") +
    geom_point(
      data = filter(plot_df, Label), shape = 21, alpha = 0.6, size = 3.1
      ) +
    geom_text_repel(
      data = filter(plot_df, Label),
      aes(x = Group, y = Mean_bottomk, label = Symbol, fontface = "italic"),
      max.overlaps = 30,
      force = 0.5,
      # nudge_x = 500,
      hjust = 0,
      direction = "y",
      size = 5,
      segment.size = 0.1, 
      segment.color = "grey50") +
    ggtitle(plot_title) +
    ylab(expr("Mean Bottom"[!!k])) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 30),
          plot.margin = margin(c(10, 20, 10, 10)))
}


p10a <- plot_divergent(divergent_l = divergent_hg,
                       null_df = summ_null_hg$Bottomk,
                       k = k,
                       plot_title = "Human")


p10b <- plot_divergent(divergent_l = divergent_mm,
                       null_df = summ_null_mm$Bottomk,
                       k = k,
                       plot_title = "Mouse")


p10 <- plot_grid(p10a, p10b, nrow = 1)


ggsave(p10, height = 10, width = 20, device = "png", dpi = 600,
       filename = file.path(paste0(plot_dir, "divergent_similarity_topk=", k, ".png")))



# Scatter plot of ratio of topk to null and bottomk to null. Colour TFs that 
# have a topk ratio < 1 (mean topk is less than the typical null topk) and a
# bottomk ratio > 1 (mean bottomk is greater than the typical null bottomk) 

overlap_ratio_scatter <- function(ratio_df, title) {
  
  pdf <- mutate(ratio_df, Group = Topk_ratio < 1 & Btmk_ratio > 1)
  
  ggplot(pdf) +
    geom_point(aes(x = Topk_ratio, y = Btmk_ratio), 
               data = filter(pdf, !Group),
               shape = 21, size = 1, alpha = 0.2) +
    geom_point(aes(x = Topk_ratio, y = Btmk_ratio), 
               data = filter(pdf, Group),
               shape = 21, size = 3, colour = "firebrick") +
    geom_hline(yintercept = 1) +
    geom_vline(xintercept = 1) +
    ggtitle(title) +
    xlab(paste0("Top K=", k, " ratio observed to null")) +
    ylab(paste0("Bottom K=", k, " ratio observed to null")) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20))
  
}


p11a <- overlap_ratio_scatter(ratio_hg, title = "Human")
p11b <- overlap_ratio_scatter(ratio_mm, title = "Mouse")
p11 <- plot_grid(p11a, p11b, nrow = 1)

ggsave(p11, height = 7, width = 14, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "ratio_similarity_topk=", k, ".png")))



# Scatterplot of mean top/bottom overlap of ortho TFs (each TF's similarity
# is calculated within species)

p12a <- qplot(ortho_rank_sim$Topk, xvar = "Mean_human", yvar = "Mean_mouse") +
  geom_smooth(method = "lm", colour = "red") +
  xlab(expr("Human Mean Top"[!!k])) +
  ylab(expr("Mouse Mean Top"[!!k])) +
  ggtitle("Consistency of positive profiles")

p12b <- qplot(ortho_rank_sim$Bottomk, xvar = "Mean_human", yvar = "Mean_mouse") +
  geom_smooth(method = "lm", colour = "red") +
  xlab(expr("Human Mean Bottom"[!!k])) +
  ylab(expr("Mouse Mean Bottom"[!!k])) +
  ggtitle("Consistency of negative profiles")

p12 <- plot_grid(p12a, p12b, nrow = 1)


ggsave(p12, height = 6, width = 12, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("ortho_cor_consistencty_k=", k, ".png")))


# Boxplot of average similarities by TF family

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


p13a <- plot_family_similarity(summ_sub_tf_hg$Topk)
p13b <- plot_family_similarity(summ_sub_tf_mm$Topk)
