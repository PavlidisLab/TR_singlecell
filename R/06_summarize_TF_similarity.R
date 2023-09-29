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
sim_tf_hg <- readRDS(sim_tf_hg_path)
sim_tf_mm <- readRDS(sim_tf_mm_path)

# L/S ribo topk overlap
sim_ribo_hg <- readRDS(sim_ribo_hg_path)
sim_ribo_mm <- readRDS(sim_ribo_mm_path)

# Null topk overlap
null_topk_hg <- readRDS(null_topk_hg_path)
null_topk_mm <- readRDS(null_topk_mm_path)



# Functions
# ------------------------------------------------------------------------------


# Returns a list of dataframes of summary stats for each TF's similarity

get_summary_df <- function(sim_l, msr_mat) {
  
  stats <- c("Scor", "Topk", "Bottomk", "Jaccard")
  
  df_l <- lapply(stats, function(stat) {
    
    lapply(sim_l, function(x) summary(x[[stat]])) %>% 
      do.call(rbind, .) %>%
      as.data.frame() %>% 
      rownames_to_column(var = "Symbol") %>% 
      mutate(N_exp = rowSums(msr_mat[names(sim_l), ])) %>%
      arrange(desc(Median)) %>% 
      mutate(Symbol = factor(Symbol, levels = unique(Symbol)))
  })
  
  names(df_l) <- stats
  return(df_l)
}


# Organizing summary of topk intersects
# ------------------------------------------------------------------------------


# Summarize each TF's similarity and organize into a df
summ_tf_hg <- get_summary_df(sim_tf_hg, msr_hg)
summ_tf_mm <- get_summary_df(sim_tf_mm, msr_mm)


# Require a minimum count of measured experiments for global trends
summ_sub_tf_hg <- lapply(summ_tf_hg, filter, N_exp >= 5)
summ_sub_tf_mm <- lapply(summ_tf_mm, filter, N_exp >= 5)


# Summarize ribo similarity and organize into a df
summ_ribo_hg <- get_summary_df(sim_ribo_hg, msr_hg)
summ_ribo_mm <- get_summary_df(sim_ribo_mm, msr_mm)


# Summary of null topk overlap
summ_null_hg <- do.call(rbind, lapply(null_topk_hg, function(x) summary(x$Topk)))
summ_null_mm <- do.call(rbind, lapply(null_topk_mm, function(x) summary(x$Topk)))


# Relationship between median topk intersect and count of measured experiments
topk_na_cor_hg <- suppressWarnings(cor.test(summ_tf_hg$Topk$Median, summ_tf_hg$Topk$N_exp, method = "spearman"))
topk_na_cor_mm <- suppressWarnings(cor.test(summ_tf_mm$Topk$Median, summ_tf_mm$Topk$N_exp, method = "spearman"))


# Count of TFs whose median topk is greater than the null
topk_gt_null_hg <- sum(summ_sub_tf_hg$Topk$Median > summ_null_hg["Median"]) / nrow(summ_sub_tf_hg)
topk_gt_null_mm <- sum(summ_sub_tf_mm$Median > summ_null_mm["Median"]) / nrow(summ_sub_tf_mm)


# TFs whose median topk is greater than the median topk for ribosomal
topk_gt_ribo_hg <- filter(summ_sub_tf_hg$Topk, Median > median(summ_ribo_hg$Topk[["Median"]]))
topk_gt_ribo_mm <- filter(summ_sub_tf_mm$Topk, Median > median(summ_ribo_mm$Topk[["Median"]]))



# Examples of best/worst pairs for a given TF
# ------------------------------------------------------------------------------



# ASCL1 "GSE156793" "GSE137537" 0 count overlap Scor = 0.3130 is this because low 
# non-tied k? Looks it - "GSE156793" k=148, although still none of these overlap
# "Elmentaite2021"  "Han2022" Scor 0.58 and Topk 2

# FOXF1 "Stewart2019" "Young2018" suspicious k=998


tf <- "FOXF1"
mat <- gene_vec_to_mat(agg_tf_hg, tf)
mat <- subset_to_measured(mat, msr_hg, tf)
id1 <- "Stewart2019"  # "Elmentaite2021"
id2 <- "Young2018"   # "Han2022"

vec1 <- mat[rownames(mat) != tf, id1]
vec2 <- mat[rownames(mat) != tf, id2]

vec1_sort <- sort(vec1, decreasing = TRUE)
vec2_sort <- sort(vec2, decreasing = TRUE)

check_k(vec1_sort, 1000)
check_k(vec2_sort, 1000)

topk_intersect(topk_sort(vec1, k = 1000), topk_sort(vec2, k = 1000))

plot(mat[, id1], mat[, id2])
cor(mat[, id1], mat[, id2], method = "spearman")
view(data.frame(Symbol = names(vec1), Vec1 = vec1, Vec2 = vec2))




# Example of TFs with good measurement coverage but low similarity. 
# ------------------------------------------------------------------------------


arrange(summ_tf_hg, Median, desc(N_exp)) %>% head
arrange(summ_tf_mm, Median, desc(N_exp)) %>% head





# Look at similarity by TF family
# TODO: repetition of factor ordering after join
# ------------------------------------------------------------------------------


summ_tf_hg <- summ_tf_hg %>% 
  left_join(., tfs_hg[, c("Symbol", "Family")], by = "Symbol") %>% 
  arrange(Median) %>% 
  mutate(Symbol = factor(Symbol, levels = unique(Symbol)))


family_hg <- summ_tf_hg %>% 
  group_by(Family) %>% 
  summarise(
    Med_topk = median(Median), 
    N_TFs = n(),
    Med_exp_msrd = median(N_exp)
    ) %>% 
  arrange(desc(Med_topk))


boxplot(summ_tf_hg$Med ~ summ_tf_hg$Family)





# Plotting
# ------------------------------------------------------------------------------


gene_hg <- "ASCL1"

# Relationship between similarity stats

p1 <- ggplot(sim_tf_hg[[gene_hg]], aes(x = Scor, y = Topk)) +
  geom_point(shape = 21, size = 2.1) +
  ggtitle(gene_hg) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



# Scatterplot of TF median topk intersect with text overlay for top n TFs


scatter_topk_median <- function(summary_df, 
                                summ_nullary, 
                                topn_label = 30,
                                title) {
  
  summary_df <- summary_df %>%  
    # mutate(Group = Symbol %in% slice_max(summary_df, Median, n = topn_label)$Symbol)
    mutate(Group = Symbol %in% slice_max(summary_df, Mean, n = topn_label)$Symbol)
  
  # ggplot(summary_df, aes(y = Median, x = Symbol)) +
  ggplot(summary_df, aes(y = Mean, x = Symbol)) +
    geom_point() +
    geom_text_repel(data = filter(summary_df, Group),
                    # aes(x = Symbol, y = Median, label = Symbol, fontface = "italic"),
                    aes(x = Symbol, y = Mean, label = Symbol, fontface = "italic"),
                    max.overlaps = 30, 
                    force = 0.5,
                    nudge_x = -0.25,
                    hjust = 0,
                    size = 5,
                    segment.size = 0.2) +
    # geom_hline(yintercept = summ_nullary["Median"], colour = "firebrick") +
    geom_hline(yintercept = summ_nullary["Mean"], colour = "firebrick") +
    geom_hline(yintercept = summ_nullary["1st Qu."], colour = "grey85") +
    geom_hline(yintercept = summ_nullary["3rd Qu."], colour = "grey85") +
    ggtitle(title) +
    # ylab(paste0("Median Top k (k=", k, ")")) +
    ylab(paste0("Mean Top k (k=", k, ")")) +
    expand_limits(x = nrow(summ_tf_hg) + 50) +  # prevent point cut off
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 10, 10, 10)))
}



p2a <- scatter_topk_median(summ_tf_hg, summ_null_hg, title = "Human")
p2b <- scatter_topk_median(summ_tf_mm, summ_null_mm, title = "Mouse")
p2 <- plot_grid(p2a, p2b)

ggsave(p2, height = 9, width = 18, device = "png", dpi = 300,
       filename = file.path(plot_dir, "scatter_topk_median.png"))



# Boxplot (show 1st and 3rd IQR) of TF median topk intersect



boxplot_topk_median <- function(summary_df, summ_nullary) {
  
  ggplot(summary_df, aes(x = Symbol)) +
    geom_boxplot(
      # aes(ymin = Min., lower = `1st Qu.`, middle = Median, upper = `3rd Qu.`, ymax = Max.),
      aes(ymin = `1st Qu.`, lower = `1st Qu.`, middle = Median, upper = `3rd Qu.`, ymax = `3rd Qu.`),
      stat = "identity") +
    geom_hline(yintercept = summ_nullary["Median"], colour = "firebrick", linewidth = 1.6) +
    geom_hline(yintercept = summ_nullary["1st Qu."], colour = "grey85") +
    geom_hline(yintercept = summ_nullary["3rd Qu."], colour = "grey40") +
    ylab(paste0("Median Top k (k=", k, ")")) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}


p3a <- boxplot_topk_median(summ_tf_hg, summ_null_hg)
p3b <- boxplot_topk_median(summ_tf_mm, summ_null_mm)


# Density plot of topk intersect for select genes + null


density_topk <- function(plot_df) {
  
  ggplot(plot_df, aes(x = Topk, fill = Group)) +
    geom_density(alpha = 0.6) +
    theme_classic() +
    ylab("Density") +
    xlab(paste0("Topk intersect (k=", k, ")")) +
    scale_fill_manual(values = c("#1b9e77", "#d95f02", "lightgrey", "#7570b3")) +
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

max_tf_hg <- as.character(slice_max(summ_tf_hg, Median)$Symbol)
max_tf_mm <- as.character(slice_max(summ_tf_mm, Median)$Symbol)

example_ribo_hg <- "RPL7A"
example_ribo_mm <- "Rpl3"

max_tf_df_hg <- sim_tf_hg[[max_tf_hg]]
max_tf_df_mm <- sim_tf_mm[[max_tf_mm]]

tf_df_hg <- sim_tf_hg[[example_tf_hg]]
tf_df_mm <- sim_tf_mm[[example_tf_mm]]

ribo_df_hg <- sim_ribo_hg[[example_ribo_hg]]
ribo_df_mm <- sim_ribo_mm[[example_ribo_mm]]


set.seed(154)
rep_null_hg <- null_topk_hg[[sample(1:length(null_topk_hg), 1)]]
rep_null_mm <- null_topk_mm[[sample(1:length(null_topk_mm), 1)]]



plot_df_hg <- data.frame(
  Group = c(rep(max_tf_hg, nrow(max_tf_df_hg)), 
            rep(example_tf_hg, nrow(tf_df_hg)), 
            rep(example_ribo_hg, nrow(ribo_df_hg)), 
            rep("Null", nrow(rep_null_hg))),
  Topk = c(max_tf_df_hg$Topk, tf_df_hg$Topk, ribo_df_hg$Topk, rep_null_hg$Topk)
)



plot_df_mm <- data.frame(
  Group = c(rep(max_tf_mm, nrow(max_tf_df_mm)), 
            rep(example_tf_mm, nrow(tf_df_mm)), 
            rep(example_ribo_mm, nrow(ribo_df_mm)), 
            rep("Null", nrow(rep_null_mm))),
  Topk = c(max_tf_df_mm$Topk, tf_df_mm$Topk, ribo_df_mm$Topk, rep_null_mm$Topk)
)



p4a <- density_topk(plot_df_hg)
p4b <- density_topk(plot_df_mm)
p4 <- plot_grid(p4a, p4b)

ggsave(p4, height = 9, width = 18, device = "png", dpi = 300,
       filename = file.path(plot_dir, "density_topk_example.png"))



boxplot(log10(plot_df_hg$Topk+1) ~ plot_df_hg$Group)
boxplot(log10(plot_df_mm$Topk+1) ~ plot_df_mm$Group)



# Look at distribution of expected values for TF versus null
plot(density(summ_sub_tf_hg$Topk[, "Mean"]), ylim = c(0, 4))
lines(density(summ_null_hg[, "Mean"]), col = "red")




# Relationship between topk median and count of experiments

# plot(summ_tf_hg$N_exp, y = summ_tf_hg$Median)
# plot(summ_tf_mm$N_exp, y = summ_tf_mm$Median)

n_break <- length(ids_mm)

p5a <- summ_tf_hg %>% 
  mutate(Group_nexp = cut(N_exp, n_break, include.lowest = TRUE)) %>% 
  ggplot(aes(x = Group_nexp, y = Median)) +
  geom_boxplot() +
  theme_classic()


p5b <- summ_tf_mm %>% 
  mutate(Group_nexp = cut(N_exp, n_break, include.lowest = TRUE)) %>% 
  ggplot(aes(x = Group_nexp, y = Median)) +
  geom_boxplot() +
  theme_classic()


