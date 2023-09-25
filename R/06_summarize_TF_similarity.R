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

k <- 1000 

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

# Query/subject all rank matrices
# allrank_dir <- "/space/scratch/amorin/R_objects/04-07-2023/"
# allrank_l <- lapply(list.files(allrank_dir, full.names = TRUE), function(x) readRDS(x))
# names(allrank_l) <- str_replace(list.files(allrank_dir), "\\.RDS", "")


# Functions
# ------------------------------------------------------------------------------


# Returns a dataframe of summary stats for topk for each TF in sim_l

get_summary_df <- function(sim_l, msr_mat) {
  
  lapply(sim_l, function(x) summary(x$Topk)) %>% 
    do.call(rbind, .) %>%
    as.data.frame() %>% 
    rownames_to_column(var = "Symbol") %>% 
    mutate(N_exp = rowSums(msr_mat[names(sim_l), ])) %>%
    arrange(Median) %>% 
    mutate(Symbol = factor(Symbol, levels = unique(Symbol)))
  
}


# Organizing summary of topk intersects
# ------------------------------------------------------------------------------


# Summarize each TF's similarity and organize into a df
summ_tf_hg <- get_summary_df(sim_tf_hg, msr_hg)
summ_tf_mm <- get_summary_df(sim_tf_mm, msr_mm)


# Require a minimum count of measured experiments for global trends
summ_sub_tf_hg <- filter(summ_tf_hg, N_exp >= 5)
summ_sub_tf_mm <- filter(summ_tf_mm, N_exp >= 5)


# Summarize ribo similarity and organize into a df
summ_ribo_hg <- get_summary_df(sim_ribo_hg, msr_hg)
summ_ribo_mm <- get_summary_df(sim_ribo_mm, msr_mm)


# Summary of null topk overlap
summ_null_hg <- summary(unlist(lapply(null_topk_hg, function(x) median(x$Topk))))
summ_null_mm <- summary(unlist(lapply(null_topk_mm, function(x) median(x$Topk))))


# Relationship between median topk intersect and count of measured experiments
topk_na_cor_hg <- suppressWarnings(cor.test(summ_tf_hg$Median, summ_tf_hg$N_exp, method = "spearman"))
topk_na_cor_mm <- suppressWarnings(cor.test(summ_tf_mm$Median, summ_tf_mm$N_exp, method = "spearman"))


# Count of TFs whose median topk is greater than the null
topk_gt_null_hg <- sum(summ_sub_tf_hg$Median > summ_null_hg["Median"]) / nrow(summ_sub_tf_hg)
topk_gt_null_mm <- sum(summ_sub_tf_mm$Median > summ_null_mm["Median"]) / nrow(summ_sub_tf_mm)


# TFs whose median topk is greater than the median topk for ribosomal
topk_gt_ribo_hg <- filter(summ_sub_tf_hg, Median > median(summ_ribo_hg[["Median"]]))
topk_gt_ribo_mm <- filter(summ_sub_tf_mm, Median > median(summ_ribo_mm[["Median"]]))



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
    mutate(Group = Symbol %in% slice_max(summary_df, Median, n = topn_label)$Symbol)
  
  ggplot(summary_df, aes(y = Median, x = Symbol)) +
    geom_point() +
    geom_text_repel(data = filter(summary_df, Group),
                    aes(x = Symbol, y = Median, label = Symbol, fontface = "italic"),
                    max.overlaps = 30, 
                    force = 0.5,
                    nudge_x = -0.25,
                    hjust = 0,
                    size = 5,
                    segment.size = 0.2) +
    geom_hline(yintercept = summ_nullary["Median"], colour = "firebrick") +
    geom_hline(yintercept = summ_nullary["1st Qu."], colour = "grey85") +
    geom_hline(yintercept = summ_nullary["3rd Qu."], colour = "grey85") +
    ggtitle(title) +
    ylab(paste0("Median Top k (k=", k, ")")) +
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

example_tf_hg <- "RUNX1"
example_tf_mm <- "Runx1"

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


# Heatmap of query/subject topk all ranks


p6 <- pheatmap(allrank_mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_col = "black",
         color = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
         display_numbers = TRUE,
         number_format = "%1.0f",
         number_color = "black",
         na_col = "black",
         fontsize = 22,
         cellwidth = 50,
         cellheight = 50,
         height = 18,
         width = 18)



# Histogram of query/subject all ranks: Pile up at the left of the histogram 
# (lower/better ranks) suggests a signal, as random ranking would be uniform


p7 <- mat_to_df(allrank_mat, value_name = "Rank") %>% 
  ggplot(., aes(x = Rank)) +
  geom_histogram(bins = 30) +
  ylab("Count") +
  ggtitle(gene_hg) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


# hist(replicate(length(allrank_mat), sample(1:length(pc_hg$Symbol), 1)), breaks = 20)


# Histogram of the example query/subject topk rank, with overlay of example gene


p8 <- data.frame(Topk = topk_rank$Topk) %>% 
  ggplot(aes(x = Topk)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = topk_rank$Topk[gene_hg], col = "red") +
  ggtitle(gene_hg) +
  ylab("Count of genes") +
  xlab("Topk intersect") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




# For inspecting all rank query/subject topk
# ------------------------------------------------------------------------------


# gene_hg <- "ASCL1"
# 
# 
# # Keep only measured genes and set diag to NA to remove self-ranking
# 
# keep_exp <- which(msr_hg[gene_hg, ] == 1)
# allrank_mat <- allrank_l[[gene_hg]][keep_exp, keep_exp, drop = FALSE]
# diag(allrank_mat) <- NA
# 
# 
# # Inspecting top pairs: by topk magnitude, and by rank. The highest value/most
# # similar pair across all experiment pairs by magnitude may not actually be the
# # top ranked gene within the subject dataset
# 
# sim_df <- sim_tf_hg[[gene_hg]]$Sim_df
# best_value_id <- slice_max(sim_df, Topk, n = 1)
# best_rank_ix <- which(allrank_mat == min(allrank_mat, na.rm = TRUE), arr.ind = TRUE)
# 
# query_vec <- load_agg_mat_list(ids = best_value_id$Row,
#                                genes = pc_hg$Symbol,
#                                sub_genes = gene_hg)[[1]][, gene_hg]
# 
# subject_mat <- load_agg_mat_list(ids = best_value_id$Col, genes = pc_hg$Symbol)[[1]]
# 
# 
# topk_rank <- query_gene_rank_topk(
#   query_vec = query_vec,
#   subject_mat = subject_mat,
#   gene = gene_hg,
#   ncores = ncore)
