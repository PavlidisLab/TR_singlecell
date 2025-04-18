## Examining how similar ortho TF aggregate profiles are between species.
## This explicitly considers only ortho genes, with each ranking generated only
## from data from the respective species.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggrepel)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 200
min_exp <- 5

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Saved list RDS of the ortho ranks
rank_tf_ortho <- readRDS(rank_tf_ortho_path)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Within species cross-dataset similarity lists
topk_l <- readRDS(paste0(output_dir, "human_mouse_topk=", k, "_similarity_df.RDS"))



# Functions
# ------------------------------------------------------------------------------


# Return a df of the spearman cor between human and mouse RSR rankings

calc_ortho_cor <- function(ortho_rank_l, ncores = 1) {
  
  symbols <- names(ortho_rank_l)
  
  cor_l <- mclapply(ortho_rank_l, function(x) {
    cor(x$Avg_aggr_coexpr_hg, x$Avg_aggr_coexpr_mm, 
        method = "spearman",
        use = "pairwise.complete.obs")
  }, mc.cores = ncores)
  
  cor_df <- data.frame(Symbol = symbols, Cor = unlist(cor_l))
  
  return(cor_df)
}


# Flip sign of average RSR (rank metric) to prioritize bottom of ranks

reverse_ranks <- function(df) {
  mutate(df, 
         Avg_aggr_coexpr_hg = -Avg_aggr_coexpr_hg, 
         Avg_aggr_coexpr_mm = -Avg_aggr_coexpr_mm)
}


# Generate a gene x gene matrix whose elements represent the count of top K
# elements shared between the ortho ranks.
# Rows are human in mouse, cols are mouse in human
# NOTE: This is taking the top k AFTER filtering for ortho, so not necessarily
# within the top k ranks of the unfiltered (non-ortho) data

calc_overlap_mat <- function(ortho_l, k, reverse = FALSE, ncores = 1) {
  
  # Reverse rankings to prioritize bottom of ranks (neg coexpression)?
  if (reverse) ortho_l <- lapply(ortho_l, reverse_ranks)
  
  # First get list of the top K elements for human and mouse for each gene rank
  genes <- names(ortho_l)
  
  # Get list of topK genes (char vec) for each gene and species
  topk_l <- mclapply(ortho_l, function(x) {
    list(
      Human = slice_max(x, Avg_aggr_coexpr_hg, n = k)$Symbol_hg,
      Mouse = slice_max(x, Avg_aggr_coexpr_mm, n = k)$Symbol_hg
    )
  }, mc.cores = ncores)
  names(topk_l) <- genes
  
  # Generate overlap matrix
  overlap_mat <- matrix(0, nrow = length(topk_l), ncol = length(topk_l))
  rownames(overlap_mat) <- colnames(overlap_mat) <- genes
  
  for (i in 1:nrow(overlap_mat)) {
    for (j in 1:ncol(overlap_mat)) {
      overlap_mat[i, j] <- length_intersect(topk_l[[i]]$Human, topk_l[[j]]$Mouse)
    }
  }
  
  return(overlap_mat)
}


# Given gene x gene overlap_mat, return a data.frame containing the count of
# overlap for the orthologous match, as well as the quantile of the match
# relative to all other TFs for each species

gen_quantile_df <- function(overlap_mat, ncores = 1) {
  
  genes <- intersect(rownames(overlap_mat), colnames(overlap_mat))
  
  quant_l <- mclapply(genes, function(x) {
    
    count_overlap <- overlap_mat[x, x]
    hg_in_mm <- overlap_mat[x, setdiff(colnames(overlap_mat), x)]
    mm_in_hg <- overlap_mat[setdiff(rownames(overlap_mat), x), x]
    
    data.frame(
      Symbol = x,
      Count = count_overlap,
      Quant_hg_in_mm = ecdf(hg_in_mm)(count_overlap),
      Quant_mm_in_hg = ecdf(mm_in_hg)(count_overlap))
    
  }, mc.cores = ncores)
  
  quant_df <- do.call(rbind, quant_l)
  quant_df$Quant_ortho <- rowMeans(quant_df[, c("Quant_hg_in_mm", "Quant_mm_in_hg")])
  
  return(quant_df)
}



# Calculate null overlap by shuffling pairs of TFs

calc_null_overlap <- function(ortho_l, k, reverse = FALSE, ncores = 1) {
  
  # Reverse rankings to prioritize bottom of ranks (neg coexpression)?
  if (reverse) ortho_l <- lapply(ortho_l, reverse_ranks)
  
  # Helper to get size of overlap between mouse and human for sampled TFs
  calc_topk <- function(k, sample_tfs) {
    length_intersect(
      slice_max(ortho_l[[sample_tfs[1]]], Avg_aggr_coexpr_hg, n = k)$Symbol_hg,
      slice_max(ortho_l[[sample_tfs[2]]], Avg_aggr_coexpr_mm, n = k)$Symbol_hg
    )
  }
  
  # Iteratively sample TFs and calculate overlap between mouse and human
  null_topk <- mclapply(1:1000, function(x) {
    sample_tfs <- sample(names(ortho_l), 2, replace = FALSE)
    calc_topk(k, sample_tfs)
  }, mc.cores = ncores)
  
  return(unlist(null_topk))
}



# Calculate null Spearman's cor by shuffling pairs of TFs

calc_null_cor <- function(ortho_l, ncores = 1) {
  
  # Keep same gene order
  genes <- ortho_l[[1]]$Symbol_hg
  
  null_cor <- mclapply(1:1000, function(x) {
    
    # Need to remove sampled TFs (which are self NAs)
    sample_tfs <- sample(names(ortho_l), 2, replace = FALSE)
    genes <- setdiff(genes, sample_tfs)
    
    tf1 <- ortho_l[[sample_tfs[1]]] %>% 
      filter(Symbol_hg %in% genes) %>% 
      arrange(match(Symbol_hg, genes))
    
    tf2 <- ortho_l[[sample_tfs[2]]] %>% 
      filter(Symbol_hg %in% genes) %>% 
      arrange(match(Symbol_hg, genes))
    
    cor(tf1$Rank_aggr_coexpr_hg, tf2$Rank_aggr_coexpr_mm, 
        method = "spearman",
        use = "pairwise.complete.obs")
    
  }, mc.cores = ncore)
  
  return(unlist(null_cor))
}



# Only keeping ortho TFs measured in a minimum amount of experiments Note that
# I also tried filtering within each TF for the same minimal amount of 
# co-measurement and it made a negligible difference, so skipping.
# ------------------------------------------------------------------------------


# All ortho protein coding genes
pc_ortho <- filter(pc_ortho,
                   Symbol_hg %in% pc_hg$Symbol &
                     Symbol_mm %in% pc_mm$Symbol)

# Ortho TFs where a ranking was generated
tf_ortho <- filter(pc_ortho, Symbol_hg %in% names(rank_tf_ortho))

# Measurement counts for each TF
msr_df <- data.frame(
  Symbol = tf_ortho$Symbol_hg,
  N_hg = rowSums(msr_hg[tf_ortho$Symbol_hg, ]),
  N_mm = rowSums(msr_mm[tf_ortho$Symbol_mm, ]))


# Remove TFs that were rarely measured (problems of small n)
rm_tfs <- filter(msr_df, N_hg < min_exp | N_mm < min_exp)
rank_tf_ortho <- rank_tf_ortho[setdiff(names(rank_tf_ortho), rm_tfs$Symbol)]
tf_ortho <- filter(tf_ortho, Symbol_hg %!in% rm_tfs$Symbol)



# Generate the similarity objects. Human symbols used for representing ortho TFs 
# ------------------------------------------------------------------------------


# Get the Spearman's correlation of rankings between every ortho TF
cor_df <- calc_ortho_cor(rank_tf_ortho, ncores = ncore)


# Top and bottom K matrices 
topk_mat <- calc_overlap_mat(rank_tf_ortho, k = k, ncores = ncore)
bottomk_mat <- calc_overlap_mat(rank_tf_ortho, k = k, reverse = TRUE, ncores = ncore)


# Nulls
set.seed(5)
null_topk <- calc_null_overlap(rank_tf_ortho, k = k, ncores = ncore)
null_bottomk <- calc_null_overlap(rank_tf_ortho, k = k, reverse = TRUE, ncores = ncore)
null_cor <- calc_null_cor(rank_tf_ortho, ncores = ncore)


# Prepare and join all dfs
# ------------------------------------------------------------------------------


# Quantiles dfs of top/bottom k
topk_df <- gen_quantile_df(topk_mat, ncores = ncore)
topk_df <- rename_with(topk_df, ~paste0("Topk_", .), -c("Symbol"))

bottomk_df <- gen_quantile_df(bottomk_mat, ncores = ncore)
bottomk_df <- rename_with(bottomk_df, ~paste0("Bottomk_", .), -c("Symbol"))


# Mean Top k within species (pairs of TF profiles across datasets)
topk_win_hg <- topk_l$Human %>%
  left_join(., pc_ortho, by = c("Symbol" = "Symbol_hg")) %>% 
  dplyr::select(Symbol, Mean) %>%
  dplyr::rename(Topk_win_human = Mean)


topk_win_mm <- topk_l$Mouse %>%
  left_join(., pc_ortho, by = c("Symbol" = "Symbol_mm")) %>% 
  dplyr::select(Symbol_hg, Mean) %>%
  dplyr::rename(Symbol = Symbol_hg, Topk_win_mouse = Mean)


sim_df <- plyr::join_all(
  list(cor_df, topk_df, bottomk_df, topk_win_hg, topk_win_mm, msr_df),
  by = "Symbol")



# Sandbox for exploring ranks between species for a given TF. Look at difference
# in ranks, adding the average RSR values, and rank product
# ------------------------------------------------------------------------------


check_tf <- "ASCL1"


rank_df <- rank_tf_ortho[[check_tf]] %>% 
  mutate(Diff_rank = Rank_aggr_coexpr_hg - Rank_aggr_coexpr_mm,
         Add = Avg_aggr_coexpr_hg + Avg_aggr_coexpr_mm)


# Inspect check TFs top and bottom k overlap across all TFs
overlap_df <- data.frame(
  Symbol = rownames(topk_mat),
  Topk_hg_in_mm = topk_mat[check_tf, ],
  Topk_mm_in_hg = topk_mat[, check_tf],
  Bottomk_hg_in_mm = bottomk_mat[check_tf, ],
  Bottomk_mm_in_hg = bottomk_mat[, check_tf]
)


# Genes in the topk for both species for check TF
overlap_genes <- intersect(
  slice_min(rank_tf_ortho[[check_tf]], Rank_aggr_coexpr_hg, n = k)$Symbol_hg,
  slice_min(rank_tf_ortho[[check_tf]], Rank_aggr_coexpr_mm, n = k)$Symbol_hg
)


# Correlation across metrics. See a greater relationship of Scor with Bottomk 
# than Topk. Still see that paired TFs tend to have greater Bottomk than null,
# but am more convinced by Top of the rankings than Bottom
# ------------------------------------------------------------------------------


sim_cor <- cor(select_if(sim_df, is.numeric), method = "spearman")



# Exploring Spearman's correlations
# High median Scor (0.71). Few TFs that have low correlation, with HESX1 lowest.
# Cautious about interpreting this as species differential... different coverage
# of cell types and inspection doesn't provide immediately striking examples. 
# ------------------------------------------------------------------------------


summ_cor <- summary(sim_df$Cor)
cor_low <- filter(sim_df, Cor < -0.2) %>% arrange(Cor)
cor_high <- filter(sim_df, Cor > 0.8) %>% arrange(-Cor)


# Topk overlap
# ------------------------------------------------------------------------------


summ_topk <- summary(Filter(is.numeric, topk_df))


extremes_topk <- list(
  Count0 = sum(topk_df$Topk_Count == 0),
  Quant_hg1 = sum(topk_df$Topk_Quant_hg_in_mm == 1),
  Quant_mm1 = sum(topk_df$Topk_Quant_mm_in_hg == 1),
  Quant_ortho1 = sum(topk_df$Topk_Quant_ortho == 1)
)


# Mean across rows (human) and cols (mouse) to get most common/generic 
common_topk <- data.frame(
  Symbol = rownames(topk_mat),
  Hg_in_mm = rowMeans(topk_mat),
  Mm_in_hg = colMeans(topk_mat)
)



# Inspect TFs with highest generic human in mouse overlap. Find that they still
# generally are specific to their matched ortholog, despite having more common
# overlap in general. FOXO4 exception: high Mm_in_hg, but Quant_mm_in_hg ~ 0.68
most_common_topk <- common_topk %>% 
  filter(Symbol %in% slice_max(., Hg_in_mm, n = 5)$Symbol | 
         Symbol %in% slice_max(., Mm_in_hg, n = 5)$Symbol) %>% 
  left_join(., sim_df, by = "Symbol")


# TFs with no overlap
no_topk <- filter(sim_df, Topk_Count == 0)



# Bottomk overlap 
# ------------------------------------------------------------------------------


summ_bottomk <- summary(Filter(is.numeric, bottomk_df))


extremes_bottomk <- list(
  Count0 = sum(bottomk_df$Bottomk_Count == 0),
  Quant_hg1 = sum(bottomk_df$Bottomk_Quant_hg_in_mm == 1),
  Quant_mm1 = sum(bottomk_df$Bottomk_Quant_mm_in_hg == 1),
  Quant_ortho1 = sum(bottomk_df$Bottomk_Quant_ortho == 1)
)


# Mean across rows (human) and cols (mouse) to get most common/generic 
common_bottomk <- data.frame(
  Symbol = rownames(bottomk_mat),
  Hg_in_mm = rowMeans(bottomk_mat),
  Mm_in_hg = colMeans(bottomk_mat)
)


# Inspect TFs with highest generic human in mouse overlap. Find that they still
# generally are specific to their matched ortholog, perhaps a touch less than 
# the common topk comparison
most_common_bottomk <- common_bottomk %>% 
  filter(Symbol %in% slice_max(., Hg_in_mm, n = 5)$Symbol | 
         Symbol %in% slice_max(., Mm_in_hg, n = 5)$Symbol) %>% 
  left_join(., sim_df, by = "Symbol")


# TFs with minimal or no overlap
no_bottomk <- filter(sim_df, Bottomk_Count == 0)
minimal_overlap <- filter(sim_df, Topk_Quant_ortho < 0.5 & Bottomk_Quant_ortho < 0.5)


# TFs preserved in bottom but not top
bottom_not_top <- filter(sim_df, Topk_Quant_ortho < 0.5 & Bottomk_Quant_ortho > 0.9)


# TFs showing top specificity at top and bottom
most_specific <- list(
  Top1 = filter(sim_df, Topk_Quant_ortho == 1 & Bottomk_Quant_ortho == 1),
  Top99 = filter(sim_df, Topk_Quant_ortho > 0.99 & Bottomk_Quant_ortho > 0.99)
)



# Plots
# ------------------------------------------------------------------------------


# Histogram of similarity overlaid with median of null overlap

hist_ortho_similarity <- function(sim_df, 
                                  stat, 
                                  null_vec, 
                                  title = NULL,
                                  xlab, 
                                  colour) {
  
  
  ggplot(sim_df, aes(x = !!sym(stat))) +
    geom_histogram(bins = 100, col = colour, fill = colour) +
    geom_vline(xintercept = median(sim_df[[stat]]), linewidth = 1.4, col = "black") +
    geom_vline(xintercept = median(null_vec), linewidth = 1.4, col = "grey") +
    ggtitle(title) +
    xlab(xlab) +
    ylab("Frequency") +
    theme_classic() +
    theme(text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5))
  
}



p1a <- hist_ortho_similarity(sim_df, 
                             stat = "Cor", 
                             null_vec = null_cor, 
                             title = "Cross-species coexpression agreement",
                             xlab = "Spearman's correlation", 
                             # xlim = c(-0.4, 1), 
                             colour = "slategrey") 
p1a <- p1a + xlim(c(-0.4, 1))


p1b <- hist_ortho_similarity(sim_df, 
                             stat = "Topk_Count", 
                             null_vec = null_topk, 
                             xlab = expr("Top"[!!k]),
                             colour = "#d53e4f")


p1c <- hist_ortho_similarity(sim_df, 
                             stat = "Bottomk_Count", 
                             null_vec = null_bottomk, 
                             xlab = expr("Bottom"[!!k]),
                             colour = "#3288bd")


p1 <- plot_grid(p1a, p1b, p1c, ncol = 1)

ggsave(p1, height = 12, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("ortho_similarity_hist_k=", k, ".png")))



# Scatter of mouse versus human RSR for a given TF

p2 <- qplot(rank_tf_ortho[[check_tf]],
            xvar = "Avg_aggr_coexpr_hg",
            yvar = "Avg_aggr_coexpr_mm",
            title = check_tf) +
  xlab("Aggregate coexpression (human)") +
  ylab("Aggregate coexpression (mouse)")



# Scatter of overlap quantiles in mouse/human

p3a <- qplot(sim_df, xvar = "Topk_Quant_hg_in_mm", yvar = "Topk_Quant_mm_in_hg") +
  xlab("Ortholog retrieval score (human in mouse)") + 
  ylab("Ortholog retrieval score (mouse in human)")


p3b <- qplot(sim_df, xvar = "Bottomk_Quant_hg_in_mm", yvar = "Bottomk_Quant_mm_in_hg") +
  xlab(expr("Ortholog retrieval score (human in mouse)")) +
  ylab(expr("Ortholog retrieval score (mouse in human)")) +
  ggtitle(expr("Bottom"[!!k] ~ "(negative coexpresson)"))


p3c <- qplot(sim_df, xvar = "Bottomk_Quant_ortho", yvar = "Topk_Quant_ortho") +
  xlab(expr("Bottom"[!!k] ~ "average ortholog retrieval score")) +
  ylab(expr("Top"[!!k] ~ "average ortholog retrieval score"))


ggsave(p3a, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("top", k, "_quantile_between_species.png")))


ggsave(p3b, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("bottom", k, "_quantile_between_species.png")))


ggsave(p3c, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("top_vs_bottom", k, "_quantile_between_species.png")))



# Scatter of topk vs bottom, and raw overlap vs quantile
p5a <- qplot(sim_df, xvar = "Topk_Count", yvar = "Bottomk_Count")
p5b <- qplot(sim_df, xvar = "Topk_Count", yvar = "Topk_Quant_ortho")
p5c <- qplot(sim_df, xvar = "Bottomk_Count", yvar = "Bottomk_Quant_ortho")


# Scatter plot of topk counts between species. Label the matched TF in red,
# and label TFs that were high in both species as black

label_genes <- overlap_df %>% 
  filter(Symbol %in% slice_max(overlap_df, Topk_hg_in_mm, n = 100)$Symbol &
         Symbol %in% slice_max(overlap_df, Topk_mm_in_hg, n = 100)$Symbol &
         Symbol != check_tf) %>% 
  pull(Symbol)


pdf1 <- mutate(overlap_df,
               tf_label = Symbol == check_tf,
               gene_label = Symbol %in% label_genes)

p6 <- 
  ggplot(pdf1, aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg)) +
  geom_jitter(shape = 21, size = 2.4, width = 0.5, height = 0.5) +
  geom_text_repel(
    data = filter(pdf1, tf_label),
    aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg, label = Symbol, fontface = "italic"),
    size = 6,
    nudge_y = 2,
    segment.size = 0.1,
    segment.color = "grey50",
    color = "red") +
  geom_text_repel(
    data = filter(pdf1, gene_label),
    aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg, label = Symbol, fontface = "italic"),
    size = 4,
    segment.size = 0.1,
    segment.color = "grey50") +
  ggtitle(paste0(check_tf, " ortholog retrieval")) +
  xlab(expr("Top"[!!k] ~ "human in mouse")) +
  ylab(expr("Top"[!!k] ~ "mouse in human")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


ggsave(p6, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0(check_tf, "_top", k, "_count_between_species.png")))


# Hist of topk counts for each species for check TF

p7a <- plot_hist(pdf1, 
                 stat_col = "Topk_hg_in_mm", 
                 title = paste("Human", check_tf), 
                 nbins = 30) + 
  geom_vline(xintercept = filter(pdf1, Symbol == check_tf)$Topk_hg_in_mm,
             linewidth = 1.6,
             col = "royalblue") +
  xlab(expr("Top"[!!k] ~ "human in mouse"))


p7b <- plot_hist(pdf1, 
                 stat_col = "Topk_mm_in_hg", 
                 title = paste("Mouse", str_to_title(check_tf)), 
                 nbins = 30) +
  geom_vline(xintercept = filter(pdf1, Symbol == check_tf)$Topk_mm_in_hg, 
             linewidth = 1.6,
             col = "goldenrod") +
  xlab(expr("Top"[!!k] ~ "mouse in human"))


p7 <- plot_grid(p7a, p7b, nrow = 2)


ggsave(p7, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0(check_tf, "_top", k, "_count_hist.png")))


# Minimal density plot of species overlap to use in schematic

p8a <- ggplot(overlap_df, aes(x = Topk_hg_in_mm)) + 
  geom_density(linewidth = 2.5) + 
  geom_vline(xintercept = filter(overlap_df, Symbol == check_tf)$Topk_hg_in_mm, 
             colour = "royalblue", 
             linewidth = 3) +
  theme_nothing() + 
  theme(text = element_blank())


p8b <- ggplot(overlap_df, aes(x = Topk_mm_in_hg)) + 
  geom_density(linewidth = 2.5) + 
  geom_vline(xintercept = filter(overlap_df, Symbol == check_tf)$Topk_mm_in_hg, 
             colour = "goldenrod", 
             linewidth = 3) +
  theme_nothing() + 
  theme(text = element_blank())


ggsave(p8a, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "demo_topk_distn_ortho_overlap_hg.png"))


ggsave(p8b, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "demo_topk_distn_ortho_overlap_mm.png"))


# Plots of common/generic

p9a <- plot_grid(
  plot_hist(common_topk, stat_col = "Hg_in_mm") + xlab("Human in mouse"),
  plot_hist(common_topk, stat_col = "Hg_in_mm") + xlab("Mouse in human"),
  nrow = 2)


p9b <- plot_grid(
  plot_hist(common_bottomk, stat_col = "Hg_in_mm") + xlab("Human in mouse"),
  plot_hist(common_bottomk, stat_col = "Hg_in_mm") + xlab("Mouse in human"),
  nrow = 2)

p9c <- qplot(common_topk, xvar = "Hg_in_mm", yvar = "Mm_in_hg")
p9d <- qplot(common_bottomk, xvar = "Hg_in_mm", yvar = "Mm_in_hg")


# Scatter of Scor versus overlap

p10a <- qplot(sim_df, xvar = "Bottomk_Count", yvar = "Cor")
p10b <- qplot(sim_df, xvar = "Topk_Count", yvar = "Cor")


# Relationship between the mean top K within species (paired datasets) and the
# top k between species of the aggregate profiles

p11a <- qplot(sim_df, xvar = "Topk_Count", yvar = "Topk_win_human") +
  geom_smooth(method = "lm", colour = "red") +
  xlab(expr("Cross-species Top"[!!k])) +
  ylab(expr("Within species mean Top"[!!k])) +
  ggtitle("Human")

p11b <- qplot(sim_df, xvar = "Topk_Count", yvar = "Topk_win_mouse") +
  geom_smooth(method = "lm", colour = "red") +
  xlab(expr("Cross-species Top"[!!k])) +
  ylab(expr("Within species mean Top"[!!k])) +
  ggtitle("Mouse")


p11 <- plot_grid(p11a, p11b, nrow = 1)


ggsave(p11, height = 6, width = 12, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("within_vs_across_species_topk_k=", k, ".png")))
