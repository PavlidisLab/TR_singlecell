## Examining aggregate coexpression profiles between mouse and human
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

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Saved list RDS of the ranks
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)


# Only keep ortho genes measured in both species

pc_ortho <- filter(pc_ortho,
                   Symbol_hg %in% pc_hg$Symbol &
                   Symbol_mm %in% pc_mm$Symbol)

tfs_ortho <- filter(pc_ortho, 
                    Symbol_hg %in% names(rank_tf_hg) & 
                    Symbol_mm %in% names(rank_tf_mm))


# Functions
# ------------------------------------------------------------------------------


# Gene is assumed to be an ortho gene found in pc_df. Retrieve and join the 
# corresponding gene rank dfs from the mouse and human rank lists. Re-rank
# average RSR just using ortho genes

join_ortho_ranks <- function(gene,
                             pc_df = pc_ortho,
                             rank_hg = rank_tf_hg,
                             rank_mm = rank_tf_mm) {
  
  gene_ortho <- filter(pc_df, Symbol_hg == gene | Symbol_mm == gene)
  
  if (nrow(gene_ortho) != 1) stop("A one to one match was not found")
  
  df_hg <- 
    left_join(rank_hg[[gene_ortho$Symbol_hg]], 
              pc_ortho[, c("Symbol_hg", "ID")],
              by = c("Symbol" = "Symbol_hg")) %>% 
    filter(!is.na(ID)) %>% 
    mutate(Rank_RSR = rank(-Avg_RSR, ties.method = "min"))
  
  df_mm <- 
    left_join(rank_mm[[gene_ortho$Symbol_mm]], 
              pc_ortho[, c("Symbol_mm", "ID")], 
              by = c("Symbol" = "Symbol_mm")) %>% 
    filter(!is.na(ID)) %>% 
    mutate(Rank_RSR = rank(-Avg_RSR, ties.method = "min"))
  
  sim_df <- 
    left_join(df_hg, df_mm,
              by = "ID",
              suffix = c("_hg", "_mm")) %>% 
    filter(!is.na(Avg_RSR_hg) & !is.na(Avg_RSR_mm))
  
  return(sim_df)
}


# Return a df of the spearman cor between human and mouse RSR rankings

calc_ortho_cor <- function(ortho_rank_l, ncores = 1) {
  
  symbols <- names(ortho_rank_l)
  
  cor_l <- mclapply(ortho_rank_l, function(x) {
    cor(x$Avg_RSR_hg, x$Avg_RSR_mm, method = "spearman")
  }, mc.cores = ncores)
  
  cor_df <- data.frame(Cor = unlist(cor_l), Symbol = symbols)
  
  return(cor_df)
}


# Flip sign of average RSR (rank metric) to prioritize bottom of ranks

reverse_ranks <- function(df) {
  mutate(df, Avg_RSR_hg = -Avg_RSR_hg, Avg_RSR_mm = -Avg_RSR_mm)
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
      Human = slice_max(x, Avg_RSR_hg, n = k)$Symbol_hg,
      Mouse = slice_max(x, Avg_RSR_mm, n = k)$Symbol_hg
    )
  }, mc.cores = ncores)
  names(topk_l) <- genes
  
  # Generate overlap matrix
  overlap_mat <- matrix(0, nrow = length(topk_l), ncol = length(topk_l))
  rownames(overlap_mat) <- colnames(overlap_mat) <- genes
  
  for (i in 1:nrow(overlap_mat)) {
    for (j in 1:ncol(overlap_mat)) {
      overlap_mat[i, j] <- topk_intersect(topk_l[[i]]$Human, topk_l[[j]]$Mouse)
    }
  }
  
  return(overlap_mat)
}


# Given gene x gene overlap_mat, return a data.frame containing the count of
# overlap for the orthologous match, as well as the percentile of the match
# relative to all other TFs for each species

gen_percentile_df <- function(overlap_mat, ncores = 1) {
  
  genes <- intersect(rownames(overlap_mat), colnames(overlap_mat))
  
  perc <- mclapply(genes, function(x) {
    
    count_overlap <- overlap_mat[x, x]
    hg_in_mm <- overlap_mat[x, setdiff(tfs_ortho$Symbol_hg, x)]
    mm_in_hg <- overlap_mat[setdiff(tfs_ortho$Symbol_hg, x), x]
    
    data.frame(
      Symbol = x,
      Count = count_overlap,
      Perc_hg_in_mm = ecdf(hg_in_mm)(count_overlap),
      Perc_mm_in_hg = ecdf(mm_in_hg)(count_overlap))
    
  }, mc.cores = ncores)
  
  perc_df <- do.call(rbind, perc)
  perc_df$Perc_ortho <- rowMeans(perc_df[, c("Perc_hg_in_mm", "Perc_mm_in_hg")])
  
  return(perc_df)
}


# Join the human and mouse rankings for only orthologous genes, then generate
# the similarity objects. The human symbols used for representing ortho TFs. 
# -------------------------------------------------------------------------------


# List of joined human and mouse TF rankings
rank_tf_ortho <- mclapply(tfs_ortho$Symbol_hg, join_ortho_ranks, mc.cores = ncore)
names(rank_tf_ortho) <- tfs_ortho$Symbol_hg


# Get the Spearman's correlation of rankings between every ortho TF
cor_df <- calc_ortho_cor(rank_tf_ortho, ncores = ncore)


# Top and bottom K matrices 
topk_mat <- calc_overlap_mat(rank_tf_ortho, k = k, ncores = ncore)
bottomk_mat <- calc_overlap_mat(rank_tf_ortho, k = k, reverse = TRUE, ncores = ncore)


# Percentiles dfs of top/bottom k
topk_df <- gen_percentile_df(topk_mat, ncores = ncore)
bottomk_df <- gen_percentile_df(bottomk_mat, ncores = ncore)


# Measurement counts for each TF
msr_df <- data.frame(
  Symbol = tfs_ortho$Symbol_hg,
  N_hg = rowSums(msr_hg[tfs_ortho$Symbol_hg, ]),
  N_mm = rowSums(msr_mm[tfs_ortho$Symbol_mm, ])
)


# Prepare and join all dfs
topk_df <- rename_with(topk_df, ~paste0("Topk_", .), -c("Symbol"))
bottomk_df <- rename_with(bottomk_df, ~paste0("Bottomk_", .), -c("Symbol"))
sim_df <- plyr::join_all(list(cor_df, topk_df, bottomk_df, msr_df), by = "Symbol")



# Sandbox for exploring ranks between species for a given TF. Look at difference
# in ranks, adding the average RSR values, and rank product
# ------------------------------------------------------------------------------


check_tf <- "ASCL1"


rank_common <- rank_tf_ortho[[check_tf]] %>% 
  mutate(Diff_rank = Rank_RSR_hg - Rank_RSR_mm,
         Add = Avg_RSR_hg + Avg_RSR_mm,
         RP = rank(Rank_RSR_hg * Rank_RSR_mm)) %>% 
  relocate(Symbol_hg, RP, Rank_RSR_hg, Rank_RSR_mm, Add, Diff_rank, Avg_RSR_hg, Avg_RSR_mm)


# Inspect TF overlap

overlap_df <- data.frame(
  Symbol = rownames(topk_mat),
  Topk_hg_in_mm = topk_mat[check_tf, ],
  Topk_mm_in_hg = topk_mat[, check_tf],
  Bottomk_hg_in_mm = bottomk_mat[check_tf, ],
  Bottomk_mm_in_hg = bottomk_mat[, check_tf]
)


overlap_genes <- intersect(
  slice_min(rank_tf_ortho$ASCL1, Rank_RSR_hg, n = k)$Symbol_hg,
  slice_min(rank_tf_ortho$ASCL1, Rank_RSR_mm, n = k)$Symbol_hg
)



# TODO: remove/replace, but schematic of topk
pjj1 <- ggplot(overlap_df, aes(x = Topk_hg_in_mm)) + 
  geom_density(linewidth = 2.5) + 
  geom_vline(xintercept = filter(overlap_df, Symbol == check_tf)$Topk_hg_in_mm, colour = "royalblue", linewidth = 3) +
  theme_nothing() + 
  theme(text = element_blank())


pjj2 <- ggplot(overlap_df, aes(x = Topk_mm_in_hg)) + 
  geom_density(linewidth = 2.5) + 
  geom_vline(xintercept = filter(overlap_df, Symbol == check_tf)$Topk_mm_in_hg, colour = "goldenrod", linewidth = 3) +
  theme_nothing() + 
  theme(text = element_blank())


ggsave(pjj1, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "demo_topk_distn_ortho_overlap_hg.png"))


ggsave(pjj2, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "demo_topk_distn_ortho_overlap_mm.png"))


# Correlation across metrics. The elevated correlation of bottom K with
# Scor makes me concerned that the high Scor values are being elevated by the 
# bottom of the rankings, which are less reliable. Also see strong cor between
# measurement and Scor. 
# ------------------------------------------------------------------------------


sim_cor <- cor(select_if(sim_df, is.numeric), method = "spearman")



# Exploring Spearman's correlations
# ------------------------------------------------------------------------------


# High median Scor (0.71). Few TFs that have low correlation, with HESX1 lowest.
# Cautious about interpreting this as species differential... different coverage
# of cell types and inspection doesn't provide immediately striking examples. 

summ_cor <- summary(sim_df$Cor)
cor_low <- filter(sim_df, Cor < -0.2) %>% arrange(Cor)
cor_high <- filter(sim_df, Cor > 0.8) %>% arrange(-Cor)


# Topk overlap
# ------------------------------------------------------------------------------


# Find that 174 (14.0%) of human and 171 (13.8%) of mouse TFs have best match
# with their ortholog.

summ_topk <- summary(Filter(is.numeric, topk_df))

extremes_topk <- list(
  Count0 = sum(topk_df$Topk_Count == 0),
  Perc_hg1 = sum(topk_df$Topk_Perc_hg_in_mm == 1),
  Perc_mm1 = sum(topk_df$Topk_Perc_mm_in_hg == 1),
  Perc_ortho1 = sum(topk_df$Topk_Perc_ortho == 1)
)


# Mean across rows (human) and cols (mouse) to get most common/generic 

common_topk <- data.frame(
  Symbol = rownames(topk_mat),
  Hg_in_mm = rowMeans(topk_mat),
  Mm_in_hg = colMeans(topk_mat)
)



# Inspect TFs with highest generic human in mouse overlap. Find that they still
# generally are specific to their matched ortholog, despite having more common
# overlap in general. FOXO4 exception: high Mm_in_hg, but Perc_mm_in_hg ~ 0.68

most_common_topk <- common_topk %>% 
  filter(Symbol %in% slice_max(., Hg_in_mm, n = 5)$Symbol | 
         Symbol %in% slice_max(., Mm_in_hg, n = 5)$Symbol) %>% 
  left_join(., sim_df, by = "Symbol")


# TFs with no overlap

no_topk <- filter(sim_df, Topk_Count == 0)



# Bottomk overlap 
# ------------------------------------------------------------------------------


# Find that 53 (4.3%) of human and 70 (5.6%) of mouse TFs have best match
# with their ortholog.

summ_bottomk <- summary(Filter(is.numeric, bottomk_df))

extremes_bottomk <- list(
  Count0 = sum(bottomk_df$Bottomk_Count == 0),
  Perc_hg1 = sum(bottomk_df$Bottomk_Perc_hg_in_mm == 1),
  Perc_mm1 = sum(bottomk_df$Bottomk_Perc_mm_in_hg == 1),
  Perc_ortho1 = sum(bottomk_df$Bottomk_Perc_ortho == 1)
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
minimal_overlap <- filter(sim_df, Topk_Perc_ortho < 0.5 & Bottomk_Perc_ortho < 0.5)

# TFs preserved in bottom but not top

bottom_not_top <- filter(sim_df, Topk_Perc_ortho < 0.5 & Bottomk_Perc_ortho > 0.9)


# TFs showing top specificity at top and bottom

most_specific <- list(
  Top1 = filter(sim_df, Topk_Perc_ortho == 1 & Bottomk_Perc_ortho == 1),
  Top99 = filter(sim_df, Topk_Perc_ortho > 0.99 & Bottomk_Perc_ortho > 0.99)
)


# Looking at the correlation between topk across and within species

topk_l <- readRDS(paste0("/space/scratch/amorin/R_objects/human_mouse_topk=", k, "_similarity_df.RDS"))
topk_l$Mouse <- left_join(topk_l$Mouse, pc_ortho, by = c("Symbol" = "Symbol_mm"))

sim_df2 <- sim_df %>% 
  left_join(., topk_l$Human, by = "Symbol") %>% 
  left_join(., topk_l$Mouse[, c("Symbol_hg", "Mean")], 
            by = c("Symbol" = "Symbol_hg"),
            suffix = c("_human", "_mouse"))


ggplot(sim_df2, 
         aes(x = fct_reorder(Family, desc(Topk_Count), .fun = mean), 
             y = Topk_Count)) +
    geom_boxplot(fill = "slategrey", outlier.shape = NA) +
    geom_jitter(width = 0.25, height = 0.1, shape = 21) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.margin = margin(c(10, 20, 10, 10)))


p10a <- qplot(sim_df2, xvar = "Topk_Count", yvar = "Mean_human") +
  geom_smooth(method = "lm", colour = "red") +
  xlab(expr("Cross-species Top"[!!k])) +
  ylab(expr("Within species mean Top"[!!k])) +
  ggtitle("Human")

p10b <- qplot(sim_df2, xvar = "Topk_Count", yvar = "Mean_mouse") +
  geom_smooth(method = "lm", colour = "red") +
  xlab(expr("Cross-species Top"[!!k])) +
  ylab(expr("Within species mean Top"[!!k])) +
  ggtitle("Mouse")


p10 <- plot_grid(p10a, p10b, nrow = 1)


ggsave(p10, height = 6, width = 12, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("within_vs_across_species_topk_k=", k, ".png")))


# Plots
# ------------------------------------------------------------------------------


# Histogram of ortho TF spearman cors

px1 <- plot_hist(cor_df, 
                stat_col = "Cor", 
                xlab = "Spearman's correlation", 
                title = "n=1241 orthologous TFs") + 
  xlim(c(-0.4, 1))



ggsave(px1, height = 5, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "ortho_rank_similarity.png"))


# Scatter of mouse versus human RSR for a given TF

px2 <- qplot(rank_tf_ortho[[check_tf]], 
             xvar = "Avg_RSR_hg", 
             yvar = "Avg_RSR_mm", 
             title = check_tf)


# Scatter of overlap percentiles in mouse/human

px3a <- qplot(sim_df, xvar = "Topk_Perc_hg_in_mm", yvar = "Topk_Perc_mm_in_hg") +
  xlab(expr("Top"[!!k] ~ "quantile human in mouse")) +
  ylab(expr("Top"[!!k] ~ "quantile mouse in human"))


px3b <- qplot(sim_df, xvar = "Bottomk_Perc_mm_in_hg", yvar = "Bottomk_Perc_hg_in_mm") +
  xlab(expr("Bottom"[!!k] ~ "quantile mouse in human")) +
  ylab(expr("Bottom"[!!k] ~ "quantile human in mouse"))


px3c <- qplot(sim_df, xvar = "Bottomk_Perc_ortho", yvar = "Topk_Perc_ortho") +
  xlab(expr("Bottom"[!!k] ~ "quantile orthologous")) +
  ylab(expr("Top"[!!k] ~ "quantile orthologous"))


ggsave(px3a, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "topk_percentile_between_species.png"))


ggsave(px3b, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "bottomk_percentile_between_species.png"))


ggsave(px3c, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "topk_vs_bottomk_percentile_between_species.png"))



# Plots of the raw overlap counts

px4a <- plot_hist(sim_df, nbins = 50, stat_col = "Topk_Count")
px4b <- plot_hist(sim_df, nbins = 50, stat_col = "Bottomk_Count")
px4c <- qplot(sim_df, xvar = "Topk_Count", yvar = "Bottomk_Count")


# Scatter of perc ortho versus top k

px5a <- qplot(sim_df, xvar = "Topk_Count", yvar = "Topk_Perc_ortho")
px5b <- qplot(sim_df, xvar = "Bottomk_Count", yvar = "Bottomk_Perc_ortho")


# Scatter plot of topk counts between species. Label the matched TF in red,
# and label TFs that were high in both species as black


label_genes <- overlap_df %>% 
  filter(Symbol %in% slice_max(overlap_df, Topk_hg_in_mm, n = 100)$Symbol &
         Symbol %in% slice_max(overlap_df, Topk_mm_in_hg, n = 100)$Symbol &
         Symbol != check_tf) %>% 
  pull(Symbol)


plot_df6 <- mutate(overlap_df, 
                   tf_label = Symbol == check_tf,
                   gene_label = Symbol %in% label_genes)


px6a <- 
  ggplot(plot_df6, aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg)) +
  geom_jitter(shape = 21, size = 2.4, width = 0.5, height = 0.5) +
  geom_text_repel(
    data = filter(plot_df6, tf_label),
    aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg, label = Symbol, fontface = "italic"),
    size = 5,
    nudge_y = 2,
    segment.size = 0.1,
    segment.color = "grey50",
    color = "red") +
  geom_text_repel(
    data = filter(plot_df6, gene_label),
    aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg, label = Symbol, fontface = "italic"),
    size = 5,
    segment.size = 0.1,
    segment.color = "grey50") +
  xlab(expr("Top"[!!k] ~ "human in mouse")) +
  ylab(expr("Top"[!!k] ~ "mouse in human")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


ggsave(px6a, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0(check_tf, "_topk_count_between_species.png")))

  
# Hist of topk counts for each species
# TODO: Consider specific use of gghist; split calls


px6b <- plot_grid(
  
  plot_hist(plot_df6, stat_col = "Topk_hg_in_mm", title = paste("Human", check_tf), nbins = 30) + 
    geom_vline(xintercept = filter(plot_df6, Symbol == check_tf)$Topk_hg_in_mm,
               size = 1.6,
               col = "royalblue") +
    xlab(expr("Top"[!!k] ~ "human in mouse")),
  
  plot_hist(plot_df6, stat_col = "Topk_mm_in_hg", title = paste("Mouse", str_to_title(check_tf)), nbins = 30) +
    geom_vline(xintercept = filter(plot_df6, Symbol == check_tf)$Topk_mm_in_hg, 
               size = 1.6,
               col = "goldenrod") +
    xlab(expr("Top"[!!k] ~ "mouse in human")),
  
  nrow = 2)


ggsave(px6b, height = 7, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0(check_tf, "_topk_count_hist.png")))


# Plots of common/generic


px7a <- plot_grid(
  plot_hist(common_topk, stat_col = "Hg_in_mm") + xlab("Human in mouse"),
  plot_hist(common_topk, stat_col = "Hg_in_mm") + xlab("Mouse in human"),
  nrow = 2)


px7b <- plot_grid(
  plot_hist(common_bottomk, stat_col = "Hg_in_mm") + xlab("Human in mouse"),
  plot_hist(common_bottomk, stat_col = "Hg_in_mm") + xlab("Mouse in human"),
  nrow = 2)


px7c <- qplot(common_topk, xvar = "Hg_in_mm", yvar = "Mm_in_hg")
px7d <- qplot(common_bottomk, xvar = "Hg_in_mm", yvar = "Mm_in_hg")


# Scatter of Scor versus overlap

px8a <- qplot(sim_df, xvar = "Bottomk_Count", yvar = "Cor")
px8b <- qplot(sim_df, xvar = "Topk_Count", yvar = "Cor")



# Demoing same stacked barchart idea, but just for ortho coexpression at different
# cut-offs

ortho_coexpr_counts <- mclapply(rank_tf_ortho, function(x) {
  
  k200 = filter(x, Rank_RSR_hg <= 200 & Rank_RSR_mm <= 200) 
  k500 = filter(x, Rank_RSR_hg <= 500 & Rank_RSR_mm <= 500 & (!Symbol_hg %in% k200$Symbol_hg)) 
  k1000 = filter(x, Rank_RSR_hg <= 1000 & Rank_RSR_mm <= 1000 & (!Symbol_hg %in% k500$Symbol_hg))
  
  data.frame(K200 = nrow(k200),
             K500 = nrow(k500),
             K1000 = nrow(k1000))
}, mc.cores = ncore)


n_topk_ortho <- data.frame(
  Symbol = names(rank_tf_ortho),
  do.call(rbind, ortho_coexpr_counts)
)


plot_df9 <- pivot_longer(n_topk_ortho,
                         cols = c("K200", "K500", "K1000"),
                         names_to = "Cutoff",
                         values_to = "Count") %>% 
  mutate(Cutoff = factor(Cutoff, levels = c("K1000", "K500", "K200")))


# p9_cols <- c("#f2f0f7","#cbc9e2", "#9e9ac8", "#6a51a3")
# p9_cols <- c("#9e9ac8","#6a51a3", "#54278f", "#3f007d")
# p9_cols <- c("#9e9ac8","#6a51a3", "#3f007d", "black")
p9_cols <- c("#9e9ac8","#6a51a3", "#3f007d")


# plot_df9 <- filter(plot_df9, Symbol %in% sample(names(rank_tf_ortho), 100))


p9a <- ggplot(plot_df9, aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Cutoff, colour = Cutoff)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of orthologous interactions") +
  xlab("Transcription regulator") +
  ggtitle(paste0("N=", nrow(n_topk_ortho))) +
  scale_fill_manual(values = p9_cols) +
  scale_colour_manual(values = p9_cols) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 25),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        plot.margin = margin(c(10, 20, 10, 10)))


ggsave(p9a, height = 8, width = 12, device = "png", dpi = 300,
       filename = file.path(plot_dir, paste0("ortho_topk_counts_all_tf.png")))



# Null distn of overlap


calc_null_topk <- function(rank_tf_ortho, ncores = 1) {
  
  # Helper to get size of overlap between mouse and human for sampled TFs
  
  calc_topk <- function(k, sample_tfs) {
    topk_intersect(
      slice_max(rank_tf_ortho[[sample_tfs[1]]], Avg_RSR_hg, n = k)$Symbol_hg,
      slice_max(rank_tf_ortho[[sample_tfs[2]]], Avg_RSR_mm, n = k)$Symbol_hg
    )
  }
  
  # Iteratively sample TFs and calculate overlap between mouse and human
  
  null_topk <- mclapply(1:1000, function(x) {
    
    sample_tfs <- sample(names(rank_tf_ortho), 2, replace = FALSE)
    
    data.frame(
      K200 = calc_topk(k = 200, sample_tfs),
      K500 = calc_topk(k = 500, sample_tfs),
      K1000 = calc_topk(k = 1000, sample_tfs))
    
  }, mc.cores = ncores)
  
  null_df <- do.call(rbind, null_topk)
  return(null_df)
}


set.seed(5)

n_topk_null <- calc_null_topk(rank_tf_ortho, ncores = ncore)


p9b1 <- plot_hist(n_topk_ortho, stat_col = "K200") + geom_vline(xintercept = median(n_topk_null$K200))
p9b2 <- plot_hist(n_topk_null, stat_col = "K200")

p9b3 <- plot_hist(n_topk_ortho, stat_col = "K500") + geom_vline(xintercept = median(n_topk_null$K500))
p9b4 <- plot_hist(n_topk_null, stat_col = "K500")

p9b5 <- plot_hist(n_topk_ortho, stat_col = "K1000") + geom_vline(xintercept = median(n_topk_null$K1000))
p9b6 <- plot_hist(n_topk_null, stat_col = "K1000")

plot_grid(p9b1, p9b3, p9b5, ncol = 1)



# TODO: why are K200 and K500 shaped identically

# Not because of reorder...
plot_df9 <- pivot_longer(n_topk_ortho,
                         cols = c("K200", "K500", "K1000"),
                         names_to = "Cutoff",
                         values_to = "Count") %>% 
  mutate(Cutoff = factor(Cutoff, levels = c("K1000", "K500", "K200")))


tf_order <- plot_df9 %>% 
  group_by(Symbol) %>% 
  summarise(Med = median(Count)) %>% 
  arrange(Med)


plot_df9 <- mutate(plot_df9, Symbol = factor(Symbol, levels = unique(tf_order$Symbol)))


ggplot(plot_df9, aes(x = Symbol, y = Count, fill = Cutoff, colour = Cutoff)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of orthologous interactions") +
  xlab("Transcription factor") +
  ggtitle(paste0("N=", nrow(n_topk_ortho))) +
  scale_fill_manual(values = p9_cols) +
  scale_colour_manual(values = p9_cols) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


# Not because of fill and colour...
ggplot(plot_df9, aes(x = Symbol, y = Count, fill = Cutoff)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of orthologous interactions") +
  xlab("Transcription factor") +
  ggtitle(paste0("N=", nrow(n_topk_ortho))) +
  scale_fill_manual(values = p9_cols) +
  # scale_colour_manual(values = p9_cols) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


# Looks like because of factor order of cutoff?

plot_df9 <- pivot_longer(n_topk_ortho,
                         cols = c("K200", "K500", "K1000"),
                         names_to = "Cutoff",
                         values_to = "Count") %>% 
  mutate(Cutoff = factor(Cutoff, levels = c("K1000", "K500", "K200")))


plot_df9 <- filter(plot_df9, Symbol %in% sample(names(rank_tf_ortho), 100))


ggplot(plot_df9, aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Cutoff, colour = Cutoff)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of orthologous interactions") +
  xlab("Transcription factor") +
  ggtitle(paste0("N=", nrow(n_topk_ortho))) +
  scale_fill_manual(values = p9_cols) +
  scale_colour_manual(values = p9_cols) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = c(0.85, 0.85),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))
