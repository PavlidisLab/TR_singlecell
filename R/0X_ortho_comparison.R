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

# Ribosomal genes
sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))

# Saved list RDS of the ranks
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)
rank_ribo_hg <- readRDS(rank_ribo_hg_path)
rank_ribo_mm <- readRDS(rank_ribo_mm_path)

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
# corresponding gene ranks from the mouse and human rank lists.

join_ortho_ranks <- function(gene,
                             pc_df = pc_ortho,
                             rank_hg = rank_tf_hg,
                             rank_mm = rank_tf_mm) {
  
  gene_ortho <- filter(pc_df, Symbol_hg == gene | Symbol_mm == gene)
  
  if (nrow(gene_ortho) != 1) stop("A one to one match was not found")
  
  df_hg <- left_join(rank_hg[[gene_ortho$Symbol_hg]], 
                     pc_ortho[, c("Symbol_hg", "ID")],
                     by = c("Symbol" = "Symbol_hg")) %>% 
    filter(!is.na(ID))
  
  
  df_mm <- left_join(rank_mm[[gene_ortho$Symbol_mm]], 
                     pc_ortho[, c("Symbol_mm", "ID")], 
                     by = c("Symbol" = "Symbol_mm")) %>% 
    filter(!is.na(ID))
  
  
  sim_df <- left_join(df_hg, df_mm, 
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


# TFs with no overlap

no_bottomk <- filter(sim_df, Bottomk_Count == 0)


# TFs showing top specificity at top and bottom

most_specific <- list(
  Top1 = filter(sim_df, Topk_Perc_ortho == 1 & Bottomk_Perc_ortho == 1),
  Top99 = filter(sim_df, Topk_Perc_ortho > 0.99 & Bottomk_Perc_ortho > 0.99)
)


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

px3a <- qplot(sim_df, xvar = "Topk_Perc_mm_in_hg", yvar = "Topk_Perc_hg_in_mm") +
  xlab(paste0("Top K=", k, " percentile mouse in human")) +
  ylab(paste0("Top K=", k, " percentile human in mouse"))


px3b <- qplot(sim_df, xvar = "Bottomk_Perc_mm_in_hg", yvar = "Bottomk_Perc_hg_in_mm") +
  xlab(paste0("Bottom K=", k, " percentile mouse in human")) +
  ylab(paste0("Bottom K=", k, " percentile human in mouse"))


px3c <- qplot(sim_df, xvar = "Bottomk_Perc_ortho", yvar = "Topk_Perc_ortho") +
  xlab(paste0("Bottom K=", k, " percentile ortho")) +
  ylab(paste0("Top K=", k, " percentile ortho"))


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
  xlab("Human in mouse") +
  ylab("Mouse in human") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



# Hist of topk counts for each species


px6b <- plot_grid(
  
  plot_hist(plot_df6, stat_col = "Topk_hg_in_mm") + 
    geom_vline(xintercept = filter(plot_df6, Symbol == check_tf)$Topk_hg_in_mm,
               # width = 1.6,
               col = "royalblue") +
    xlab("Human in mouse"),
  
  plot_hist(plot_df6, stat_col = "Topk_mm_in_hg") + 
    geom_vline(xintercept = filter(plot_df6, Symbol == check_tf)$Topk_mm_in_hg, 
               # width = 1.6,
               col = "goldenrod") +
    xlab("Mouse in human"),
  
  nrow = 2)


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
