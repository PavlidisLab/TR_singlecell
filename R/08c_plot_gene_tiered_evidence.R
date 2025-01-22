## Representing the tiers of evidence (+literature evidence) for a TF as heatmap
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggrepel)
library(pheatmap)
library(grid)
library(gridExtra)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 500

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Integrated rankings
rank_hg <- readRDS(rank_int_hg_path)
rank_mm <- readRDS(rank_int_mm_path)

# Tiered list (grouped by TFs)
tier_l <- readRDS(tiered_evidence_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)


# Setting up rank info
# ------------------------------------------------------------------------------


plot_tf <- "ASCL1"

# Isolating evidence for the given TF
rank_l <- tier_l[[plot_tf]]
tf_hg <- rank_hg[[plot_tf]][, c("Symbol", "Rank_bind", "Rank_aggr_coexpr")]
tf_mm <- rank_mm[[str_to_title(plot_tf)]][, c("Symbol", "Rank_bind", "Rank_aggr_coexpr")]


# For species-specific interactions, join the ranking for the other species 
# (which would have not made the cut, but want to visualize as well)

rank_l$Human_specific <- rank_l$Human_specific %>%
  dplyr::rename(Rank_aggr_coexpr_hg = Rank_aggr_coexpr, 
                Rank_bind_hg = Rank_bind) %>%
  mutate(Symbol_hg = Symbol) %>% 
  left_join(pc_ortho, by = "Symbol_hg") %>% 
  left_join(tf_mm, by = c("Symbol_mm" = "Symbol")) %>%
  dplyr::rename(Rank_aggr_coexpr_mm = Rank_aggr_coexpr, Rank_bind_mm = Rank_bind)


# For mouse, if human symbol doesn't exist, replace with mouse symbol
rank_l$Mouse_specific <- rank_l$Mouse_specific %>%
  dplyr::rename(Rank_aggr_coexpr_mm = Rank_aggr_coexpr, 
                Rank_bind_mm = Rank_bind) %>%
  mutate(Symbol_mm = Symbol) %>% 
  left_join(pc_ortho, by = "Symbol_mm") %>% 
  left_join(tf_hg,by = c("Symbol_hg" = "Symbol")) %>%
  dplyr::rename(Rank_aggr_coexpr_hg = Rank_aggr_coexpr,
                Rank_bind_hg = Rank_bind) %>% 
  mutate(Symbol_hg = ifelse(is.na(Symbol_hg), Symbol_mm, Symbol_hg))


# Remove Mouse/Human tiers, which contain interactions also found in elevated/stringent
rank_l <- rank_l[setdiff(names(rank_l), c("Human", "Mouse"))]



# Making a dataframe tracking discrete status of k-cutoff. Also including a 
# group of 'near' misses (501-1000)

build_plotdf <- function(rank_l, k) {
  
  stat_cols <- c("Rank_aggr_coexpr_hg", "Rank_aggr_coexpr_mm", 
                 "Rank_bind_hg", "Rank_bind_mm")
  
  pdf <- do.call(rbind, lapply(rank_l, `[`, stat_cols))
  
  pdf <- mutate(
    
    pdf,
    
    Rank_aggr_coexpr_hg = case_when(
      Rank_aggr_coexpr_hg <= k ~ 0,
      Rank_aggr_coexpr_hg > k & Rank_aggr_coexpr_hg <= k + 500 ~ 1,
      is.na(Rank_aggr_coexpr_hg) ~ NA_real_,
      TRUE ~ 2
    ),
    
    Rank_aggr_coexpr_mm = case_when(
      Rank_aggr_coexpr_mm <= k ~ 0,
      Rank_aggr_coexpr_mm > k & Rank_aggr_coexpr_mm <= k + 500 ~ 1,
      is.na(Rank_aggr_coexpr_mm) ~ NA_real_,
      TRUE ~ 2), 
    
    Rank_bind_hg = case_when(
      Rank_bind_hg <= k ~ 0,
      Rank_bind_hg > k & Rank_bind_hg <= k + 500 ~ 1,
      is.na(Rank_bind_hg) ~ NA_real_,
      TRUE ~ 2),
    
    Rank_bind_mm = case_when(
      Rank_bind_mm <= k ~ 0,
      Rank_bind_mm > k & Rank_bind_mm <= k + 500 ~ 1,
      is.na(Rank_bind_mm) ~ NA_real_,
      TRUE ~ 2)
    )
  
  
  rownames(pdf) <- unlist(lapply(rank_l, pluck, "Symbol_hg"))
  return(pdf)
}


pdf <- build_plotdf(rank_l, k)




# Padding between species and genes in evidence tiers 

gap_genes <- rep(
  head(cumsum(unlist(lapply(rank_l, nrow))), -1),
  each = 4)

# gap_evidence <- rep(1:4, each = 4)
gap_evidence <- c(rep(1, 4), rep(2, 8), rep(3, 4), rep(4, 2))


# Heatmap of k-cutoff status, with data types/species labels
pheatmap(t(pdf),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         na_col = "grey",
         gaps_col = gap_genes,
         gaps_row = gap_evidence,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90,
         width = 14,
         filename = file.path(plot_dir, paste0(plot_tf, "_tiered_evidence_targets_heatmap_label.png"))
)


# Without data types/species labels (for matching curated heatmap in illustrator)
pheatmap(t(pdf),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         na_col = "grey",
         gaps_col = gap_genes,
         gaps_row = gap_evidence,
         show_rownames = FALSE,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90,
         width = 14,
         filename = file.path(plot_dir, paste0(plot_tf, "_tiered_evidence_targets_heatmap_nolabel.png"))
)


# Binary heatmap of curation status

labels_curated <- get_curated_labels(tf = plot_tf, 
                                     curated_df = curated, 
                                     ortho_df = pc_ortho,
                                     pc_df = pc_hg, 
                                     species = "Human", 
                                     remove_self = TRUE)


curated_vec <- setNames(as.integer(rownames(pdf) %in% labels_curated), rownames(pdf))


# With gene symbols
pheatmap(t(curated_vec),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         border_color = "black",
         gaps_col = gap_genes,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90,
         width = 14,
         filename = file.path(plot_dir, paste0(plot_tf, "_curation_heatmap_label.png"))
)


# Without gene symbols (for matching tiered evidence heatmap in illustrator)
pheatmap(t(curated_vec),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         border_color = "black",
         gaps_col = gap_genes,
         show_colnames = FALSE,
         cellheight = 10,
         legend = FALSE,
         angle_col = 90,
         width = 14,
         filename = file.path(plot_dir, paste0(plot_tf, "_curation_heatmap_nolabel.png"))
)



# Saving a histogram of ASCL1 aggr coexpr values for schematic in Supp Fig1

p_hist <- ggplot(rank_hg[[plot_tf]], aes(x = Avg_aggr_coexpr)) + 
  geom_histogram(bins = 50, fill = "black", col = "black") +
  ylab("Count of genes") +
  xlab("Average coexpression profile") +
  theme_classic() +
  theme(text = element_text(size = 25))


ggsave(p_hist, height = 6, width = 10, device = "png", dpi = 300, bg = "white",
       filename = file.path(plot_dir, "demo_ASCL1_average_coexpr_profile.png"))
