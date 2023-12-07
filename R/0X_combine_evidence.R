## TODO
# TODO: ortho join
# TODO: winner takes all binding evidence?
# TODO: get distn of all average bind scores? look for outliers

## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggrepel)
library(pheatmap)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 500

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Saved list of the aggregate coexpression ranks
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Average bind scores
collection <- "Permissive"
stopifnot(collection %in% c("Robust", "Permissive"))
bind_summary_path <- paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_bindscore_summary.RDS")
bind_summary <- readRDS(bind_summary_path)
bind_dat <- readRDS("/space/scratch/amorin/R_objects/processed_unibind_data.RDS")

# TODO: revisit mouse casing (Unibind mouse upper case)
colnames(bind_summary$Mouse_TF) <- str_to_title(colnames(bind_summary$Mouse_TF))
# stopifnot(all(colnames(bind_summary$Mouse_TF) %in% pc_mm$Symbol))
setdiff(colnames(bind_summary$Mouse_TF), pc_mm$Symbol)
common_hg <- intersect(names(rank_tf_hg), colnames(bind_summary$Human_TF))
common_mm <- intersect(names(rank_tf_mm), colnames(bind_summary$Mouse_TF))

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)


# Functions
# ------------------------------------------------------------------------------


# TODO: consider separate object or single list with reference to common TFs

join_coexpr_bind <- function(coexpr_l, bind_mat) {
  
  common_genes <- intersect(names(coexpr_l), colnames(bind_mat))
  
  join_l <- lapply(common_genes, function(x) {
    
    bind_df <- data.frame(Symbol = rownames(bind_mat), 
                          Bind_score = bind_mat[, x])
    
    join_df <- left_join(coexpr_l[[x]], bind_df, by = "Symbol") %>%
      filter(Symbol != x) %>% 
      mutate(
        Rank_bind = rank(-Bind_score, ties.method = "min"),
        RP = rank(Rank_RSR * Rank_bind)
      ) %>% 
      relocate(Rank_RSR, .before = RP) %>% 
      arrange(RP)
  })
  
  names(join_l) <- common_genes
  return(join_l)
}


# Subset each rank df in rank_l to genes that are in the top K of both coexpr 
# and binding. Check k in case of ties, using the lowest value of k that ensures
# there are no ties in either rank

subset_topk <- function(rank_l, k, ncores = 1) {
  
  topk_l <- mclapply(rank_l, function(x) {
    
    # k <- min(
    #   check_k(sort(x$Bind_score, decreasing = TRUE), k = k),
    #   check_k(sort(x$Avg_RSR, decreasing = TRUE), k = k))
    
    k_bind <- check_k(sort(x$Bind_score, decreasing = TRUE), k = k)
    k_rsr <- check_k(sort(x$Avg_RSR, decreasing = TRUE), k = k)
    
    filter(x, Rank_RSR <= k_rsr & Rank_bind <= k_bind)
    
  }, mc.cores = ncore)
  
  names(topk_l) <- names(rank_l)
  return(topk_l)
}


# TODO:

subset_bottomk <- function(rank_l, k, ncores = 1) {
  
  topk_l <- mclapply(rank_l, function(x) {
    
    # k <- min(
    #   check_k(sort(x$Bind_score, decreasing = TRUE), k = k),
    #   check_k(sort(x$Avg_RSR, decreasing = TRUE), k = k))
    
    k_bind <- check_k(sort(x$Bind_score, decreasing = TRUE), k = k)
    k_rsr <- check_k(sort(x$Avg_RSR, decreasing = TRUE), k = k)
    
    filter(x, Rank_RSR > (max(x$Rank_RSR) - k_rsr) & Rank_bind <= k_bind)
    
  }, mc.cores = ncore)
  
  names(topk_l) <- names(rank_l)
  return(topk_l)
}




# Join aggregate coexpression and binding
# ------------------------------------------------------------------------------


rank_hg <- join_coexpr_bind(rank_tf_hg, bind_summary$Human_TF)
rank_mm <- join_coexpr_bind(rank_tf_mm, bind_summary$Mouse_TF)


# Rank product (RP) can prioritize instances where one rank is near maximal and
# the other is low. Select top K in both data types

topk_hg <- subset_topk(rank_hg, k = k, ncores = ncore)
topk_mm <- subset_topk(rank_mm, k = k, ncores = ncore)

n_topk_hg <- sort(unlist(lapply(topk_hg, nrow)), decreasing = TRUE)
n_topk_mm <- sort(unlist(lapply(topk_mm, nrow)), decreasing = TRUE)


# Bottom of coexpr ranks to look for repression

bottomk_hg <- subset_bottomk(rank_hg, k = k, ncores = ncore)
bottomk_mm <- subset_bottomk(rank_mm, k = k, ncores = ncore)

n_bottomk_hg <- sort(unlist(lapply(bottomk_hg, nrow)), decreasing = TRUE)
n_bottomk_mm <- sort(unlist(lapply(bottomk_mm, nrow)), decreasing = TRUE)


# Look for orthologous topk interactions
# TODO: best to create single ortho list upfront, or just overlap topk objects?
# ------------------------------------------------------------------------------


tf_ortho <- filter(pc_ortho,
                   Symbol_hg %in% names(rank_hg) &
                   Symbol_mm %in% names(rank_mm))


# TODO: this is largely copied from ortho_comparison

sub_cols <- c("Symbol", "Avg_RSR", "Bind_score")


rank_ortho <- mclapply(1:nrow(tf_ortho), function(x) {
  
  df_hg <- left_join(rank_hg[[tf_ortho$Symbol_hg[x]]][, sub_cols], 
                     pc_ortho[, c("Symbol_hg", "ID")],
                     by = c("Symbol" = "Symbol_hg")) %>% 
    filter(!is.na(ID))
  
  
  df_mm <- left_join(rank_mm[[tf_ortho$Symbol_mm[x]]][, sub_cols], 
                     pc_ortho[, c("Symbol_mm", "ID")], 
                     by = c("Symbol" = "Symbol_mm")) %>% 
    filter(!is.na(ID))
  
  
  df_ortho <- left_join(df_hg, df_mm,
                        by = "ID",
                        suffix = c("_hg", "_mm")) %>%
    filter(!is.na(Avg_RSR_hg) & !is.na(Avg_RSR_mm)) %>%
    mutate(
      Rank_RSR_hg = rank(-Avg_RSR_hg, ties.method = "min"),
      Rank_RSR_mm = rank(-Avg_RSR_mm, ties.method = "min"),
      Rank_bind_hg = rank(-Bind_score_hg, ties.method = "min"),
      Rank_bind_mm = rank(-Bind_score_mm, ties.method = "min"),
      RP = rank(log(Rank_RSR_hg + Rank_RSR_mm + Rank_bind_hg + Rank_bind_mm))
    )
}, mc.cores = ncore)

names(rank_ortho) <- tf_ortho$Symbol_hg



# Require top K in both species for each data type
# TODO: topk check

topk_ortho_stringent <- mclapply(rank_ortho, function(x) {
  
  k_bind_hg <- check_k(sort(x$Bind_score_hg, decreasing = TRUE), k = k)
  k_bind_mm <- check_k(sort(x$Bind_score_mm, decreasing = TRUE), k = k)
  k_rsr_hg <- check_k(sort(x$Avg_RSR_hg, decreasing = TRUE), k = k)
  k_rsr_mm <- check_k(sort(x$Avg_RSR_mm, decreasing = TRUE), k = k)
  
  # k <- min(
  #   check_k(sort(x$Bind_score_hg, decreasing = TRUE), k = k),
  #   check_k(sort(x$Bind_score_mm, decreasing = TRUE), k = k),
  #   check_k(sort(x$Avg_RSR_hg, decreasing = TRUE), k = k),
  #   check_k(sort(x$Avg_RSR_mm, decreasing = TRUE), k = k))
  
  filter(x, 
        (Rank_RSR_hg <= k_rsr_hg & Rank_RSR_mm <= k_rsr_mm) & 
        (Rank_bind_hg <= k_bind_hg & Rank_bind_mm <= k_bind_mm))
  
}, mc.cores = ncore)


# Require top K in at least one species for each data type

topk_ortho_relaxed <- mclapply(rank_ortho, function(x) {
  
  k_bind_hg <- check_k(sort(x$Bind_score_hg, decreasing = TRUE), k = k)
  k_bind_mm <- check_k(sort(x$Bind_score_mm, decreasing = TRUE), k = k)
  k_rsr_hg <- check_k(sort(x$Avg_RSR_hg, decreasing = TRUE), k = k)
  k_rsr_mm <- check_k(sort(x$Avg_RSR_mm, decreasing = TRUE), k = k)
  
  filter(x, 
        (Rank_RSR_hg <= k_rsr_hg | Rank_RSR_mm <= k_rsr_mm) & 
        (Rank_bind_hg <= k_bind_hg | Rank_bind_mm <= k_bind_mm))
  
}, mc.cores = ncore)


# Require top K in both for species for one data type, but allow top K in only
# one species for the other data type

topk_ortho_middle <- mclapply(rank_ortho, function(x) {
  
  k_bind_hg <- check_k(sort(x$Bind_score_hg, decreasing = TRUE), k = k)
  k_bind_mm <- check_k(sort(x$Bind_score_mm, decreasing = TRUE), k = k)
  k_rsr_hg <- check_k(sort(x$Avg_RSR_hg, decreasing = TRUE), k = k)
  k_rsr_mm <- check_k(sort(x$Avg_RSR_mm, decreasing = TRUE), k = k)
  
  filter(x,
         ((Rank_RSR_hg <= k_rsr_hg | Rank_RSR_mm <= k_rsr_mm) & 
          (Rank_bind_hg <= k_bind_hg & Rank_bind_mm <= k_bind_mm)
         ) |
         ((Rank_bind_hg <= k_bind_hg | Rank_bind_mm <= k_bind_mm) &
          (Rank_RSR_hg <= k_rsr_hg & Rank_RSR_mm <= k_rsr_mm)
         ))
  
}, mc.cores = ncore)



# Total counts (can have duplicates)

# n_topk_ortho <- data.frame(
#   Symbol = names(rank_ortho),
#   Stringent = unlist(lapply(topk_ortho_stringent, nrow)),
#   Middle = unlist(lapply(topk_ortho_middle, nrow)),
#   Relaxed = unlist(lapply(topk_ortho_relaxed, nrow))
# )

# Counts without duplicates

topk_ortho_middle_dedup <- lapply(names(topk_ortho_middle), function(x) {
  setdiff(topk_ortho_middle[[x]], topk_ortho_stringent[[x]])
})
names(topk_ortho_middle_dedup) <- names(topk_ortho_middle)


topk_ortho_relaxed_dedup <- lapply(names(topk_ortho_relaxed), function(x) {
  setdiff(topk_ortho_relaxed[[x]], topk_ortho_middle[[x]])
})
names(topk_ortho_relaxed_dedup) <- names(topk_ortho_relaxed)


n_topk_ortho <- data.frame(
  Symbol = names(rank_ortho),
  Stringent = unlist(lapply(topk_ortho_stringent, nrow)),
  Middle = unlist(lapply(topk_ortho_middle_dedup, nrow)),
  Relaxed = unlist(lapply(topk_ortho_relaxed_dedup, nrow))
)


plot_df1 <- pivot_longer(n_topk_ortho,
                         cols = c("Stringent", "Middle", "Relaxed"),
                         names_to = "Scheme",
                         values_to = "Count") %>% 
  mutate(Scheme = factor(Scheme, levels = c("Relaxed", "Middle", "Stringent")))
  


p1 <- ggplot(plot_df1, aes(x = reorder(Symbol, Count, FUN = median), y = Count, fill = Scheme, colour = Scheme)) +
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of reproducible interactions") +
  xlab("Transcription factor") +
  ggtitle(paste0("N=", nrow(n_topk_ortho))) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


# Heatmap of tiered evidence for a given TF

plot_tf <- "ASCL1"

scale01 <- function(x) (x - min(x)) / (max(x) - min(x))


genes <- c(arrange(topk_ortho_stringent[[plot_tf]], RP)$Symbol_hg,
           arrange(topk_ortho_middle_dedup[[plot_tf]], RP)$Symbol_hg,
           arrange(topk_ortho_relaxed_dedup[[plot_tf]], RP)$Symbol_hg)


plot_df2 <- rank_ortho[[plot_tf]] %>% 
  mutate(Avg_RSR_hg = scale01(Avg_RSR_hg),
         Avg_RSR_mm = scale01(Avg_RSR_mm),
         Bind_score_hg = scale01(Bind_score_hg),
         Bind_score_mm = scale01(Bind_score_mm),
         Rank_RSR_hg = as.integer(Rank_RSR_hg <= k), 
         Rank_RSR_mm = as.integer(Rank_RSR_mm <= k), 
         Rank_bind_hg = as.integer(Rank_bind_hg <= k), 
         Rank_bind_mm = as.integer(Rank_bind_mm <= k)) %>%  
  filter(Symbol_hg %in% genes) %>% 
  arrange(match(Symbol_hg, genes))

rownames(plot_df2) <- plot_df2$Symbol_hg


# Padding between species (columns) and genes in evidence tiers (rows)

gaps_row <- rep(c(nrow(topk_ortho_stringent[[plot_tf]]),
                  nrow(topk_ortho_middle_dedup[[plot_tf]])),
                each = 4)

gaps_col <- rep(1, 4)


pheatmap(plot_df2[, c("Avg_RSR_hg", "Avg_RSR_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000'),
         breaks = seq(0, 1, length.out = 9),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_coexpr_heatmap.png")
)


pheatmap(plot_df2[, c("Bind_score_hg", "Bind_score_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c('#fff7fb','#ece7f2','#d0d1e6','#a6bddb','#74a9cf','#3690c0','#0570b0','#045a8d','#023858'),
         breaks = seq(0, 1, length.out = 9),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_binding_heatmap.png")
)


# Binarizing status


pheatmap(plot_df2[, c("Rank_RSR_hg", "Rank_RSR_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c('white', '#7f0000'),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_coexpr_binary_heatmap.png")
)


pheatmap(plot_df2[, c("Rank_bind_hg", "Rank_bind_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c('white', '#045a8d'),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_binding_binary_heatmap.png")
)


# Bin status


plot_df2 <- rank_ortho[[plot_tf]] %>% 
  
  mutate(
    
    Avg_RSR_hg = scale01(Avg_RSR_hg),
    Avg_RSR_mm = scale01(Avg_RSR_mm),
    Bind_score_hg = scale01(Bind_score_hg),
    Bind_score_mm = scale01(Bind_score_mm), 
    
    Rank_RSR_hg = case_when(Rank_RSR_hg <= k ~ 0,
                            Rank_RSR_hg > k & Rank_RSR_hg <= k + 500 ~ 1,
                            TRUE ~ 2), 
    
    Rank_RSR_mm = case_when(Rank_RSR_mm <= k ~ 0,
                            Rank_RSR_mm > k & Rank_RSR_mm <= k + 500 ~ 1,
                            TRUE ~ 2),
    
    Rank_bind_hg = case_when(Rank_bind_hg <= k ~ 0,
                             Rank_bind_hg > k & Rank_bind_hg <= k + 500 ~ 1,
                             TRUE ~ 2),
    
    Rank_bind_mm = case_when(Rank_bind_mm <= k ~ 0,
                             Rank_bind_mm > k & Rank_bind_mm <= k + 500 ~ 1,
                             TRUE ~ 2)) %>% 
 
  filter(Symbol_hg %in% genes) %>% 
  arrange(match(Symbol_hg, genes))

rownames(plot_df2) <- plot_df2$Symbol_hg



pheatmap(plot_df2[, c("Rank_RSR_hg", "Rank_RSR_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#fc8d59", '#7f0000')),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_coexpr_step_heatmap.png")
)


pheatmap(plot_df2[, c("Rank_bind_hg", "Rank_bind_mm")],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = rev(c('white', "#a6bddb",'#045a8d')),
         border_color = "black",
         gaps_row = gaps_row,
         gaps_col = gaps_col,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_stringent_binding_step_heatmap.png")
)


# Binary curated status


labels_curated <- get_curated_labels(tf = plot_tf, 
                                     curated_df = curated, 
                                     pc_df = pc_hg, 
                                     species = "Human", 
                                     remove_self = TRUE)

curated_vec <- setNames(as.integer(genes %in% labels_curated), genes)

pheatmap(curated_vec,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         border_color = "black",
         gaps_row = gaps_row,
         cellwidth = 10,
         cellheight = 10,
         filename = file.path(plot_dir, "demo_curated_heatmap.png")
)

