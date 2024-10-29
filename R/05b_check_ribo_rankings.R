## Checking that the top ranked genes for each ribo gene ranking are other ribos
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Ribosomal genes
ribo_genes <- read.delim(ribo_path, stringsAsFactors = FALSE)

# Saved list RDS of the ranks
rank_ribo_hg <- readRDS(rank_ribo_hg_path)
rank_ribo_mm <- readRDS(rank_ribo_mm_path)



# Ribosomal rankings. Looking at top 82 (because there are 82 ortho L/S ribo 
# genes) coexpressed genes with each ribo gene
# ------------------------------------------------------------------------------


# Get the top ribo count for each ribo gene

count_top_ribo <- function(ribo_l, topn = 82) {
  
  ribo_genes <- names(ribo_l)
  
  sum_l <- lapply(ribo_genes, function(x) {
    sum(arrange(ribo_l[[x]], Rank_aggr_coexpr)$Symbol[1:topn] %in% ribo_genes)
  })
  
  return(data.frame(Symbol = ribo_genes, Count = unlist(sum_l)))
}



# Across all ribo genes, tally the top non-ribo genes

top_non_ribo <- function(rank_l, ribo_genes, topn = 82) {
  
  non_ribo <- lapply(rank_l, function(x) {
    setdiff(arrange(x, Rank_aggr_coexpr)$Symbol[1:topn], ribo_genes)
  })
  
  non_ribo_counts <- sort(table(unlist(non_ribo)))
  return(non_ribo_counts)
}



# Find a median of 71/82 genes are other L/S ortho ribo genes

count_ribo_hg <- count_top_ribo(rank_ribo_hg)
count_ribo_mm <- count_top_ribo(rank_ribo_mm)


summ_hg <- summary(count_ribo_hg$Count)
summ_mm <- summary(count_ribo_mm$Count)


# Most common "non-ribo" genes are actually also ribo, just outside the set
# of 82 L/S ortho ribo genes

non_ribo_hg <- top_non_ribo(rank_ribo_hg, ribo_genes$Symbol_hg)
non_ribo_mm <- top_non_ribo(rank_ribo_mm, ribo_genes$Symbol_mm)



# Plots 
# ------------------------------------------------------------------------------


# Barchart of count of ribo top-ranked partners that are also ribo

ribo_barchart <- function(ribo_df, title) {
  
  ggplot(ribo_df, aes(x = reorder(Symbol, Count), y = Count)) +
    geom_bar(stat = "identity", colour = "black", fill = "slategrey") +
    ylab("Count of top 82 coexpression partners that are also ribosomal") +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title = element_text(size = 25),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 25),
          plot.margin = margin(c(10, 10, 10, 10)))
  
}


p1a <- ribo_barchart(count_ribo_hg, "Human")

p1b <- ribo_barchart(count_ribo_mm, "Mouse")

p1 <- plot_grid(p1a, p1b, ncol = 1)

ggsave(p1a, height = 12, width = 18, device = "png", dpi = 300,
       filename = file.path(plot_dir, "ribosomal_top_rank_overlap_human.png"))


ggsave(p1b, height = 12, width = 18, device = "png", dpi = 300,
       filename = file.path(plot_dir, "ribosomal_top_rank_overlap_mouse.png"))



# Heatmap of binarized topk status

stopifnot(identical(names(rank_ribo_hg), ribo_genes$Symbol_hg))

ribo_mat_hg <- 
  lapply(rank_ribo_hg, function(x) x[ribo_genes$Symbol_hg, "Rank_aggr_coexpr"]) %>% 
  do.call(cbind, .)

colnames(ribo_mat_hg) <- rownames(ribo_mat_hg) <- ribo_genes$Symbol_hg

ribo_mat_hg <- (ribo_mat_hg <= 200) * 1

pheatmap(ribo_mat_hg,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = c("white", "black"),
         na_col = "red",
         border_color = NA)
