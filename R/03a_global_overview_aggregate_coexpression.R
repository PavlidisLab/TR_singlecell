## This script averages and sorts all aggregate coexpression networks for the
## purpose of identifying the max pos and neg correlation gene pairs.
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

force_resave <- FALSE

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)



# Get the average aggregate coexpression across all datasets, returned as a
# data frame. Slow due to having to load every gene x gene matrix.
# ------------------------------------------------------------------------------


gen_avg_coexpr <- function(ids, genes) {
  
  avg_mat <- matrix(0, nrow = length(genes), ncol = length(genes))
  colnames(avg_mat) <- rownames(avg_mat) <- genes
  
  for (id in ids) {
    mat <- load_agg_mat_list(id, genes = genes)[[1]]
    avg_mat <- avg_mat + mat
  }
  
  avg_mat <- avg_mat / length(ids)
  
  avg_df <- mat_to_df(avg_mat, symmetric = TRUE, value_name = "Avg_coexpr")
  avg_df <- arrange(avg_df, desc(Avg_coexpr))
  
  return(avg_df)
  
}


# Note that even loading is quite slow.


if (!file.exists(avg_coexpr_hg_path) || force_resave) {
  avg_coexpr_hg <- gen_avg_coexpr(ids_hg, genes = pc_hg$Symbol)
  saveRDS(avg_coexpr_hg, avg_coexpr_hg_path)
} else {
  avg_coexpr_hg <- readRDS(avg_coexpr_hg_path)
}



if (!file.exists(avg_coexpr_mm_path) || force_resave) {
  avg_coexpr_mm <- gen_avg_coexpr(ids_mm, genes = pc_mm$Symbol)
  saveRDS(avg_coexpr_mm, avg_coexpr_mm_path)
} else {
  avg_coexpr_mm <- readRDS(avg_coexpr_mm_path)
}



topn <- 20


extremes <- list(
  Top_human = head(avg_coexpr_hg, topn),
  Bottom_human = tail(avg_coexpr_hg, topn),
  Top_mouse = head(avg_coexpr_mm, topn),
  Bottom_mouse = tail(avg_coexpr_mm, topn)
)
