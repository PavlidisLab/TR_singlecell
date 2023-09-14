library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
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


agg_tf_cpm_hg_path <- "/space/scratch/amorin/R_objects/TF_agg_mat_list_human_CPM.RDS"
agg_tf_cpm_mm_path <- "/space/scratch/amorin/R_objects/TF_agg_mat_list_mouse_CPM.RDS"
rank_tf_ln_hg_path <- "/space/scratch/amorin/R_objects/ranking_TF_human_lognorm.RDS"
rank_tf_cpm_hg_path <- "/space/scratch/amorin/R_objects/ranking_TF_human_CPM.RDS"
rank_tf_ln_mm_path <- "/space/scratch/amorin/R_objects/ranking_TF_mouse_lognorm.RDS"
rank_tf_cpm_mm_path <- "/space/scratch/amorin/R_objects/ranking_TF_mouse_CPM.RDS"


# Saved list RDS of the summaries
rank_tf_ln_hg <- readRDS(rank_tf_ln_hg_path)
rank_tf_ln_mm <- readRDS(rank_tf_ln_mm_path)

rank_tf_cpm_hg <- readRDS(rank_tf_cpm_hg_path)
rank_tf_cpm_mm <- readRDS(rank_tf_cpm_mm_path)



stopifnot(identical(rank_tf_ln_hg$PAX6$Symbol, rank_tf_cpm_hg$PAX6$Symbol))
stopifnot(identical(rank_tf_ln_mm$Pax6$Symbol, rank_tf_cpm_mm$Pax6$Symbol))

stopifnot(identical(names(rank_tf_ln_hg), names(rank_tf_cpm_hg)))
stopifnot(identical(names(rank_tf_ln_mm), names(rank_tf_cpm_mm)))

tfs_hg <- names(rank_tf_ln_hg)
tfs_mm <- names(rank_tf_ln_mm)



# For each TF ranking, compare Spearman cor and topk between lognorm and CPM
# ------------------------------------------------------------------------------





top_comparison <- function(rank_l1, 
                           rank_l2, 
                           stat_col = "Avg_RSR", 
                           k = 1000,
                           ncores = 1) {
  
  tfs <- intersect(names(rank_l1), names(rank_l2))
  
  top_l <- mclapply(tfs, function(tf) {

    message(tf)
    
    vec1 <- sort(setNames(rank_l1[[tf]][[stat_col]], rank_l1[[tf]]$Symbol), decreasing = TRUE)
    vec2 <- sort(setNames(rank_l2[[tf]][[stat_col]], rank_l2[[tf]]$Symbol), decreasing = TRUE)
    
    k1 <- check_k(vec1, k)
    k2 <- check_k(vec2, k)
    k_min <- min(k1, k2)
    
    
    topk <- topk_intersect(
      topk_sort(vec1, k = k1, check_k_arg = FALSE),
      topk_sort(vec2, k = k2, check_k_arg = FALSE)
    ) 
    
    scor <- cor(rank_l1[[tf]][[stat_col]],
                rank_l2[[tf]][[stat_col]],
                method = "spearman",
                use = "pairwise.complete.obs")
    
    data.frame(Symbol = tf, 
               Topk = topk, Scor = scor,
               K_min = k_min,
               Topk_prop = topk / k_min)
    
  }, mc.cores = ncores)

  top_df <- as.data.frame(do.call(rbind, top_l))
  return(top_df)
}




top_hg <- top_comparison(rank_l1 = rank_tf_cpm_hg, rank_l2 = rank_tf_ln_hg, ncores = 8)
top_mm <- top_comparison(rank_l1 = rank_tf_cpm_mm, rank_l2 = rank_tf_ln_mm, ncores = 8)



# Inspecting the best/worst agreements


arrange(top_hg, desc(Scor)) %>% head
arrange(top_hg, desc(Topk)) %>% head
arrange(top_mm, desc(Scor)) %>% head
arrange(top_mm, desc(Topk)) %>% head

summary(Filter(is.numeric, top_hg))
summary(Filter(is.numeric, top_mm))



qplot <- function(df, xvar, yvar) {
  
  ggplot(df, aes(x = !!sym(xvar), y = !!sym(yvar))) +
    geom_point(shape = 21) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
}



qplot(top_hg, "Scor", "Topk")
qplot(top_hg, "Scor", "Topk_prop")
qplot(top_hg, "Topk_prop", "K_min")
qplot(top_hg, "Topk", "K_min")


qplot(top_mm, "Scor", "Topk")
qplot(top_mm, "Scor", "Topk_prop")
qplot(top_mm, "Topk_prop", "K_min")
qplot(top_mm, "Topk", "K_min")



gene_hg <- "PAX6"
gene_mm <- "Pax6"


plot_df_hg <- data.frame(
  Symbol = rank_tf_ln_hg[[gene_hg]]$Symbol,
  LN = rank_tf_ln_hg[[gene_hg]]$Avg_RSR, 
  CPM = rank_tf_cpm_hg[[gene_hg]]$Avg_RSR)


plot_df_mm <- data.frame(
  Symbol = rank_tf_ln_mm[[gene_mm]]$Symbol,
  LN = rank_tf_ln_mm[[gene_mm]]$Avg_RSR, 
  CPM = rank_tf_cpm_mm[[gene_mm]]$Avg_RSR)


qplot(plot_df_hg, xvar = "LN", yvar = "CPM")
qplot(plot_df_mm, xvar = "LN", yvar = "CPM")




# Load best and worst tfs of note


sub_hg <- c("PAX6", "EGR1")
sub_mm <- c("Pax6", "Egr1")

agg_ln_hg <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = sub_hg)
agg_ln_mm <- load_agg_mat_list(ids = ids_mm, genes = pc_hg$Symbol, sub_genes = sub_mm)

agg_cpm_hg <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = sub_hg, pattern = "_RSR_allrank_CPM.tsv")
agg_cpm_mm <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = sub_hg, pattern = "_RSR_allrank_CPM.tsv")

