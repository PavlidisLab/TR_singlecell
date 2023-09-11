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

# TODO: replace when done
# Saved list RDS of the summaries
# rank_tf_ln_hg <- readRDS(rank_tf_ln_hg_path)  
rank_tf_ln_hg <- readRDS(rank_tf_hg_path)
# rank_tf_ln_mm <- readRDS(rank_tf_ln_mm_path)

# rank_tf_cpm_hg <- readRDS(rank_tf_cpm_hg_path)
rank_tf_cpm_hg <- rank_tf_ln_hg
# rank_tf_cpm_mm <- readRDS(rank_tf_cpm_mm_path)



stopifnot(identical(rank_tf_ln_hg$PAX6$Symbol, rank_tf_cpm_hg$PAX6$Symbol))
stopifnot(identical(rank_tf_ln_mm$Pax6$Symbol, rank_tf_cpm_hg$Pax6$Symbol))

stopifnot(identical(names(rank_tf_ln_hg), names(rank_tf_cpm_hg)))
stopifnot(identical(names(rank_tf_ln_mm), names(rank_tf_cpm_mm)))

tfs_hg <- names(rank_tf_ln_hg)
tfs_mm <- names(rank_tf_ln_mm)



# For each TF ranking, compare Spearman cor and topk between lognorm and CPM
# ------------------------------------------------------------------------------


tf <- "ZNF845"
stat_col <- "Avg_RSR"

rank_l1 <- rank_tf_cpm_hg
rank_l2 <- rank_tf_ln_hg
vec1 <- sort(setNames(rank_l1[[tf]][[stat_col]], rank_l1[[tf]]$Symbol), decreasing = TRUE)
vec2 <- sort(setNames(rank_l2[[tf]][[stat_col]], rank_l2[[tf]]$Symbol), decreasing = TRUE)





top_comparison <- function(rank_l1, rank_l2, stat_col = "Avg_RSR", ncores = 1) {
  
  tfs <- intersect(names(rank_l1), names(rank_l2))
  
  top_l <- mclapply(tfs, function(tf) {

    message(tf)
    
    vec1 <- setNames(rank_l1[[tf]][[stat_col]], rank_l1[[tf]]$Symbol)
    vec2 <- setNames(rank_l2[[tf]][[stat_col]], rank_l2[[tf]]$Symbol)
    
    topk <- topk_intersect(
      topk_sort(vec1, k = k, check_k_arg = TRUE),
      topk_sort(vec2, k = k, check_k_arg = TRUE)
    ) 
    
    scor <- cor(rank_l1[[tf]][[stat_col]],
                rank_l2[[tf]][[stat_col]],
                method = "spearman",
                use = "pairwise.complete.obs")
    
    data.frame(Symbol = tf, Topk = topk, Scor = scor)
    
  }, mc.cores = ncores)

  top_df <- as.data.frame(do.call(rbind, top_l))
  return(top_df)
}




top_hg <- top_comparison(rank_l1 = rank_tf_cpm_hg, rank_l2 = rank_tf_ln_hg, ncores = 8)
top_mm <- top_comparison(rank_l1 = rank_tf_cpm_mm, rank_l2 = rank_tf_ln_mm, ncores = 8)



# Inspecting the best/worst agreements


arrange(top_hg, Scor) %>% head
arrange(top_mm, Scor) %>% head

summary(Filter(is.numeric, top_hg))
summary(Filter(is.numeric, top_mm))

ggplot(top_hg, aes(x = Scor, y = Topk)) +
  geom_point(shape = 21) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))
  

# Load best and worst tfs of note


sub_hg <- c("PAX6")
sub_mm <- c("Pax6")

agg_ln_hg <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = sub_hg)
agg_ln_mm <- load_agg_mat_list(ids = ids_mm, genes = pc_hg$Symbol, sub_genes = sub_mm)

agg_cpm_hg <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = sub_hg, pattern = "_RSR_allrank_CPM.tsv")
agg_cpm_mm <- load_agg_mat_list(ids = ids_hg, genes = pc_hg$Symbol, sub_genes = sub_hg, pattern = "_RSR_allrank_CPM.tsv")

