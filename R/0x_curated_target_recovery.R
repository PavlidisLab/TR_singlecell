## TODO:
## -----------------------------------------------------------------------------

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

n_samps <- 1000

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

# For each dataset, load the subset gene x TF aggregation matrix 
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

# Saved list RDS of the summaries
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)


# Curated TFs and targets

tfs_curated_hg <- intersect(tfs_hg$Symbol, str_to_upper(curated$TF_Symbol))
tfs_curated_mm <- intersect(tfs_mm$Symbol, str_to_title(curated$TF_Symbol))

targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))



# tf <- "SPI1"
# curated_df <- curated
# rank_l <- rank_tf_hg
# pc_df <- pc_hg
# measure <- "AUPRC"
# species <- "Human"
# labels_all <- targets_curated_hg
# 


get_null_performance <- function(labels_all,
                                 measure,
                                 vec_scores,
                                 n_target, 
                                 n_samps) {
  
  null_perf <- vapply(1:n_samps, function(x) {
    null_labels <- sample(labels_all, n_target, replace = FALSE)
    vec_auc(vec_scores = vec_scores, vec_labels = null_labels, measure = measure)
  }, numeric(1))
  
  return(null_perf)
}



get_performance <- function(tf, 
                            curated_df, 
                            labels_all,
                            rank_l, 
                            pc_df,
                            measure,
                            species,
                            n_samps) {
  
 
  curated_tf <- curated_df %>%
    filter(str_to_upper(TF_Symbol) == str_to_upper(tf) &
             !(str_to_upper(Target_Symbol) == str_to_upper(tf)))
  
  
  if (species == "Human") {
    labels <- unique(str_to_upper(curated_tf$Target_Symbol))
  } else if (species == "Mouse") {
    labels <- unique(str_to_title(curated_tf$Target_Symbol))
  }
    
  labels <- labels[labels %in% pc_df$Symbol]
  n_target <- length(labels)
  
  rank <- rank_l[[tf]] %>% 
    filter(Symbol != tf) %>% 
    arrange(desc(Avg_RSR))
  
  scores <- setNames(rank$Avg_RSR, rank$Symbol)
  
  auc <- vec_auc(vec_scores = scores, vec_labels = labels, measure = measure)
  
  null <- get_null_performance(labels_all = labels_all,
                               measure = measure,
                               vec_scores = scores,
                               n_target = n_target, 
                               n_samps = n_samps)
  
  df <- data.frame(
    Symbol = tf,
    N_targets = n_target,
    AUC = auc,
    Percentile_observed = ecdf(null)(auc),
    Ratio_observed = auc / median(null),
    Diff_observed = auc - median(null)
  )
  
  return(list(Perf_df = df, Labels = labels, Null = null)) 
  
}



save_all_performance <- function(path,
                                 tfs, 
                                 curated_df, 
                                 labels_all,
                                 rank_l, 
                                 pc_df,
                                 measure,
                                 species,
                                 n_samps,
                                 force_resave = FALSE)  {
  
  if (!file.exists(path) || force_resave) {
    
    perf_l <- lapply(tfs, function(tf) {
      
      message(paste(tf, Sys.time()))
      
      try(get_performance(tf = tf,
                          curated_df = curated_df,
                          labels_all = labels_all,
                          rank_l = rank_l,
                          pc_df = pc_df,
                          measure = measure,
                          species = species,
                          n_samps = n_samps))
      
    })  

    names(perf_l) <- tfs
    saveRDS(perf_l, path)
    return(invisible(NULL))
  }
}



curated_auprc_hg_path <- "/space/scratch/amorin/R_objects/curated_auprc_hg.RDS"
curated_auprc_mm_path <- "/space/scratch/amorin/R_objects/curated_auprc_mm.RDS"
curated_auroc_hg_path <- "/space/scratch/amorin/R_objects/curated_auroc_hg.RDS"
curated_auroc_mm_path <- "/space/scratch/amorin/R_objects/curated_auroc_mm.RDS"



# Human AURPC

save_all_performance(path = curated_auprc_hg_path,
                     tfs = tfs_curated_hg, 
                     curated_df = curated, 
                     labels_all = targets_curated_hg,
                     rank_l = rank_tf_hg, 
                     pc_df = pc_hg,
                     measure = "AUPRC",
                     species = "Human",
                     n_samps = n_samps,
                     force_resave = FALSE)


# Mouse AUPRC

save_all_performance(path = curated_auprc_mm_path,
                     tfs = tfs_curated_mm, 
                     curated_df = curated, 
                     labels_all = targets_curated_mm,
                     rank_l = rank_tf_mm, 
                     pc_df = pc_mm,
                     measure = "AUPRC",
                     species = "Mouse",
                     n_samps = n_samps,
                     force_resave = FALSE)



# Human AUROC

save_all_performance(path = curated_auroc_hg_path,
                     tfs = tfs_curated_hg, 
                     curated_df = curated, 
                     labels_all = targets_curated_hg,
                     rank_l = rank_tf_hg, 
                     pc_df = pc_hg,
                     measure = "AUROC",
                     species = "Human",
                     n_samps = n_samps,
                     force_resave = FALSE)


# Mouse AROC

save_all_performance(path = curated_auroc_mm_path,
                     tfs = tfs_curated_mm, 
                     curated_df = curated, 
                     labels_all = targets_curated_mm,
                     rank_l = rank_tf_mm, 
                     pc_df = pc_mm,
                     measure = "AUROC",
                     species = "Mouse",
                     n_samps = n_samps,
                     force_resave = FALSE)



# Demo a single TF


# tf_hg <- "EGR1"
# tf_mm <- "Egr1"
# 
# 
# perf_hg <- get_performance(tf = tf_hg,
#                            curated_df = curated,
#                            labels_all = targets_curated_hg,
#                            rank_l = rank_tf_hg,
#                            pc_df = pc_hg,
#                            measure = "AUPRC",
#                            species = "Human",
#                            n_samps = n_samps)
# 
# 
# perf_mm <- get_performance(tf = tf_mm,
#                            curated_df = curated,
#                            labels_all = targets_curated_mm,
#                            rank_l = rank_tf_mm,
#                            pc_df = pc_mm,
#                            measure = "AUPRC",
#                            species = "Mouse",
#                            n_samps = n_samps)
# 
# 
# 
# summary(perf_hg$Null)
# plot(density(perf_hg$Null), xlim = c(min(perf_hg$Null) * 0.9, perf_hg$Perf_df$AUC * 1.1))
# abline(v = perf_hg$Perf_df$AUC, col = "red")
# 
# 
# summary(perf_mm$Null)
# plot(density(perf_mm$Null), xlim = c(min(perf_mm$Null) * 0.9, perf_mm$Perf_df$AUC * 1.1))
# abline(v = perf_mm$Perf_df$AUC, col = "red")
