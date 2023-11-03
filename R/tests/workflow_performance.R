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

# Curated TFs with ChIP-seq and all targets for null
tfs_curated_hg <- intersect(tfs_hg$Symbol, str_to_upper(curated$TF_Symbol))
tfs_curated_mm <- intersect(tfs_mm$Symbol, str_to_title(curated$TF_Symbol))
targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))

# Rankings from Morin 2023
evidence_l <- readRDS(evidence_path)


# TODO: formalize unibind loading
# Processed list of meta and matrices
# bind_dat_path <- "/space/scratch/amorin/R_objects/processed_unibind_data.RDS"
# bind_dat <- readRDS(bind_dat_path)

# Average bind scores and output of binding specificity model
# bind_summary_path <- "/space/scratch/amorin/R_objects/unibind_bindscore_summary.RDS"
# bind_model_path <- "/space/scratch/amorin/R_objects/unibind_bindscore_modelfit.RDS"
# bind_summary <- readRDS(bind_summary_path)
# bind_model <- readRDS(bind_model_path)




# Inspecting retrieved ranks for a given TF


agg_l <- agg_tf_hg
msr_mat <- msr_hg
score_mat <- gene_vec_to_mat(agg_l, tf)
score_mat <- subset_to_measured(score_mat, msr_mat = msr_mat, gene = tf)
score_mat <- cbind(score_mat, Average = rowMeans(score_mat))
labels_curated <- get_curated_labels(tf = tf, curated_df = curated, pc_df = pc_hg, species = "Human", remove_self = TRUE)


rank_l <- lapply(1:ncol(score_mat), function(x) {
  score_rank <- rank(-score_mat[, x], ties.method = "min")
  # score_rank[labels_curated]
  # median(score_rank[labels_curated])
  # mean(score_rank[labels_curated])
  # sort(score_rank[labels_curated])[1:10]
  median(sort(score_rank[labels_curated])[1:10])
})
names(rank_l) <- colnames(score_mat)

rank_l[auc_df$ID]