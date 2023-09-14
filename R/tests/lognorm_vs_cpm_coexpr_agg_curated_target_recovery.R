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

# Saved list RDS of the summaries

rank_tf_ln_hg_path <- "/space/scratch/amorin/R_objects/ranking_TF_human_lognorm.RDS"
rank_tf_cpm_hg_path <- "/space/scratch/amorin/R_objects/ranking_TF_human_CPM.RDS"
rank_tf_ln_mm_path <- "/space/scratch/amorin/R_objects/ranking_TF_mouse_lognorm.RDS"
rank_tf_cpm_mm_path <- "/space/scratch/amorin/R_objects/ranking_TF_mouse_CPM.RDS"



rank_tf_ln_hg <- readRDS(rank_tf_ln_hg_path)
rank_tf_ln_mm <- readRDS(rank_tf_ln_mm_path)
rank_tf_cpm_hg <- readRDS(rank_tf_cpm_hg_path)
rank_tf_cpm_mm <- readRDS(rank_tf_cpm_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Curated TFs with ChIP-seq and all targets for null
tfs_curated_hg <- intersect(tfs_hg$Symbol, str_to_upper(curated$TF_Symbol))
tfs_curated_mm <- intersect(tfs_mm$Symbol, str_to_title(curated$TF_Symbol))
targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))



# TODO:
# ------------------------------------------------------------------------------



coexpr_recover_curated_ln_hg_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_ln_hg.RDS"
coexpr_recover_curated_ln_mm_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_ln_mm.RDS"

coexpr_recover_curated_cpm_hg_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_cpm_hg.RDS"
coexpr_recover_curated_cpm_mm_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_cpm_mm.RDS"



set.seed(5)


# Human 

save_curated_auc_list(path = coexpr_recover_curated_ln_hg_path,
                      tfs = tfs_curated_hg,
                      rank_l = rank_tf_ln_hg,
                      score_col = "Avg_RSR",
                      curated_df = curated,
                      label_all = targets_curated_hg,
                      pc_df = pc_hg,
                      species = "Human",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)



save_curated_auc_list(path = coexpr_recover_curated_cpm_hg_path,
                      tfs = tfs_curated_hg,
                      rank_l = rank_tf_cpm_hg,
                      score_col = "Avg_RSR",
                      curated_df = curated,
                      label_all = targets_curated_hg,
                      pc_df = pc_hg,
                      species = "Human",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)


# Mouse 

save_curated_auc_list(path = coexpr_recover_curated_ln_mm_path,
                      tfs = tfs_curated_mm,
                      rank_l = rank_tf_ln_mm,
                      score_col = "Avg_RSR",
                      curated_df = curated,
                      label_all = targets_curated_mm,
                      pc_df = pc_mm,
                      species = "Mouse",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)



save_curated_auc_list(path = coexpr_recover_curated_cpm_mm_path,
                      tfs = tfs_curated_mm,
                      rank_l = rank_tf_cpm_mm,
                      score_col = "Avg_RSR",
                      curated_df = curated,
                      label_all = targets_curated_mm,
                      pc_df = pc_mm,
                      species = "Mouse",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)
