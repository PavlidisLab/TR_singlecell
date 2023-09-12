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

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# TODO: formalize unibind loading
# Processed list of meta and matrices
bind_dat_path <- "/space/scratch/amorin/R_objects/processed_unibind_data.RDS"
bind_dat <- readRDS(bind_dat_path)

# Average bind scores and output of binding specificity model
bind_summary_path <- "/space/scratch/amorin/R_objects/unibind_bindscore_summary.RDS"
bind_model_path <- "/space/scratch/amorin/R_objects/unibind_bindscore_modelfit.RDS"
bind_summary <- readRDS(bind_summary_path)
bind_model <- readRDS(bind_model_path)

# Curated TFs with ChIP-seq and all targets for null
tfs_curated_hg <- intersect(colnames(bind_summary$Human_TF), str_to_upper(curated$TF_Symbol))
tfs_curated_mm <- intersect(colnames(bind_summary$Mouse_TF), str_to_upper(curated$TF_Symbol))
targets_curated_hg <- intersect(pc_hg$Symbol, str_to_upper(curated$Target_Symbol))
targets_curated_mm <- intersect(pc_mm$Symbol, str_to_title(curated$Target_Symbol))


# TODO:
# ------------------------------------------------------------------------------


# Split matrix of summarized bind scores col/TF-wise into a list


split_to_list <- function(mat) {
  
  bind_l <- asplit(mat, 2)
  
  bind_l <- lapply(bind_l, function(x) {
    data.frame(Bind_score = x) %>% 
      rownames_to_column(var = "Symbol") %>% 
      arrange(desc(Bind_score))
  })
}


bind_hg <- split_to_list(bind_summary$Human_TF)
bind_mm <- split_to_list(bind_summary$Mouse_TF)


# TODO: review if necessary/can be avoided
names(bind_mm) <- str_to_title(names(bind_mm))
tfs_curated_mm <- str_to_title(tfs_curated_mm)


set.seed(5)


# Human 

save_curated_auc_list(path = unibind_recover_curated_hg_path,
                      tfs = tfs_curated_hg,
                      rank_l = bind_hg,
                      score_col = "Bind_score",
                      curated_df = curated,
                      label_all = targets_curated_hg,
                      pc_df = pc_hg,
                      species = "Human",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)


# Mouse 

save_curated_auc_list(path = unibind_recover_curated_mm_path,
                      tfs = tfs_curated_mm,
                      rank_l = bind_mm,
                      score_col = "Bind_score",
                      curated_df = curated,
                      label_all = targets_curated_mm,
                      pc_df = pc_mm,
                      species = "Mouse",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)
