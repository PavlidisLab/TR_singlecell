## Save out a list of the AUC performances of aggregated ChIP-seq's ability
## to recover curated targets
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
collection <- "Permissive"

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

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Average bind scores
stopifnot(collection %in% c("Robust", "Permissive"))
bind_summary_path <- paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_bindscore_summary.RDS")
bind_summary <- readRDS(bind_summary_path)

# TODO: finalize
unibind_auc_hg_path <- paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_recover_curated_hg.RDS")
unibind_auc_mm_path <- paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_recover_curated_mm.RDS")


# Functions
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



# Get ortho-matched symbols of TFs with available data, as well as all targets
# which are used for null
# ------------------------------------------------------------------------------


# Human

ortho_tf_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_hg %in% colnames(bind_summary$Human_TF)) %>% 
  pull(Symbol_hg)

ortho_target_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_hg %in% rownames(bind_summary$Human_TF)) %>% 
  pull(Symbol_hg)

tf_hg <- union(
  intersect(colnames(bind_summary$Human_TF), str_to_upper(curated$TF_Symbol)),
  ortho_tf_hg)

target_hg <- union(
  intersect(rownames(bind_summary$Human_TF), str_to_upper(curated$Target_Symbol)),
  ortho_target_hg)


# Mouse

ortho_tf_mm <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_mm %in% colnames(bind_summary$Mouse_TF)) %>% 
  pull(Symbol_mm)

ortho_target_mm <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_mm %in% rownames(bind_summary$Mouse_TF)) %>% 
  pull(Symbol_mm)

tf_mm <- union(
  intersect(colnames(bind_summary$Mouse_TF), str_to_title(curated$TF_Symbol)),
  ortho_tf_mm)

target_mm <- union(
  intersect(rownames(bind_summary$Mouse_TF), str_to_title(curated$Target_Symbol)),
  ortho_target_mm)



# The following generates and saves a list containing, for each TF, a summary
# dataframe of AUC performance, as well as another list of the null AUCs
# ------------------------------------------------------------------------------


set.seed(5)


# Human 

save_curated_auc_list(path = unibind_auc_hg_path,
                      tfs = tf_hg,
                      rank_l = bind_hg,
                      score_col = "Bind_score",
                      curated_df = curated,
                      label_all = target_hg,
                      ortho_df = pc_ortho,
                      pc_df = pc_hg,
                      species = "Human",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)


# Mouse 

save_curated_auc_list(path = unibind_auc_mm_path,
                      tfs = tf_mm,
                      rank_l = bind_mm,
                      score_col = "Bind_score",
                      curated_df = curated,
                      label_all = target_mm,
                      ortho_df = pc_ortho,
                      pc_df = pc_mm,
                      species = "Mouse",
                      n_samps = 1000,
                      ncores = 8,
                      verbose = TRUE,
                      force_resave = TRUE)
