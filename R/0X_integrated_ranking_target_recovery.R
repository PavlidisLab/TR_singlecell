## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/00_config.R")

n_samps <- 1000
set.seed(5)
force_resave <- TRUE

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
rank_tf_hg <- readRDS(rank_int_hg_path)
rank_tf_mm <- readRDS(rank_int_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)


int_auc_hg_path <- "/space/scratch/amorin/R_objects/integrated_recover_curated_hg.RDS"
int_auc_mm_path <- "/space/scratch/amorin/R_objects/integrated_recover_curated_mm.RDS"



# Get ortho-matched symbols of TFs with available data, as well as all ortho
# targets, which are used for null
# ------------------------------------------------------------------------------


# Human

ortho_tf_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_hg %in% names(rank_tf_hg)) %>% 
  pull(Symbol_hg)

ortho_target_hg <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_hg %in% rownames(rank_tf_hg[[1]])) %>% 
  pull(Symbol_hg)

tf_hg <- union(
  intersect(names(rank_tf_hg), str_to_upper(curated$TF_Symbol)),
  ortho_tf_hg)

target_hg <- union(
  intersect(rank_tf_hg[[1]]$Symbol, str_to_upper(curated$Target_Symbol)),
  ortho_target_hg)


# Mouse

ortho_tf_mm <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$TF_Symbol | Symbol_hg %in% curated$TF_Symbol) %>% 
  filter(Symbol_mm %in% names(rank_tf_mm)) %>% 
  pull(Symbol_mm)

ortho_target_mm <- pc_ortho %>% 
  filter(Symbol_mm %in% curated$Target_Symbol | Symbol_hg %in% curated$Target_Symbol) %>% 
  filter(Symbol_mm %in% rownames(rank_tf_mm[[1]])) %>% 
  pull(Symbol_mm)

tf_mm <- union(
  intersect(names(rank_tf_mm), str_to_title(curated$TF_Symbol)),
  ortho_tf_mm)

target_mm <- union(
  intersect(rank_tf_mm[[1]]$Symbol, str_to_title(curated$Target_Symbol)),
  ortho_target_mm)



# Integrated ranking (coexpr and binding)
# ------------------------------------------------------------------------------


# Reverse integrated rank column: process assumes higher values are more important

rank_tf_hg <- lapply(rank_tf_hg, function(x) mutate(x, Rank_integrated = -Rank_integrated))
rank_tf_mm <- lapply(rank_tf_mm, function(x) mutate(x, Rank_integrated = -Rank_integrated))


# Human
save_function_results(
  path = int_auc_hg_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
    tfs = tf_hg,
    rank_l = rank_tf_hg,
    score_col = "Rank_integrated",
    curated_df = curated,
    label_all = target_hg,
    ortho_df = pc_ortho,
    pc_df = pc_hg,
    species = "Human",
    n_samps = n_samps,
    ncores = ncore,
    verbose = TRUE
  ),
  force_resave = force_resave
)



# Mouse
save_function_results(
  path = int_auc_mm_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
    tfs = tf_mm,
    rank_l = rank_tf_mm,
    score_col = "Rank_integrated",
    curated_df = curated,
    label_all = target_mm,
    ortho_df = pc_ortho,
    pc_df = pc_mm,
    species = "Mouse",
    n_samps = n_samps,
    ncores = ncore,
    verbose = TRUE
  ),
  force_resave = force_resave
)
