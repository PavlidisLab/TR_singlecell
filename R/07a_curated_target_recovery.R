## Save out lists summarizing the ability of TF-target rankings to recover
## literature curated targets. 4 comparisons are done: individual coexpr profiles
## relative to aggregate coexpr; aggregate positive coexpr profiles relative to
## a sampled null; ditto but for negative coexpr; and Unibind aggregate ChIP-seq
## bind scores relative to a null
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

# For each dataset, load the subset gene x TF aggregation matrix 
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

# Saved list RDS of the summaries
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Average bind scores
bind_summary <- readRDS(bind_summary_path)



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



# Aggregate coexpr rankings relative to single dataset profiles
# ------------------------------------------------------------------------------


# Human
save_function_results(
  path = avg_vs_ind_auc_hg_path,
  fun = get_colwise_curated_auc_list,
  args = list(
    tfs = tf_hg,
    agg_l = agg_tf_hg,
    msr_mat = msr_hg,
    curated_df = curated,
    ortho_df = pc_ortho,
    pc_df = pc_hg,
    species = "Human",
    ncores = ncore,
    verbose = TRUE
  ),
  force_resave = force_resave
)



# Mouse
save_function_results(
  path = avg_vs_ind_auc_mm_path,
  fun = get_colwise_curated_auc_list,
  args = list(
    tfs = tf_mm,
    agg_l = agg_tf_mm,
    msr_mat = msr_mm,
    curated_df = curated,
    ortho_df = pc_ortho,
    pc_df = pc_mm,
    species = "Mouse",
    ncores = ncore,
    verbose = TRUE
  ),
  force_resave = force_resave
)



# Positive coexpr (top of aggr coexpr rankings)
# ------------------------------------------------------------------------------


# Human
save_function_results(
  path = coexpr_auc_hg_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
    tfs = tf_hg,
    rank_l = rank_tf_hg,
    score_col = "Avg_aggr_coexpr",
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



# # Mouse
save_function_results(
  path = coexpr_auc_mm_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
    tfs = tf_mm,
    rank_l = rank_tf_mm,
    score_col = "Avg_aggr_coexpr",
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



# Negative coexpr (bottom of aggr coexpr rankings)
# ------------------------------------------------------------------------------


# Flip sign of Average_RSR (used as importance score) to make the bottom of
# the rankings (most consistent negative coexpr) more important

flip_avgrsr <- function(rank_df) {
  rank_df$Avg_aggr_coexpr <- -rank_df$Avg_aggr_coexpr
  return(rank_df)
}


rank_tf_hg <- lapply(rank_tf_hg, flip_avgrsr)
rank_tf_mm <- lapply(rank_tf_mm, flip_avgrsr)


# Human
save_function_results(
  path = rev_coexpr_auc_hg_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
    tfs = tf_hg,
    rank_l = rank_tf_hg,
    score_col = "Avg_aggr_coexpr",
    curated_df = curated,
    label_all = target_hg,
    ortho_df = pc_ortho,
    pc_df = pc_hg,
    species = "Human",
    n_samps = n_samps,
    ncores = 8,
    verbose = TRUE
  ),
  force_resave = force_resave
)



# Mouse
save_function_results(
  path = rev_coexpr_auc_mm_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
    tfs = tf_mm,
    rank_l = rank_tf_mm,
    score_col = "Avg_aggr_coexpr",
    curated_df = curated,
    label_all = target_mm,
    ortho_df = pc_ortho,
    pc_df = pc_mm,
    species = "Mouse",
    n_samps = n_samps,
    ncores = 8,
    verbose = TRUE
  ),
  force_resave = force_resave
)



# Binding scores
# ------------------------------------------------------------------------------


# Split matrix of summarized bind scores col/TF-wise into a list of dataframes

split_bindmat_to_list <- function(mat) {
  
  bind_l <- asplit(mat, 2)
  
  bind_l <- lapply(bind_l, function(x) {
    data.frame(Bind_score = x) %>% 
      rownames_to_column(var = "Symbol") %>% 
      arrange(desc(Bind_score))
  })
}


# List of TF binding scores
bind_hg <- split_bindmat_to_list(bind_summary$Human_TF)
bind_mm <- split_bindmat_to_list(bind_summary$Mouse_TF)


# Smaller subset of TFs with binding evidence
tf_hg <- intersect(tf_hg, names(bind_hg))
tf_mm <- intersect(tf_mm, names(bind_mm))



# Human
save_function_results(
  path = unibind_auc_hg_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
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
    verbose = TRUE
  ),
  force_resave = force_resave
)



# Mouse
save_function_results(
  path = unibind_auc_mm_path,
  fun = curated_obs_and_null_auc_list,
  args = list(
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
    verbose = TRUE
  ),
  force_resave = force_resave
)
