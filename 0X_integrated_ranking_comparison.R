## TODO:
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# The minimum count of curated targets for a TF for reporting
min_targets <- 5

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Aggregate coexpr profiles versus null 
coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)

# Integrated ranking versus null 
int_auc_hg <- readRDS(int_auc_hg_path)
int_auc_mm <- readRDS(int_auc_mm_path)

# Aggregate binding profiles versus null
bind_auc_hg <- readRDS(unibind_auc_hg_path)
bind_auc_mm <- readRDS(unibind_auc_mm_path)


# Functions
# ------------------------------------------------------------------------------


# Join the summary dataframes from the list of coexpression and bind AUCs

join_auc_df <- function(coexpr_l, bind_l, int_l, min_targets) {
  
  coexpr_df <- do.call(rbind, lapply(coexpr_l, `[[`, "Perf_df"))
  bind_df <- do.call(rbind, lapply(bind_l, `[[`, "Perf_df"))
  int_df <- do.call(rbind, lapply(int_l, `[[`, "Perf_df"))
  
  join_df <- 
    left_join(coexpr_df, bind_df,
              by = c("Symbol", "N_targets"),
              suffix = c("_coexpr", "_bind")) %>% 
    left_join(., int_df, 
              suffix = c("", "_integrated")) %>% 
    filter(N_targets >= min_targets)
  
  return(join_df)
}



#
# ------------------------------------------------------------------------------

agg_df_hg <- join_auc_df(coexpr_l = coexpr_auc_hg, 
                         bind_l = bind_auc_hg, 
                         int_l = int_auc_hg,
                         min_targets = min_targets)


agg_df_mm <- join_auc_df(coexpr_l = coexpr_auc_mm, 
                         bind_l = bind_auc_mm,
                         int_auc_mm,
                         min_targets = min_targets)


agg_df_hg$N_coexpr_msr <- rowSums(msr_hg[agg_df_hg$Symbol, ])
agg_df_mm$N_coexpr_msr <- rowSums(msr_mm[agg_df_mm$Symbol, ])




