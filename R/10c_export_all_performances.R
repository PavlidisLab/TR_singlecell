## Organize aggregate performance metrics from 07a_ and 10a_ to export as tables
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")

# Aggregate coexpr profiles versus null 
coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)

# Aggregate reversed coexpr (prioritize negative cor) profiles versus null 
rev_coexpr_auc_hg <- readRDS(rev_coexpr_auc_hg_path)
rev_coexpr_auc_mm <- readRDS(rev_coexpr_auc_mm_path)

# Aggregate binding profiles versus null
bind_auc_hg <- readRDS(unibind_auc_hg_path)
bind_auc_mm <- readRDS(unibind_auc_mm_path)

# Integrated ranking versus null 
int_auc_hg <- readRDS(int_auc_hg_path)
int_auc_mm <- readRDS(int_auc_mm_path)


# Join the summary dataframes from the list of coexpression, binding, and 
# integrated performa
# ------------------------------------------------------------------------------

join_auc_df <- function(coexpr_l, negcoexpr_l, bind_l, int_l) {
  
  coexpr_df <- do.call(rbind, lapply(coexpr_l, `[[`, "Perf_df"))
  negcoexpr_df <- do.call(rbind, lapply(negcoexpr_l, `[[`, "Perf_df"))
  bind_df <- do.call(rbind, lapply(bind_l, `[[`, "Perf_df"))
  int_df <- do.call(rbind, lapply(int_l, `[[`, "Perf_df"))
  
  colnames(coexpr_df)[3:ncol(coexpr_df)] <- paste0(colnames(coexpr_df)[3:ncol(coexpr_df)], "_coexpr")
  colnames(negcoexpr_df)[3:ncol(negcoexpr_df)] <- paste0(colnames(negcoexpr_df)[3:ncol(negcoexpr_df)], "_negcoexpr")
  colnames(bind_df)[3:ncol(bind_df)] <- paste0(colnames(bind_df)[3:ncol(bind_df)], "_bind")
  colnames(int_df)[3:ncol(int_df)] <- paste0(colnames(int_df)[3:ncol(int_df)], "_integrated")
  
  join_df <- plyr::join_all(list(coexpr_df, negcoexpr_df, bind_df, int_df),
                            by = c("Symbol", "N_targets"))

  return(join_df)
}



agg_df_hg <- join_auc_df(coexpr_l = coexpr_auc_hg, 
                         negcoexpr_l = rev_coexpr_auc_hg,
                         bind_l = bind_auc_hg, 
                         int_l = int_auc_hg)


agg_df_mm <- join_auc_df(coexpr_l = coexpr_auc_mm, 
                         negcoexpr_l = rev_coexpr_auc_mm,
                         bind_l = bind_auc_mm, 
                         int_l = int_auc_mm)


write.table(agg_df_hg, sep = "\t", row.names = FALSE, quote = FALSE,
            file = auc_table_hg_path)


write.table(agg_df_mm, sep = "\t", row.names = FALSE, quote = FALSE,
            file = auc_table_mm_path)
