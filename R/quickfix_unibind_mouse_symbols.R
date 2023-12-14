# NOTE: This is a hacky quick fix to coerce mouse symbol names for unibind data
# that are not using the official gene name (namely TP53 -> Trp53 and TP63 -> Trp63)
# Need to incorporate this in the unibind processing scripts when I get around
# to integrating that workflow with this project

library(tidyverse)
source("R/00_config.R")

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Average bind scores
collection <- "Permissive"
stopifnot(collection %in% c("Robust", "Permissive"))
bind_summary_path <- paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_bindscore_summary.RDS")
bind_summary <- readRDS(bind_summary_path)


setdiff(colnames(bind_summary$Mouse_TF), pc_ortho$Symbol_hg)

old <- colnames(bind_summary$Mouse_TF)

replace_df <- data.frame(
  Old = old,
  Match = match(old, pc_ortho$Symbol_hg))

replace_df$New <- pc_ortho$Symbol_mm[replace_df$Match]
replace_df$New <- ifelse(is.na(replace_df$New), replace_df$Old, replace_df$New)

colnames(bind_summary$Mouse_TF) <- replace_df$New

saveRDS(bind_summary, bind_summary_path)
