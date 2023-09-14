library(googlesheets4)
library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")

# Gsheets review sheet
review <- read_sheet(gsheets_id, sheet = "Review_celltypes", trim_ws = TRUE)

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

ct_l <- readRDS("~/sc_output/celltype_list.RDS")


# Dendritic cells that had integer delim removed (may be actual identity rather
# than cluster identifier)

dcs <- lapply(ct_l, function(x) {
  res <- NA
  ct <- str_to_lower(x$Ct_count$Cell_type)
  which_dc <- str_detect(ct, "dendritic|dc")
  if (any(which_dc)) res <- ct[which_dc]
  return(res)
})

dcs <- dcs[!is.na(dcs)]


dcs_collapse <- intersect(names(dcs), 
                          filter(review, as.logical(Collapsed_celltype))$ID)

dcs[dcs_collapse]



## Blank cell type labels

blank <- lapply(sc_meta$ID, function(x) {
  
  dat <- load_dat_list(x)[[1]]
  
  meta <- dat$Meta %>% 
  filter(is.na(Cell_type) | Cell_type == "")
  
  if (nrow(meta) == 0) {
    return(NA)
  }
  
  return(meta)

})

names(blank) <- sc_meta$ID
blank <- blank[!is.na(blank)]

saveRDS(blank, "~/scratch/R_objects/blank.RDS")


###

# agg_l = agg_tf_mm
# msr_mat = msr_mm
# genes = tfs_mm$Symbol
# 
# x <- "Tgif2lx1"
# id <- "GSE160523"
# 
# agg <- load_agg_mat_list(id, genes = pc_mm$Symbol)
# dat <- load_dat_list(id)
# mat <- dat[[id]]$Mat
# meta <- dat[[id]]$Meta
# na_mat <- load_agg_mat_list(id, genes = pc_mm$Symbol, pattern = "_NA_mat.tsv")
# 
# 

# ct_l$GSE160523
# filter(sc_meta, ID == id)$N_celltypes
# max(na_mat$GSE160523)
# unique(meta$Cell_type)
# sum(is.na(meta$Cell_type))
# sum(meta$Cell_type == "")
# 
# 
# table(meta$Cell_type)
# table(meta$Cell_type, useNA = "ifany")
# 
# 
# 
# for (ct in unique(meta$Cell_type)) {
#   print(ct)
#   print(sum(mat[x, filter(meta, Cell_type == ct)$ID] != 0))
# }
# 
# 
# sum(mat[x, ] != 0)
# 
# sort(table(filter(meta, ID %in% names(which(mat[x, ] != 0)))$Cell_type))
# 
# 
# mat[, filter(meta, is.na(Cell_type))$ID]
# 
# 
# ###