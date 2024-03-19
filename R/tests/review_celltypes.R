## Hacky script for examining the saved cell types

library(googlesheets4)
library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")

# Gsheets review sheet
review <- read_sheet(gsheets_id, sheet = "Review_celltypes", trim_ws = TRUE)

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

ct_l <- readRDS(celltype_list_path)


# Inpsect presence of specific cell type


ct <- "oligo"  # "microglia|mcg"


check_ct <- lapply(ct_l, function(x) {
  ct_vec <- str_to_lower(x$Ct_count$Cell_type)
  ct_which <- str_detect(ct_vec, ct)
  
  if (sum(ct_which) == 0) {
    return(NA)
  }
  
  x$Ct_count[ct_which, ]
  
})

check_ct <- check_ct[!is.na(check_ct)]


tally_ct <- table(filter(sc_meta, ID %in% names(check_ct))$Species)


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

