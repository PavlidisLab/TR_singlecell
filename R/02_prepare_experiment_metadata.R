## Load a google sheets metadata table tracking scRNA-seq experiments, then
## format and save a local copy for available experiments
## TODO: finalize loading structure
## TODO: counts should be only considering new data (very slow otherwise)
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")

meta <- read_sheet(gsheets_id, sheet = "Main", trim_ws = TRUE)


stopifnot(all(meta$Species %in% c("Human", "Mouse")))


loaded <- lapply(meta$ID, function(x) {
  file.exists(file.path(amat_dir, x, paste0(x, "_RSR_allrank.tsv")))
})


meta_loaded <- meta[unlist(loaded), ]



# For inspecting those that are not loaded/failed
setdiff(meta$ID, meta_loaded$ID)


loaded <- unlist(loaded[!is.na(loaded)])


all_paths <- data.frame(ID = names(loaded), Path = loaded)


meta_loaded <- filter(meta, ID %in% all_paths$ID) %>% 
  left_join(all_paths, by = "ID")
  

# TODO: consider accounting for number of cells by loading NA mat and checking NA
# counts against count of cell types

# TODO: figure out how to consolidate this with gsheets and only check missing
# don't have to rerun over everything each time

add_meta_cols <- function(id) {
  
  dat <- load_dat_list(id)[[1]]
  
  data.frame(
    ID = id,
    N_cells = ncol(dat$Mat),
    N_celltypes = n_distinct(dat$Meta$Cell_type)
  )
}


meta_cols <- do.call(rbind, lapply(meta_loaded$ID, add_meta_cols))


stopifnot(identical(meta_cols$ID, meta_loaded$ID))
meta_loaded$N_cells <- meta_cols$N_cells
meta_loaded$N_celltypes <- meta_cols$N_celltypes


# meta2 <- left_join(
#   dplyr::select(meta, -c(N_cells, N_celltypes)),
#   meta_cols, 
#   by = "ID")


write.table(meta_loaded,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = meta_path)
