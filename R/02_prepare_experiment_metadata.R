## Load a google sheets metadata table tracking scRNA-seq experiments, then
## format and save a local copy for available experiments
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")

sc_meta <- read_sheet(gsheets_id, sheet = "Main", trim_ws = TRUE)


stopifnot(all(sc_meta$Species %in% c("Human", "Mouse")))


loaded <- lapply(sc_meta$ID, function(x) {
  f <- list.files(file.path(amat_dir, x), pattern = "RSR_allrank.tsv", full.names = TRUE)
  if (length(f) == 0) f <- NA
  return(f)
})

names(loaded) <- sc_meta$ID

loaded <- unlist(loaded[!is.na(loaded)])

all_paths <- data.frame(ID = names(loaded), Path = loaded)


sc_meta <- filter(sc_meta, ID %in% all_paths$ID) %>% 
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


meta_cols <- do.call(rbind, lapply(sc_meta$ID, add_meta_cols))


stopifnot(identical(meta_cols$ID, sc_meta$ID))
sc_meta$N_cells <- meta_cols$N_cells
sc_meta$N_celltypes <- meta_cols$N_celltypes


# sc_meta2 <- left_join(
#   dplyr::select(sc_meta, -c(N_cells, N_celltypes)),
#   meta_cols, 
#   by = "ID")


write.table(sc_meta,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = sc_meta_path)
