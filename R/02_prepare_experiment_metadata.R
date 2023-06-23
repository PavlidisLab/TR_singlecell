## Load a google sheets metadata table tracking scRNA-seq experiments, then
## format and save a local copy for available experiments
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
source("R/00_config.R")

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
  
  
write.table(sc_meta,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = sc_meta_path)
