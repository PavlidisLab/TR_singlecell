## TODO: 
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
source("R/00_config.R")

sc_meta <- read_sheet(gsheets_id, sheet = "Main", trim_ws = TRUE)

stopifnot(all(sc_meta$Species %in% c("Human", "Mouse")))


# TODO: This hacky crap is because the aggregate matrices are currently split
# between an old ranking run (no seed) and an update with set seed. So need
# to collect


old_dir <- file.path(amat_dir, "old")

old_ids <- list.files(old_dir) %>% 
  str_split(pattern = "_", simplify = TRUE) %>% 
  .[, 1]
  
# eugh
old_ids[old_ids == "Tabula"] <- "Tabula_Sapiens"
  

current <- lapply(sc_meta$ID, function(x) {
  f <- list.files(file.path(amat_dir, x), pattern = "RSR_allrank.RDS", full.names = TRUE)
  if (length(f) == 0) f <- NA
  return(f)
})

names(current) <- sc_meta$ID

current <- unlist(current[!is.na(current)])


# remove old that have a current generated
old_ids <- old_ids[!old_ids %in% names(current)]
old_paths <- file.path(old_dir, paste0(old_ids, "_RSR_allrank.RDS"))
names(old_paths) <- old_ids

all_paths <- data.frame(ID = c(names(current), old_ids),
                        Path = c(current, old_paths))


sc_meta <- filter(sc_meta, ID %in% all_paths$ID) %>% 
  left_join(all_paths, by = "ID")
  
  
write.table(sc_meta,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = sc_meta_path)
