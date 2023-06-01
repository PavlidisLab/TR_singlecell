##
## -----------------------------------------------------------------------------

library(googlesheets4)
source("R/00_config.R")

sc_meta <- read_sheet(gsheets_id, sheet = "Main", trim_ws = TRUE)

stopifnot(all(sc_meta$Species %in% c("Human", "Mouse")))

write.table(sc_meta,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = sc_meta_path)
