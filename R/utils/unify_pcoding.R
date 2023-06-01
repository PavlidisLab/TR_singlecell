## This script ensures that ENSEMBL protein coding tables and Refseq tables have
## the same gene symbols represented, to facilitate the same dimensionality when
## working with a count matrix using either IDs.
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")


ref_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
ref_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
ens_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
ens_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)


ens_hg <- ens_hg %>% 
  filter(Symbol %in% ref_hg$Symbol & Symbol != "") %>% 
  right_join(ref_hg[, c("Symbol", "Refseq_ID")], by = "Symbol")


ens_mm <- ens_mm %>% 
  filter(Symbol %in% ref_mm$Symbol & Symbol != "") %>% 
  right_join(ref_mm[, c("Symbol", "Refseq_ID")], by = "Symbol")



write.table(ens_hg, sep = "\t", quote = FALSE, file = ens_hg_path)
write.table(ens_mm, sep = "\t", quote = FALSE, file = ens_mm_path)
