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


# Keep only X copies of duplicated genes in human

dupl_genes <- ref_hg %>% group_by(Symbol) %>% filter(n() > 1)

ref_hg <- filter(ref_hg, !(Symbol %in% dupl_genes$Symbol & Chromosome == "Y"))

  
# Reduce Ensembl to Refseq select, keeping only a single entry per symbol, and
# only the X copy of duplicated counts. As using gene-summarized counts, the 
# exact transcript is unimportant


ens_hg <- ens_hg %>%
  filter(Symbol %in% ref_hg$Symbol & Symbol != "") %>%
  right_join(ref_hg[, c("Symbol", "Refseq_ID")], by = "Symbol") %>%
  filter(!(Symbol %in% dupl_genes$Symbol & Chromosome == "Y")) %>%
  distinct(Symbol, .keep_all = TRUE) %>%
  arrange(match(Symbol, ref_hg$Symbol))


ens_mm <- ens_mm %>% 
  filter(Symbol %in% ref_mm$Symbol & Symbol != "") %>% 
  right_join(ref_mm[, c("Symbol", "Refseq_ID")], by = "Symbol") %>% 
  distinct(Symbol, .keep_all = TRUE) %>%
  arrange(match(Symbol, ref_mm$Symbol))



stopifnot(identical(ens_hg$Symbol, ref_hg$Symbol))
stopifnot(identical(ens_mm$Symbol, ref_mm$Symbol))


write.table(ref_hg, sep = "\t", quote = FALSE, file = ref_hg_path)
write.table(ref_mm, sep = "\t", quote = FALSE, file = ref_mm_path)
write.table(ens_hg, sep = "\t", quote = FALSE, file = ens_hg_path)
write.table(ens_mm, sep = "\t", quote = FALSE, file = ens_mm_path)
