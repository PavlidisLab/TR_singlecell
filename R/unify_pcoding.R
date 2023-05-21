library(tidyverse)
source("R/00_config.R")

ref_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
ref_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
ens_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
ens_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

ens_hg <- filter(ens_hg, Symbol %in% ref_hg$Symbol & Symbol != "")
ens_mm <- filter(ens_mm, Symbol %in% ref_mm$Symbol & Symbol != "")


write.table(ens_hg, sep = "\t", quote = FALSE, file = ens_hg_path)
write.table(ens_mm, sep = "\t", quote = FALSE, file = ens_mm_path)

