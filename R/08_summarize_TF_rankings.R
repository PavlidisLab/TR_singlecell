## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Saved list RDS of the summaries
summ_hg <- readRDS(tf_summ_hg_path)
summ_mm <- readRDS(tf_summ_mm_path)

evidence_l <- readRDS(evidence_path)

arrange(summ_hg$E2F8, desc(Topk_count)) %>% head(30)
arrange(summ_mm$E2f8, desc(Topk_count)) %>% head(30)


# Create a gene x TF matrix of summarized ranks


rank_mat_hg <- as.matrix(do.call(cbind, lapply(summ_hg, `[`, "Rank_RSR")))
colnames(rank_mat_hg) <- names(summ_hg)


rank_mat_mm <- as.matrix(do.call(cbind, lapply(summ_mm, `[`, "Rank_RSR")))
colnames(rank_mat_mm) <- names(summ_mm)


# Most commonly highly ranked genes
# TODO: are these genes high because of technical aspect?

rank_order_hg <- sort(rowMeans(rank_mat_hg, na.rm = TRUE))
head(rank_order_hg)
head(sort(rank_order_hg, decreasing = TRUE))
head(sort(rank_mat_hg["RPL3",]), 50)
head(sort(rank_mat_hg["DLL3",]), 50)
head(sort(rank_mat_hg["DIXDC1",], decreasing = TRUE), 50)



rank_order_mm <- sort(rowMeans(rank_mat_mm, na.rm = TRUE))
head(sort(rank_mat_mm["Bnc1",]), 50)
head(sort(rank_mat_mm["Esrrg",]), 50)
head(sort(rank_mat_mm["Nbl1",], decreasing = TRUE), 50)
