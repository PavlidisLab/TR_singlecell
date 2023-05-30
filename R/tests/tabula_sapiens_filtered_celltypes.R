library(Seurat)
library(UpSetR)
library(cowplot)
library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

id <- "Tabula_Sapiens"
sc_dir <- file.path("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datasets", id)
dat_path <- file.path(sc_dir, paste0(id, "_cellxgene_seurat.RDS"))

pc <- read.delim(ens_hg_path, stringsAsFactors = FALSE)

dat <- readRDS(dat_path)

mat <- GetAssayData(dat, slot = "counts")

meta <- dat@meta.data %>% 
  dplyr::rename(Cell_type = cell_type) %>% 
  rownames_to_column(var = "ID") %>% 
  add_count_info(mat = mat)

# Separate the four QC metrics used for filtering

rm_list <- list(
  UMI_counts = filter(meta, UMI_counts < 500)$ID,
  Gene_counts = filter(meta, Gene_counts < 250)$ID,
  RNA_novelty = filter(meta, RNA_novelty < 0.8)$ID,
  MT_ratio = filter(meta, MT_ratio >= 0.2)$ID
)


p1 <- upset(fromList(rm_list), empty.intersections = "on", text.scale = 2.5)


# Group all removed, and by gene/UMI counts versus mito ratio

rm_all <- Reduce(union, rm_list)


# Proportion of filtered cell types

n_ct_all <- table(meta$Cell_type)
n_ct_rm <- lapply(rm_list, function(x) table(filter(meta, ID %in% x)$Cell_type))
prop_ct_rm <- lapply(n_ct_rm, function(x) x / n_ct_all)


stopifnot(identical(names(prop_ct_rm$RNA_novelty), names(n_ct_all)))


ct_df <- data.frame(
  Cell_type = names(n_ct_all),
  N = as.integer(n_ct_all),
  do.call(cbind, prop_ct_rm)
)


p2 <- ggplot(ct_df, aes(x = log10(N), y = RNA_novelty)) +
  geom_point(size = 3, shape = 19) + 
  ylab("Proportion of cells removed by RNA novelty filter") +
  xlab("Log10 Count of cells in cell type") +
  theme_classic() + 
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))


p3 <- ggplot(ct_df, aes(x = log10(N), y = MT_ratio)) +
  geom_point(size = 3, shape = 19) + 
  ylab("Proportion of cells removed by MT ratio filter") +
  xlab("Log10 Count of cells in cell type") +
  theme_classic() + 
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))


p4 <- ggplot(ct_df, aes(x = RNA_novelty, y = MT_ratio)) +
  geom_point(size = 3, shape = 19) + 
  ylab("Proportion of cells removed by MT ratio filter") +
  xlab("Proportion of cells removed by RNA novelty filter") +
  theme_classic() + 
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))


plot_grid(p2, p4)


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5126291/
# Notes that pancreatic acinar cell are expected to have low RNA complexity




