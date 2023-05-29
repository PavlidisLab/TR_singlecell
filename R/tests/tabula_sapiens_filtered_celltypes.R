library(WGCNA)
library(tidyverse)
library(Seurat)
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


mat_filt <- rm_low_qc_cells(mat, meta) %>%
  ensembl_to_symbol(ensembl_df = pc) %>% 
  get_pcoding_only(pcoding_df = pc) %>% 
  Seurat::LogNormalize(., verbose = FALSE)


meta_filt <- filter(meta, ID %in% colnames(mat_filt))


meta_rm <- filter(meta, ID %in% setdiff(meta$ID, meta_filt$ID))


n_rm <- table(meta_rm$Cell_type)
n_all <- table(meta$Cell_type)


rm_df <- data.frame(
  Cell_type = names(n_rm),
  N_rm = as.integer(n_rm),
  N_all = as.integer(n_all)
)
rm_df$Ratio <- rm_df$N_rm / rm_df$N_all


rm_df <- rm_df %>% 
  arrange(Ratio) %>% 
  mutate(Cell_type = factor(Cell_type, levels = unique(Cell_type)))


ggplot(rm_df, aes(x = log10(N_all+1), y = Ratio)) +
  geom_point() +
  theme_classic()


# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5126291/
# Notes that pancreatic acinar cell are expected to have low RNA complexity




