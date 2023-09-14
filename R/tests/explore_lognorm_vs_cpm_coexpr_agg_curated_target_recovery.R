## TODO:
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

n_samps <- 1000

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Saved list of performances

coexpr_recover_curated_ln_hg_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_ln_hg.RDS"
coexpr_recover_curated_ln_mm_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_ln_mm.RDS"

coexpr_recover_curated_cpm_hg_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_cpm_hg.RDS"
coexpr_recover_curated_cpm_mm_path <- "/space/scratch/amorin/R_objects/coexpr_recover_curated_cpm_mm.RDS"


perf_ln_hg <- readRDS(coexpr_recover_curated_ln_hg_path)
perf_cpm_hg <- readRDS(coexpr_recover_curated_cpm_hg_path)




# TODO: can keep be done upstream?

join_perf_df <- function(rank_l1, rank_l2) {
  
  keep1 <- names(which(unlist(lapply(rank_l1, function(x) "Perf_df" %in% names(x)))))
  df1 <- do.call(rbind, lapply(rank_l1[keep1], `[[`, "Perf_df"))
  
  keep2 <- names(which(unlist(lapply(rank_l2, function(x) "Perf_df" %in% names(x)))))
  df2 <- do.call(rbind, lapply(rank_l2[keep2], `[[`, "Perf_df"))
  
  left_join(df1, df2,
            by = c("Symbol", "N_targets"),
            suffix = c("_LN", "_CPM")
            ) %>% 
    mutate(
      Diff_AUPRC = AUPRC_percentile_observed_LN - AUPRC_percentile_observed_CPM,
      Diff_AUROC = AUROC_percentile_observed_LN - AUROC_percentile_observed_CPM
    )
}




perf_hg <- join_perf_df(perf_ln_hg, perf_cpm_hg) %>% filter(N_targets >= 5)

qplot(perf_hg, "AUPRC_percentile_observed_LN", "AUPRC_percentile_observed_CPM")
qplot(perf_hg, "AUROC_percentile_observed_LN", "AUROC_percentile_observed_CPM")

hist(perf_hg$Diff_AUPRC, breaks = 100, xlim = c(-1, 1), main = "AUPRC")
hist(perf_hg$Diff_AUROC, breaks = 100, xlim = c(-1, 1), main = "AUROC")
