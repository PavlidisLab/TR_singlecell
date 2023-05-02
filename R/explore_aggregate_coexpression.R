library(tidyverse)
library(Seurat)
library(pheatmap)
library(RColorBrewer)
source("R/utils/agg_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

pc_ortho <- read.delim(pc_ortho_path)

# Ranked targets from paper
evidence_l <- readRDS(evidence_path)
tfs_mm <- names(rank_l$Mouse)
tfs_hg <- names(rank_l$Human)


# IDs for scRNA-seq datasets
ids_hg <- c("Velmeshev", "ROSMAP", "GSE180928", "GSE216019")
ids_mm <- c("Posner2022", "MKA")

# Load aggregate matrix into list

rsr1_hg <- lapply(ids_hg, function(x) readRDS(paste0(amat_dir, x, "_RSR1.RDS")))
names(rsr1_hg) <- ids_hg

rsr1_mm <- lapply(ids_mm, function(x) readRDS(paste0(amat_dir, x, "_RSR1.RDS")))
z1_mm <- lapply(ids_mm, function(x) readRDS(paste0(amat_dir, x, "_Z1.RDS")))
names(rsr1_mm) <- names(z1_mm) <- ids_mm

# Load mat and meta

dat_mm <- lapply(ids_mm, function(x) readRDS(paste0(amat_dir, x, "_mat_and_meta.RDS")))
names(dat_mm) <- ids_mm

dat_hg <- lapply(ids_hg, function(x) readRDS(paste0(amat_dir, x, "_mat_and_meta.RDS")))
names(dat_hg) <- ids_hg


# TODO: won't be necessary when start with fixed pcoding genes
common_mm <- Reduce(intersect, lapply(rsr1_mm, rownames))

gene <- "Mecp2"

stopifnot(gene %in% common_mm)



# get_df <- function(rsr1_l, z1_l, rank_df, sub_gene) {
#   
#   ids <- intersect(names(rsr1_l), names(z1_l))
#   
#   
#   tt <- lapply(ids, function(x) {
#     
#     rsr1_l[[x]]
#     
#   })
# }


df <- data.frame(Symbol = common_mm,
                 Col_rank_Posner2022 = rsr1_mm$Posner2022[common_mm, gene],
                 Zscore_all_Posner2022 = z1_mm$Posner2022$Avg_all[common_mm, gene],
                 Zscore_nonNA_Posner2022 = z1_mm$Posner2022$Avg_nonNA[common_mm, gene],
                 Count_NA_cor_Posner2022 = z1_mm$Posner2022$NA_mat[common_mm, gene],
                 Col_rank_MKA = rsr1_mm$MKA[common_mm, gene],
                 Zscore_all_MKA = z1_mm$MKA$Avg_all[common_mm, gene],
                 Zscore_nonNA_MKA = z1_mm$MKA$Avg_nonNA[common_mm, gene],
                 Count_NA_cor_MKA = z1_mm$MKA$NA_mat[common_mm, gene])



df <- left_join(evidence_l$Mouse[[gene]], df, by = "Symbol") %>%
  filter(Symbol != gene)


# cor(select_if(df, is.numeric), use = "pairwise.complete.obs", method = "spearman")


df$Posner2022_RS <- rank(-df$Col_rank_Posner2022)/length(common_mm)
df$MKA_RS <- rank(-df$Col_rank_MKA)/length(common_mm)

plot(df$Zscore_all_Posner2022, df$Zscore_all_MKA)
plot(df$Zscore_nonNA_Posner2022, df$Zscore_nonNA_MKA)
plot(df$Col_rank_Posner2022, df$Col_rank_MKA)
plot(df$Posner2022_RS, df$MKA_RS)



df %>% 
  filter(Rank_integrated <= 500) %>% 
  slice_max(Zscore_all, n = 10)


mat <- dat_mm$Posner2022[[1]]
meta <- dat_mm$Posner2022[[2]]

mat <- dat_hg$Velmeshev[[1]]
meta <- dat_hg$Velmeshev[[2]]

gene <- "RUNX1"
boxplot(mat[gene, ] ~ meta$Cell_type)
all_celltype_cor(mat = mat, meta = meta, gene1 = gene, gene2 = "Fubp1")
