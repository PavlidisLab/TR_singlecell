## Look at co-measurement of TF targets and the TF itself across scRNA data
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
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)

# Ranked targets from Morin2023
evidence_l <- readRDS(evidence_path)



# Examining instances of genes with strong regulation evidence but are never
# measured with the TF itself
# TODO: compare using mismatched targets
# ------------------------------------------------------------------------------


# Loading

na_hg <- load_agg_mat_list(ids_hg, 
                           pattern = "_NA_mat.tsv",
                           genes = pc_hg$Symbol, 
                           sub_genes = names(evidence_l$Human))


na_mm <- load_agg_mat_list(ids_mm, 
                           pattern = "_NA_mat.tsv",
                           genes = pc_mm$Symbol, 
                           sub_genes = names(evidence_l$Mouse))




# Convert count matrices to proportions

n_ct_hg <- filter(sc_meta, ID %in% ids_hg)$N_celltype

n_ct_mm <- filter(sc_meta, ID %in% ids_mm)$N_celltype

na_hg <- lapply(1:length(na_hg), function(x) sweep(na_hg[[x]], 2, n_ct_hg[x], `/`))
names(na_hg) <- ids_hg

na_mm <- lapply(1:length(na_mm), function(x) sweep(na_mm[[x]], 2, n_ct_mm[x], `/`))
names(na_mm) <- ids_mm

# Extract a single TF's vector across datasets into a single matrix

gene_hg <- "RUNX1"
gene_mm <- str_to_title(gene_hg)
stopifnot(gene_mm %in% pc_mm$Symbol)

gene_mat_hg <- gene_vec_to_mat(na_hg, gene_hg)
gene_mat_mm <- gene_vec_to_mat(na_mm, gene_mm)

# Order genes by their average measurement with the gene of interest. 
# A value of 1 means that the gene pair were never measured together (possibly
# because the gene itself isn't expressed). Lower values mean that a gene pair
# are more commonly measured in the same cell types

gene_order_hg <- sort(rowMeans(gene_mat_hg))
gene_order_mm <- sort(rowMeans(gene_mat_mm))

hist(gene_order_hg, breaks = 100)
hist(gene_order_mm, breaks = 100)

# Compare this order to the regulation evidence

gene_evidence_hg <- left_join(
  evidence_l$Human[[gene_hg]],
  rownames_to_column(data.frame(NA_prop = gene_order_hg), var = "Symbol"),
  by = "Symbol")


gene_evidence_mm <- left_join(
  evidence_l$Mouse[[gene_mm]],
  rownames_to_column(data.frame(NA_prop = gene_order_mm), var = "Symbol"),
  by = "Symbol")


evidence_col <- "Rank_integrated"


cor(gene_evidence_hg[[evidence_col]], gene_evidence_hg$NA_prop, method = "spearman")
cor(gene_evidence_mm[[evidence_col]], gene_evidence_mm$NA_prop, method = "spearman")


# Cut evidence into discrete groups for visualizing relationship with NA proportion.


p3a <- gene_evidence_hg %>% 
  mutate(Group_evidence = cut(gene_evidence_hg[[evidence_col]], breaks = 50)) %>% 
  ggplot(aes(x = Group_evidence, y = NA_prop)) + 
  geom_boxplot() +
  xlab("Binned regulation evidence") + 
  ylab("Proportion of NA measurements") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 20))



p3b <- gene_evidence_mm %>% 
  mutate(Group_evidence = cut(gene_evidence_mm[[evidence_col]], breaks = 50)) %>% 
  ggplot(aes(x = Group_evidence, y = NA_prop)) + 
  geom_boxplot() +
  xlab("Binned regulation evidence") + 
  ylab("Proportion of NA measurements") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 20))


# Inspecting highest ranked gene that is not ever detected with the TF, and 
# whether they are actually measured at all


top_hg <- gene_evidence_hg %>% 
  filter(NA_prop > 0.9) %>% 
  slice_min(Rank_integrated, n = 20) %>% 
  select(Symbol, Count_DE, Count_NA, Rank_binding, Rank_perturbation, Rank_integrated)


top_mm <- gene_evidence_mm %>% 
  filter(NA_prop == 1) %>% 
  slice_min(Rank_integrated, n = 20) %>% 
  select(Symbol, Count_DE, Count_NA, Rank_binding, Rank_perturbation, Rank_integrated)


rowSums(msr_mat_hg[top_hg$Symbol, ])
rowSums(msr_mat_mm[top_mm$Symbol, ])