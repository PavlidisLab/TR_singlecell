## Examine gene measurement/coverage across experiments. A gene is considered
## measured if has at least 1 count in at least 20 cells in at least one cell
## type. Exports a binary gene x experiment matrix that tracks measurement.
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



# Create a binary gene by experiment matrix that tracks whether a given gene
# was measured (at least one count in at least 20 cells in at least one cell
# type) for each dataset/id.
# ------------------------------------------------------------------------------


get_gene_msr_mat <- function(ids, meta, genes) {
  
  msr_mat <- matrix(0, nrow = length(genes), ncol = length(ids))
  colnames(msr_mat) <- ids
  rownames(msr_mat) <- genes
  
  # For each dataset, load the matrix tracking NA counts of gene msrmt. pairs, 
  # using the diagonal (identity) to binarize if a gene was measured at least
  # once (1) or not (0)
  
  for (id in ids) {
    
    na_mat <- load_agg_mat_list(id, genes = genes, pattern = "_NA_mat.tsv")[[1]]
    n_celltype <- filter(meta, ID == id)$N_celltype
    binary_msr <- as.integer(diag(na_mat) != n_celltype)
    msr_mat[, id] <- binary_msr
    rm(na_mat)
    gc(verbose = FALSE)
    
  }
  
  return(msr_mat)
}



if (!file.exists(msr_mat_hg_path)) {
  msr_mat_hg <- get_gene_msr_mat(ids_hg, sc_meta, pc_hg$Symbol)
  saveRDS(msr_mat_hg, msr_mat_hg_path)
} else {
  msr_mat_hg <- readRDS(msr_mat_hg_path)
}


if (!file.exists(msr_mat_mm_path)) {
  msr_mat_mm <- get_gene_msr_mat(ids_mm, sc_meta, pc_mm$Symbol)
  saveRDS(msr_mat_mm, msr_mat_mm_path)
} else {
  msr_mat_mm <- readRDS(msr_mat_mm_path)
}


# Get the average/proportion of gene measurement across experiments
# Human: all genes are measured at least twice and 5,140 genes always measured
# Mouse: all genes are measured at least 4 times and 864 genes are always measured
# Ortho: 649 genes measured in every mouse and human dataset
# TFs: In both species TF genes show small trend of being measured more than non-TF genes 
# ------------------------------------------------------------------------------


count_msr_hg <- rowSums(msr_mat_hg)
prop_msr_hg <- count_msr_hg / ncol(msr_mat_hg)

count_msr_mm <- rowSums(msr_mat_mm)
prop_msr_mm <- count_msr_mm / ncol(msr_mat_mm)


gene_msr_df <- data.frame(
  Count_msr = c(count_msr_hg, count_msr_mm),
  Proportion_msr = c(prop_msr_hg, prop_msr_mm),
  Symbol = c(names(prop_msr_hg), names(prop_msr_mm)),
  Species = c(rep("Human", length(prop_msr_hg)), rep("Mouse", length(prop_msr_mm))),
  TF = c(names(prop_msr_hg) %in% tfs_hg$Symbol, names(prop_msr_mm) %in% tfs_mm$Symbol)
)


never_msr <- gene_msr_df %>% 
  filter(Proportion_msr == 0) %>%
  split(.$Species)


rarely_msr <- gene_msr_df %>% 
  filter(Proportion_msr <= 0.1) %>%
  arrange(Count_msr) %>% 
  split(.$Species)


always_msr <- gene_msr_df %>% 
  filter(Proportion_msr == 1) %>%
  split(.$Species)


mostly_msr <- gene_msr_df %>% 
  filter(Proportion_msr >= 0.9) %>%
  split(.$Species)


always_ortho <- filter(pc_ortho, 
                       Symbol_hg %in% always_msr$Human$Symbol &
                       Symbol_mm %in% always_msr$Mouse$Symbol)


tf_med_count <- gene_msr_df %>% 
  group_by(Species, TF) %>% 
  summarise(Median_count = median(Count_msr), .groups = "keep")


tf_wilx_count <- gene_msr_df %>% 
  group_by(Species) %>% 
  do(W = wilcox.test(Count_msr ~ TF, data = ., paired = FALSE)) %>% 
  summarise(Species, Wilcox = W$p.value)



# Look at experiment-wise gene coverage
# Human: TabulaSapiens and Olah2020 measure every gene, GSE85241 min at 8953
# Mouse: No experiment measures all genes, but TabulaMuris closest at 20339.
# GSE160193 min at 2046 genes.
# ------------------------------------------------------------------------------


exp_hg <- sort(colSums(msr_mat_hg), decreasing = TRUE)
exp_mm <- sort(colSums(msr_mat_mm), decreasing = TRUE)


exp_df <- data.frame(
  Gene_count = c(exp_hg, exp_mm),
  Symbol = c(names(exp_hg), names(exp_mm)),
  Species = c(rep("Human", length(exp_hg)), rep("Mouse", length(exp_mm)))
)


exp_summ <- lapply(split(exp_df, exp_df$Species), function(x) summary(x$Gene_count))


# Plots
# ------------------------------------------------------------------------------


# Binary heatmap of gene measurement: gene x experiment where colour denotes
# that the gene was measured in the given experiment


pheatmap(
  msr_mat_hg,
  col = c("black", "royalblue"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  height = 12,
  width = 12,
  treeheight_row = 0,
  treeheight_col = 0,
  filename = file.path(plot_dir, "measurement_heatmap_human.png")
)


pheatmap(
  msr_mat_mm,
  col = c("black", "royalblue"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  height = 12,
  width = 12,
  treeheight_row = 0,
  treeheight_col = 0,
  filename = file.path(plot_dir, "measurement_heatmap_mouse.png")
)


# Look at experiment counts of gene measurement 

p1 <- ggplot(gene_msr_df, aes(x = Proportion_msr)) +
  facet_wrap(~Species, nrow = 2) +
  geom_histogram(bins = 100) +
  ylab("Count of genes") +
  xlab("Proportion of experiments with measurement") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 25))
        

ggsave(p1, height = 9, width = 12, device = "png", dpi = 300,
       filename = file.path(plot_dir, "gene_proportion_measurement_hist.png"))


# Looking at the proportion of gene measurement split by TF status

p2 <- ggplot(gene_msr_df, aes(y = Proportion_msr, x = TF)) +
  facet_wrap(~Species) +
  geom_boxplot() +
  ylab("Proportion of experiments with measurement") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 25))


ggsave(p2, height = 9, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, "TF_proportion_measurement_boxplot.png"))


# Look at experiment counts of gene measurement 

p3 <- ggplot(exp_df, aes(x = Gene_count)) +
  facet_wrap(~Species, nrow = 2) +
  geom_histogram(bins = 20) +
  ylab("Count of experiments") +
  xlab("Count of genes measured") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 25))

ggsave(p3, height = 9, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, "experiment_gene_measurement_hist.png"))
