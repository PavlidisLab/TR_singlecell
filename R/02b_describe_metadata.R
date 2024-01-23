## Examine gene measurement/coverage across experiments. A gene is considered
## measured if has at least 1 count in at least 20 cells in at least one cell
## type. Exports a binary gene x experiment matrix that tracks measurement.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

force_resave <- FALSE

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



# Get the average/proportion of gene measurement across experiments
# Human: 450 genes are never measured and 5,155 genes always measured
# Mouse: 887 genes are never measured and 838 genes are always measured
# Ortho: 630 genes measured in every mouse and human dataset
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
# Human: HPA max coverage at 18299 genes, GSE85241 min at 8952
# Mouse: TabulaMuris max coverage at 17977 genes, GSE160193 min at 1988 genes.
# ------------------------------------------------------------------------------


exp_hg <- sort(colSums(msr_mat_hg), decreasing = TRUE)
exp_mm <- sort(colSums(msr_mat_mm), decreasing = TRUE)


exp_df <- data.frame(
  Gene_count = c(exp_hg, exp_mm),
  Symbol = c(names(exp_hg), names(exp_mm)),
  Species = c(rep("Human", length(exp_hg)), rep("Mouse", length(exp_mm)))
)


exp_summ <- lapply(split(exp_df, exp_df$Species), function(x) summary(x$Gene_count))

summ_cells <- lapply(split(meta_loaded, meta_loaded$Species), function(x) summary(x$N_cells))
summ_cts <- lapply(split(meta_loaded, meta_loaded$Species), function(x) summary(x$N_celltypes))



# Plot 
# ------------------------------------------------------------------------------


# Barchart of data set counts by species

p1 <- count(meta_loaded, Species) %>% 
  ggplot(., aes(x = Species, y = n)) +
  # geom_bar(stat = "identity", fill = "#1c9099", col = "black", width = 0.6) +
  geom_bar(stat = "identity", fill = "slategrey", col = "black", width = 0.6) +
  ylab("Count of data sets") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank(),
        plot.margin = margin(10, 5, 5, 5))



ggsave(p1, height = 6, width = 4.5, device = "png", dpi = 300,
       filename = file.path(plot_dir, "dataset_counts_by_species.png"))


# Histogram of log10 cell counts by species

p2 <- ggplot(meta_loaded, aes(x = log10(N_cells), fill = Species, colour = Species)) +
  geom_histogram(bins = 30) +
  ylab("Count of data sets") +
  xlab("Log10 count of cells") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  scale_colour_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = c(0.85, 0.75),
        plot.margin = margin(10, 5, 5, 5))


ggsave(p2, height = 6, width = 9, device = "png", dpi = 600,
       filename = file.path(plot_dir, "cell_counts_by_species.png"))


# Histogram of cell type counts by species

p3 <- ggplot(meta_loaded, aes(x = N_celltypes, fill = Species)) +
  geom_histogram(bins = 50) +
  ylab("Count of data sets") +
  xlab("Count of cell types") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = c(0.85, 0.75))


ggsave(p3, height = 6, width = 9, device = "png", dpi = 600,
       filename = file.path(plot_dir, "celltype_counts_by_species.png"))



# Binary heatmap of gene measurement: gene x experiment where colour denotes
# that the gene was measured in the given experiment


pheatmap(
  msr_mat_hg,
  col = c("black", "royalblue"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,
  legend = FALSE,
  filename = file.path(plot_dir, "measurement_heatmap_human.png")
)


pheatmap(
  msr_mat_mm,
  col = c("black", "royalblue"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  treeheight_col = 0,
  legend = FALSE,
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
