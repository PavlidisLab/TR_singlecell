## Summarizing and plotting details about gene measurement and cell counts
## -----------------------------------------------------------------------------

library(tidyverse)
library(ComplexHeatmap)
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

# Binary measurement matrices
msr_mat_hg <- readRDS(msr_mat_hg_path)
msr_mat_mm <- readRDS(msr_mat_mm_path)


# Describing metadata
# ------------------------------------------------------------------------------


# Gene-wise summaries of measurement across experiments

count_msr_hg <- rowSums(msr_mat_hg)
count_msr_mm <- rowSums(msr_mat_mm)

prop_msr_hg <- count_msr_hg / ncol(msr_mat_hg)
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


# Organize counts

n_msr <- do.call(rbind, list(
  Never = lapply(never_msr, nrow),
  Rarely = lapply(rarely_msr, nrow),
  Always = lapply(always_msr, nrow),
  Mostly = lapply(mostly_msr, nrow)
))



# Median count of datasets that a TF is measured in

tf_med_count <- gene_msr_df %>% 
  group_by(Species, TF) %>% 
  summarise(Median_count = median(Count_msr), .groups = "keep")


# Do TFs tend to be measured differentially compared to non-TF genes?

tf_wilx_count <- gene_msr_df %>% 
  group_by(Species) %>% 
  do(W = wilcox.test(Count_msr ~ TF, data = ., paired = FALSE)) %>% 
  summarise(Species, Wilcox = W$p.value)


# Summarize count of cells and cell types across experiments
summ_cells <- lapply(split(sc_meta, sc_meta$Species), function(x) summary(x$N_cells))
summ_cts <- lapply(split(sc_meta, sc_meta$Species), function(x) summary(x$N_celltypes))


# Summarize gene measurement across experiments
summ_gene_msr <- lapply(split(sc_meta, sc_meta$Species), function(x) summary(x$N_genes))


# Plot 
# ------------------------------------------------------------------------------


# Barchart of data set counts by species

p1 <- count(sc_meta, Species) %>% 
  ggplot(., aes(x = Species, y = n)) +
  geom_bar(stat = "identity", fill = "slategrey", col = "black", width = 0.6) +
  ylab("Count of datasets") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank(),
        plot.margin = margin(10, 5, 5, 5))

ggsave(p1, height = 6, width = 4.5, device = "png", dpi = 300,
       filename = file.path(plot_dir, "dataset_counts_by_species.png"))


# Histogram of log10 cell counts by species

p2 <- ggplot(sc_meta, aes(x = log10(N_cells), fill = Species, colour = Species)) +
  geom_histogram(bins = 30) +
  ylab("Count of datasets") +
  xlab(bquote(~Log[10]~ "count of cells")) +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  scale_colour_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.position = c(0.85, 0.85),
        plot.margin = margin(10, 10, 5, 5))

ggsave(p2, height = 6, width = 9, device = "png", dpi = 600,
       filename = file.path(plot_dir, "cell_counts_by_species.png"))


# Histogram of cell type counts by species

p3 <- ggplot(sc_meta, aes(x = N_celltypes, fill = Species)) +
  geom_histogram(bins = 50) +
  ylab("Count of datasets") +
  xlab("Count of cell types") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.position = c(0.85, 0.85),
        plot.margin = margin(10, 10, 5, 5))

ggsave(p3, height = 6, width = 9, device = "png", dpi = 600,
       filename = file.path(plot_dir, "celltype_counts_by_species.png"))


# Look at experiment counts of gene measurement 

p4 <- ggplot(gene_msr_df, aes(x = Proportion_msr)) +
  facet_wrap(~Species, nrow = 2) +
  geom_histogram(bins = 100) +
  ylab("Count of genes") +
  xlab("Proportion of experiments with measurement") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 25))
        

ggsave(p4, height = 9, width = 12, device = "png", dpi = 300,
       filename = file.path(plot_dir, "gene_proportion_measurement_hist.png"))


# Looking at the proportion of gene measurement split by TF status

p5 <- ggplot(gene_msr_df, aes(y = Proportion_msr, x = TF)) +
  facet_wrap(~Species) +
  geom_boxplot() +
  ylab("Proportion of experiments with measurement") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 25))


ggsave(p5, height = 9, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, "TF_proportion_measurement_boxplot.png"))


# Look at experiment counts of gene measurement 

p6 <- ggplot(sc_meta, aes(x = N_genes)) +
  facet_wrap(~Species, nrow = 2) +
  geom_histogram(bins = 20) +
  ylab("Count of experiments") +
  xlab("Count of genes measured") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 25))

ggsave(p6, height = 9, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, "experiment_gene_measurement_hist.png"))



# Binary heatmap of gene measurement: gene x experiment where colour denotes
# that the gene was measured in the given experiment


p7a <- Heatmap(msr_mat_hg,
               col = c("black", "royalblue"),
               heatmap_legend_param = list(
                 title = "Measured",
                 labels = c("TRUE", "FALSE")
               ),
               cluster_rows = TRUE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               show_column_names = FALSE,
               show_row_dend = FALSE,
               show_column_dend = FALSE,
               use_raster = FALSE)



png(file = file.path(plot_dir, "measurement_heatmap_human.png"), 
    res = 1200, 
    height = 6,
    width = 6,
    units = "in")

draw(p7a)
dev.off()



p7b <- Heatmap(msr_mat_mm,
               col = c("black", "royalblue"),
               heatmap_legend_param = list(
                 title = "Measured",
                 labels = c("TRUE", "FALSE")
               ),
               cluster_rows = TRUE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               show_column_names = FALSE,
               show_row_dend = FALSE,
               show_column_dend = FALSE,
               use_raster = FALSE)


png(file = file.path(plot_dir, "measurement_heatmap_mouse.png"), 
    res = 1200, 
    height = 6,
    width = 6,
    units = "in")

draw(p7b)
dev.off()
