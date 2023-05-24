## Process and save data as a Seurat object.
## https://www.biorxiv.org/content/biorxiv/early/2022/10/25/2022.10.25.513733.full.pdf
## -----------------------------------------------------------------------------

library(tidyverse)
library(Seurat)
library(parallel)
library(future)
library(cowplot)
source("R/utils/functions.R")
source("R/00_config.R")

plot_dir <- file.path(plot_dir, "Hochgerner2022")

# Where the output should be saved
seurat_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_seurat.RDS"
marker_path <- "/space/scratch/amorin/R_objects/Hochgerner2022_markers.RDS"

# Number of PCs to use 
n_pc <- 50

# Loading the text file containing the counts and metadata 
dat_path <- "/space/scratch/Hochgerner2022_amygdala_biorvx/Amy_FC_allcells_with_metadata_31-Jul-2022.txt"
dat <- read.delim(dat_path, stringsAsFactors = FALSE)


# Restructuring/isolating metadata and count matrix
# ------------------------------------------------------------------------------


meta <- data.frame(
  Cell_ID = colnames(dat[1, 2:ncol(dat)]),
  Cell_type = unname(unlist(dat[1, 2:ncol(dat), drop = TRUE])),
  Sample = unname(unlist(dat[2, 2:ncol(dat), drop = TRUE])),
  FC_time = unname(unlist(dat[3, 2:ncol(dat), drop = TRUE])),
  Batch = unname(unlist(dat[4, 2:ncol(dat), drop = TRUE]))
)


genes <- dat[5:nrow(dat), 1]

mat <- apply(as.matrix(dat[5:nrow(dat), 2:ncol(dat)]), 2, as.numeric)

rownames(mat) <- genes


stopifnot(identical(
  mat["Dnpep", "GCCAGGTCATGGGAAC.1_76.1"],
  as.numeric(filter(dat, cellID == "Dnpep")$GCCAGGTCATGGGAAC.1_76.1)
))


stopifnot(all(meta$Cell_ID %in% colnames(mat)))

rownames(meta) <- meta$Cell_ID


# Remove bulky unprocessed data
rm(dat)


# Creating Seurat object: Provided counts but uncertain what (if any) processing
# was involved. The dataset has 3 batches that have already been merged into one
# count matrix. Going to just use default processing noted in Seurat tutorial.
# ------------------------------------------------------------------------------


sdat <- CreateSeuratObject(counts = mat, 
                           meta.data = meta, 
                           min.cells = 3,
                           min.features = 200,
                           project = "Hochgerner2022")

genes <- rownames(sdat)


# Count summaries. Note bimodal nature of cell UMI counts


cell_counts <- colSums(sdat@assays$RNA@counts)
gene_counts <- log10(rowSums(sdat@assays$RNA@counts) + 1)


p1a <- ggplot(data.frame(Cell_counts = cell_counts), aes(x = Cell_counts)) +
  geom_histogram(bins = 50) +
  xlab("Count of UMIs within cells") +
  ylab("Count of cells") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))
  
p1b <- ggplot(data.frame(Gene_counts = gene_counts), aes(x = Gene_counts)) +
  geom_histogram(bins = 50) +
  xlab("log10 Count of UMIs across cells") +
  ylab("Count of genes") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

p1 <- plot_grid(p1a, p1b)

ggsave(p1, device = "png", dpi = 300, height = 9, width = 15,
       filename = file.path(plot_dir, "umi_count_histograms.png"))


# Mt genes: see extensive amount of mitochondrial transcripts. Seurat tutorial
# removes cells with >5% mt counts... that would remove 73% of cells here!

sdat <- PercentageFeatureSet(sdat, pattern = "^mt-", col.name = "percent_mt")

message(paste("Proportion mitochondrial:",
              round(sum(sdat$percent_mt > 5) / ncol(sdat), 3)))

p2 <- VlnPlot(sdat,
              features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
              ncol = 3)

ggsave(p2, device = "png", dpi = 300, height = 9, width = 15,
       filename = file.path(plot_dir, "qc_feature_vplot.png"))


# Default normalization
sdat <- NormalizeData(sdat, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)


# Highly variable genes
sdat <- FindVariableFeatures(sdat, 
                             selection.method = "vst", 
                             nfeatures = 2000)

top10 <- head(VariableFeatures(sdat), 10)

p3 <- LabelPoints(plot = VariableFeaturePlot(sdat),
                  points = top10,
                  repel = TRUE)


ggsave(p3, device = "png", dpi = 300, height = 9, width = 15,
       filename = file.path(plot_dir, "var_genes_scatter.png"))


# Scaling
sdat <- ScaleData(sdat, features = genes)


# PCA
sdat <- RunPCA(sdat, npcs = n_pc)

p4 <- ElbowPlot(sdat, ndims = n_pc)

ggsave(p4, device = "png", dpi = 300, height = 6, width = 6,
       filename = file.path(plot_dir, paste0("PC_elbowplot_npc=", n_pc, ".png")))


# Cluster + UMAP
sdat <- FindNeighbors(sdat, dims = 1:n_pc)
sdat <- FindClusters(sdat, resolution = 0.5)
sdat <- RunUMAP(sdat, dims = 1:n_pc)


# Diff expr, setting identities to provided labels
Idents(sdat) <- sdat$Cell_type
marker_all <- FindAllMarkers(sdat)

# Save out
saveRDS(sdat, seurat_path)
saveRDS(marker_all, marker_path)
