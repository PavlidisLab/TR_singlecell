library(WGCNA)
library(tidyverse)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")

id <- "GSE180928"
sc_dir <- paste0("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/alex_sc_requests/human/has_celltype_metadata/", id)
dat_path <- file.path(sc_dir, paste0(id, "_filtered_cell_counts.csv"))
meta_path <- file.path(sc_dir, paste0(id, "_metadata.csv"))
out_path <- paste0("/space/scratch/amorin/R_objects/", id, "_mat_and_meta.RDS")
pc <- read.delim("/home/amorin/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)


# This is for cellxgene data
# meta <- dat@meta.data
# mat <- GetAssayData(dat, slot = "counts")
# plot(meta$nCount_RNA, meta$nFeature_RNA)
# meta["GCGCGATCATACGCCG-alexsc", "nFeature_RNA"]
# sum(mat[, "GCGCGATCATACGCCG-alexsc"] != 0)
# rna_complexity <- log10(meta$nFeature_RNA) / log10(meta$nCount_RNA)
# plot(density(rna_complexity))
# filter_ix_numi <- which(meta$nCount_RNA < 500)
# filter_ix_ngene <- which(meta$nFeature_RNA < 250)
# filter_ix_complexity <- which(rna_complexity < 0.8)
#



# Add metadata columns corresponding to the count of UMIs per cell, the count
# of non-zero genes, and the RNA complexity/novelty score
# https://hbctraining.github.io/scRNA-seq_online/lessons/04_SC_quality_control.html

umi_counts <- colSums(mat)
gene_counts <- colSums(mat > 0)
complexity <- log10(gene_counts) / log10(umi_counts)


mat[1:5, 1:10]
hist(umi_counts, breaks = 100)
hist(gene_counts, breaks = 100)
plot(density(complexity))


# TODO: filter logic by focusing on the lowest counts?
hi_counts <- sort(cell_counts, decreasing = TRUE)[1:40e3]
hist(hi_counts, breaks = 100)
# low_counts <- sort(gene_counts, decreasing = FALSE)[1:(length(keep_genes) - 40)]
low_counts <- sort(cell_counts, decreasing = FALSE)[1:20e3]
hist(low_counts, breaks = 100)



p1a <- ggplot(data.frame(Cell_counts = cell_counts), aes(x = Cell_counts)) +
  geom_histogram(bins = 50) +
  xlab("Column (cell)-wise counts") +
  ylab("nCells") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

p1b <- ggplot(data.frame(Gene_counts = gene_counts), aes(x = Gene_counts)) +
  geom_histogram(bins = 50) +
  xlab("Row (gene)-wise counts") +
  ylab("nGenes") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))

p1 <- cowplot::plot_grid(p1a, p1b)

# ggsave(p1, device = "png", dpi = 300, height = 9, width = 15,
#        filename = file.path(plot_dir, "umi_count_histograms.png"))


#


if (!file.exists(out_path)) {
  
  dat <- read.delim(dat_path, sep = ",")
  meta <- read.delim(meta_path, sep = ",")
  
  # Preparing count matrix
  rownames(dat) <- dat$X
  dat$X <- NULL
  colnames(dat) <- str_replace_all(colnames(dat), "\\.", "_")
  
  mat <- as.matrix(dat)
  mat <- get_pcoding_only(mat, pc)
  
  # Ready metadata
  meta <- meta %>% 
    dplyr::rename(ID = X, Cell_type = Cluster) %>% 
    mutate(ID = str_replace_all(ID, "-", "_"))
  
  stopifnot(all(colnames(mat) %in% meta$ID))
  
  
  
  saveRDS(list(mat, meta), file = out_path)
  
  rm(dat)
  gc()
  
} else {
  
  dat <- readRDS(out_path)
  meta <- dat[[2]]
  mat <- dat[[1]]
  
}


stopifnot(identical(colnames(mat), meta$ID))


rsr1 <- all_RSR_aggregate1(mat, meta)
saveRDS(rsr1, file = "/space/scratch/amorin/R_objects/GSE180928_RSR1.RDS")

z1 <- all_zscore_aggregate(mat, meta)
saveRDS(z1, file = "/space/scratch/amorin/R_objects/GSE180928_Z1.RDS")

rsr2 <- all_RSR_aggregate2(mat, meta)
saveRDS(rsr2, file = "/space/scratch/amorin/R_objects/GSE180928_RSR2.RDS")

rsr3 <- all_RSR_aggregate3(mat, meta)
saveRDS(rsr3, file = "/space/scratch/amorin/R_objects/GSE180928_RSR3.RDS")

rsr4 <- all_RSR_aggregate4(mat, meta)
saveRDS(rsr4, file = "/space/scratch/amorin/R_objects/GSE180928_RSR4.RDS")

rsr5 <- all_RSR_aggregate5(mat, meta)
saveRDS(rsr5, file = "/space/scratch/amorin/R_objects/GSE180928_RSR5.RDS")

rsr6 <- all_RSR_aggregate6(mat, meta)
saveRDS(rsr6, file = "/space/scratch/amorin/R_objects/GSE180928_RSR6.RDS")



# rsr1 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR1.RDS")
# rsr2 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR2.RDS")
# rsr3 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_RSR3.RDS")
# z1 <- readRDS("/space/scratch/amorin/R_objects/GSE180928_Z1.RDS")

