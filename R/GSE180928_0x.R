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





# TODO: a warning if 0 is not min
# TODO: this is costly - rowsums. If 0 not min, then use apply with min (with flag)
zero_genes <- names(which(apply(mat, 1, function(x) sum(x != 0)) == 0))
zero_genes <- which(rowSums(mat) == 0)
keep_genes <- unique(setdiff(rownames(mat), zero_genes))
keep_genes <- keep_genes[keep_genes != ""]


cell_counts <- colSums(mat)
cell_log_counts <- colSums(log10(mat + 1))
gene_counts <- rowSums(mat)
gene_log_counts <- rowSums(log10(mat + 1))


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
  
  # First column is gene symbols, move to rownames and clean col names
  rownames(dat) <- dat$X
  dat$X <- NULL
  colnames(dat) <- str_replace_all(colnames(dat), "\\.", "_")
  mat <- as.matrix(dat)
  
  # Ready metadata
  meta <- meta %>% 
    dplyr::rename(ID = X, Cell_type = Cluster) %>% 
    mutate(ID = str_replace_all(ID, "-", "_"))
  
  stopifnot(all(colnames(mat) %in% meta$ID))
  
  # Only protein coding genes
  common_genes <- intersect(pc$Symbol, rownames(mat))
  mat <- mat[common_genes, meta$ID]
  
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

