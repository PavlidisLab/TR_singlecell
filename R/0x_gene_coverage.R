## Inspect NA pairs of TFs and target genes
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

# pcoding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Ranked targets from paper
evidence_l <- readRDS(evidence_path)
tfs_hg <- names(evidence_l$Human)
tfs_mm <- names(evidence_l$Mouse)


# TODO

load_na_mat_list <- function(ids,
                             dir = "/space/scratch/amorin/TR_singlecell/",
                             sub_genes = NULL) {
  
  if (!is.null(sub_genes)) sub_genes <- c("V1", sub_genes)
  
  mat_l <- lapply(ids, function(x) {
    path <- file.path(amat_dir, x, paste0(x, "_NA_mat.tsv"))
    dat <- fread(path, sep = "\t", select = sub_genes)
    mat <- as.matrix(dat[, -1, drop = FALSE])
    rownames(mat) <- dat$V1
    mat
  })
  
  names(mat_l) <- ids
  gc(verbose = FALSE)
  
  return(mat_l)
}



get_gene_msr_mat <- function(ids, meta, genes) {
  
  msr_mat <- matrix(0, nrow = length(genes), ncol = length(ids))
  colnames(msr_mat) <- ids
  rownames(msr_mat) <- genes
  
  for (id in ids) {
    
    na_mat <- load_na_mat_list(id)[[1]][genes, ]
    n_celltype <- filter(meta, ID == id)$N_celltype
    msr_mat[, id] <- as.integer(diag(na_mat) != n_celltype)
    rm(na_mat)
    gc(verbose = FALSE)
    
  }
  
  return(msr_mat)
}




if (!file.exists(msr_mat_hg_path)) {
  msr_mat_hg <- get_gene_msr_mat(ids_hg, sc_meta, pc_hg$Symbol)
} else {
  msr_mat_hg <- readRDS(msr_mat_hg_path)
}


if (!file.exists(msr_mat_mm_path)) {
  msr_mat_mm <- get_gene_msr_mat(ids_mm, sc_meta, pc_mm$Symbol)
} else {
  msr_mat_mm <- readRDS(msr_mat_mm_path)
}





na_hg <- load_na_mat_list(ids_hg, sub_genes = tfs_hg)
na_mm <- load_na_mat_list(ids_mm, sub_genes = tfs_mm)


genes_hg <- rownames(na_hg[[1]])
genes_mm <- rownames(na_mm[[1]])


# TODO: this needs to be replaced with upstream ordering
na_hg <- lapply(na_hg, function(x) x[genes_hg, ])
na_mm <- lapply(na_mm, function(x) x[genes_mm, ])


stopifnot(all(unlist(lapply(na_hg, function(x) identical(rownames(x), genes_hg)))))
stopifnot(all(unlist(lapply(na_mm, function(x) identical(rownames(x), genes_mm)))))


# Convert count matrices to proportions

n_ct_hg <- filter(sc_meta, ID %in% ids_hg)$N_celltype
n_ct_mm <- filter(sc_meta, ID %in% ids_mm)$N_celltype


na_hg <- lapply(1:length(na_hg), function(x) sweep(na_hg[[x]], 2, n_ct_hg[x], `/`))
names(na_hg) <- ids_hg

na_mm <- lapply(1:length(na_mm), function(x) sweep(na_mm[[x]], 2, n_ct_mm[x], `/`))
names(na_mm) <- ids_mm


# Example of describing a single gene's coverage over all experiments
# TODO: check that ensures naming consistency

gene_hg <- "ASCL1"
gene_mm <- str_to_title(gene_hg)
stopifnot(gene_mm %in% genes_mm)

gene_mat_hg <- gene_vec_to_mat(na_hg, gene_hg)
gene_mat_mm <- gene_vec_to_mat(na_mm, gene_mm)

# If the gene itself equals 1, then it is not expressed at all in the data set

names(which(gene_mat_hg[gene_hg, ] == 1))
names(which(gene_mat_mm[gene_mm, ] == 1))

hist(gene_mat_hg[gene_hg, ], breaks = 30)
hist(gene_mat_mm[gene_mm, ], breaks = 30)

# Order genes by their average measurement with the gene of interest. 
# A value of 1 means that the gene pair were never measured together (possibly
# because the gene itself isn't expressed). Lower values mean that a gene pair
# are measured in the same cell types

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
plot(gene_evidence_hg[[evidence_col]], gene_evidence_hg$NA_prop)


cor(gene_evidence_mm[[evidence_col]], gene_evidence_mm$NA_prop, method = "spearman")
plot(gene_evidence_mm[[evidence_col]], gene_evidence_mm$NA_prop)


gene_evidence_hg$Group_evidence <- cut(gene_evidence_hg[[evidence_col]], breaks = 50)

ggplot(gene_evidence_hg, aes(x = Group_evidence, y = NA_prop)) + 
  geom_boxplot() +
  xlab("Binned regulation evidence") + 
  ylab("Proportion of NA measurements") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 20))




gene_evidence_mm$Group_evidence <- cut(gene_evidence_mm[[evidence_col]], breaks = 50)

ggplot(gene_evidence_mm, aes(x = Group_evidence, y = NA_prop)) + 
  geom_boxplot() +
  xlab("Binned regulation evidence") + 
  ylab("Proportion of NA measurements") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 20))


# What is the highest ranked gene that is not ever detected with the gene?


gene_evidence_hg %>% 
  filter(NA_prop == 1) %>% 
  slice_min(Rank_integrated, n = 20) %>% 
  select(Symbol, Count_DE, Count_NA, Rank_binding, Rank_perturbation, Rank_integrated)


# And is this gene actually expressed at all? Requires loading the gene of interest
tt1 <- load_na_mat_list(ids_hg, sub_genes = "UBL4B")
tt2 <- gene_vec_to_mat(tt1, "UBL4B")
tt3 <- sweep(tt2, 2, n_ct_hg, `/`)
tt4 <- rowMeans(tt3)
tt4["UBL4B"]





# checking paired nature of NA mat

na_mat <- load_na_mat_list("GSE180928")[[1]]
# colnames(na_mat) <- rownames(na_mat)
na_mat[1:10, 1:10]
isSymmetric(na_mat, check.attributes = FALSE)





