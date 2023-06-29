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

# Protein coding genes
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Ranked targets from paper
evidence_l <- readRDS(evidence_path)



# Load matrices that track counts of NA pairs (lack of mutual measurement)
# across datasets/ids into a list

load_na_mat_list <- function(ids,
                             dir = "/space/scratch/amorin/TR_singlecell/",
                             genes,
                             sub_genes = NULL) {
  
  mat_l <- lapply(ids, function(x) {
    
    path <- file.path(amat_dir, x, paste0(x, "_NA_mat.tsv"))
    
    if (!is.null(sub_genes)) {
      
      dat <- fread(path, sep = "\t", select = c("V1", sub_genes))
      mat <- as.matrix(dat[, -1, drop = FALSE])
      rownames(mat) <- dat$V1
      colnames(mat) <- sub_genes
      mat <- mat[genes, sub_genes]
   
     } else {
       
       dat <- fread(path, sep = "\t")
       mat <- as.matrix(dat[, -1, drop = FALSE])
       rownames(mat) <- colnames(mat) <- dat$V1
       mat <- mat[genes, genes]
     }
    
    return(mat)
  })
  
  names(mat_l) <- ids
  gc(verbose = FALSE)
  
  return(mat_l)
}



# Create a binary gene by experiment matrix that tracks whether a given gene
# was measured (at least one count in at least 20 cells in at least one cell
# type) for each dataset/id.

get_gene_msr_mat <- function(ids, meta, genes) {
  
  msr_mat <- matrix(0, nrow = length(genes), ncol = length(ids))
  colnames(msr_mat) <- ids
  rownames(msr_mat) <- genes
  
  for (id in ids) {
    
    na_mat <- load_na_mat_list(id, genes = genes)[[1]]
    n_celltype <- filter(meta, ID == id)$N_celltype
    msr_mat[, id] <- as.integer(diag(na_mat) != n_celltype)
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



# Get the average/proportion of measurement across experiments
# ------------------------------------------------------------------------------


avg_msr_hg <- sort(rowMeans(msr_mat_hg), decreasing = TRUE)
avg_msr_mm <- sort(rowMeans(msr_mat_mm), decreasing = TRUE)


never_msr_hg <- avg_msr_hg[avg_msr_hg == 0]
never_msr_mm <- avg_msr_mm[avg_msr_mm == 0]


always_msr_hg <- avg_msr_hg[avg_msr_hg == 1]
always_msr_mm <- avg_msr_mm[avg_msr_mm == 1]



# Focus on TFs


tf_avg_msr_hg <- avg_msr_hg[tfs_hg$Symbol]
tf_avg_msr_mm <- avg_msr_mm[tfs_mm$Symbol]


avg_msr_df <- data.frame(
  Prop_msr = c(avg_msr_hg, avg_msr_mm),
  Symbol = c(names(avg_msr_hg), names(avg_msr_mm)),
  Species = c(rep("Human", length(avg_msr_hg)), rep("Mouse", length(avg_msr_mm))),
  TF = c(names(avg_msr_hg) %in% tfs_hg$Symbol, names(avg_msr_mm) %in% tfs_mm$Symbol)
)


p1a <- ggplot(avg_msr_df, aes(x = Prop_msr, fill = TF)) +
  facet_wrap(~Species) +
  geom_density(alpha = 0.4, position = "stack") +
  theme_classic()
  

p1b <- ggplot(avg_msr_df, aes(y = Prop_msr, x = TF)) +
  facet_wrap(~Species) +
  geom_boxplot(alpha = 0.4) +
  theme_classic()



# Look at experiment-wise gene coverage
# ------------------------------------------------------------------------------


exp_hg <- sort(colSums(msr_mat_hg), decreasing = TRUE)
exp_mm <- sort(colSums(msr_mat_mm), decreasing = TRUE)


exp_df <- data.frame(
  Gene_count = c(exp_hg, exp_mm),
  Symbol = c(names(exp_hg), names(exp_mm)),
  Species = c(rep("Human", length(exp_hg)), rep("Mouse", length(exp_mm)))
)


p2 <- ggplot(exp_df, aes(x = Gene_count)) +
  facet_wrap(~Species) +
  geom_histogram(bins = 20) +
  theme_classic()




# Examining instances of genes with strong regulation evidence but are never
# measured with the TF itself
# ------------------------------------------------------------------------------


# Loading

na_hg <- load_na_mat_list(ids_hg, genes = pc_hg$Symbol, sub_genes = names(evidence_l$Human))
na_mm <- load_na_mat_list(ids_mm, genes = pc_mm$Symbol, sub_genes = names(evidence_l$Mouse))

stopifnot(all(unlist(lapply(na_hg, function(x) identical(rownames(x), pc_hg$Symbol)))))
stopifnot(all(unlist(lapply(na_mm, function(x) identical(rownames(x), pc_mm$Symbol)))))

# Convert count matrices to proportions

n_ct_hg <- filter(sc_meta, ID %in% ids_hg)$N_celltype
n_ct_mm <- filter(sc_meta, ID %in% ids_mm)$N_celltype

na_hg <- lapply(1:length(na_hg), function(x) sweep(na_hg[[x]], 2, n_ct_hg[x], `/`))
names(na_hg) <- ids_hg
na_mm <- lapply(1:length(na_mm), function(x) sweep(na_mm[[x]], 2, n_ct_mm[x], `/`))
names(na_mm) <- ids_mm

# Extract a single TF's vector across datasets ino a single matrix

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

# hist(gene_order_hg, breaks = 100)
# hist(gene_order_mm, breaks = 100)

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


# cor(gene_evidence_hg[[evidence_col]], gene_evidence_hg$NA_prop, method = "spearman")
# cor(gene_evidence_mm[[evidence_col]], gene_evidence_mm$NA_prop, method = "spearman")


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
  filter(NA_prop == 1) %>% 
  slice_min(Rank_integrated, n = 20) %>% 
  select(Symbol, Count_DE, Count_NA, Rank_binding, Rank_perturbation, Rank_integrated)


top_mm <- gene_evidence_mm %>% 
  filter(NA_prop == 1) %>% 
  slice_min(Rank_integrated, n = 20) %>% 
  select(Symbol, Count_DE, Count_NA, Rank_binding, Rank_perturbation, Rank_integrated)


rowSums(msr_mat_hg[top_hg$Symbol, ])
rowSums(msr_mat_mm[top_mm$Symbol, ])
