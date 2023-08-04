## Simulate scRNA-seq counts and explore the relationship of coexpression with
## mean expression level across different processing stragies.
## -----------------------------------------------------------------------------


# TODO: scatter plot that shows average expression ~ average cor by group
# TODO: show distn of representative gene from select quantiles
# TODO: show density of all correlations for a representative Q99 gene verus mid
# 5: Set aside raw counts, CPM counts, and lognorm counts
# 6: Generate coexpr matrix for each of the three count matrices
# 7: Apply coexpr correction to each
# 8: For the 6 matrices (raw, CPM, lognorm and +/- correction) plot cor within each quantile



library(tidyverse)
library(WGCNA)
library(Seurat)
library(parallel)
library(spqn)
source("/home/amorin/Projects/TR_singlecell/R/utils/functions.R")

set.seed(5)


# 1: Load data and focus on a single cell type
# Cheng2018: showed good correspondence between lognorm and CPM aggregate coexpr
# HPAP: show poor correspondence
# ------------------------------------------------------------------------------


dat <- readRDS("/space/scratch/amorin/coexpr_sim/cheng2018_rawcounts_meta.RDS")
# dat <- readRDS("/space/scratch/amorin/coexpr_sim/HPAP_rawcounts_meta.RDS")

meta <- dat$Meta
ct <- as.character(unique(meta$Cell_type))[5]
ref_mat <- as.matrix(dat$Mat[, filter(meta, Cell_type == ct)$ID])


# 2: Get average counts per gene and # of 0s. Bin average into quantiles
# ------------------------------------------------------------------------------

gene_mean <- sort(rowMeans(ref_mat))
gene_n0 <- apply(ref_mat, 1, function(x) sum(x == 0))[names(gene_mean)]
sum(gene_n0 == ncol(ref_mat))


# List of genes split in increasing order expression  https://stackoverflow.com/questions/38547998/
nqtls <- 100
gene_qtls <- split(gene_mean, floor(nqtls * seq.int(0, length(gene_mean) - 1) / length(gene_mean)))
names(gene_qtls) <- paste0("Qtl", 1:nqtls)

# What is the first bin to have non-zeros
first_non0 <- which(!unlist(lapply(gene_qtls, function(x) all(x == 0))))[1]


# 3: Simulate expression, drawing from each quantile
# ------------------------------------------------------------------------------


# Given a gene x cell count matrix (ref_mat) and sampled genes, return a 
# sample_genes x n_cells matrix where each gene/row consists of a re-sampling of
# counts from the corresponding sampled gene in ref_mat

sample_counts_to_mat <- function(ref_mat, sample_genes, n_cells) {
  
  sample_vec <- lapply(sample_genes, function(gene) {
    sample(ref_mat[gene, ], n_cells, replace = TRUE)
  })
  
  sample_mat <- do.call(rbind, sample_vec)
  return(sample_mat)
}


# Given a list of grouped genes whose average expression are in ordered quantiles,
# generate a sampled count matrix for each quantile of expression. Bind these
# into a single final matrix

generate_sample_counts <- function(gene_qtls, ref_mat, ncore = 8) {
  
  n_cells <- ncol(ref_mat)
  n_genes <- nrow(ref_mat)
  n_genes_per_qtl <- floor(n_genes / length(gene_qtls))
  
  sample_l <- mclapply(gene_qtls, function(x) {
    sample_genes <- sample(names(x), n_genes_per_qtl, replace = TRUE)
    sample_mat <- sample_counts_to_mat(ref_mat, sample_genes, n_cells)
  }, mc.cores = ncore)
  
  sample_mat <- do.call(rbind, sample_l)
  
  # Names used for tracking are of form "Qtl:[quantile_bin]_[gene_number]
  
  rownames(sample_mat) <- paste0("Qtl", 
                                 rep(1:nqtls, each = n_genes_per_qtl), 
                                 paste0("_", (1:n_genes_per_qtl)))
  
  return(sample_mat)
}


sample_mat <- generate_sample_counts(gene_qtls, ref_mat)


# 4: Normalize both reference and sampled mat, and generate gene-gene cor for all
# ------------------------------------------------------------------------------


get_cor_mat <- function(count_mat) {
  
  t(count_mat) %>% 
    under_min_count_to_na() %>% 
    WGCNA::cor(., nThreads = 8, use = "pairwise.complete.obs") %>%
    na_to_zero()
  
}


ref_mat_cpm <- Seurat::NormalizeData(ref_mat, normalization.method = "RC", scale.factor = 1e6)
sample_mat_cpm <- Seurat::NormalizeData(sample_mat, normalization.method = "RC", scale.factor = 1e6)

ref_mat_lognorm <- Seurat::NormalizeData(ref_mat, normalization.method = "LogNormalize")
sample_mat_lognorm <- Seurat::NormalizeData(sample_mat, normalization.method = "LogNormalize")


mat_l <- list(
  Ref = ref_mat,
  Ref_CPM = ref_mat_cpm,
  Ref_lognorm = ref_mat_lognorm,
  Sample = sample_mat,
  Sample_CPM = sample_mat_cpm,
  Sample_lognorm = sample_mat_lognorm
)

mat_l <- lapply(mat_l, as.matrix)


cor_l <- mclapply(mat_l, get_cor_mat, mc.cores = 8)



# Hacky: map the reference matrix gene names to the same Qtl-x_x name structure
# as the sampled matrices for downstream convenience


map_names <- lapply(1:length(gene_qtls), function(x) {
  
  data.frame(Qtl = paste0("Qtl", x, "_", 1:length(gene_qtls[[x]])),
             Symbol = names(gene_qtls[[x]]))
  
})


map_names <- do.call(rbind, map_names)


rownames(cor_l$Ref) <- colnames(cor_l$Ref) <- map_names$Qtl[match(rownames(cor_l$Ref), map_names$Symbol)]
rownames(cor_l$Ref_CPM) <- colnames(cor_l$Ref_CPM) <- map_names$Qtl[match(rownames(cor_l$Ref_CPM), map_names$Symbol)]
rownames(cor_l$Ref_lognorm) <- colnames(cor_l$Ref_CPM) <- map_names$Qtl[match(rownames(cor_l$Ref_lognorm), map_names$Symbol)]




# 5: Plotting
# ------------------------------------------------------------------------------


get_plot_df <- function(cor_mat, nqtls, ncore = 8) {
  
  df_l <- mclapply(1:nqtls, function(x) {
    ix <- which(str_detect(rownames(cor_mat), paste0("Qtl", x, "_")))
    data.frame(Group = paste0("Qtl", x), mat_to_df(cor_mat[ix, ix]))
  }, mc.cores = ncore)
  
  return(do.call(rbind, df_l))
}



cor_boxplot <- function(plot_df, title) {
  
  ggplot(plot_df, aes(x = Group, y = Value)) +
    geom_boxplot() +
    ggtitle(title) +
    ylab("Pearson's correlation") +
    xlab("Binned genes") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}





plot_df_l <- lapply(cor_l, get_plot_df, nqtls = nqtls)




boxplot_l <- lapply(1:length(plot_df_l), function(x) {
  cor_boxplot(plot_df_l[[x]], title = names(plot_df_l)[x])
}) 
                    
                    



# TODO: why is there is a low quantile group with nonzero cor? Qtl3


# Inspect the averaged correlations for top quantile/bin across all genes, not
# just within the bin. Do same for lower bins and compare


mat_l$Sample

filter(pdf1, str_detect(Row, "Qtl:99") | str_detect(col, "Qtl:99"))
filter(pdf1, str_detect(Row, "Qtl:75") | str_detect(col, "Qtl:75"))






# Inspect a gene's count distribution across different quantiles


plot_df1 <- rbind(
  data.frame(Group = "Q26", Counts = ref_mat[sample(names(gene_qtls$`26`), 1), ]),
  data.frame(Group = "Q50", Counts = ref_mat[sample(names(gene_qtls$`50`), 1), ]),
  data.frame(Group = "Q75", Counts = ref_mat[sample(names(gene_qtls$`75`), 1), ]),
  data.frame(Group = "Q99", Counts = ref_mat[sample(names(gene_qtls$`99`), 1), ])
)


boxplot(plot_df1$Counts ~ plot_df1$Group)




# Inspecting distribution of sorted gene averages
sample_mean <- rowMeans(sample_mat)
plot(sample_mean)
plot(gene_mean)

# Inspect raw versus norm
plot(sample_mat["Qtl:67_67", ],sample_mat_cpm["Qtl:67_67", ])
plot(sample_mat["Qtl:67_67", ],sample_mat_lognorm["Qtl:67_67", ])





