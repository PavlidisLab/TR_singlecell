## PCA of all experiments for a given set of TFs, and their aggregate profile
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggrepel)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# IDs for scRNA-seq datasets
ids_hg <- filter(sc_meta, Species == "Human")$ID
ids_mm <- filter(sc_meta, Species == "Mouse")$ID

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Ribosomal genes
sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))

# Saved list RDS of the ranks
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)
rank_ribo_hg <- readRDS(rank_ribo_hg_path)
rank_ribo_mm <- readRDS(rank_ribo_mm_path)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Loading the TF aggregate matrices
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)
agg_ribo_hg <- load_or_generate_agg(path = agg_ribo_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = ribo_genes$Symbol_hg)
agg_ribo_mm <- load_or_generate_agg(path = agg_ribo_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = ribo_genes$Symbol_mm)


# Only keep ortho genes measured in both species

pc_ortho <- filter(pc_ortho,
                   Symbol_hg %in% pc_hg$Symbol &
                     Symbol_mm %in% pc_mm$Symbol)

tfs_ortho <- filter(pc_ortho, 
                    Symbol_hg %in% names(rank_tf_hg) & 
                      Symbol_mm %in% names(rank_tf_mm))


# TODO:

bind_ortho_matrix <- function(gene, 
                              agg_l_hg, 
                              agg_l_mm,
                              msr_mat_hg = msr_hg,
                              msr_mat_mm = msr_mm,
                              pc_df = pc_ortho) {
  
  gene_ortho <- filter(pc_df, Symbol_hg == gene | Symbol_mm == gene)
  if (nrow(gene_ortho) != 1) stop("A one to one match was not found")
  gene_hg <- gene_ortho$Symbol_hg
  gene_mm <- gene_ortho$Symbol_mm
  
  mat_hg <- subset_to_measured(gene_vec_to_mat(agg_l_hg, gene_hg), msr_mat_hg, gene_hg)
  mat_mm <- subset_to_measured(gene_vec_to_mat(agg_l_mm, gene_mm), msr_mat_mm, gene_mm)
  
  mat_hg <- mat_hg[pc_df$Symbol_hg, ]
  mat_mm <- mat_mm[pc_df$Symbol_mm, ]
  mat_ortho <- cbind(mat_hg, mat_mm)
  rownames(mat_ortho) <- pc_df$Symbol_hg
  
  return(mat_ortho)  
}



# Performs PCA with prcomp and returns list of the resulting 
# object as well as the variance explained
# TODO: expected shape of mat

pca_and_var <- function(mat, scale_arg = TRUE) {
  
  # prcomp expects samples as rows, features (genes) as columns so transpose
  pcmat <- prcomp(t(mat), scale = scale_arg)
  
  # variance explained by the PCs
  prc_var <- pcmat$sdev ^ 2
  var_explained <- round(prc_var / sum(prc_var) * 100, 2)
  cumvar_explained <- cumsum(var_explained)/sum(var_explained)
  
  return(list(PC = pcmat, 
              Var_explained = var_explained, 
              Cumvar_explained = cumvar_explained))
}



pc_scatter <- function(df,
                       pc_list,
                       pc_x,
                       pc_y,
                       title = NULL) {
  
  ggplot(df, aes(x = !!sym(paste0("PC", pc_x)),
                 y = !!sym(paste0("PC", pc_y)))) +
    geom_point(aes(fill = Symbol, shape = Species), size = 6) +
    xlab(paste0("PC", pc_x, "(", pc_list$Var_explained[pc_x], "%)")) +
    ylab(paste0("PC", pc_y, "(", pc_list$Var_explained[pc_y], "%)")) +
    ggtitle(title) +
    scale_shape_manual(values = c(22, 24)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    theme_classic() +
    theme(axis.title = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
}



gene1 <- "PAX6"
gene2 <- "RUNX1"
gene3 <- "ASCL1"


gene_mat1 <- bind_ortho_matrix(gene = gene1,
                               agg_l_hg = agg_tf_hg,
                               agg_l_mm = agg_tf_mm)


gene_mat2 <- bind_ortho_matrix(gene = gene2,
                               agg_l_hg = agg_tf_hg,
                               agg_l_mm = agg_tf_mm)


gene_mat3 <- bind_ortho_matrix(gene = gene3,
                               agg_l_hg = agg_tf_hg,
                               agg_l_mm = agg_tf_mm)


# gene_mat3 <- bind_ortho_matrix(gene = gene3,
#                                agg_l_hg = agg_ribo_hg,
#                                agg_l_mm = agg_ribo_mm)


ortho_mat <- cbind(gene_mat1, gene_mat2, gene_mat3)


ortho_pca <- pca_and_var(ortho_mat, scale_arg = FALSE)


npcs <- 10

# top genes
# sort(abs(ortho_pca$PC$rotation[, 1]), decreasing = TRUE)[1:50]
# sort(abs(ortho_pca$PC$rotation[, 2]), decreasing = TRUE)[1:50]


pca_df <- data.frame(
  ortho_pca$PC$x[, 1:npcs],
  Symbol = c(rep(gene1, ncol(gene_mat1)), 
             rep(gene2, ncol(gene_mat2)), 
             rep(gene3, ncol(gene_mat3))),
  ID = colnames(ortho_mat)
) %>% 
  left_join(sc_meta[, c("ID", "Species")], by = "ID")



p1a <- pc_scatter(df = pca_df, pc_list = ortho_pca, pc_x = 1, pc_y = 2)
p1b <- pc_scatter(df = pca_df, pc_list = ortho_pca, pc_x = 2, pc_y = 3)
p1c <- pc_scatter(df = pca_df, pc_list = ortho_pca, pc_x = 3, pc_y = 4) 
p1d <- pc_scatter(df = pca_df, pc_list = ortho_pca, pc_x = 4, pc_y = 5) 
p1_leg <- get_legend(p1d)
p1d <- p1d + theme(legend.position = "none")
p1 <- plot_grid(p1a, p1b, p1c, p1d, nrow = 2)
p1 <- plot_grid(p1, p1_leg, rel_widths = c(2, 0.3))



# genes <- c("ASCL1", "RUNX1", "PAX6", "HES1", "MEF2C", "NEUROD1", "MECP2", "TCF4")
# genes <- c(gene1, gene2, gene3)
# genes <- slice_max(ortho_df, Cor, n = 8)$Symbol
genes <- slice_max(ortho_df, Topk_count, n = 8)$Symbol


# Rank has self-gene removed. Need to add back in

b1 <- lapply(genes, function(x) {
  
  df <- rank_tf_ortho[[x]][, c("Symbol_hg", "Avg_RSR_hg", "Avg_RSR_mm")]
  colnames(df) <- c("Symbol", paste0(x, c("_Human", "_Mouse")))
  
  df <- rbind(df, c(x, 1, 1))
  
  df <- df %>% 
    arrange(match(Symbol, pc_ortho$Symbol_hg)) %>% 
    select(-Symbol) %>% 
    mutate_if(is.character, as.numeric)
  
  
  return(df)
})

b2 <- do.call(cbind, b1)
rownames(b2) <- pc_ortho$Symbol_hg


agg_pca <- pca_and_var(b2, scale_arg = FALSE)


agg_pca_df <- data.frame(
  agg_pca$PC$x[, 1:length(genes)],
  Symbol = rep(genes, each = 2),
  Species = rep(c("Human", "Mouse"), times = length(genes))
) 


pc_scatter(df = agg_pca_df, pc_list = agg_pca, pc_x = 1, pc_y = 2)
pc_scatter(df = agg_pca_df, pc_list = agg_pca, pc_x = 3, pc_y = 4)
pc_scatter(df = agg_pca_df, pc_list = agg_pca, pc_x = 1, pc_y = 3)


tt1 <- umap::umap(t(b2))
tt2 <- data.frame(tt1$layout)
colnames(tt2) <- c("UMAP1", "UMAP2")
tt3 <- cbind(tt2, agg_pca_df)


ggplot(tt3, aes(x = UMAP1, y = UMAP2)) +
  geom_point(aes(fill = Symbol, shape = Species), size = 6) +
  scale_shape_manual(values = c(22, 24)) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
