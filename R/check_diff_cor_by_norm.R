## Looking at the difference in gene-gene correlation when generated from raw
## UMI counts vs the cell x gene pre-processing (rankit tranformation on genes)

library(WGCNA)
library(tidyverse)
library(Matrix)
source("R/utils/functions.R")
source("R/utils/agg_functions.R")

id <- "GSE216019"
sc_dir <- paste0("/cosmos/data/downloaded-data/sc_datasets_w_supplementary_files/lab_projects_datasets/amorin_sc_datases/Human/", id)
dat_path <- file.path(sc_dir, paste0(id, ".RDS"))
pc <- read.delim("/home/amorin/Data/Metadata/ensembl_human_protein_coding_105.tsv", stringsAsFactors = FALSE)

# These are the rank sum rank aggregate matrices using all rank
rsr_norm <- readRDS("/space/scratch/amorin/R_objects/GSE216019_RSR2.RDS")
rsr_counts <- readRDS("/space/scratch/amorin/R_objects/GSE216019_RAW_RSR2.RDS")

# These are the gene-gene cor matrices generated across all cells
# Pearson/Spearman cor, and using raw counts or cxg (rankit) normalization
all_pcor_counts <- readRDS("/space/scratch/amorin/R_objects/GSE216019_pcor_all_counts.RDS")
all_pcor_norm <- readRDS("/space/scratch/amorin/R_objects/GSE216019_pcor_all_norm.RDS")
# all_scor_counts <- readRDS("/space/scratch/amorin/R_objects/GSE216019_scor_all_counts.RDS")
# all_scor_norm <- readRDS("/space/scratch/amorin/R_objects/GSE216019_scor_all_norm.RDS")

# Load cxg dat, and extract meta, raw counts, and normalized 

dat <- readRDS(dat_path)


meta <- dat@meta.data %>% 
  dplyr::rename(Cell_type = cell_type) %>% 
  rownames_to_column(var = "ID")


mat_norm <- GetAssayData(dat, slot = "data") %>% 
  as.matrix() %>% 
  ensembl_to_symbol(pc) %>% 
  get_pcoding_only(pc)
  

mat_counts <- GetAssayData(dat, slot = "counts") %>% 
  as.matrix() %>% 
  ensembl_to_symbol(pc) %>% 
  get_pcoding_only(pc)


# demo norm vs counts
gene <- "ACTB"
hist(mat_counts[gene, ], breaks = 100)
hist(mat_norm[gene, ], breaks = 100)
plot(mat_counts[gene, ], mat_norm[gene, ])
cor(mat_counts[gene, ], mat_norm[gene, ], method = "spearman")

hist(mat_counts[, which.max(meta$nCount_RNA)], breaks = 100)
hist(mat_norm[, which.max(meta$nCount_RNA)], breaks = 100)
hist(mat_counts[, which.min(meta$nCount_RNA)], breaks = 100)
hist(mat_norm[, which.min(meta$nCount_RNA)], breaks = 100)
#


# Comparing all cor

genes <- colnames(all_pcor_counts)
genes <- genes[genes != ""]


diag(all_pcor_counts) <- NA
diag(all_pcor_norm) <- NA


all_pcor_l <- lapply(genes, function(x) {
  cor(all_pcor_counts[, x], all_pcor_norm[, x], method = "spearman", use = "pairwise.complete.obs")
})
names(all_pcor_l) <- genes


summary(unlist(all_pcor_l))
plot(density(unlist(all_pcor_l), na.rm = TRUE))
which.min(all_pcor_l)


tfs <- c("ASCL1", "HES1", "MECP2", "MEF2C", "NEUROD1", "PAX6", "RUNX1", "TCF4")
unlist(all_pcor_l)[tfs]


gene <- "RPL41"
all_pcor_l[[gene]]

gene_counts_df <- data.frame(Counts = mat_counts[gene, ], Norm = mat_norm[gene, ])
plot(gene_counts_df$Counts, gene_counts_df$Norm)

gene_cor_df <- data.frame(Symbol = rownames(all_pcor_counts),
                          Counts = all_pcor_counts[, gene],
                          Norm = all_pcor_norm[, gene])

gene_cor_df$Diff <- abs(gene_cor_df$Counts - gene_cor_df$Norm)

plot(gene_cor_df$Counts, gene_cor_df$Norm)



plot(mat_counts["MECP2", ], mat_counts["RPL41", ])
cor(mat_counts["MECP2", ], mat_counts["RPL41", ])

plot(mat_norm["MECP2", ], mat_norm["RPL41", ])
cor(mat_norm["MECP2", ], mat_norm["RPL41", ])



# Comparing RSR


rsr_counts <- lowertri_to_symm(rsr_counts)
rsr_norm <- lowertri_to_symm(rsr_norm)


genes <- colnames(rsr_counts)
genes <- genes[genes != ""]


all_rsr_l <- lapply(genes, function(x) {
  cor(rsr_counts[, x], rsr_norm[, x], method = "spearman", use = "pairwise.complete.obs")
})



gene <- "ASCL1"
plot(rsr_counts[, gene], rsr_norm[, gene])
cor(rsr_counts[, gene], rsr_norm[, gene], method = "spearman", use = "pairwise.complete.obs")


summary(unlist(all_rsr_l))



# Average exprs across all cells

avg_df <- data.frame(Symbol = rownames(mat_norm),
                     Avg_count = rowMeans(mat_counts),
                     Avg_norm = rowMeans(mat_norm))

avg_df$Avg_count_rank <- rank(avg_df$Avg_count) / nrow(avg_df)
avg_df$Avg_norm_rank <- rank(avg_df$Avg_norm) / nrow(avg_df)


cor(avg_df$Avg_count, avg_df$Avg_norm, method = "spearman")
plot(avg_df$Avg_count, avg_df$Avg_norm)
plot(avg_df$Avg_count_rank, avg_df$Avg_norm_rank)


avg_df <- data.frame(Count_norm_cor = unlist(all_pcor_l)) %>% 
  rownames_to_column(var = "Symbol") %>% 
  left_join(avg_df, by = "Symbol")


plot(avg_df$Count_norm_cor, avg_df$Avg_count)
plot(avg_df$Count_norm_cor, avg_df$Avg_norm)

plot(avg_df$Count_norm_cor, avg_df$Avg_count)



# Showing relationship of gene cor using raw counts or normalized. Want to show
# over a range of low, average, and highly expressed genes

set.seed(45)

avg_lo <- avg_df %>% filter(Count_rank > 0 & Count_rank < 0.2) %>% slice_sample(n = 100)
avg_med <- avg_df %>% filter(Count_rank > 0.4 & Count_rank < 0.6) %>% slice_sample(n = 100)
avg_hi <- avg_df %>% filter(Count_rank > 0.8) %>% slice_sample(n = 30)



get_cor_df_list <- function(mat_counts, mat_norm, subset_genes, ncores = 8) {
  
  cor_l <- lapply(subset_genes, function(z) {
    
    cor_count <- WGCNA::cor(x = mat_counts[, z, drop = FALSE], 
                            y = mat_counts,
                            nThreads = ncores)
    
    cor_norm <- WGCNA::cor(x = mat_norm[, z, drop = FALSE], 
                           y = mat_norm,
                           nThreads = ncores)
    
    df <- data.frame(Symbol = colnames(mat_counts), 
                     Cor_count = cor_count[1, ], 
                     Cor_norm = cor_norm[1, ],
                     Diff = abs(cor_count[1, ] - cor_norm[1, ]))
    
    df <- filter(df, Symbol != z)
  
  })
  names(cor_l) <- subset_genes
  
  return(cor_l)
}



cor_df_list <- get_cor_df_list(mat_counts, 
                               mat_norm, 
                               c(avg_lo$Symbol, avg_med$Symbol, avg_hi$Symbol))


cor_of_cor <- unlist(lapply(cor_df_list, function(x) {
  cor(x$Cor_count, x$Cor_norm, use = "pairwise.complete.obs")
}))


summary_cor_df <- data.frame(
  Symbol = c(avg_lo$Symbol, avg_med$Symbol, avg_hi$Symbol),
  Cor = cor_of_cor,
  Group = c(
    rep("Low_exprs", nrow(avg_lo)),
    rep("Med_exprs", nrow(avg_med)),
    rep("High_exprs", nrow(avg_hi))
  )
)

summary_cor_df$Group <- factor(summary_cor_df$Group, 
                               levels = unique(summary_cor_df$Group))


ggplot(summary_cor_df, aes(x = Group, y = Cor)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25))



# Example of low agreement between gene correlation vectors generated from
# raw count matrix versus norm count


min_cor <- summary_cor_df %>% 
  group_by(Group) %>% 
  filter(Cor == min(Cor))


p_list <- lapply(min_cor$Symbol, function(x) {
  
  ggplot(cor_df_list[[x]], aes(x = Cor_norm, Cor_count)) +
  geom_point(shape = 21) + 
  ggtitle(x) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 20))
  
})
names(p_list) <- min_cor$Symbol


# TODO: it looks like a lack of agreement may be caused by genes that have a
# small range of expression for the integer counts, which get streteched out 
# for rankit norm
# TODO: do all cor and cor_of_cor ~ standard_rank + max_count + nonzero_cells
plot(mat_counts[, "AP5Z1"], mat_norm[, "AP5Z1"])
plot(mat_counts[, "RUNX1"], mat_norm[, "RUNX1"])

gene1 <- "MECP2"
head(sort(rsr_norm[gene1, ]), 10)
gene2 <- "ZBTB20"
all_celltype_cor(mat = t(mat_norm), meta = meta, gene1 = gene1, gene2 = gene2)
ct <- "endothelial cell"
plot(mat_norm[filter(meta, Cell_type == ct)$ID, gene1], mat_norm[filter(meta, Cell_type == ct)$ID, gene2])
plot(mat_counts[filter(meta, Cell_type == ct)$ID, gene1], mat_counts[filter(meta, Cell_type == ct)$ID, gene2])
