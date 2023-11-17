## Examining aggregate coexpression between mouse and human
## Most analysis would focus on the 1:1 DIOPT orthologs.
## But it would be interesting to look for gene family members presence

## For each TR, get the ability of one species to recover the other species's TR.
## AUPRC, ROC. Have to decide on a cut-off. Could show plot of recovery at different k
## Null would be... sample of matched genes?
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
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

# Genomic evidence from Morin 2023 
evidence_l <- readRDS(evidence_path)

# Measurement matrices used for filtering when a gene was never expressed
# msr_hg <- readRDS(msr_mat_hg_path)
# msr_mm <- readRDS(msr_mat_mm_path)

# Loading the TF aggregate matrices
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)
# agg_ribo_hg <- load_or_generate_agg(path = agg_ribo_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = ribo_genes$Symbol_hg)
# agg_ribo_mm <- load_or_generate_agg(path = agg_ribo_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = ribo_genes$Symbol_mm)


# Only keep ortho genes measured in both species

pc_ortho <- filter(pc_ortho,
                   Symbol_hg %in% pc_hg$Symbol &
                   Symbol_mm %in% pc_mm$Symbol)

tfs_ortho <- filter(pc_ortho, 
                    Symbol_hg %in% names(rank_tf_hg) & 
                      Symbol_mm %in% names(rank_tf_mm))




# 
# ------------------------------------------------------------------------------




gene_hg <- "ASCL1"
gene_mm <- "Ascl1"


df_hg <- left_join(rank_tf_hg[[gene_hg]], pc_ortho,
                   by = c("Symbol" = "Symbol_hg")) %>%
  filter(!is.na(ID))


df_mm <- left_join(rank_tf_mm[[gene_mm]], pc_ortho, 
                   by = c("Symbol" = "Symbol_mm")) %>% 
  filter(!is.na(ID))


ortho_df <- left_join(df_hg, df_mm, 
                      by = "ID", 
                      suffix = c("_Human", "_Mouse")) %>% 
  filter(!is.na(Avg_RSR_Human) & !is.na(Avg_RSR_Mouse))


plot(ortho_df$Avg_RSR_Human, ortho_df$Avg_RSR_Mouse)
cor(ortho_df$Avg_RSR_Human, ortho_df$Avg_RSR_Mouse, method = "spearman")
filter(ortho_df, Rank_RSR_Human <= 100 & Rank_RSR_Mouse <= 100)



ortho_cor <- mclapply(1:nrow(tfs_ortho), function(x) {
  
  df_hg <- left_join(rank_tf_hg[[tfs_ortho$Symbol_hg[x]]], pc_ortho, by = c("Symbol" = "Symbol_hg")) %>% filter(!is.na(ID))
  df_mm <- left_join(rank_tf_mm[[tfs_ortho$Symbol_mm[x]]], pc_ortho, by = c("Symbol" = "Symbol_mm")) %>% filter(!is.na(ID))
  ortho_df <- left_join(df_hg, df_mm, by = "ID", suffix = c("_Human", "_Mouse")) %>% filter(!is.na(Avg_RSR_Human) & !is.na(Avg_RSR_Mouse))
  cor(ortho_df$Avg_RSR_Human, ortho_df$Avg_RSR_Mouse, method = "spearman")
  
}, mc.cores = ncore)


ortho_cor <- data.frame(Cor = unlist(ortho_cor), Symbol = tfs_ortho$Symbol_hg)
summary(ortho_cor$Cor)


px <- plot_hist(ortho_cor, 
                stat_col = "Cor", 
                xlab = "Spearman's correlation", 
                title = "n=1241 orthologous TFs") + 
  xlim(c(-0.4, 1))



ggsave(px, height = 5, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "ortho_rank_similarity.png"))



# Topk overlap of targets and relative to other TFs



# For each ortho TF, gets its top k overlap between species
# NOTE: this is taking the top k AFTER filtering for ortho, so not necessarily
# within the top k ranks of the unfiltered (non-ortho) data
# TODO: no reason to do so many joins, just use once and reference main
# TODO: just get list of top 100 for each species and perform overlap after.
# TODO: also want to comment on generic overlaps


topk_relative <- function(gene_query, 
                          tfs_ortho, 
                          pc_ortho, 
                          rank_tf_hg, 
                          rank_tf_mm, 
                          k = 100) {
  
  gene_ortho <- filter(pc_ortho, 
                       Symbol_hg == gene_query | Symbol_mm == gene_query)
  
  
  topk_query_hg <- left_join(rank_tf_hg[[gene_ortho$Symbol_hg]], pc_ortho,
                             by = c("Symbol" = "Symbol_hg")) %>%
    filter(!is.na(ID)) %>% 
    slice_min(Rank_RSR, n = k) 
  
  
  topk_query_mm <- left_join(rank_tf_mm[[gene_ortho$Symbol_mm]], pc_ortho,
                             by = c("Symbol" = "Symbol_mm")) %>%
    filter(!is.na(ID)) %>% 
    slice_min(Rank_RSR, n = k) 
  
  
  # length(intersect(topk_query_hg$ID, topk_query_mm$ID))
  
  
  topk_hg_in_mm <- lapply(tfs_ortho$Symbol_mm, function(x) {
    
    topk <- filter(rank_tf_mm[[x]], Symbol %in% pc_ortho$Symbol_mm) %>%
      slice_min(Rank_RSR, n = k) %>% 
      pull(Symbol)
    
    length(intersect(topk_query_hg$Symbol_mm, topk))
    
  })
  
  names(topk_hg_in_mm) <- tfs_ortho$Symbol_mm
  topk_hg_in_mm <- unlist(topk_hg_in_mm)
  
  
  
  topk_mm_in_hg <- lapply(tfs_ortho$Symbol_hg, function(x) {
    
    topk <- filter(rank_tf_hg[[x]], Symbol %in% pc_ortho$Symbol_hg) %>%
      slice_min(Rank_RSR, n = k) %>% 
      pull(Symbol)
    
    length(intersect(topk_query_mm$Symbol_hg, topk))
    
  })
  
  names(topk_mm_in_hg) <- tfs_ortho$Symbol_hg
  topk_mm_in_hg <- unlist(topk_mm_in_hg)
  
  
  topk_df <- data.frame(
    Symbol = tfs_ortho$Symbol_hg, 
    Topk_hg_in_mm = topk_hg_in_mm, 
    Topk_mm_in_hg = topk_mm_in_hg
  )
  
  return(topk_df)
  
}



topk_all <- mclapply(tfs_ortho$Symbol_hg, function(x) {
  
  
  topk_relative(gene_query = x,
                tfs_ortho,
                pc_ortho,
                rank_tf_hg,
                rank_tf_mm,
                k = 100)
  
}, mc.cores = 8)
names(topk_all) <- tfs_ortho$Symbol_hg
# saveRDS(topk_all, "/space/scratch/amorin/R_objects/ortho_tf_topk=100.RDS")
topk_all <- readRDS("/space/scratch/amorin/R_objects/ortho_tf_topk=100.RDS")


x <- "ASCL1"



# For each species, calculate the percentile of the matched ortho TF relative to
# all others


ortho_perc <- lapply(tfs_ortho$Symbol_hg, function(x) {
  
  df <- topk_all[[x]]
  tf_in <- filter(df, Symbol == x)
  tf_out <- filter(df, Symbol != x)
  
  data.frame(
    Symbol = x,
    Count_mm_in_hg = tf_in$Topk_mm_in_hg,
    Perc_mm_in_hg = ecdf(tf_out$Topk_mm_in_hg)(tf_in$Topk_mm_in_hg),
    Count_hg_in_mm = tf_in$Topk_hg_in_mm,
    Perc_hg_in_mm = ecdf(tf_out$Topk_hg_in_mm)(tf_in$Topk_hg_in_mm)
  )
})


ortho_perc_df <- do.call(rbind, ortho_perc)


summary(Filter(is.numeric, ortho_perc_df))
sum(ortho_perc_df$Perc_mm_in_hg == 1)
sum(ortho_perc_df$Perc_hg_in_mm == 1)



qplot(ortho_perc_df, xvar = "Perc_mm_in_hg", yvar = "Perc_hg_in_mm") +
  xlab("Percentile mouse in human") +
  ylab("Percentile human in mouse")


gene <- "ASCL1"

plot_df <- mutate(topk_all[[gene]], Label = Symbol == gene)


ggplot(plot_df, aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg)) +
  geom_jitter(shape = 21, size = 2.4, width = 0.5, height = 0.5) +
  geom_text_repel(
    data = filter(plot_df, Label),
    aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg, label = Symbol, fontface = "italic"),
    size = 5,
    segment.size = 0.1,
    segment.color = "grey50") +
  xlab("Human in mouse") +
  ylab("Mouse in human") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



plot_grid(
  
  plot_hist(plot_df, stat_col = "Topk_hg_in_mm") + 
    geom_vline(xintercept = filter(plot_df, Symbol == gene)$Topk_hg_in_mm,
               linewidth = 1.6,
               col = "royalblue") +
    xlab("Human in mouse"),
  
  plot_hist(plot_df, stat_col = "Topk_mm_in_hg") + 
    geom_vline(xintercept = filter(plot_df, Symbol == gene)$Topk_mm_in_hg, 
               linewidth = 1.6,
               col = "goldenrod") +
    xlab("Mouse in human"),
  
  nrow = 1)





### TODO: OLD below, to be removed/updated



# Functions
# ------------------------------------------------------------------------------



get_ortho_mat <- function(genes, ids, meta, pc_df, ncores = 1) {
  
  stopifnot(all(ids %in% meta$ID))
  meta_hg <- filter(meta, ID %in% ids & Species == "Human")
  meta_mm <- filter(meta, ID %in% ids & Species == "Mouse")
  ortho_genes <- filter(pc_df, Symbol_hg %in% genes | Symbol_mm %in% genes)
  
  # load data and only keep requested gene 
  
  agg_hg <- mclapply(1:nrow(meta_hg), function(x) {
    load_agg_mat_list(ids = meta_hg$ID[x], paths = meta_hg$Path[x])[[1]][pc_df$Symbol_hg, ortho_genes$Symbol_hg]
  }, mc.cores = ncores)
  
  
  agg_mm <- mclapply(1:nrow(meta_mm), function(x) {
    load_agg_mat_list(ids = meta_mm$ID[x], paths = meta_mm$Path[x])[[1]][pc_df$Symbol_mm, ortho_genes$Symbol_mm]
  }, mc.cores = ncores)
  
  
  ortho_l <- mclapply(1:nrow(ortho_genes), function(i) {
    
    hg_l <- lapply(agg_hg, function(x) x[, ortho_genes$Symbol_hg[i]])
    mm_l <- lapply(agg_mm, function(x) x[, ortho_genes$Symbol_mm[i]])
    
    ortho_mat <- do.call(cbind, c(hg_l, mm_l))
    rownames(ortho_mat) <- pc_df$ID
    colnames(ortho_mat) <- paste0(c(meta_hg$ID, meta_mm$ID), "_", ortho_genes$Symbol_hg[i])
    
    return(ortho_mat)
  }, mc.cores = ncores)
  
  names(ortho_l) <- ortho_genes$ID
  gc(verbose = FALSE)
  
  return(ortho_l)
}



# Binarize matrix such that top k and bottom k is 1, mid is 0

make_discrete_k_mat <- function(mat, k = 1000) {
  
  bin_mat <- apply(mat, 2, function(x) {
    sort_values <- sort(x, decreasing = TRUE)
    topk <- sort_values[k]
    btmk <- sort_values[length(x) - k + 1]
    ifelse(x >= topk | x <= btmk, 1, 0)
  })
  
  return(bin_mat)
}



# https://stackoverflow.com/a/66594545 Jaccard faster than nested loop 

binary_jaccard <- function(x, y) {
  sum(x & y) / sum(x | y)  
}


get_jaccard_matrix <- function(mat) {
  tmp <- asplit(mat, 2)
  jacc_mat <- outer(tmp, tmp, Vectorize(binary_jaccard))
  return(jacc_mat)
}



# Performs PCA with prcomp and returns list of the resulting 
# object as well as the variance explained

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



# Generate ortho mat for given gene
# ------------------------------------------------------------------------------


tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")

genes <- pc_df %>% 
  filter(Symbol_hg %in% c(tfs, ribo_genes$Symbol_hg) | 
         Symbol_mm %in% c(tfs, ribo_genes$Symbol_mm)) %>% 
  pull(Symbol_hg)
  
ids <- c(meta_hg$ID, meta_mm$ID)



# outfile <- "/space/scratch/amorin/R_objects/16-06-2023_ortho_aggcoexp.RDS"
outfile <- "/space/scratch/amorin/R_objects/19-06-2023_ortho_aggcoexp.RDS"



if (!file.exists(outfile)) {
  
  message("Starting to generate list of ortho matrices: ", Sys.time())
  
  ortho_l <- get_ortho_mat(genes = genes, 
                           ids = ids, 
                           meta = sc_meta, 
                           pc_df = pc_df,
                           ncores = ncore)
  
  saveRDS(ortho_l, outfile)
  
  message("Done: ", Sys.time())
  
} else {
  
  ortho_l <- readRDS(outfile)

}




# Extract a single matrix of interest

# sub_genes <- filter(pc_df, Symbol_hg %in% c("PAX6", "RPL3", "RPL13"))$ID
# sub_genes <- filter(pc_df, Symbol_mm %in% tfs)$ID
# sub_genes <- c(filter(pc_df, Symbol_mm %in% tfs)$ID, filter(pc_df, Symbol_hg %in% c("RPL3", "RPL13"))$ID)
# sub_genes <- filter(pc_df, Symbol_hg %in% c("RPL3", "RPL13"))$ID
sub_genes <- names(ortho_l)
# sub_genes <- filter(pc_df, Symbol_hg %in% c("RUNX1", "HES1"))$ID


ortho_mat <- do.call(cbind, ortho_l[sub_genes])
ortho_mat_bin <- make_discrete_k_mat(ortho_mat, k = 1000)




# PCA 
# ------------------------------------------------------------------------------

# ortho_pca <- pca_and_var(ortho_mat, scale_arg = TRUE)
ortho_pca <- pca_and_var(ortho_mat, scale_arg = FALSE)
# ortho_pca <- pca_and_var(ortho_mat_bin, scale_arg = FALSE)



# top genes
npcs <- 10
# sort(abs(ortho_pca$PC$rotation[, 1]), decreasing = TRUE)[1:50]
# sort(abs(ortho_pca$PC$rotation[, 2]), decreasing = TRUE)[1:50]


df <- data.frame(
  ortho_pca$PC$x[, 1:npcs],
  Matrix_ID = colnames(ortho_mat),
  ID = rep(ids, length(sub_genes)),
  Symbol = rep(sub_genes, each = length(ids))
) %>%
  left_join(sc_meta[, c("ID", "Species")], by = "ID")



p1a <- pc_scatter(df = df, pc_list = ortho_pca, pc_x = 1, pc_y = 2) + theme(legend.position = "none")
p1b <- pc_scatter(df = df, pc_list = ortho_pca, pc_x = 2, pc_y = 3) + theme(legend.position = "none")
p1c <- pc_scatter(df = df, pc_list = ortho_pca, pc_x = 3, pc_y = 4) + theme(legend.position = "none")
p1d <- pc_scatter(df = df, pc_list = ortho_pca, pc_x = 4, pc_y = 5) 
p1_leg <- get_legend(p1d)
p1d <- p1d + theme(legend.position = "none")
p1 <- plot_grid(p1a, p1b, p1c, p1d, nrow = 2)
p1 <- plot_grid(p1, p1_leg, rel_widths = c(2, 0.3))



# Get average per gene
ortho_avg <- do.call(
  cbind,
  lapply(sub_genes, function(x) rowMeans(ortho_mat[, filter(df, Symbol == x)$Matrix_ID]))
)
colnames(ortho_avg) <- sub_genes


lapply(1:ncol(ortho_avg), function(x) head(sort(ortho_avg[, x], decreasing = TRUE), 10))

# ortho_mat_rank <- colrank_mat(ortho_mat)
# head(sort(rowMeans(ortho_mat), decreasing = TRUE), 20)
# head(sort(rowMeans(ortho_mat_rank), decreasing = FALSE), 20)




 
# Similarity by TF status. Use Jaccard of top/bottomk and cor

jacc_mat <- get_jaccard_matrix(ortho_mat_bin)
cor_mat <- WGCNA::cor(ortho_mat, method = "pearson")


# Organize into df of unique pairs with group status


# hacky name fix for _ delim in IDs
split_ids <- str_split(colnames(jacc_mat), pattern = "_")

fixed_names <- lapply(split_ids, function(x) {
  if (length(x) > 2) {
    x <- paste0(x[1], x[2], "_", x[3])
  } else {
    x <- paste0(x[1], "_", x[2])
  }
  return(x)
})


colnames(jacc_mat) <- rownames(jacc_mat) <- unlist(fixed_names) 


format_pair_df <- function(mat, meta, symmetric = TRUE) {
  
  df <- mat_to_df(mat, symmetric = symmetric)
  
  row_split <- str_split(df$Row, "_", simplify = TRUE)
  col_split <- str_split(df$Col, "_", simplify = TRUE)
  
  id1 <- row_split[, 1]
  id2 <- col_split[, 1]
  
  tf1 <- str_to_upper(row_split[, 2])
  tf2 <- str_to_upper(col_split[, 2])

  species1 <- meta$Species[match(id1, meta$ID)]
  species2 <- meta$Species[match(id2, meta$ID)]

  
  group <- vapply(1:nrow(row_split), function(i) {
    if (tf1[i] == tf2[i] & species1[i] == species2[i]) {
      return("In_gene_in_species")
    } else if (tf1[i] == tf2[i] & species1[i] != species2[i]) {
      return("In_gene_out_species")
    } else {
      return("Out")
    }
  }, FUN.VALUE = character(1))
  

  group <- factor(group,
                  levels = c("In_gene_in_species", "In_gene_out_species", "Out"))
  
  pair_df <- data.frame(df, 
                        TF1 = tf1, 
                        TF2 = tf2, 
                        Group = group)
  
  return(pair_df)
}



pair_df <- format_pair_df(jacc_mat, meta = sc_meta, symmetric = TRUE) %>% 
  dplyr::rename(Jaccard = Value) %>% 
  mutate(Cor = mat_to_df(cor_mat, symmetric = TRUE)$Value)



stat <- "Cor"
boxplot(pair_df[[stat]] ~ pair_df$Group)
plot(density(filter(pair_df, Group == "In_gene_in_species")[[stat]]))
lines(density(filter(pair_df, Group == "In_gene_out_species")[[stat]]), col = "red")
lines(density(filter(pair_df, Group == "Out")[[stat]]), col = "blue")





# Why are there numerous examples of perfect similarities between datasets...
# Because they are not expressed 

top <- filter(pair_df, Jaccard == 1)

id1 <- "GSE115469"
id2 <- "GSE145928"
gene1 <- "PAX6"


plot(ortho_mat[, paste0(id1, "_", gene1)], ortho_mat[, paste0(id2, "_", gene1)])
cor(ortho_mat[, paste0(id1, "_", gene1)], ortho_mat[, paste0(id2, "_", gene1)], method = "spearman")

sub_meta <- filter(sc_meta, ID %in% c(id1, id2))

agg_l <- load_agg_mat_list(ids = sub_meta$ID, paths = sub_meta$Path)
dat_l <- load_dat_list(ids = sub_meta$ID)

sub_mat <- do.call(cbind, lapply(agg_l, function(x) x[, gene1]))
plot(sub_mat[, 1], sub_mat[, 2])
cor(sub_mat, method = "spearman")

gene2 <- names(head(sort(sub_mat[, 2], decreasing = TRUE)))[2]

all_celltype_cor(mat = dat_l[[id1]]$Mat, meta = dat_l[[id1]]$Meta, gene1 = gene1, gene2 = gene2)
all_celltype_cor(mat = dat_l[[id1]]$Mat, meta = dat_l[[id1]]$Meta, gene1 = "RPL3", gene2 = "RPL11")


# Train a classifier to check model performance at group discrimination




