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
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Loading the TF aggregate matrices
# agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
# agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)
# agg_ribo_hg <- load_or_generate_agg(path = agg_ribo_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = ribo_genes$Symbol_hg)
# agg_ribo_mm <- load_or_generate_agg(path = agg_ribo_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = ribo_genes$Symbol_mm)


# Only keep ortho genes measured in both species

pc_ortho <- filter(pc_ortho,
                   Symbol_hg %in% pc_hg$Symbol &
                   Symbol_mm %in% pc_mm$Symbol)

tfs_ortho <- filter(pc_ortho, 
                    Symbol_hg %in% names(rank_tf_hg) & 
                      Symbol_mm %in% names(rank_tf_mm))




# Functions
# ------------------------------------------------------------------------------


# Given a gene character, return a df of the gene rank information for human and
# mouse for only ortho genes

subset_ortho_data <- function(gene, 
                              pc_df = pc_ortho, 
                              rank_hg = rank_tf_hg,
                              rank_mm = rank_tf_mm) {
  
  gene_ortho <- filter(pc_ortho, Symbol_hg == gene | Symbol_mm == gene)
  
  if (nrow(gene_ortho) != 1) stop("A one to one match was not found")
  
  df_hg <- left_join(rank_hg[[gene_ortho$Symbol_hg]], 
                     pc_ortho[, c("Symbol_hg", "ID")],
                     by = c("Symbol" = "Symbol_hg")) %>% 
    filter(!is.na(ID))
  
  
  df_mm <- left_join(rank_mm[[gene_ortho$Symbol_mm]], 
                     pc_ortho[, c("Symbol_mm", "ID")], 
                     by = c("Symbol" = "Symbol_mm")) %>% 
    filter(!is.na(ID))
  
  ortho_df <- left_join(df_hg, df_mm, 
                        by = "ID", 
                        suffix = c("_hg", "_mm")) %>% 
    filter(!is.na(Avg_RSR_hg) & !is.na(Avg_RSR_mm))
  
  return(ortho_df)
}



# Return a df of the spearman cor between human and mouse RSR rankings

calc_ortho_cor <- function(ortho_rank_l, ncores = 1) {
  
  symbols <- names(ortho_rank_l)
  
  cor_l <- mclapply(ortho_rank_l, function(x) {
    cor(x$Avg_RSR_hg, x$Avg_RSR_mm, method = "spearman")
  }, mc.cores = ncores)
  
  cor_df <- data.frame(Cor = unlist(cor_l), Symbol = symbols)
  
  return(cor_df)
}



# Create a unified ranking for each ortho TF and explore mouse/human correlation
# -------------------------------------------------------------------------------


# List of joined human and mouse TF rankings

rank_tf_ortho <- mclapply(tfs_ortho$Symbol_hg, subset_ortho_data, mc.cores = ncore)
names(rank_tf_ortho) <- tfs_ortho$Symbol_hg


# Get the Spearman correlation of rankings between every ortho TF

ortho_cor <- calc_ortho_cor(rank_tf_ortho, ncores = ncore)


# Exploring correlations

cor_summ <- summary(ortho_cor$Cor)

cor_low <- filter(ortho_cor, Cor < -0.2) %>% arrange(Cor)
cor_high <- filter(ortho_cor, Cor > 0.8) %>% arrange(-Cor)


check_tf <- "HNF4A"


rank_common <- rank_tf_ortho[[check_tf]] %>% 
  mutate(Diff = Avg_RSR_hg - Avg_RSR_mm,
         Add = Avg_RSR_hg + Avg_RSR_mm,
         RP = rank(Rank_RSR_hg * Rank_RSR_mm)) %>% 
  relocate(Symbol_hg, RP, Rank_RSR_hg, Rank_RSR_mm, Add, Diff)



# Check data coverage of extreme TFs, find they are generally well-measured

ext_counts_hg <- rowSums(msr_hg[c(cor_low$Symbol, cor_high$Symbol), ])
ext_mm <- filter(pc_ortho, Symbol_hg %in% c(cor_low$Symbol, cor_high$Symbol))$Symbol_mm
ext_counts_mm <- rowSums(msr_mm[ext_mm, ])





# Topk overlap of matched ortho TFs
# ------------------------------------------------------------------------------



# For each ortho TF, gets its top k overlap between species
# NOTE: this is taking the top k AFTER filtering for ortho, so not necessarily
# within the top k ranks of the unfiltered (non-ortho) data
# TODO: no reason to do so many joins, just use once and reference main
# TODO: just get list of top 100 for each species and perform overlap after.
# TODO: also want to comment on generic overlaps



k <- 100


topk_l <- mclapply(rank_tf_ortho, function(x) {
  list(
    Human = slice_max(x, Avg_RSR_hg, n = k)$Symbol_hg,
    Mouse = slice_max(x, Avg_RSR_mm, n = k)$Symbol_hg
  )
}, mc.cores = ncore)
names(topk_l) <- tfs_ortho$Symbol_hg



# Rows are human in mouse, cols are mouse in human

overlap_mat <- matrix(0, nrow = length(topk_l), ncol = length(topk_l))
rownames(overlap_mat) <- colnames(overlap_mat) <- tfs_ortho$Symbol_hg

for (i in 1:nrow(overlap_mat)) {
  for (j in 1:ncol(overlap_mat)) {
    overlap_mat[i, j] <- length(intersect(topk_l[[i]]$Human, topk_l[[j]]$Mouse))
  }
}



# Sums across rows (human) and cols (mouse) to get most common/generic 

overlap_row <- rowSums(overlap_mat)
overlap_col <- colSums(overlap_mat)


# Inspect the human in mouse with the highest overlap
head(sort(-overlap_row))
overlap_mat["MEF2C", "MEF2C"]
head(sort(-overlap_mat["MEF2C", ]), 20)
sum(msr_hg["MEF2C", ])
intersect(topk_l$MEF2C$Human, topk_l$FLI1$Mouse)

# Inspect the mouse in human with the highest overlap
head(sort(-overlap_col))
overlap_mat["RXRA", "RXRA"]
head(sort(-overlap_mat[, "RXRA"]), 20)
sum(msr_mm["Rxra", ])
intersect(topk_l$MEF2C$Human, topk_l$RXRA$Mouse)

# Inspect the human in mouse with lowest overlap
head(sort(overlap_row))
overlap_mat["ZSCAN5B", "ZSCAN5B"]
head(sort(-overlap_mat["ZSCAN5B", ]), 20)
sum(msr_hg["ZSCAN5B", ])
intersect(topk_l$ZSCAN5B$Human, topk_l$ASCL2$Mouse)

# Inspect mouse in human with lowest overlap
head(sort(overlap_col))
overlap_mat["DMRTC2", "DMRTC2"]
head(sort(-overlap_mat[, "DMRTC2"]), 20)
sum(msr_mm["Dmrtc2", ])
intersect(topk_l$KLF17$Human, topk_l$DMRTC2$Mouse)


# TODO: relative positioning of CENPA



# Df of the counts and relative percentiles of the ortho top k overlap


ortho_perc <- mclapply(tfs_ortho$Symbol_hg, function(x) {
  
  count_overlap <- overlap_mat[x, x]
  hg_in_mm <- overlap_mat[setdiff(tfs_ortho$Symbol_hg, x), x]
  mm_in_hg <- overlap_mat[x, setdiff(tfs_ortho$Symbol_hg, x)]
  
  data.frame(
    Symbol = x,
    Topk_count = count_overlap,
    Perc_hg_in_mm = ecdf(hg_in_mm)(count_overlap),
    Perc_mm_in_hg = ecdf(mm_in_hg)(count_overlap))
  
}, mc.cores = ncore)


ortho_perc_df <- do.call(rbind, ortho_perc)
ortho_perc_df$Perc_ortho <- rowMeans(ortho_perc_df[, c("Perc_hg_in_mm", "Perc_mm_in_hg")])


summary(Filter(is.numeric, ortho_perc_df))
sum(ortho_perc_df$Perc_mm_in_hg == 1)
sum(ortho_perc_df$Perc_hg_in_mm == 1)
sum(ortho_perc_df$Perc_ortho == 1)


# NEUROG3 example of no overlap but elevated percentiles

inspect_df <- data.frame(
  Symbol = rownames(overlap_mat),
  Hg_in_mm = overlap_mat["NEUROG3", ],
  Mm_in_hg = overlap_mat[, "NEUROG3"]
)



# Inspect bottom k
# ------------------------------------------------------------------------------


bottomk_l <- mclapply(rank_tf_ortho, function(x) {
  list(
    Human = slice_min(x, Avg_RSR_hg, n = k)$Symbol_hg,
    Mouse = slice_min(x, Avg_RSR_mm, n = k)$Symbol_hg
  )
}, mc.cores = ncore)
names(bottomk_l) <- tfs_ortho$Symbol_hg



# Rows are human in mouse, cols are mouse in human

bottomk_mat <- matrix(0, nrow = length(bottomk_l), ncol = length(bottomk_l))
rownames(bottomk_mat) <- colnames(bottomk_mat) <- tfs_ortho$Symbol_hg

for (i in 1:nrow(bottomk_mat)) {
  for (j in 1:ncol(bottomk_mat)) {
    bottomk_mat[i, j] <- length(intersect(bottomk_l[[i]]$Human, bottomk_l[[j]]$Mouse))
  }
}


bottomk_perc <- mclapply(tfs_ortho$Symbol_hg, function(x) {
  
  count_overlap <- bottomk_mat[x, x]
  hg_in_mm <- bottomk_mat[setdiff(tfs_ortho$Symbol_hg, x), x]
  mm_in_hg <- bottomk_mat[x, setdiff(tfs_ortho$Symbol_hg, x)]
  
  data.frame(
    Symbol = x,
    Bottomk_count = count_overlap,
    Perc_hg_in_mm = ecdf(hg_in_mm)(count_overlap),
    Perc_mm_in_hg = ecdf(mm_in_hg)(count_overlap))
  
}, mc.cores = ncore)


bottomk_perc_df <- do.call(rbind, bottomk_perc)
bottomk_perc_df$Perc_ortho <- rowMeans(bottomk_perc_df[, c("Perc_hg_in_mm", "Perc_mm_in_hg")])


###



# Plots
# ------------------------------------------------------------------------------


# Histogram of ortho TF spearman cors

px <- plot_hist(ortho_cor, 
                stat_col = "Cor", 
                xlab = "Spearman's correlation", 
                title = "n=1241 orthologous TFs") + 
  xlim(c(-0.4, 1))



ggsave(px, height = 5, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "ortho_rank_similarity.png"))


# Scatter of mouse versus human RSR for a given TF

qplot(rank_tf_ortho[[check_tf]], xvar = "Avg_RSR_hg", yvar = "Avg_RSR_mm", title = check_tf)



# Scatterplot of percentile in mouse/human

qplot(ortho_perc_df, xvar = "Perc_mm_in_hg", yvar = "Perc_hg_in_mm") +
  xlab("Percentile mouse in human") +
  ylab("Percentile human in mouse")


# Hist of actual Top K counts

plot_hist(ortho_perc_df, nbins = 50, stat_col = "Topk_count")


# Scatter of perc ortho versus top k

qplot(ortho_perc_df, xvar = "Topk_count", yvar = "Perc_ortho") +
  xlab("Top K count") +
  ylab("Percentile ortho")



# Scatter plot of topk counts between species

plot_gene <- "NEUROG3"

plot_dfx <- data.frame(
  Symbol = rownames(overlap_mat),
  Topk_hg_in_mm = overlap_mat[plot_gene, ],
  Topk_mm_in_hg = overlap_mat[, plot_gene],
  Label = rownames(overlap_mat) == plot_gene)


ggplot(plot_dfx, aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg)) +
  geom_jitter(shape = 21, size = 2.4, width = 0.5, height = 0.5) +
  geom_text_repel(
    data = filter(plot_dfx, Label),
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



# Hist of topk counts for each species


plot_grid(
  
  plot_hist(plot_dfx, stat_col = "Topk_hg_in_mm") + 
    geom_vline(xintercept = filter(plot_dfx, Symbol == plot_gene)$Topk_hg_in_mm,
               # width = 1.6,
               col = "royalblue") +
    xlab("Human in mouse"),
  
  plot_hist(plot_dfx, stat_col = "Topk_mm_in_hg") + 
    geom_vline(xintercept = filter(plot_dfx, Symbol == plot_gene)$Topk_mm_in_hg, 
               # width = 1.6,
               col = "goldenrod") +
    xlab("Mouse in human"),
  
  nrow = 2)




### TODO: OLD below, to be removed/updated



# Functions
# ------------------------------------------------------------------------------


# TODO: replace - get a matrix of aggregate vectors for one gene both species

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



# TODO: check overlap with similarity functions

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




### TODO: for loop vs mclapply vs outer


# TODO: the existance of the nested list object is somehow greatly slowing
# down R session... speeds back up when removed

# a1 <- mclapply(1:length(topk_l), function(i) {
#   lapply(1:length(topk_l), function(j) {
#     topk_intersect(topk_l[[i]]$Human, topk_l[[j]]$Mouse)
#   })
# }, mc.cores = ncore)


a1 <- lapply(1:length(topk_l), function(i) {
  lapply(1:length(topk_l), function(j) {
    topk_intersect(topk_l[[i]]$Human, topk_l[[j]]$Mouse)
  })
})


a2 <- as.numeric(unlist(a1))
a3 <- matrix(a2, byrow = TRUE, nrow = length(topk_l), ncol = length(topk_l))
rownames(a3) <- colnames(a3) <- names(topk_l)
identical(overlap_mat, a3)



# https://stackoverflow.com/a/66594545 Jaccard faster than nested loop 

binary_jaccard <- function(x, y) {
  sum(x & y) / sum(x | y)  
}


get_jaccard_matrix <- function(mat) {
  tmp <- asplit(mat, 2)
  jacc_mat <- outer(tmp, tmp, Vectorize(binary_jaccard))
  return(jacc_mat)
}


###