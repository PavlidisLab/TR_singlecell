## Examining aggregate coexpression between mouse and human
## Most analysis would focus on the 1:1 DIOPT orthologs.
## But it would be interesting to look for gene family members presence

## For each TR, get the ability of one species to recover the other species's TR.
## AUPRC, ROC. Have to decide on a cut-off. Could show plot of recovery at different k
## Null would be... sample of matched genes?
## -----------------------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(randomForest)
library(caret)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

# Protein coding tables
pc_ortho <- read.delim(pc_ortho_path, stringsAsFactors = FALSE)
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# Ribo genes for comparison
sribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Sribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
lribo_hg <- read.table("/space/grp/amorin/Metadata/HGNC_human_Lribosomal_genes.csv", stringsAsFactors = FALSE, skip = 1, sep = ",", header = TRUE)
ribo_genes <- filter(pc_ortho, Symbol_hg %in% c(sribo_hg$Approved.symbol, lribo_hg$Approved.symbol))


# Metadata and IDs
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# Hacky merge until I get the IDs better sorted. "_" delim in ID names problematic
sc_meta <- mutate(sc_meta, ID = str_replace(ID, "_", ""))

meta_hg <- filter(sc_meta, Species == "Human")
meta_mm <- filter(sc_meta, Species == "Mouse")



# Loading aggregate coexpression matrices
# agg_hg <- load_agg_mat_list(ids_hg[1:2])
# agg_mm <- load_agg_mat_list(ids_mm[1:2])
# agg_hg <- load_agg_mat_list(ids = meta_hg$ID, paths = meta_hg$Path)
# agg_mm <- load_agg_mat_list(ids = meta_mm$ID, paths = meta_mm$Path)


# Only keep ortho genes
pc_df <- filter(pc_ortho, 
                 Symbol_hg %in% pc_hg$Symbol, 
                 Symbol_mm %in% pc_mm$Symbol)


# agg_hg <- lapply(agg_hg, function(x) x[pc_df$Symbol_hg, pc_df$Symbol_hg])
# agg_mm <- lapply(agg_mm, function(x) x[pc_df$Symbol_mm, pc_df$Symbol_mm])




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




