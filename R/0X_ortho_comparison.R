## Examining aggregate coexpression profiles between mouse and human
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggrepel)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

k <- 200

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


# Only keep ortho genes measured in both species

pc_ortho <- filter(pc_ortho,
                   Symbol_hg %in% pc_hg$Symbol &
                   Symbol_mm %in% pc_mm$Symbol)

tfs_ortho <- filter(pc_ortho, 
                    Symbol_hg %in% names(rank_tf_hg) & 
                    Symbol_mm %in% names(rank_tf_mm))


# Functions
# ------------------------------------------------------------------------------


# Gene is assumed to be an ortho gene found in pc_df. Retrieve and join the 
# corresponding gene ranks from the mouse and human rank lists.

join_ortho_ranks <- function(gene,
                             pc_df = pc_ortho,
                             rank_hg = rank_tf_hg,
                             rank_mm = rank_tf_mm) {
  
  gene_ortho <- filter(pc_df, Symbol_hg == gene | Symbol_mm == gene)
  
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


# Flip sign of average RSR (rank metric) to prioritize bottom of ranks

reverse_ranks <- function(df) {
  mutate(df, Avg_RSR_hg = -Avg_RSR_hg, Avg_RSR_mm = -Avg_RSR_mm)
}


# Generate a gene x gene matrix whose elements represent the count of top K
# elements shared between the ortho ranks.
# Rows are human in mouse, cols are mouse in human
# NOTE: This is taking the top k AFTER filtering for ortho, so not necessarily
# within the top k ranks of the unfiltered (non-ortho) data

calc_overlap_mat <- function(ortho_l, k, reverse = FALSE, ncores = 1) {
  
  # Reverse rankings to prioritize bottom of ranks (neg coexpression)?
  if (reverse) ortho_l <- lapply(ortho_l, reverse_ranks)
  
  # First get list of the top K elements for human and mouse for each gene rank
  genes <- names(ortho_l)
  
  # Get list of topK genes (char vec) for each gene and species
  topk_l <- mclapply(ortho_l, function(x) {
    list(
      Human = slice_max(x, Avg_RSR_hg, n = k)$Symbol_hg,
      Mouse = slice_max(x, Avg_RSR_mm, n = k)$Symbol_hg
    )
  }, mc.cores = ncores)
  names(topk_l) <- genes
  
  # Generate overlap matrix
  overlap_mat <- matrix(0, nrow = length(topk_l), ncol = length(topk_l))
  rownames(overlap_mat) <- colnames(overlap_mat) <- genes
  
  for (i in 1:nrow(overlap_mat)) {
    for (j in 1:ncol(overlap_mat)) {
      overlap_mat[i, j] <- topk_intersect(topk_l[[i]]$Human, topk_l[[j]]$Mouse)
    }
  }
  
  return(overlap_mat)
}


# Given gene x gene overlap_mat, return a data.frame containing the count of
# overlap for the orthologous match, as well as the percentile of the match
# relative to all other TFs for each species

gen_percentile_df <- function(overlap_mat, ncores = 1) {
  
  genes <- intersect(rownames(overlap_mat), colnames(overlap_mat))
  
  perc <- mclapply(genes, function(x) {
    
    count_overlap <- overlap_mat[x, x]
    hg_in_mm <- overlap_mat[x, setdiff(tfs_ortho$Symbol_hg, x)]
    mm_in_hg <- overlap_mat[setdiff(tfs_ortho$Symbol_hg, x), x]
    
    data.frame(
      Symbol = x,
      Count = count_overlap,
      Perc_hg_in_mm = ecdf(hg_in_mm)(count_overlap),
      Perc_mm_in_hg = ecdf(mm_in_hg)(count_overlap))
    
  }, mc.cores = ncores)
  
  perc_df <- do.call(rbind, perc)
  perc_df$Perc_ortho <- rowMeans(perc_df[, c("Perc_hg_in_mm", "Perc_mm_in_hg")])
  
  return(perc_df)
}


# Join the human and mouse rankings for only orthologous genes, then generate
# the similarity objects. The human symbols used for representing ortho TFs. 
# -------------------------------------------------------------------------------


# List of joined human and mouse TF rankings
rank_tf_ortho <- mclapply(tfs_ortho$Symbol_hg, join_ortho_ranks, mc.cores = ncore)
names(rank_tf_ortho) <- tfs_ortho$Symbol_hg


# Get the Spearman's correlation of rankings between every ortho TF
cor_df <- calc_ortho_cor(rank_tf_ortho, ncores = ncore)


# Top and bottom K matrices 
topk_mat <- calc_overlap_mat(rank_tf_ortho, k = k, ncores = ncore)
bottomk_mat <- calc_overlap_mat(rank_tf_ortho, k = k, reverse = TRUE, ncores = ncore)


# Percentiles dfs of top/bottom k
topk_df <- gen_percentile_df(topk_mat, ncores = ncore)
bottomk_df <- gen_percentile_df(bottomk_mat, ncores = ncore)


# Prepare and join all similarity dfs
topk_df <- rename_with(topk_df, ~paste0("Topk_", .), -c("Symbol"))
bottomk_df <- rename_with(bottomk_df, ~paste0("Bottomk_", .), -c("Symbol"))
sim_df <- plyr::join_all(list(cor_df, topk_df, bottomk_df), by = "Symbol")



# Exploring the Spearman's correlations
# ------------------------------------------------------------------------------


cor_summ <- summary(cor_df$Cor)

cor_low <- filter(cor_df, Cor < -0.2) %>% arrange(Cor)
cor_high <- filter(cor_df, Cor > 0.8) %>% arrange(-Cor)


check_tf <- "ASCL1"


rank_common <- rank_tf_ortho[[check_tf]] %>% 
  mutate(Diff_rank = Rank_RSR_hg - Rank_RSR_mm,
         Add = Avg_RSR_hg + Avg_RSR_mm,
         RP = rank(Rank_RSR_hg * Rank_RSR_mm)) %>% 
  relocate(Symbol_hg, RP, Rank_RSR_hg, Rank_RSR_mm, Add, Diff_rank, Avg_RSR_hg, Avg_RSR_mm)



# Check data coverage of extreme TFs, find they are generally well-measured

ext_counts_hg <- rowSums(msr_hg[c(cor_low$Symbol, cor_high$Symbol), ])
ext_mm <- filter(pc_ortho, Symbol_hg %in% c(cor_low$Symbol, cor_high$Symbol))$Symbol_mm
ext_counts_mm <- rowSums(msr_mm[ext_mm, ])





# Topk overlap of matched ortho TFs
# ------------------------------------------------------------------------------




tt <- calc_overlap_mat(ortho_l = rank_tf_ortho, 
                       k = 200, 
                       ncores = ncore)



tt2 <- calc_overlap_mat(ortho_l = rank_tf_ortho, 
                        reverse = TRUE,
                        k = 200, 
                        ncores = ncore)



tt3 <- gen_percentile_df(tt1, ncores = ncore)
tt4 <- gen_percentile_df(tt2, ncores = ncore)





# Mean across rows (human) and cols (mouse) to get most common/generic 

common_overlap <- data.frame(
  Symbol = rownames(overlap_mat),
  Hg_in_mm = rowMeans(overlap_mat),
  Mm_in_hg = colMeans(overlap_mat)
)



# Inspect the human in mouse with the highest overlap
slice_max(common_overlap, Hg_in_mm)
overlap_mat["MEF2C", "MEF2C"]
head(sort(-overlap_mat["MEF2C", ]), 20)
sum(msr_hg["MEF2C", ])
intersect(topk_l$MEF2C$Human, topk_l$FLI1$Mouse)

# Inspect the mouse in human with the highest overlap
slice_max(common_overlap, Mm_in_hg)
overlap_mat["RXRA", "RXRA"]
head(sort(-overlap_mat[, "RXRA"]), 20)
sum(msr_mm["Rxra", ])
intersect(topk_l$MEF2C$Human, topk_l$RXRA$Mouse)

# Inspect the human in mouse with lowest overlap
slice_min(common_overlap, Hg_in_mm)
overlap_mat["ZSCAN5B", "ZSCAN5B"]
head(sort(-overlap_mat["ZSCAN5B", ]), 20)
sum(msr_hg["ZSCAN5B", ])
intersect(topk_l$ZSCAN5B$Human, topk_l$ASCL2$Mouse)

# Inspect mouse in human with lowest overlap
slice_min(common_overlap, Mm_in_hg)
overlap_mat["DMRTC2", "DMRTC2"]
head(sort(-overlap_mat[, "DMRTC2"]), 20)
sum(msr_mm["Dmrtc2", ])
intersect(topk_l$KLF17$Human, topk_l$DMRTC2$Mouse)


# Relative positioning of CENPA
which(arrange(common_overlap, desc(Hg_in_mm))$Symbol == "CENPA")





summary(Filter(is.numeric, ortho_perc_df))
summary(Filter(is.numeric, filter(ortho_perc_df, Topk_count > 0)))

sum(ortho_perc_df$Topk_count == 0)
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



summary(Filter(is.numeric, bottomk_perc_df))


# Combining similarities

a1 <- rename_with(ortho_perc_df, ~paste0("Topk_", .), -c("Symbol", "Topk_count"))
a2 <- rename_with(bottomk_perc_df, ~paste0("Bottomk_", .), -c("Symbol", "Bottomk_count"))
ortho_df <- left_join(ortho_cor, a1, by = "Symbol") %>% left_join(a2, by = "Symbol")


cor(select_if(ortho_df, is.numeric), method = "spearman")

qplot(ortho_df, xvar = "Bottomk_Perc_mm_in_hg", yvar = "Bottomk_Perc_hg_in_mm")

qplot(ortho_df, xvar = "Bottomk_Perc_ortho", yvar = "Topk_Perc_ortho")

qplot(ortho_df, xvar = "Topk_count", yvar = "Bottomk_count")

qplot(ortho_df, xvar = "Bottomk_count", yvar = "Bottomk_Perc_ortho")


qplot(ortho_df, xvar = "Bottomk_count", yvar = "Cor")
qplot(ortho_df, xvar = "Topk_count", yvar = "Cor")



filter(ortho_df, Topk_Perc_ortho == 1 & Bottomk_Perc_ortho == 1)
filter(ortho_df, Topk_Perc_ortho > 0.99 & Bottomk_Perc_ortho > 0.99)






# Plots
# ------------------------------------------------------------------------------


# Histogram of ortho TF spearman cors

px <- plot_hist(cor_df, 
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

plot_gene <- "ASCL1"

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


# Plots of common


plot_grid(
  plot_hist(common_overlap, stat_col = "Hg_in_mm") + xlab("Human in mouse"),
  plot_hist(common_overlap, stat_col = "Hg_in_mm") + xlab("Mouse in human"),
  nrow = 2)


qplot(common_overlap, xvar = "Hg_in_mm", yvar = "Mm_in_hg")




### TODO: OLD below, to be removed/updated



# Functions
# ------------------------------------------------------------------------------






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