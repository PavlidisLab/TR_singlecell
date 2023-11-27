## TODO
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
# agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
# agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)
# agg_ribo_hg <- load_or_generate_agg(path = agg_ribo_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = ribo_genes$Symbol_hg)
# agg_ribo_mm <- load_or_generate_agg(path = agg_ribo_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = ribo_genes$Symbol_mm)



# Create a gene x TF matrix of summarized ranks for general inspection
# ------------------------------------------------------------------------------


gen_rank_matrix <- function(rank_l) {
  
  rank_mat <- as.matrix(do.call(cbind, lapply(rank_l, `[`, "Rank_RSR")))
  colnames(rank_mat) <- names(rank_l)
  return(rank_mat)
}


rank_mat_hg <- gen_rank_matrix(rank_tf_hg)
rank_mat_mm <- gen_rank_matrix(rank_tf_mm)



# Most commonly highly ranked genes
# TODO: are these genes high because of technical aspect?

rank_order_hg <- sort(rowMeans(rank_mat_hg, na.rm = TRUE))
rank_order_mm <- sort(rowMeans(rank_mat_mm, na.rm = TRUE))


# head(rank_order_hg)
# tail(rank_order_hg)
# head(rank_order_mm)
# tail(rank_order_mm)

head(sort(rank_mat_hg["RPL18A",]), 50)
tail(sort(rank_mat_hg["RPL18A",]), 100)
hist(rank_mat_hg["RPL18A",])



# microglia
head(sort(rank_mat_hg["TRIM28",]), 50)
head(sort(rank_mat_hg["CSF1R",]), 50)
head(sort(rank_mat_hg["SPI1",]), 50)

head(sort(rank_mat_mm["Trim28",]), 50)
head(sort(rank_mat_mm["Csf1r",]), 50)
head(sort(rank_mat_mm["Spi1",]), 50)



# Ribosomal rankings
# ------------------------------------------------------------------------------


count_top_ribo <- function(ribo_l, topn = 82) {
  
  ribo_genes <- names(ribo_l)
  
  sum_l <- lapply(ribo_genes, function(x) {
    sum(arrange(ribo_l[[x]], Rank_RSR)$Symbol[1:topn] %in% ribo_genes)
  })
  
  return(data.frame(Symbol = ribo_genes, Count = unlist(sum_l)))
}


count_ribo_hg <- count_top_ribo(rank_ribo_hg)
count_ribo_mm <- count_top_ribo(rank_ribo_mm)

summary(count_ribo_hg$Count)
summary(count_ribo_mm$Count)



# Top coexpr partner for each TF
# ------------------------------------------------------------------------------


top_hg <- lapply(rank_tf_hg, slice_min, Rank_RSR) %>% 
  do.call(rbind, .,) %>% 
  rownames_to_column(var = "TF") %>% 
  arrange(desc(Topk_count))


top_mm <- lapply(rank_tf_mm, slice_min, Rank_RSR) %>% 
  do.call(rbind, .,) %>% 
  rownames_to_column(var = "TF") %>% 
  arrange(desc(Topk_count))


# Largest delta between top ranked by Avg_RSR and second highest... evidence
# of neighbouring genes?

calc_top_delta <- function(rank_l) {
  
  tfs <- names(rank_l)
  
  delta_l <- lapply(tfs, function(x) {
    top <- slice_max(rank_l[[x]], Avg_RSR, n = 2) %>% arrange(desc(Avg_RSR))
    data.frame(TF = x, Symbol = top$Symbol[1], Delta = abs(diff(top$Avg_RSR)))
  })
  
  delta_df <- data.frame(do.call(rbind, delta_l)) %>% arrange(desc(Delta))
  
  return(delta_df)
}


delta_hg <- calc_top_delta(rank_tf_hg)
delta_mm <- calc_top_delta(rank_tf_mm)



source("~/TRagg/R/utils/range_table_functions.R")
pc_gr_hg <- pc_to_gr(filter(pc_hg, !is.na(Start)))
pc_gr_mm <- pc_to_gr(filter(pc_mm, !is.na(Start)))


# Removing NA start removes some TFs with rankings
# TODO: check why ZNF883 removed
delta_hg <- filter(delta_hg, TF %in% pc_gr_hg$Symbol)
delta_mm <- filter(delta_hg, TF %in% pc_gr_mm$Symbol)


# Distance
# TODO: NA necessary?

# dist_hg <- lapply(1:nrow(delta_hg), function(x) {
#   
#   tf_gr <- pc_gr_hg[pc_gr_hg$Symbol %in% c(delta_hg$TF[x], delta_hg$Symbol[x])]
#   chroms <- as.character(seqnames(tf_gr))
#   
#   if (!identical(chroms[1], chroms[2])) {
#     return(NA)
#   } else {
#     abs(distance(tf_gr[1], tf_gr[2], ignore.strand = TRUE))
#   }
# })


# delta_hg$Distance <- unlist(dist_hg)


# TODO: nearest gene

dist_hg <- lapply(1:nrow(delta_hg), function(x) {
  
  tf_gr <- pc_gr_hg[pc_gr_hg$Symbol == delta_hg$TF[x]]
  out_gr <- pc_gr_hg[pc_gr_hg$Symbol != delta_hg$TF[x]]
  
  out_gr$Distance <- distance(tf_gr, out_gr, ignore.strand = TRUE)
  out_gr <- out_gr[order(out_gr$Distance)]
  
  top_dist <- out_gr[which(out_gr$Symbol == delta_hg$Symbol[x])]$Distance
  top_rank <- which(out_gr$Symbol == delta_hg$Symbol[x])
  
  data.frame(Distance = top_dist, Rank = top_rank)
  
})


delta_hg2 <- cbind(delta_hg, do.call(rbind, dist_hg))




# Inspecting targets from the most similar TFs
arrange(rank_tf_hg$E2F8, desc(Avg_RSR)) %>% head(30)
arrange(rank_tf_mm$E2f8, desc(Avg_RSR)) %>% head(30)



# Correlate the final rankings between TF to find similar/dissimilar
# ------------------------------------------------------------------------------


rank_similarity <- function(rank_mat) {
  
  rank_cor <- mat_to_df(colwise_cor(rank_mat), symmetric = TRUE)
  
  rank_topk <-
    mat_to_df(colwise_topk_intersect(rank_mat, k = 1000), symmetric = TRUE)
  
  rank_bottomk <-
    mat_to_df(colwise_topk_intersect(rank_mat, k = 1000, decreasing = FALSE),
              symmetric = TRUE)
  
  rank_df <- data.frame(Row = rank_cor$Row, 
                        Col = rank_cor$Col,
                        Cor = rank_cor$Value,
                        Topk = rank_topk$Value,
                        Bottomk = rank_bottomk$Value)
  
  return(rank_df)
}



rank_sim_hg <- rank_similarity(rank_mat_hg)
rank_sim_mm <- rank_similarity(rank_mat_mm)





# TODO: Rank_RSR versus stat var

slice_topk_genes <- function(df, topk) {
  stopifnot(c("Symbol", "Rank_RSR") %in% colnames(df))
  slice_min(df, Rank_RSR, n = topk)[["Symbol"]]
}


# get_stat <- function() {
#   lapply(rank_l, function(y) filter(y, Symbol == x)$Rank_RSR)
# }



get_stat_by_genes <- function(rank_l, genes) {
  
  # stopifnot(stat %in% colnames(rank_l[[1]]), genes %in% rownames(rank_l[[1]]))
  stopifnot(genes %in% rownames(rank_l[[1]]))
  
  lapply(genes, function(gene) {
    unlist(
      lapply(rank_l, function(rank_df) filter(rank_df, Symbol == gene)[["Rank_RSR"]])
      )
    })
}





get_topk_rank_mat <- function(rank_l, topk = 5) {
  
  topk_each_l <- lapply(rank_l, slice_topk_genes, topk = topk)
  
  all_topk_genes <- unlist(topk_each_l, use.names = FALSE)
  
  if (any(duplicated(all_topk_genes))) message("Duplicate topk genes between TFs")
  
  topk_all_l <- get_stat_by_genes(rank_l, genes = all_topk_genes)
  topk_all_mat <- as.matrix(do.call(rbind, topk_all_l))
  rownames(topk_all_mat) <- all_topk_genes
  
  return(topk_all_mat)
}



# TODO: brush up stringr str extract+NA or str detect logical

extract_pc_pattern <- function(pattern, pc_hg, pc_mm, pc_ortho) {
  
  hg <- str_detect(str_to_lower(pc_hg$Symbol), pattern)
  hg <- sort(pc_hg$Symbol[hg])
  
  mm <- str_detect(str_to_lower(pc_mm$Symbol), pattern)
  mm <- sort(pc_mm$Symbol[mm])
  
  ortho <- filter(pc_ortho, Symbol_hg %in% hg | Symbol_mm %in% mm)
  
  return(list(Human = hg, Mouse = mm, Ortho = ortho))
  
}


pax <- extract_pc_pattern("^pax[:digit:]+$", pc_hg, pc_mm, pc_ortho)
pax_rank_hg <- get_topk_rank_mat(rank_tf_hg[pax$Human], topk = 8)
pax_rank_mm <- get_topk_rank_mat(rank_tf_mm[pax$Mouse], topk = 8)

ascl <- extract_pc_pattern("^ascl[:digit:]+$", pc_hg, pc_mm, pc_ortho)
ascl_rank_hg <- get_topk_rank_mat(rank_tf_hg[ascl$Human], topk = 8)
ascl_rank_mm <- get_topk_rank_mat(rank_tf_mm[ascl$Mouse], topk = 8)

sample <- sample(names(rank_tf_hg), 9)  # testing sample
sample_rank <- get_topk_rank_mat(rank_tf_hg[sample])


# TODO: rank range is not informative. Consider breaks like 1-500, 501-1000,
# 2000-3000...

rank_heatmap <- function(rank_mat) {
  
  pheatmap(rank_mat, 
           color = rev(c('#fcfbfd','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d')),
           cluster_rows = FALSE, 
           cluster_cols = FALSE)
}



# TODO: better way of binarizing to mat? logical->integer dropped str

rank_binary_heatmap <- function(rank_mat, rank_cutoff) {
  
  rank_mat <- rank_mat < rank_cutoff
  rank_mat <- ifelse(rank_mat, 1, 0)
  
  pheatmap(rank_mat, 
           color = c("white", "black"),
           cluster_rows = FALSE, 
           cluster_cols = FALSE)
  
}



rank_binary_heatmap(ascl_rank_mm, rank_cutoff = 100)
rank_heatmap(ascl_rank_mm)


# Plots 
# ------------------------------------------------------------------------------


# Barchart of count of ribo top-ranked partners that are also ribo

ggplot(count_ribo_hg, aes(x = reorder(Symbol, Count), y = Count)) +
  geom_bar(stat = "identity", colour = "black", fill = "slategrey") +
  ylab("Count of top coexpression partners that are also ribosomal") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))



plot_hist(count_ribo_hg, stat_col = "Count")



# Distn of top k for a given TF

ggplot(rank_tf_hg$ASCL1, aes(x = Topk_count)) +
  geom_histogram(bins = 20) +
  geom_vline(xintercept = filter(rank_tf_hg$ASCL1, Symbol == "HES6")$Topk_count) +
  geom_vline(xintercept = filter(rank_tf_hg$ASCL1, Symbol == "DLL1")$Topk_count) +
  geom_vline(xintercept = filter(rank_tf_hg$ASCL1, Symbol == "DLL3")$Topk_count) +
  geom_vline(xintercept = filter(rank_tf_hg$ASCL1, Symbol == "DLL4")$Topk_count) +
  ggtitle("Human ASCL1") +
  xlab("Top k=1000 count") +
  ylab("Gene count") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))



# Inspecting the topk count/proportion of the top coexpr partner for each TF
plot_hist(top_hg, stat_col = "Topk_count")
plot_hist(top_mm, stat_col = "Topk_count")

plot_hist(top_hg, stat_col = "Topk_proportion")
plot_hist(top_mm, stat_col = "Topk_proportion")



# RSR heatmap
# ------------------------------------------------------------------------------


# TODO: remove when finalized

# agg_tf_cpm <- agg_tf_hg
# rank_cpm <- rank_tf_hg
# msr_cpm <- msr_hg

# agg_tf_ln <- readRDS("/space/scratch/amorin/R_objects/agg_mat_TF_list_hg_lognorm.RDS")
# rank_ln <- readRDS("/space/scratch/amorin/R_objects/ranking_agg_TF_hg_lognorm.RDS")
# msr_ln <- readRDS("~/Robj/binary_measurement_matrix_hg_lognorm.RDS")

# msr_mat <- msr_cpm
# agg_l <- agg_tf_cpm
# rank_l <- rank_cpm

# msr_mat <- msr_ln
# agg_l <- agg_tf_ln
# rank_l <- rank_ln


gene <- "ASCL1"

msr_mat <- msr_hg
agg_l <- agg_tf_hg  # agg_ribo_hg
rank_l <- rank_tf_hg


rank_df <- filter(rank_l[[gene]], Symbol != gene) %>% arrange(Rank_RSR)

#
gene_mat <- subset_to_measured(gene_vec_to_mat(agg_l, gene), msr_mat, gene)
gene_mat <- gene_mat[setdiff(rownames(gene_mat), gene), ]


# Logical matrix of whether TF-genes were co-measured in an experiment
comsr <- lapply(rownames(gene_mat), function(y) {
  msr_mat[gene, colnames(gene_mat)] & msr_mat[y, colnames(gene_mat)]
})
comsr <- do.call(rbind, comsr)
rownames(comsr) <- rownames(gene_mat)


# Replace non-measured ties with NAs
gene_mat_na <- gene_mat
gene_mat_na[comsr == FALSE] <- NA

# Replace non-measured ties with median of the entire matrix
med <- median(gene_mat, na.rm = TRUE)
gene_mat_med <- gene_mat
gene_mat_med[comsr == FALSE] <- med

# Replace non-measured ties with 0.5
gene_mat_05 <- gene_mat
gene_mat_05[comsr == FALSE] <- 0.5

# Replace non-measured ties with row means of measured observations
# NOTE: scrap - all NA genes will still be NA and so need a second imputation

gene_mat_rowmean <- gene_mat_na
row_means <- rowMeans(gene_mat_na, na.rm = TRUE)

for (i in 1:nrow(gene_mat_rowmean)) {
  ix <- is.na(gene_mat_rowmean[i, ])
  gene_mat_rowmean[i, ix] <- row_means[i]
}



# Organizing rank dataframe for different ways of handling NAs and aggregating.
# NAs: keep as tied RSR values, impute to median of matrix, or keep as NAs to be
# ignored in following aggregation)
# aggregation: average RSR, median RSR, average rank of RSR, median rank of RSR,
# and rank product

# TODO: Avg/Medrank_wo_imp/w_NA selects for all NAs
# TODO: rank product with NAs select for all non-NAs



rsr_avg_wo_imp <- rowMeans(gene_mat)
rsr_avg_w_med <- rowMeans(gene_mat_med)
rsr_avg_w_05 <- rowMeans(gene_mat_05)
rsr_avg_w_na <- rowMeans(gene_mat_na, na.rm = TRUE)

rsr_med_wo_imp <- apply(gene_mat, 1, median)
rsr_med_w_med <- apply(gene_mat_med, 1, median)
rsr_med_w_05 <- apply(gene_mat_05, 1, median)
rsr_med_w_na <- apply(gene_mat_na, 1, function(x) median(na.omit(x)))

rsr_avgrank_wo_imp <- rowMeans(colrank_mat(gene_mat))
rsr_avgrank_w_med <- rowMeans(colrank_mat(gene_mat_med))
rsr_avgrank_w_05 <- rowMeans(colrank_mat(gene_mat_05))
rsr_avgrank_w_na <- rowMeans(colrank_mat(gene_mat_na), na.rm = TRUE)

rsr_medrank_wo_imp <- apply(colrank_mat(gene_mat), 1, median)
rsr_medrank_w_med <- apply(colrank_mat(gene_mat_med), 1, median)
rsr_medrank_w_05 <- apply(colrank_mat(gene_mat_05), 1, median)
rsr_medrank_w_na <- apply(colrank_mat(gene_mat_na), 1, function(x) median(na.omit(x)))

# NOTE: Often run into integer overflow issue with rank product, use sum of logs
# rsr_rankprod_wo_imp <- rank(matrixStats::rowProds(colrank_mat(gene_mat)))
# rsr_rankprod_w_med <- rank(matrixStats::rowProds(colrank_mat(gene_mat_med)))
# rsr_rankprod_w_05 <- rank(matrixStats::rowProds(colrank_mat(gene_mat_05)))
# rsr_rankprod_w_na <- rank(matrixStats::rowProds(colrank_mat(gene_mat_na)))

rsr_rankprod_wo_imp <- rank(rowSums(log(colrank_mat(gene_mat))))
rsr_rankprod_w_med <- rank(rowSums(log(colrank_mat(gene_mat_med))))
rsr_rankprod_w_05 <- rank(rowSums(log(colrank_mat(gene_mat_05))))
rsr_rankprod_w_na <- rank(rowSums(log(colrank_mat(gene_mat_na))))


df <- data.frame(
  
  Symbol = names(rsr_avg_wo_imp),
  
  Avg_wo_imp = rsr_avg_wo_imp,
  Avg_w_med = rsr_avg_w_med,
  Avg_w_05 = rsr_avg_w_05,
  Avg_w_na = rsr_avg_w_na,
  
  Med_wo_imp = rsr_med_wo_imp,
  Med_w_med = rsr_med_w_med,
  Med_w_05 = rsr_med_w_05,
  Med_w_na = rsr_med_w_na,
  
  Avgrank_wo_imp = rsr_avgrank_wo_imp,
  Avgrank_w_med = rsr_avgrank_w_med,
  Avgrank_w_05 = rsr_avgrank_w_05,
  Avgrank_w_na = rsr_avgrank_w_na,
  
  Medrank_wo_imp = rsr_medrank_wo_imp,
  Medrank_w_med = rsr_medrank_w_med,
  Medrank_w_05 = rsr_medrank_w_05,
  Medrank_w_na = rsr_medrank_w_na,
  
  Rankprod_wo_imp = rsr_rankprod_wo_imp,
  Rankprod_w_med = rsr_rankprod_w_med,
  Rankprod_w_05 = rsr_rankprod_w_05,
  Rankprod_w_na = rsr_rankprod_w_na
  
)


df <- left_join(df, rank_df[, c("Symbol", "N_comeasured")], by = "Symbol")

cor_rank <- cor(select_if(df, is.numeric), method = "spearman", use = "pairwise.complete.obs")


# Heatmap colours and breaks

pal <- c('#ffffff','#f0f0f0','#d9d9d9','#bdbdbd','#969696','#737373','#525252','#252525','#000000')
na_col <- "red" 
# breaks <- seq(min(rank_df$Avg_RSR), max(rank_df$Avg_RSR), length.out = length(pal))
breaks <- seq(0, 1, length.out = length(pal))



# Look at top ranked genes

n <- 100
rank_order <- arrange(df, desc(Avg_w_med))$Symbol

# plot_mat <- gene_mat
# plot_mat <- gene_mat_med
plot_mat <- gene_mat_na


top_mat <- plot_mat[rank_order[1:n], ]
top_na_order <- sort(-apply(top_mat, 2, function(x) sum(is.na(x))))
top_mat <- top_mat[, names(top_na_order)]


pheatmap(top_mat, 
         breaks = breaks,
         color = pal,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         # show_rownames = FALSE,
         show_rownames = TRUE,
         fontsize_row = 10,
         border_color = NA,
         # cellheight = 10,
         na_col = na_col)



# Lower ranks as a comparison

pheatmap(plot_mat[rank_order[11000:(11000 + n)], ], 
         breaks = breaks,
         color = pal,
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         # show_rownames = FALSE,
         show_rownames = TRUE,
         fontsize_row = 10,
         border_color = NA,
         na_col = na_col)



# Looking at all genes/experiments with NAs coloured

na_order <- sort(-apply(gene_mat_na, 2, function(x) sum(is.na(x))))


# No gene order
pheatmap(gene_mat_na[, names(na_order)],
         color = pal,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         na_col = na_col)


# Order by average RSR with median imputation
pheatmap(gene_mat_na[arrange(df, desc(Avg_w_med))$Symbol, names(na_order)],
         color = pal,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         na_col = na_col)


# Order by average RSR without median imputation
pheatmap(gene_mat_na[arrange(df, desc(Avg_wo_imp))$Symbol, na_order],
         color = pal,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         na_col = na_col)


# Relationship between measurement and ranking


pqa <- qplot(df, yvar = "Avg_w_med", xvar = "N_comeasured")


pqb <- df %>%
  mutate(Bin = cut_width(N_comeasured, width = 10, boundary = 0)) %>%
  ggplot(aes(x = Bin, y = Avg_w_med)) +
  geom_boxplot(width = 0.2, fill = "#69b3a2") +
  xlab("N comeasured") +
  ylab("Average RSR (+imputation)") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))



pqc <- qplot(df, yvar = "Avg_wo_imp", xvar = "N_comeasured")


pqd <- df %>%
  mutate(Bin = cut_width(N_comeasured, width = 10, boundary = 0)) %>%
  ggplot(aes(x = Bin, y = Avg_wo_imp)) +
  geom_boxplot(width = 0.2, fill = "#69b3a2") +
  xlab("N comeasured") +
  ylab("Average RSR (-imputation)") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))


plot_grid(
  plot_grid(pqa, pqc, ncol = 1),
  plot_grid(pqb, pqd, ncol = 1),
  rel_widths = c(0.5, 1)
)
