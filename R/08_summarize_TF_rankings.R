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


head(rank_order_hg)
tail(rank_order_hg)

head(rank_order_mm)
tail(rank_order_mm)

head(sort(rank_mat_hg["RPL27",]), 50)
tail(sort(rank_mat_hg["RPL27",]), 50)


head(sort(rank_mat_mm["Hes6",]), 50)
tail(sort(rank_mat_mm["Pax6",]), 50)


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




# Inspect ortho
# TODO: helpers to get ortho
# TODO: move to own script
# ------------------------------------------------------------------------------


tfs_ortho <- filter(pc_ortho, 
                    Symbol_hg %in% names(rank_tf_hg) & 
                    Symbol_mm %in% names(rank_tf_mm))


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
    Mm_in_hg = ecdf(tf_out$Topk_mm_in_hg)(tf_in$Topk_mm_in_hg),
    Count_hg_in_mm = tf_in$Topk_hg_in_mm,
    Hg_in_mm = ecdf(tf_out$Topk_hg_in_mm)(tf_in$Topk_hg_in_mm)
  )
})


ortho_perc_df <- data.frame(
  Symbol = tfs_ortho$Symbol_hg,
  do.call(rbind, ortho_perc)
)



qplot(ortho_perc_df, xvar = "Mm_in_hg", yvar = "Hg_in_mm")


gene <- "PAX6"

plot_df <- mutate(topk_all[[gene]], Label = Symbol == gene)


ggplot(plot_df, aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg)) +
  geom_jitter(shape = 21, size = 2.4, width = 0.5, height = 0.5) +
  geom_text_repel(
    data = filter(plot_df, Label),
    aes(x = Topk_hg_in_mm, y = Topk_mm_in_hg, label = Symbol, fontface = "italic"),
    size = 5,
    segment.size = 0.1,
    segment.color = "grey50") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))



plot_grid(
  
  plot_hist(plot_df, stat_col = "Topk_hg_in_mm") + 
    geom_vline(xintercept = filter(plot_df, Symbol == gene)$Topk_hg_in_mm, col = "royalblue"),
          
  plot_hist(plot_df, stat_col = "Topk_mm_in_hg") + 
    geom_vline(xintercept = filter(plot_df, Symbol == gene)$Topk_mm_in_hg, col = "goldenrod"),
  
  nrow = 1)






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
# agg_tf_ln <- readRDS("/space/scratch/amorin/R_objects/agg_mat_TF_list_hg_lognorm.RDS")
# rank_cpm <- rank_tf_hg
# rank_ln <- readRDS("/space/scratch/amorin/R_objects/ranking_agg_TF_hg_lognorm.RDS")
# msr_cpm <- msr_hg
# msr_ln <- readRDS("~/Robj/binary_measurement_matrix_hg_lognorm.RDS")


gene <- "ASCL1"

# msr_mat <- msr_hg
# msr_mat <- msr_cpm
msr_mat <- msr_ln

# agg_l <- agg_ribo_hg
# agg_l <- agg_tf_hg
# agg_l <- agg_tf_cpm
agg_l <- agg_tf_ln

# rank_l <- rank_tf_hg
# rank_l <- rank_ribo_hg
# rank_l <- rank_cpm
rank_l <- rank_ln


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


#
med <- median(gene_mat, na.rm = TRUE)
gene_mat_med <- gene_mat_na <- gene_mat
gene_mat_med[comsr == FALSE] <- med
gene_mat_na[comsr == FALSE] <- NA


# Organizing rank dataframe for different ways of handling NAs and aggregating.
# NAs: keep as tied RSR values, impute to median of matrix, or keep as NAs to be
# ignored in following aggregation)
# aggregation: average RSR, median RSR, average rank of RSR, median rank of RSR,
# and rank product

# TODO: Avg/Medrank_wo_imp/w_NA selects for all NAs
# TODO: rank product with NAs select for all non-NAs

rsr_avg_wo_imp <- rowMeans(gene_mat)
rsr_avg_w_imp <- rowMeans(gene_mat_med)
rsr_avg_w_na <- rowMeans(gene_mat_na, na.rm = TRUE)

rsr_med_wo_imp <- rowMeans(gene_mat)
rsr_med_w_imp <- apply(gene_mat_med, 1, median)
rsr_med_w_na <- apply(gene_mat_na, 1, function(x) median(na.omit(x)))

rsr_avgrank_wo_imp <- rowMeans(colrank_mat(gene_mat))
rsr_avgrank_w_imp <- rowMeans(colrank_mat(gene_mat_med))
rsr_avgrank_w_na <- rowMeans(colrank_mat(gene_mat_na), na.rm = TRUE)

rsr_medrank_wo_imp <- apply(colrank_mat(gene_mat), 1, median)
rsr_medrank_w_imp <- apply(colrank_mat(gene_mat_med), 1, median)
rsr_medrank_w_na <- apply(colrank_mat(gene_mat_na), 1, function(x) median(na.omit(x)))

# NOTE: Often run into integer overflow issue with rank product
rsr_rankprod_wo_imp <- rank(matrixStats::rowProds(colrank_mat(gene_mat)))
rsr_rankprod_w_imp <- rank(matrixStats::rowProds(colrank_mat(gene_mat_med)))
rsr_rankprod_w_na <- rank(matrixStats::rowProds(colrank_mat(gene_mat_na)))


df <- data.frame(
  
  Symbol = names(rsr_avg_wo_imp),
  
  Avg_wo_imp = rsr_avg_wo_imp,
  Avg_w_imp = rsr_avg_w_imp,
  Avg_w_na = rsr_avg_w_na,
  
  Med_wo_imp = rsr_med_wo_imp,
  Med_w_imp = rsr_med_w_imp,
  Med_w_na = rsr_med_w_na,
  
  Avgrank_wo_imp = rsr_avgrank_wo_imp,
  Avgrank_w_imp = rsr_avgrank_w_imp,
  Avgrank_w_na = rsr_avgrank_w_na,
  
  Medrank_wo_imp = rsr_medrank_wo_imp,
  Medrank_w_imp = rsr_medrank_w_imp,
  Medrank_w_na = rsr_medrank_w_na,
  
  Rankprod_wo_imp = rsr_rankprod_wo_imp,
  Rankprod_w_imp = rsr_rankprod_w_imp,
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
rank_order <- arrange(df, desc(Avg_w_imp))$Symbol

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
pheatmap(gene_mat_na[arrange(df, desc(Avg_w_imp))$Symbol, names(na_order)],
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


pqa <- qplot(df, yvar = "Avg_w_imp", xvar = "N_comeasured")


pqb <- df %>%
  mutate(Bin = cut_width(N_comeasured, width = 10, boundary = 0)) %>%
  ggplot(aes(x = Bin, y = Avg_w_imp)) +
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
