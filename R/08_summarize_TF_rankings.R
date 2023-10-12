## TODO
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(WGCNA)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
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



# Create a gene x TF matrix of summarized ranks
# ------------------------------------------------------------------------------



rank_mat_hg <- as.matrix(do.call(cbind, lapply(rank_tf_hg, `[`, "Rank_RSR")))
colnames(rank_mat_hg) <- names(rank_tf_hg)


rank_mat_mm <- as.matrix(do.call(cbind, lapply(rank_tf_mm, `[`, "Rank_RSR")))
colnames(rank_mat_mm) <- names(rank_tf_mm)




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



# Global top
# ------------------------------------------------------------------------------


top_hg <- do.call(rbind, lapply(rank_tf_hg, slice_min, Rank_RSR))
top_mm <- do.call(rbind, lapply(rank_tf_mm, slice_min, Rank_RSR))

hist(top_hg$Topk_count)
hist(top_hg$Topk_proportion)


# Inspecting targets from the most similar TFs
arrange(rank_tf_hg$E2F8, desc(Avg_RSR)) %>% head(30)
arrange(rank_tf_mm$E2f8, desc(Avg_RSR)) %>% head(30)




rank_cor_hg <- mat_to_df(colwise_cor(rank_mat_hg), symmetric = TRUE)
rank_topk_hg <- mat_to_df(colwise_topk_intersect(rank_mat_hg, k = 1000), symmetric = TRUE)

rank_hg <- data.frame(Row = rank_cor_hg$Row, 
                      Col = rank_cor_hg$Col,
                      Topk = rank_topk_hg$Value,
                      Cor = rank_cor_hg$Value)



rank_cor_mm <- mat_to_df(colwise_cor(rank_mat_mm), symmetric = TRUE)
rank_topk_mm <- mat_to_df(colwise_topk_intersect(rank_mat_mm, k = 1000), symmetric = TRUE)


rank_mm <- data.frame(Row = rank_cor_mm$Row, 
                      Col = rank_cor_mm$Col,
                      Topk = rank_topk_mm$Value,
                      Cor = rank_cor_mm$Value)

# Most commonly highly ranked genes
# TODO: are these genes high because of technical aspect?

rank_order_hg <- sort(rowMeans(rank_mat_hg, na.rm = TRUE))
head(rank_order_hg)
tail(rank_order_hg)

head(sort(rank_mat_hg["TRIM28",]), 50)
tail(sort(rank_mat_hg["CSF1R",]), 50)


rank_order_mm <- sort(rowMeans(rank_mat_mm, na.rm = TRUE))
head(rank_order_mm)
tail(rank_order_mm)


head(sort(rank_mat_mm["Hes6",]), 50)
tail(sort(rank_mat_mm["Pax6",]), 50)





# Inspect ortho
# ------------------------------------------------------------------------------


tfs_ortho <- filter(pc_ortho, Symbol_hg %in% names(rank_tf_hg) & Symbol_mm %in% names(rank_tf_mm))
  

df_hg <- left_join(rank_tf_hg$ASCL1, pc_ortho, by = c("Symbol" = "Symbol_hg")) %>% filter(!is.na(ID))
df_mm <- left_join(rank_tf_mm$Ascl1, pc_ortho, by = c("Symbol" = "Symbol_mm")) %>% filter(!is.na(ID))
ortho_df <- left_join(df_hg, df_mm, by = "ID", suffix = c("_Human", "_Mouse")) %>% filter(!is.na(Avg_RSR_Human) & !is.na(Avg_RSR_Mouse))
plot(ortho_df$Avg_RSR_Human, ortho_df$Avg_RSR_Mouse)
cor(ortho_df$Avg_RSR_Human, ortho_df$Avg_RSR_Mouse, method = "spearman")


ortho_cor <- mclapply(1:nrow(tfs_ortho), function(x) {
  
  df_hg <- left_join(rank_tf_hg[[tfs_ortho$Symbol_hg[x]]], pc_ortho, by = c("Symbol" = "Symbol_hg")) %>% filter(!is.na(ID))
  df_mm <- left_join(rank_tf_mm[[tfs_ortho$Symbol_mm[x]]], pc_ortho, by = c("Symbol" = "Symbol_mm")) %>% filter(!is.na(ID))
  ortho_df <- left_join(df_hg, df_mm, by = "ID", suffix = c("_Human", "_Mouse")) %>% filter(!is.na(Avg_RSR_Human) & !is.na(Avg_RSR_Mouse))
  cor(ortho_df$Avg_RSR_Human, ortho_df$Avg_RSR_Mouse, method = "spearman")
  
}, mc.cores = ncore)


ortho_cor <- data.frame(Cor = unlist(ortho_cor), Symbol = tfs_ortho$Symbol_hg)


px <- plot_hist(ortho_cor, stat_col = "Cor", xlab = "Similarity of TF coexpression partners", title = "n=1241 orthologous TFs") + 
  xlim(c(-0.4, 1))


ggsave(px, height = 5, width = 7, device = "png", dpi = 300,
       filename = file.path(plot_dir, "ortho_rank_similarity.png"))



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





# Distn of top k

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



