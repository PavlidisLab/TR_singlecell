# Query/subject all rank matrices
allrank_dir <- "/space/scratch/amorin/R_objects/04-07-2023/"
allrank_l <- lapply(list.files(allrank_dir, full.names = TRUE), function(x) readRDS(x))
names(allrank_l) <- str_replace(list.files(allrank_dir), "\\.RDS", "")



# Heatmap of query/subject topk all ranks


p6 <- pheatmap(allrank_mat,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               border_col = "black",
               color = c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b'),
               display_numbers = TRUE,
               number_format = "%1.0f",
               number_color = "black",
               na_col = "black",
               fontsize = 22,
               cellwidth = 50,
               cellheight = 50,
               height = 18,
               width = 18)



# Histogram of query/subject all ranks: Pile up at the left of the histogram 
# (lower/better ranks) suggests a signal, as random ranking would be uniform


p7 <- mat_to_df(allrank_mat, value_name = "Rank") %>% 
  ggplot(., aes(x = Rank)) +
  geom_histogram(bins = 30) +
  ylab("Count") +
  ggtitle(gene_hg) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))


# hist(replicate(length(allrank_mat), sample(1:length(pc_hg$Symbol), 1)), breaks = 20)


# Histogram of the example query/subject topk rank, with overlay of example gene


p8 <- data.frame(Topk = topk_rank$Topk) %>% 
  ggplot(aes(x = Topk)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = topk_rank$Topk[gene_hg], col = "red") +
  ggtitle(gene_hg) +
  ylab("Count of genes") +
  xlab("Topk intersect") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20))




# For inspecting all rank query/subject topk
# ------------------------------------------------------------------------------


# gene_hg <- "ASCL1"
# 
# 
# # Keep only measured genes and set diag to NA to remove self-ranking
# 
# keep_exp <- which(msr_hg[gene_hg, ] == 1)
# allrank_mat <- allrank_l[[gene_hg]][keep_exp, keep_exp, drop = FALSE]
# diag(allrank_mat) <- NA
# 
# 
# # Inspecting top pairs: by topk magnitude, and by rank. The highest value/most
# # similar pair across all experiment pairs by magnitude may not actually be the
# # top ranked gene within the subject dataset
# 
# sim_df <- sim_tf_hg[[gene_hg]]$Sim_df
# best_value_id <- slice_max(sim_df, Topk, n = 1)
# best_rank_ix <- which(allrank_mat == min(allrank_mat, na.rm = TRUE), arr.ind = TRUE)
# 
# query_vec <- load_agg_mat_list(ids = best_value_id$Row,
#                                genes = pc_hg$Symbol,
#                                sub_genes = gene_hg)[[1]][, gene_hg]
# 
# subject_mat <- load_agg_mat_list(ids = best_value_id$Col, genes = pc_hg$Symbol)[[1]]
# 
# 
# topk_rank <- query_gene_rank_topk(
#   query_vec = query_vec,
#   subject_mat = subject_mat,
#   gene = gene_hg,
#   ncores = ncore)
