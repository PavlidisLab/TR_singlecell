## Examine the topk overlap of sampled TFs
## -----------------------------------------------------------------------------

source("R/05b_save_sampled_null_topk.R")


# Null topk overlap
null_topk_hg <- readRDS(null_topk_hg_path)
null_topk_mm <- readRDS(null_topk_mm_path)

# List of paired experiment similarities for TFs
sim_tf_hg <- readRDS(sim_tf_hg_path)
sim_tf_mm <- readRDS(sim_tf_mm_path)

# Sample matched to mouse Ascl1 (n=66 exp)
null_ascl1_mm <- readRDS("~/scratch/R_objects/03-09-2023_topk_size_matched_null.RDS")


# Looking at the distribution of count of experiments in which TFs were measured

tf_msr_mm <- rowSums(msr_mm[tfs_mm$Symbol, ])
tf_msr_hg <- rowSums(msr_hg[tfs_hg$Symbol, ])

summary(tf_msr_mm)
summary(tf_msr_hg)

hist(tf_msr_mm, breaks = 100, cex.axis = 2, cex.lab = 2, main = NA)
hist(tf_msr_hg, breaks = 100, cex.axis = 2, cex.lab = 2, main = NA)


# Looking at the distribution of how many TFs were measured in each experiment

exp_msr_mm <- colSums(msr_mm[tfs_mm$Symbol, ])
exp_msr_hg <- colSums(msr_hg[tfs_hg$Symbol, ])

summary(exp_msr_mm)
summary(exp_msr_hg)

hist(exp_msr_mm, breaks = 30, cex.axis = 2, cex.lab = 2, main = NA)
hist(exp_msr_hg, breaks = 30, cex.axis = 2, cex.lab = 2, main = NA)


# Count of TFs that were always measured

sum(tf_msr_mm == ncol(msr_mm))
sum(tf_msr_hg == ncol(msr_hg))


# Count of TFs that were never measured

sum(tf_msr_mm == 0)
sum(tf_msr_hg == 0)

# Rere as maximally measured TF

which.max(sapply(sim_tf_mm, nrow))
which.max(rowSums(msr_mm[tfs_mm$Symbol, ]))
sum(msr_mm["Rere", ])

identical(
  as.integer((length(ids_mm) * ( length(ids_mm) - 1 ) / 2)),
  as.integer(nrow(sim_tf_mm$Rere))
)


# Inspecting one realization of a sample


# For reference

# sample_topk_intersect <- function(agg_l,
#                                   genes,
#                                   msr_mat,
#                                   k = 1000,
#                                   check_k_arg = TRUE) {
#   
#   ids <- names(agg_l)
#   
#   # Sample one gene that is measured in each data set
#   sample_genes <- unlist(lapply(ids, function(x) {
#     msr_gene <- names(which(msr_mat[genes, x] == 1))
#     sample(msr_gene, 1)
#   }))
#   
#   # Bind sampled genes into a matrix
#   sample_mat <- lapply(1:length(sample_genes), function(x) {
#     agg_l[[x]][, sample_genes[x]]
#   })
#   sample_mat <- do.call(cbind, sample_mat)
#   colnames(sample_mat) <- paste0(ids, "_", sample_genes)
#   
#   # Get topk overlap between sampled genes
#   sample_topk <- colwise_topk_intersect(sample_mat, k = k, check_k_arg = check_k_arg)
#   sample_df <- mat_to_df(sample_topk, symmetric = TRUE, value_name = "Topk")
#   
#   return(sample_df)
# }


set.seed(5)


agg_l <- agg_tf_mm
genes <- tfs_mm$Symbol
msr_mat <- msr_mm
k <- 1000
check_k_arg <-  TRUE


ids <- names(agg_l)

# Sample one gene that is measured in each data set
sample_genes <- unlist(lapply(ids, function(x) {
  msr_gene <- names(which(msr_mat[genes, x] == 1))
  sample(msr_gene, 1)
}))

# Bind sampled genes into a matrix
sample_mat <- lapply(1:length(sample_genes), function(x) {
  agg_l[[x]][, sample_genes[x]]
})
sample_mat <- do.call(cbind, sample_mat)
colnames(sample_mat) <- paste0(ids, "_", sample_genes)

# Get topk overlap between sampled genes
sample_topk <- colwise_topk_intersect(sample_mat, k = k, check_k_arg = check_k_arg)
sample_df <- mat_to_df(sample_topk, symmetric = TRUE, value_name = "Topk")



# Distn of top k for one sample
hist(sample_df$Topk, breaks = 30, cex.axis = 2, cex.lab = 2, main = NA)
summary(sample_df$Topk)

# Best and bottom pairs (likely ties at 0)
top_pair <- slice_max(sample_df, Topk)
bottom_pair <- slice_min(sample_df, Topk)

head(sort(sample_mat[, top_pair$Row], decreasing = TRUE), 10)
head(sort(sample_mat[, top_pair$Col], decreasing = TRUE), 10)
plot(sample_mat[, top_pair$Row], sample_mat[, top_pair$Col])

head(sort(sample_mat[, bottom_pair$Row[1]], decreasing = TRUE), 10)
head(sort(sample_mat[, bottom_pair$Col[1]], decreasing = TRUE), 10)
plot(sample_mat[, bottom_pair$Row[1]], sample_mat[, bottom_pair$Col[1]])


topk_intersect(
  topk_sort(sample_mat[, top_pair$Row], k = 1000),
  topk_sort(sample_mat[, top_pair$Col], k = 1000)
)


topk_intersect(
  topk_sort(sample_mat[, bottom_pair$Row[1]], k = 1000),
  topk_sort(sample_mat[, bottom_pair$Col[1]], k = 1000)
)



# Example where check topk is used. 0 overlap, but if check topk is not used,
# 5 genes overlap despite being all NAs in "GSE132042SmartSeq2_Tbx3"


vec1 <- agg_tf_mm[["GSE132042SmartSeq2"]][, "Tbx3"]
vec2 <- agg_tf_mm[["MKA"]][, "Hoxc6"]


topk_sort(vec1, k = 1000, check_k_arg = TRUE)
topk_sort(vec1, k = 1000, check_k_arg = FALSE)

sort(vec1, decreasing = TRUE)[363:367]
sort(vec2, decreasing = TRUE)[999:1002]

intersect(
  topk_sort(vec1, k = 1000, check_k_arg = TRUE),
  topk_sort(vec2, k = 1000, check_k_arg = TRUE)
)


intersect(
  topk_sort(vec1, k = 1000, check_k_arg = FALSE),
  topk_sort(vec2, k = 1000, check_k_arg = FALSE)
)


# Organize all samples into a single df for plotting/inspection 

null_df <- do.call(
  rbind,
  lapply(1:length(null_topk_mm), function(x) mutate(null_topk_mm[[x]], ID = paste0("Rep_", x)))
)


# Boxplot of every sample realization

ggplot(null_df, aes(x = reorder(ID, Topk, FUN = median), y = Topk)) +
  geom_boxplot(outlier.shape = NA, coef = 0) +
  xlab("Sample_ID") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())



null_median_topk <- null_df %>% 
  group_by(ID) %>% 
  summarize(med = median(Topk)) %>% 
  arrange(desc(med))

top_median <- as.integer(str_replace(null_median_topk$ID[1], "Rep_", ""))


# Look at the distn of the sample with the highest median topk
hist(null_topk_mm[[top_median]]$Topk, breaks = 30, cex.axis = 2, cex.lab = 2, main = NA)
summary(null_topk_mm[[top_median]]$Topk)
head(arrange(null_topk_mm[[top_median]], desc(Topk)))
head(arrange(null_topk_mm[[top_median]], Topk))

# Inspect the highest global sampled pair
head(arrange(null_df, desc(Topk)))


# Inpsect observed Ascl1 topk compared to a null sample that considers only
# experiments measuring Ascl1 (n=66) instead of all (n=104)

head(arrange(sim_tf_mm$Ascl1, desc(Topk)))


null_df_ascl1 <- do.call(
  rbind,
  lapply(1:length(null_ascl1_mm), function(x) mutate(null_ascl1_mm[[x]], ID = paste0("Rep_", x)))
)


null_median_topk_ascl1 <- null_df_ascl1 %>% 
  group_by(ID) %>% 
  summarize(med = median(Topk)) %>% 
  arrange(desc(med))

top_median_ascl1 <- as.integer(str_replace(null_median_topk_ascl1$ID[1], "Rep_", ""))


# Organizing observed Ascl1 similarity pairs, null using all experiments,
# and null using only Ascl1-expressing experiments. 

topk_ascl1 <- data.frame(Topk = sim_tf_mm$Ascl1$Topk, ID = "Ascl1")

topk_null_all_best <- data.frame(Topk = null_topk_mm[[top_median]]$Topk, 
                                 ID = "All_null_best")

topk_null_all_global <- data.frame(Topk = null_df$Topk, 
                                 ID = "All_null_global")

topk_null_sub_best <- data.frame(Topk = null_ascl1_mm[[top_median_ascl1]]$Topk, 
                                 ID = "Sub_null_best")

topk_null_sub_global <- data.frame(Topk = null_df_ascl1$Topk, 
                                   ID = "Sub_null_global")



summ_null_all_global <- summary(topk_null_all_global$Topk)
summ_null_all_best <- summary(topk_null_all_best$Topk)
summ_null_sub_global <- summary(topk_null_sub_global$Topk)
summ_null_sub_best <- summary(topk_null_sub_best$Topk)
summ_ascl1 <- summary(sim_tf_mm$Ascl1$Topk)


plot_df <- do.call(rbind, list(topk_ascl1,
                               topk_null_all_best,
                               topk_null_all_global,
                               topk_null_sub_best,
                               topk_null_sub_global))


groups <- c("Ascl1", "All_null_best", "All_null_global", "Sub_null_best", "Sub_null_global")
plot_df$ID <- factor(plot_df$ID, levels = unique(groups))



ggplot(plot_df, aes(x = ID, y = Topk)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.margin = margin(c(10, 10, 10, 10)))



ggplot(plot_df, aes(x = Topk, colour = ID)) +
  geom_density(lwd = 2) +
  scale_colour_manual(values = c('#fb9a99', '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c')) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)),
        legend.position = c(0.6, 0.75))




