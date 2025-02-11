## 
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
library(parallel)
library(aggtools)
library(ggrepel)
source("R/utils/vector_comparison_functions.R")
source("R/utils/functions.R")
source("R/00_config.R")

k <- 200
agg_method = "allrank"


# Table of assembled scRNA-seq datasets
sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)

# Protein coding genes and TFs
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
tfs_mm <- read.delim(tfs_mm_path, stringsAsFactors = FALSE)


file <- file.path("/space/scratch/amorin/TRsc_output/",
                  paste0("GSE132042_10X_smartseq_comparison_", agg_method, ".RDS"))

# Measurement matrices used for filtering when a gene was never expressed
# msr_mm <- readRDS(msr_mat_mm_path)


dat_ss <- load_dat_list("GSE132042SmartSeq2")[[1]]
dat_10x <- load_dat_list("GSE132042")[[1]]

# table(dat_ss$Meta$assay)
# table(dat_10x$Meta$assay)


common_counts <- inner_join(
  count(dat_ss$Meta, Cell_type),
  count(dat_10x$Meta, Cell_type),
  by = "Cell_type", suffix = c("_smartseqV2", "_10XV2")
) %>% 
  mutate(Cell_type = as.character(Cell_type)) %>% 
  filter(n_smartseqV2 >= 100 & n_10XV2 >= 100)


common_cts <- common_counts$Cell_type


meta_ss <- filter(dat_ss$Meta, Cell_type %in% common_cts)
mat_ss <- dat_ss$Mat[, meta_ss$ID]


meta_10x <- filter(dat_10x$Meta, Cell_type %in% common_cts)
mat_10x <- dat_10x$Mat[, meta_10x$ID]


if (!file.exists(file)) {
  
  agg_ss <- aggtools::aggr_coexpr_single_dataset(
    mat = mat_ss,
    meta = meta_ss,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = agg_method
  )
  
  agg_10x <- aggtools::aggr_coexpr_single_dataset(
    mat = mat_10x,
    meta = meta_10x,
    pc_df = pc_mm,
    cor_method = "pearson",
    agg_method = agg_method
  )
  
  agg_l <- list(Agg_SS = agg_ss, Agg_10x = agg_10x)
  
  saveRDS(agg_l, file)
  
} else {
  
  agg_l <- readRDS(file)
  
}



# Examining TF expression
# ------------------------------------------------------------------------------


mat_ss_log <- log2(mat_ss + 1)
mat_10x_log <- log2(mat_10x + 1)


avg_expr_mat <- function(mat, meta, cts, tfs) {
  
  avg_l <- lapply(cts, function(ct) {
    
    ct_ids <- filter(meta, Cell_type == ct)$ID
    ct_mat <- mat[tfs, ct_ids]
    rowMeans(ct_mat)
    
  })
  
  avg_mat <- as.matrix(do.call(cbind, avg_l))
  colnames(avg_mat) <- common_cts
  rownames(avg_mat) <- tfs
  
  return(avg_mat)
}


avg_ss <- avg_expr_mat(mat_ss_log, meta_ss, common_cts, tfs_mm$Symbol)
avg_10x <- avg_expr_mat(mat_10x_log, meta_10x, common_cts, tfs_mm$Symbol)
gc()


all0_ss <- which(rowSums(avg_ss) == 0)
all0_10x <- which(rowSums(avg_10x) == 0)
keep <- setdiff(tfs_mm$Symbol, intersect(names(all0_ss), names(all0_10x)))


avg_ss <- avg_ss[keep, ]
avg_10x <- avg_10x[keep, ]


avg_topk <- pair_colwise_topk(avg_ss, avg_10x, k = k, ncores = ncore)
avg_scor <- pair_colwise_cor(avg_ss, avg_10x, cor_method = "spearman", ncores = ncore)


avg_topk_shuffle <- unlist(lapply(1:10, function(iter) {
  pair_shuffle_topk(avg_ss, avg_10x, k = k, ncores = ncore)
}))


plot(density(avg_topk_shuffle))
abline(v = mean(avg_topk))


colnames(avg_ss) <- paste0(common_cts, "_SS")
colnames(avg_10x) <- paste0(common_cts, "_10x")
all_topk_mat <- colwise_topk_intersect(cbind(avg_ss, avg_10x), k = k)
all_topk_df <- mat_to_df(all_topk_mat)

all_scor_mat <- cor(cbind(avg_ss, avg_10x), method = "spearman")
all_scor_df <- mat_to_df(all_scor_mat)



# 
# ------------------------------------------------------------------------------


ggplot(common_counts, aes(x = log10(n_smartseqV2), y = log10(n_10XV2))) +
  geom_point(shape = 21, size = 3) +
  ylab("log10 Count of cells 10X V2") +
  xlab("log10 Count of cells Smartseq2") +
  theme_classic() +
  theme(text = element_text(size = 20))





agg_topk <- pair_colwise_topk(mat1 = agg_l$Agg_SS$Agg_mat, 
                              mat2 = agg_l$Agg_10x$Agg_mat, 
                              k = k, 
                              ncores = ncore)


# agg_shuffle_topk <- pair_shuffle_topk(mat1 = agg_l$Agg_SS$Agg_mat, 
#                                       mat2 = agg_l$Agg_10x$Agg_mat, 
#                                       k = k, 
#                                       ncores = ncore)



agg_shuffle_topk <- pair_shuffle_topk(mat1 = agg_l$Agg_SS$Agg_mat[, keep], 
                                      mat2 = agg_l$Agg_10x$Agg_mat[, keep], 
                                      k = k, 
                                      ncores = ncore)


plot(density(agg_shuffle_topk), col = "red")
lines(density(agg_topk[keep]), col = "black")

summary(agg_shuffle_topk)
summary(agg_topk)
summary(agg_topk[keep])

wilcox.test(x = agg_shuffle_topk, y = agg_topk[keep])


data.frame(
  Topk = c(agg_topk[keep], agg_shuffle_topk),
  Group = c(rep("TR", length(agg_topk[keep])), rep("Shuffled", length(agg_shuffle_topk)))
) %>% 
  ggplot(., aes(x = Topk, y = Group)) +
  geom_boxplot() +
  xlab("Top200") +
  theme_classic() +
  theme(text = element_text(size = 20),
        axis.title.y = element_blank())


mat1 <- agg_l$Agg_SS$Agg_mat[, keep]
mat2 <- agg_l$Agg_10x$Agg_mat[, keep]

mat1[is.na(mat1) | is.infinite(mat1)] <- 0
mat2[is.na(mat2) | is.infinite(mat2)] <- 0




agg_scor <- cor(mat1, 
                mat2,
                method = "spearman")


agg_shuffle_scor <- pair_shuffle_cor(mat1 = mat1, 
                                     mat2 = mat2, 
                                     cor_method = "spearman",
                                     ncores = ncore)

plot(density(agg_shuffle_scor, na.rm = TRUE), xlim = c(-0.4, 1))
lines(density(diag(agg_scor), na.rm = TRUE), col = "red")


check_tf <- "Hmmr"
df <- data.frame(Agg_10x = agg_l$Agg_SS$Agg_mat[, check_tf],
                 Agg_SS = agg_l$Agg_10x$Agg_mat[, check_tf])

plot(df$Agg_10x, df$Agg_SS)
