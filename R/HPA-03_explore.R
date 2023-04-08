##
## -----------------------------------------------------------------------------

library(plyr)
library(tidyverse)
library(parallel)
library(cowplot)
library(pheatmap)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

#
expr_l <- readRDS(expr_mat_l_path)
cor_l <- readRDS(cor_mat_l_path)


# Ranked targets from paper
rank_path <- "/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS"
rank_l <- readRDS(rank_path)

# TFs
tf_hg <- read.delim("/space/grp/amorin/Metadata/human_tfs.tsv", stringsAsFactors = FALSE)
tfs <- c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4")


#
# ------------------------------------------------------------------------------


# Get a data frame of gene cor values for requested gene across all data sets 

cor_df <- function(cor_list, gene) {
  
  df_list <- lapply(cor_list, function(x) {
    
    if (!(gene %in% rownames(x))) {
      return(data.frame(Cor = NA, Symbol = rownames(x)))
    }
    
    vec <- x[gene, ]
    vec <- vec[names(vec) != gene]
    data.frame(Cor = vec) %>% rownames_to_column(var = "Symbol")
  })
  
  df <- plyr::join_all(df_list, by = "Symbol")
  colnames(df)[2:ncol(df)] <- names(cor_list)
  
  return(df)
}



# Convert df of TF-gene cors to ranks

rank_cor <- function(cor_df) {
  
  rank_l <- lapply(2:ncol(cor_df), function(x) rank(-cor_df[, x]))  # TODO this doesn't need to be lapply, no?
  mat <- do.call(cbind, rank_l)
  rownames(mat) <- cor_df$Symbol
  colnames(mat) <- colnames(cor_df)[2:ncol(cor_df)]
  mat <- cbind(mat, Average_rank = rowMeans(mat, na.rm = TRUE))
  mat <- mat[order(mat[, "Average_rank"]), ]
  
  return(mat)
}



# Variability of TFs
# ------------------------------------------------------------------------------


#  Return a df of the row/gene-wise mean/sd/CV and whether the gene is a TF

get_coefvar_df <- function(df, tf_vec) {
  
  gene_mean <- rowMeans(df)
  gene_sd <- apply(df, 1, sd)
  gene_cv <- gene_sd/gene_mean
  
  data.frame(
    Symbol = rownames(df), 
    Mean = gene_mean, 
    SD = gene_sd, 
    CV = gene_cv,
    TF = rownames(df) %in% tf_vec)
  
}


# Average the mean/SD/CV across all data sets

get_avg_coefvar_df <- function(cv_l) {
  
  df <- plyr::join_all(cv_l, by = "Symbol")
  
  gene_mean <- rowMeans(df[, str_detect(colnames(df), "Mean")], na.rm = TRUE)
  gene_sd <- rowMeans(df[, str_detect(colnames(df), "SD")], na.rm = TRUE)
  gene_cv <- rowMeans(df[, str_detect(colnames(df), "CV")], na.rm = TRUE)
  n_na <- apply(df[, str_detect(colnames(df), "Mean")], 1, function(x) sum(is.na(x)))
  
  data.frame(Symbol = df$Symbol,
             Mean_mean = gene_mean,
             Mean_SD = gene_sd,
             Mean_CV = gene_cv,
             Count_NA = n_na,
             TF = df$TF)
  
}



cv_l <- lapply(expr_l, get_coefvar_df, tf_hg)



avg_cv <- get_avg_coefvar_df(cv_l)





# Example of (sqrt for plotting) CV vs mean for HPA 
cv_l$HPA %>% 
  arrange(desc(CV)) %>% 
  mutate(Symbol = factor(Symbol, levels = unique(Symbol))) %>% 
  ggplot(aes(x = Mean, y = sqrt(CV))) +
  geom_point(shape = 21) +
  xlab("Average expression") +
  ylab("Square root of coefficient of variation") +
  ggtitle("HPA = 253 tissues") +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 25))



ggplot(tt, aes(x = sqrt(Mean_CV), colour = TF)) +
  geom_density()


plot(density(cv_l$HPA$CV))
hist(cv_l$HPA$CV, breaks = 100)
boxplot(cv_l$HPA$CV ~ cv_l$HPA$TF)
boxplot(log(cv_l$HPA$CV) ~ cv_l$HPA$TF)



cv_l$HPA %>% 
  ggplot(aes(x = CV, colour = TF)) +
  geom_density()



cv_l$HPA %>% 
  arrange(desc(CV)) %>% 
  mutate(Symbol = factor(Symbol, levels = unique(Symbol))) %>% 
  ggplot(aes(x = Symbol, y = CV)) +
  geom_point()



cv_l$HPA %>% 
  arrange(desc(CV)) %>% 
  mutate(Symbol = factor(Symbol, levels = unique(Symbol))) %>% 
  ggplot(aes(x = Mean, y = sqrt(CV))) +
  geom_point()



cv_l$HPA %>% 
  filter(TF) %>% 
  arrange(CV) %>% 
  head(20)


cv_l$HPA %>% 
  filter(TF) %>% 
  arrange(desc(CV)) %>% 
  head(20)



#
# ------------------------------------------------------------------------------


tf <- "ASCL1"
tf_cor <- cor_df(cor_l, tf)
tf_cor_rank <- rank_cor(tf_cor)


# Get the "cor of cors" - how consistent TR-gene coexpr is across data sets

all_cor <- cor(dplyr::select_if(tf_cor, is.numeric), 
               use = "pairwise.complete.obs")



#
# ------------------------------------------------------------------------------


rank_df <- left_join(
  rank_l$Human[[tf]][, c("Symbol", "Rank_integrated", "Curated_target")],
  rownames_to_column(data.frame(tf_cor_rank), var = "Symbol"),
  by = "Symbol") %>% 
  filter(Symbol != tf)
  

keep_cols <- c("Rank_integrated", "Average_rank", names(expr_l))


pr_df <- all_perf_df(rank_df, keep_cols, label_col = "Curated_target", measure = "PR")
auprc <- all_au_perf(rank_df, keep_cols, label_col = "Curated_target", measure = "AUPRC")

roc_df <- all_perf_df(rank_df, keep_cols, label_col = "Curated_target", measure = "ROC")
auroc <- all_au_perf(rank_df, keep_cols, label_col = "Curated_target", measure = "AUROC")


cols <- c(rep("lightgrey", length(keep_cols)))
names(cols) <- keep_cols
cols["Rank_integrated"] <- "black"
cols["Average_rank"] <- "red"


plot_perf(df = roc_df, auc_l = auroc, measure = "ROC", cols = cols, title = tf, ncol_legend = 2)
plot_perf(df = pr_df, auc_l = auprc, measure = "PR", cols = cols, title = tf, ncol_legend = 2)



# Plotting
# ------------------------------------------------------------------------------


# PAX6 - PCSK2 perfect cor due to sparsity, but does have good evidence otherwise


# TODO: 

gene <- "DLL1"

p_list <- lapply(names(expr_l), function(x) {
  
  mat <- expr_l[[x]]
  
  if(!all(c(tf, gene) %in% rownames(mat))) {
    return(NA)
  }
  
  df <- data.frame(tf = mat[tf,], gene = mat[gene, ])
  ggplot(df, aes(x = tf, y = gene)) + 
    geom_point(size = 2, shape = 21, fill = "#756bb1") +
    geom_smooth(method = "lm", colour = "black") +
    xlab(tf) +
    ylab(gene) +
    ggtitle(paste0(x, ": ", round(cor_l[[x]][tf, gene], 3))) +
    theme_classic() +
    theme(plot.title = element_text(size = 20))
})


p_list <- p_list[!is.na(p_list)]

cowplot::plot_grid(plotlist = p_list)

head(sort(expr_l$FANTOM[gene, ], decreasing = TRUE), 10)
head(sort(expr_l$FANTOM[tf, ], decreasing = TRUE), 10)


# Heatmap of cor of cor. Pre-cluster to get order for triangular mat

hc <-  hclust(as.dist(1 - all_cor))
all_cor <- all_cor[hc$order, hc$order]
all_cor[lower.tri(all_cor)] <-  NA
diag(all_cor) <- NA
all_cor <- all_cor[1:nrow(all_cor)-1, 2:ncol(all_cor)]


pal_length <- 100
# bluered_pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(pal_length)
bluered_pal <- colorRampPalette(c("#ca0020", "white", "#0571b0"))(pal_length)
color_breaks <- seq(min(all_cor, na.rm = TRUE), max(all_cor, na.rm = TRUE), length.out = pal_length)


pheatmap(all_cor, 
         cluster_col = FALSE, 
         cluster_row = FALSE,
         color = bluered_pal,
         breaks = color_breaks,
         na_col = "white",
         border_color = NA,
         fontsize = 18)
