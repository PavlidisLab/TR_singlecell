## TODO: all cor df should also tack on rank
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


# Ranked targets from paper and data matrices
rank_l <- readRDS("/space/scratch/amorin/R_objects/ranked_target_list_Apr2022.RDS")
dat_l <- readRDS("/space/scratch/amorin/R_objects/all_data_list_Apr2022.RDS")


# TFs
tf_hg <- read.delim("/space/grp/amorin/Metadata/human_tfs.tsv", stringsAsFactors = FALSE)
tfs <- str_to_upper(c("Ascl1", "Hes1", "Mecp2", "Mef2c", "Neurod1", "Pax6", "Runx1", "Tcf4"))


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
# TODO this doesn't need to be lapply, no?

rank_cor <- function(cor_df) {

  rank_l <- lapply(2:ncol(cor_df), function(x) rank(-cor_df[, x])) 
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
  
  genes <- rownames(df)
  tf_vec <- intersect(tf_vec, genes)
  
  gene_mean <- rowMeans(df)
  gene_sd <- apply(df, 1, sd)
  gene_cv <- gene_sd/gene_mean
  
  # Percentile rank relative to other TFs versus all genes
  
  pr_sd_all <- rank(gene_sd) / length(genes)
  pr_sd_tf <- rep(NA, length(pr_sd_all))
  names(pr_sd_tf) <- genes
  pr_sd_tf[tf_vec] <- rank(gene_sd[tf_vec]) / length(tf_vec)
  
  pr_cv_all <- rank(gene_cv) / length(genes)
  pr_cv_tf <- rep(NA, length(pr_cv_all))
  names(pr_cv_tf) <- genes
  pr_cv_tf[tf_vec] <- rank(gene_cv[tf_vec]) / length(tf_vec)
  
  data.frame(
    Symbol = rownames(df), 
    Mean = gene_mean, 
    SD = gene_sd, 
    CV = gene_cv,
    Perc_rank_SD_all = pr_sd_all,
    Perc_rank_SD_TF = pr_sd_tf,
    Perc_rank_CV_all = pr_cv_all,
    Perc_rank_CV_TF = pr_cv_tf,
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


# List of CVs for each data set
cv_l <- lapply(expr_l, get_coefvar_df, tf_hg$Symbol)


# Dataframe of average CV across data sets
avg_cv <- get_avg_coefvar_df(cv_l)


# Inspecting batch 1 TFs 
cv_l_tf <- lapply(cv_l, function(x) filter(x, Symbol %in% tfs))
avg_cv_tf <- filter(avg_cv, Symbol %in% tfs)


# stats into a mat for heatmap

stat_to_mat <- function(cv_l, stat, subset_genes) {
  
  stat_l <- lapply(cv_l, function(x) {
    vec <- x[subset_genes, stat]
    names(vec) <- subset_genes
    return(vec)
  })
  mat <- do.call(cbind, stat_l)
  
  return(mat)
}



pheatmap(stat_to_mat(cv_l, "Perc_rank_CV_TF", tfs),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 30,
         cellheight = 30,
         fontsize = 20,
         na_col = "darkgrey",
         border_color = "black")



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


# Distribution of averaged, keeping only genes measured in at least 4 data sets

avg_cv_sub <- filter(avg_cv, Count_NA <= 9)
table(avg_cv_sub$TF)
wilcox.test(avg_cv$Mean_CV ~ avg_cv$TF)
wilcox.test(avg_cv_sub$Mean_CV ~ avg_cv_sub$TF)


avg_cv_sub %>%
  ggplot(aes(x = sqrt(Mean_CV), fill = TF)) +
  # ggplot(aes(x = Mean_SD, fill = TF)) +
  geom_density(alpha = 0.4, colour = "black") +
  scale_fill_manual(values = c("grey", "red")) +
  xlab("Square root of average coefficient of variation") +
  # xlab("Average SD") +
  ylab("Density") +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = c(0.85, 0.2))


avg_cv_sub %>%
  ggplot(aes(x = TF, y = sqrt(Mean_CV))) +
  geom_boxplot(width = 0.3) +
  xlab("TF") +
  ylab("Square root of average coefficient of variation") +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20))


# Most var TFs

avg_cv_sub %>% 
  filter(TF, Mean_mean > 0.1) %>%
  arrange(desc(Mean_CV)) %>%
  # arrange(desc(Mean_SD)) %>% 
  head(20)


stat_to_mat(cv_l, "SD", "ANHX")
stat_to_mat(cv_l, "Mean", "ANHX")
stat_to_mat(cv_l, "CV", "ANHX")


# Least var TFs

avg_cv_sub %>% 
  filter(TF) %>% 
  arrange(Mean_CV, desc(Mean_SD)) %>% 
  head(20)


#
# ------------------------------------------------------------------------------


tf <- "ASCL1"
tf_cor <- cor_df(cor_l, tf)
tf_cor_rank <- rank_cor(tf_cor)


# Get the "cor of cors" - how consistent TR-gene coexpr is across data sets

all_cor <- cor(dplyr::select_if(tf_cor, is.numeric), 
               use = "pairwise.complete.obs",
               method = "spearman")


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



# Looking for repression
# ------------------------------------------------------------------------------


get_negcor_rank <- function(cor_df) {
  
  num <- select_if(cor_df, is.numeric)
  rmat <- colrank_mat(-num)
  rownames(rmat) <- cor_df$Symbol
  rmat <- cbind(rmat, Negcor_rank = rank(rowMeans(rmat, na.rm = TRUE)))
  rmat <- rmat[order(rmat[, "Negcor_rank"]), ]
  
  return(rmat)
}


tf_negcor <- get_negcor_rank(tf_cor)


rank_df2 <- left_join(
  rank_l$Human[[tf]],
  rownames_to_column(data.frame(tf_negcor), var = "Symbol"),
  by = "Symbol")


# First check top binding evidence - prioritize binding strength and neg cor

rank_df2 %>% 
  filter(Rank_binding <= 500 & Negcor_rank <= 500) %>% 
  arrange(Negcor_rank) %>% 
  head(30)


filter(rank_l$Ortho$ASCL1, Symbol == "KANK2_Kank2")




# Plotting
# ------------------------------------------------------------------------------


# PAX6 - PCSK2 perfect cor due to sparsity, but does have good evidence otherwise


# TODO: 

gene <- "KANK2"

p_list <- lapply(names(expr_l), function(x) {
  
  mat <- expr_l[[x]]
  
  if (!all(c(tf, gene) %in% rownames(mat))) {
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


##


bind_hg <- dat_l$Binding$Human$QN_log[gene, ]
bind_mm <- dat_l$Binding$Mouse$QN_log[str_to_title(gene), ]

bind_df <- data.frame(
  Experiment_ID = c(names(bind_hg), names(bind_mm)),
  Binding_score = c(bind_hg, bind_mm)
) %>% 
  left_join(dat_l$Binding$Meta[, c("Experiment_ID", "Symbol", "Species")],
            by = "Experiment_ID") %>% 
  mutate(TF = str_to_upper(Symbol) == str_to_upper(tf))




fc_hg <- dat_l$Perturbation$Human$FC_mat[gene, ]
fc_mm <- dat_l$Perturbation$Mouse$FC_mat[str_to_title(gene), ]
de_hg <- dat_l$Perturbation$Human$FDR_mat[gene, ] < 0.1
de_mm <- dat_l$Perturbation$Mouse$FDR_mat[str_to_title(gene), ] < 0.1


perturb_df <- data.frame(
  Experiment_ID = c(names(fc_hg), names(fc_mm)),
  FC = c(fc_hg, fc_mm),
  Abs_FC = abs(c(fc_hg, fc_mm)),
  DE = c(de_hg, de_mm)
) %>% 
  left_join(dat_l$Perturbation$Meta[, c("Experiment_ID", "Symbol", "Species", "Perturbation")],
            by = "Experiment_ID") %>% 
  mutate(TF = str_to_upper(Symbol) == str_to_upper(tf))



# get_gene_stats <- function(gene,
#                            exp_bind_tf, 
#                            exp_bind_ct, 
#                            mat_bind_tf, 
#                            exp_perturb_tf,
#                            exp_perturb_ct, 
#                            mat_fc_tf, 
#                            mat_fdr, 
#                            fdr = 0.1) {
#   list(
#     Binding = data.frame(
#       Experiment = exp_bind_tf,
#       Binding = mat_bind_tf[gene, ],
#       TF = exp_bind_tf %in% exp_bind_ct
#     ),
#     Perturbation = data.frame(
#       Experiment = exp_perturb_tf,
#       FC = mat_fc_tf[gene, ],
#       Abs_FC = abs(mat_fc_tf[gene, ]),
#       DE = mat_fdr_tf[gene, ] < fdr,
#       TF = exp_perturb_tf %in% exp_perturb_ct
#     )
#   )
# }


gene_boxplot <- function(df, yvar, yname, title, de = FALSE) {
  
  df <- filter(df, !is.na(!!sym(yvar)))
  
  p <- 
    ggplot(df, aes(x = TF, y = !!sym(yvar))) +
    geom_boxplot(width = 0.2, fill = "slategrey", outlier.shape = NA)
  
  if (de) {
    p <- p + 
      geom_jitter(aes(colour = Species, shape = Perturbation), 
                  size = 3, width = 0.2, alpha = 0.5)
    
  } else {
    
    p <- p + 
      geom_jitter(aes(colour = Species),
                  size = 3, shape = 21, width = 0.3, alpha = 0.5)
  }
  
  p <- p +
    ylab(yname) +
    ggtitle(title) +
    scale_colour_manual(values = c("black", "red")) +
    theme_classic() +
    theme(axis.title = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          plot.title = element_text(size = 30),
          # legend.position = "top",
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))
  
  return(p)
}


pxa <- gene_boxplot(bind_df, yvar = "Binding_score", yname = "Binding score", title = paste0(tf, " -> ", gene))
# pxb <- gene_boxplot(perturb_df, yvar = "Abs_FC", yname = "Absolute FC", title = paste0(tf, " -> ", gene), de = TRUE) + theme(plot.title = element_blank())
pxb <- gene_boxplot(perturb_df, yvar = "FC", yname = "FC", title = paste0(tf, " -> ", gene), de = TRUE)
px <- plot_grid(pxa, pxb, nrow = 2)
# px <- ggdraw(add_sub(px, tf, vpadding = grid::unit(0, "lines"), y = 6, x = 0.5, vjust = 6, size = 30))


# Heatmap of cor of cor. Pre-cluster to get order for triangular mat

hc <-  hclust(as.dist(1 - all_cor))
all_cor <- all_cor[hc$order, hc$order]
all_cor[lower.tri(all_cor)] <-  NA
diag(all_cor) <- NA
all_cor <- all_cor[1:nrow(all_cor) - 1, 2:ncol(all_cor)]


pal_length <- 100
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
