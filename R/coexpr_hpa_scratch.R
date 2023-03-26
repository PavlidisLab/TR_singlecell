#
# TODO: random forest/xboost
# https://towardsdatascience.com/feature-selection-using-random-forest-26d7b747597f
# https://corysimon.github.io/articles/feature-importance-in-random-forests-when-features-are-correlated/
# https://medium.com/@raj5287/effects-of-multi-collinearity-in-logistic-regression-svm-rf-af6766d91f1b
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-307


library(tidyverse)
library(reshape2)
library(WGCNA)
library(parallel)
library(glmnet)
library(pheatmap)
library(cowplot)


#
tfs_hg <- read.delim("~/Data/Metadata/human_tfs_lambert2018.tsv", stringsAsFactors = FALSE)
pc_hg <- read.delim("~/Data/Metadata/ensembl_human_protein_coding_105.tsv", stringsAsFactors = FALSE)

#
hpa <- read.delim("~/Data/Expression_files/HPA/rna_tissue_hpa.tsv", stringsAsFactors = FALSE)
gtex <- read.delim("~/Data/Expression_files/HPA/rna_tissue_gtex.tsv", stringsAsFactors = FALSE)
fantom <- read.delim("~/Data/Expression_files/HPA/rna_tissue_fantom.tsv", stringsAsFactors = FALSE)
allen <- read.delim("~/Data/Expression_files/HPA/rna_mouse_brain_allen.tsv", stringsAsFactors = FALSE)
sc <- read.delim("~/Data/Expression_files/HPA/rna_single_cell_type.tsv", stringsAsFactors = FALSE)
sc_tissue <- read.delim("~/Data/Expression_files/HPA/rna_single_cell_type_tissue.tsv",stringsAsFactors = FALSE)
pigbrain <- read.delim("~/Data/Expression_files/HPA/rna_pig_brain_hpa.tsv", stringsAsFactors = FALSE)
mousebrain <- read.delim("~/Data/Expression_files/HPA/rna_mouse_brain_hpa.tsv", stringsAsFactors = TRUE)
blood_hpa <- read.delim("~/Data/Expression_files/HPA/rna_blood_cell.tsv", stringsAsFactors = TRUE)
blood_monaco <- read.delim("~/Data/Expression_files/HPA/rna_blood_cell_monaco.tsv", stringsAsFactors = TRUE)
blood_schmiedel <- read.delim("~/Data/Expression_files/HPA/rna_blood_cell_schmiedel.tsv", stringsAsFactors = TRUE)
cline <- read.delim("~/Data/Expression_files/HPA/rna_celline.tsv", stringsAsFactors = TRUE)


# TCGA cancer processing - collapse sample grouped by cancer type by mean gene.
# Note that FPKM - others were some form of TPM. Also ensmebl IDs so must get symbol
cancer <- read.delim("~/Data/Expression_files/HPA/rna_cancer_sample.tsv", stringsAsFactors = FALSE)
cancer <- aggregate(cancer$FPKM, list(cancer$Gene, cancer$Cancer), FUN = mean) 
colnames(cancer) <- c("Gene.name", "Cancer", "FPKM")
cancer$Gene.name <- pc_hg$Symbol[match(cancer$Gene.name, pc_hg$Gene_ID)]
cancer <- filter(cancer, Gene.name != "")

#
# pigbrain_region <- read.delim("~/Data/Expression_files/HPA/rna_pig_brain_sample_hpa.tsv", stringsAsFactors = FALSE)
# mousebrain_region <- read.delim("~/Data/Expression_files/HPA/rna_mouse_brain_sample_hpa.tsv", stringsAsFactors = TRUE)
# blood_hpa_sample <- read.delim("~/Data/Expression_files/HPA/rna_blood_cell_sample.tsv", stringsAsFactors = TRUE)
# gtex_brain <- read.delim("~/Data/Expression_files/HPA/rna_brain_gtex.tsv", stringsAsFactors = FALSE)
# fantom_brain <- read.delim("~/Data/Expression_files/HPA/rna_brain_fantom.tsv", stringsAsFactors = FALSE)


#
# prot <- read.delim("~/Data/Expression_files/HPA/normal_tissue.tsv", stringsAsFactors = FALSE)

#
allrank <- readRDS("~/scratch/R_objects/Apr2022_ranked_target_list.RDS")
diff_bind_hg <- readRDS("~/scratch/R_objects/unibind_bindscore_robust_human.RDS")
diff_bind_mm <- readRDS("~/scratch/R_objects/unibind_bindscore_robust_mouse.RDS")


#
lt <- read.delim("~/Data/Metadata/Curated_targets_all_Oct2022.tsv")

lt <- mutate(lt,
             TF_Symbol = str_to_upper(TF_Symbol),
             Target_Symbol = str_to_upper(Target_Symbol))



# Process df before casting to matrix: Only keep genes present in all tissues,
# remove genes with all 0 counts, and log transform + 1
# TODO: replace use of eval (ill advised)


process_mat <- function(expr_df, var = "nTPM", group = "Tissue") {

  dup <- split(expr_df, expr_df$Gene.name)
  dup <- unlist(lapply(dup, nrow))
  
  # hacky means of getting the most represented count
  dom_count <- as.integer(names(sort(table(dup), decreasing = TRUE)[1]))
  
  # Remove duplicated gene symbols
  dup <- dup[dup != dom_count]
  expr_df <- filter(expr_df, !(Gene.name %in% names(dup)))
  
  # To matrix and remove all 0s
  mat <- acast(expr_df, 
               Gene.name ~ eval(as.name(group)),
               # Gene.name ~ substitute(group, list(group = group)), 
               value.var = var, 
               drop = FALSE)
  
  all0 <- which(rowSums(mat) == 0)
  mat <- mat[-all0, ]
  mat <- log2(mat+1)
  return(mat)
  
}


# single collapse for unique ID
sc_tissue <- sc_tissue %>% 
  mutate(ID = paste(Tissue, Cluster, Cell.type, sep = "_"))


mat_list <- list(
  HPA = process_mat(hpa),
  GTEX = process_mat(gtex),
  # GTEX_brain = process_mat(gtex_brain, group = "Brain.region"),
  FANTOM = process_mat(fantom, var = "Normalized.tags.per.million"),
  # FANTOM_brain = process_mat(fantom_brain, group = "Brain.region"),
  Allen = process_mat(allen, var = "Expression.energy", group = "Brain.region"),
  Scell = process_mat(sc, group = "Cell.type"),
  Scell_tissue = process_mat(sc_tissue, var = "pTPM", group = "ID"),
  Pigbrain = process_mat(pigbrain, group = "Brain.region"),
  Mousebrain = process_mat(mousebrain, group = "Brain.region"),
  Blood_HPA = process_mat(blood_hpa, group = "Blood.cell"),
  Blood_Monaco = process_mat(blood_hpa, var = "pTPM", group = "Blood.cell"),
  Blood_Schmiedel = process_mat(blood_schmiedel, var = "TPM", group = "Blood.cell"),
  Cell_line = process_mat(cline, group = "Cell.line"),
  Cancer = process_mat(cancer, var = "FPKM", group = "Cancer")
)



# Correlation matrices across all data sets


cor_list <- mclapply(mat_list, function(x) {
  WGCNA::cor(t(x), use = "pairwise.complete.obs")
}, mc.cores = 8)


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



tf <- "PAX6"
tf_cor <- cor_df(cor_list, tf)

tf_cor_rank <- lapply(2:ncol(tf_cor), function(x) {
  mat <- rank(-tf_cor[, x])
  return(mat)
})

tf_cor_rank <- do.call(cbind, tf_cor_rank)
rownames(tf_cor_rank) <- tf_cor$Symbol
colnames(tf_cor_rank) <- colnames(tf_cor)[2:ncol(tf_cor)]
order_cor <- sort(rowMeans(tf_cor_rank))
head(order_cor)


gene <- "CRB1"


p_list <- lapply(names(mat_list), function(x) {
  
  mat <- mat_list[[x]]
  
  if(!all(c(tf, gene) %in% rownames(mat))) {
    return(NA)
  }
  
  df <- data.frame(tf = mat[tf,], gene = mat[gene, ])
  ggplot(df, aes(x = tf, y = gene)) + 
    geom_point(size = 2, shape = 21, fill = "#756bb1") +
    geom_smooth(method = "lm", colour = "black") +
    xlab(tf) +
    ylab(gene) +
    ggtitle(paste0(x, ": ", round(cor_list[[x]][tf, gene], 3))) +
    theme_classic() +
    theme(plot.title = element_text(size = 20))
})


p_list <- p_list[!is.na(p_list)]

cowplot::plot_grid(plotlist = p_list)


# cor of cor
all_cor <- cor(dplyr::select_if(tf_cor, is.numeric), use = "pairwise.complete.obs")
hc <-  hclust(as.dist(1 - all_cor))
all_cor <- all_cor[hc$order, hc$order]
all_cor[lower.tri(all_cor)] <-  NA
diag(all_cor) <- NA
all_cor <- all_cor[1:nrow(all_cor)-1, 2:ncol(all_cor)]


pal_length <- 100
bluered_pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(pal_length)
color_breaks <- seq(min(all_cor, na.rm = TRUE), max(all_cor, na.rm = TRUE), length.out = pal_length)


pheatmap(all_cor, 
         cluster_col = FALSE, 
         cluster_row = FALSE,
         color = bluered_pal,
         breaks = color_breaks,
         na_col = "white",
         border_color = NA,
         fontsize = 18)


#



group_cols <- c(
  "FALSE" = "grey",
  "Binding" = "#e41a1c",
  "Integrated" = "#377eb8",
  "Perturbation" = "#4daf4a"
)


tf <- "MEF2C"
tf_cor <- cor_df(cor_list, tf)
topn <- 100
top_integrated <- arrange(allrank$Human[[tf]], Rank_integrated)$Symbol[1:topn]
top_chip <- arrange(allrank$Human[[tf]], Rank_binding)$Symbol[1:topn]
top_perturb <- arrange(allrank$Human[[tf]], Rank_perturbation)$Symbol[1:topn]
# sample_symbols <- sample(tf_cor$Symbol, topn)


# Split each group into a panel


panel_l <- lapply(names(mat_list), function(x) {
  
  if (! tf %in% rownames(mat_list[[x]])) {
    return(NA)
  }
  
  df <- tf_cor[, c("Symbol", x)] %>%
    dplyr::rename(Cor = all_of(x)) %>%
    mutate(
      Integrated = factor(
        ifelse(Symbol %in% top_integrated, "Integrated", FALSE),
        levels = c("Integrated", FALSE)
      ),
      Binding = factor(
        ifelse(Symbol %in% top_chip, "Binding", FALSE),
        levels = c("Binding", FALSE)
      ),
      Perturbation = factor(
        ifelse(Symbol %in% top_perturb, "Perturbation", FALSE),
        levels = c("Perturbation", FALSE)
      )
    )
  
  
  p1 <- ggplot(df) +
    geom_density(aes(x = Cor, colour = Integrated), size = 1.3) +
    geom_boxplot(aes(x = Cor, fill = Integrated), width = 0.05) +
    scale_colour_manual(values = c(group_cols["FALSE"], group_cols["Integrated"])) +
    scale_fill_manual(values = c(group_cols["FALSE"], group_cols["Integrated"])) +
    ylab("Density") +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20))
  
  p2 <- ggplot(df) +
    geom_density(aes(x = Cor, colour = Binding), size = 1.3) +
    geom_boxplot(aes(x = Cor, fill = Binding), width = 0.05) +
    scale_colour_manual(values = c(group_cols["FALSE"], group_cols["Binding"])) +
    scale_fill_manual(values = c(group_cols["FALSE"], group_cols["Binding"])) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20))
  
  p3 <- ggplot(df) +
    geom_density(aes(x = Cor, colour = Perturbation), size = 1.3) +
    geom_boxplot(aes(x = Cor, fill = Perturbation), width = 0.05) +
    scale_colour_manual(values = c(group_cols["FALSE"], group_cols["Perturbation"])) +
    scale_fill_manual(values = c(group_cols["FALSE"], group_cols["Perturbation"])) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20))
  
  plot_grid(p1, p2, p3, nrow = 1)
  
  
})
names(panel_l) <- names(mat_list)
panel_l <- panel_l[!is.na(panel_l)]



# One group var per plot


one_l <- lapply(names(mat_list), function(x) {
  
  if (! tf %in% rownames(mat_list[[x]])) {
    return(NA)
  }
  
  df <- tf_cor[, c("Symbol", x)] %>%
    dplyr::rename(Cor = all_of(x)) %>%
    mutate(
      Integrated = factor(
        ifelse(Symbol %in% top_integrated, "Integrated", FALSE),
        levels = c("Integrated", FALSE)
      ),
      Binding = factor(
        ifelse(Symbol %in% top_chip, "Binding", FALSE),
        levels = c("Binding", FALSE)
      ),
      Perturbation = factor(
        ifelse(Symbol %in% top_perturb, "Perturbation", FALSE),
        levels = c("Perturbation", FALSE)
      )
    )
  
  
  ggplot(df) +
    geom_density(aes(x = Cor, colour = Integrated), size = 1.3) +
    geom_density(aes(x = Cor, colour = Binding), size = 1.3) +
    geom_density(aes(x = Cor, colour = Perturbation), size = 1.3) +
    scale_color_manual(values = group_cols, name = paste0("Top ", topn)) +
    ylab("Density") +
    ggtitle(x) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          legend.position = "none")
  
  
})
names(one_l) <- names(mat_list)
one_l <- one_l[!is.na(one_l)]
plot_grid(plotlist = one_l)

ggsave(plot_grid(plotlist = one_l),
       dpi = 300, device = "png", height = 20, width = 20,
       filename = paste0("~/Plots/HPA/genomic_evidence_cor_", tf, "_", topn, ".png"))



# Pcor by curation status. Require measurement of TF in at least one RNA source
# and minimum number of curated targets



tf_na <- lapply(unique(lt$TF_Symbol), function(tf) {
  all(!unlist(lapply(mat_list, function(x) tf %in% rownames(x))))
})

names(tf_na) <- unique(lt$TF_Symbol)
tf_na <- unlist(tf_na)
tf_na <- names(tf_na[tf_na])


tfs_curated <- lt %>% 
  filter(! TF_Symbol %in% tf_na) %>% 
  group_by(TF_Symbol) %>% 
  summarize(N = n_distinct(Target_Symbol)) %>% 
  filter(N >= 3)




lt_pcor <- mclapply(tfs_curated$TF_Symbol, function(tf) {
  
  tf_cor <- cor_df(cor_list, tf)
  lt_tf <- filter(lt, TF_Symbol == tf)
  tf_cor$Curated <- tf_cor$Symbol %in% lt_tf$Target_Symbol
  
  l <- lapply(names(mat_list), function(x) {
  
    if (! tf %in% rownames(mat_list[[x]])) {
      return(NA)
    }
  
    df <- tf_cor[, c("Curated", x)] %>%
      dplyr::rename(Cor = all_of(x)) 
    
    tryCatch({
      wilcox.test(df$Cor ~ df$Curated)$p.value
    }, error = function(e) NA)
    
  })
  
  names(l) <- names(mat_list)
  data.frame(t(unlist(l)))
  
}, mc.cores = 8)


names(lt_pcor) <- tfs_curated$TF_Symbol
lt_pcor <- lt_pcor[!is.na(lt_pcor)]
lt_pcor_df <- do.call(rbind, lt_pcor)


summ_lt <- data.frame(
  N_sig = rowSums(lt_pcor_df < 0.05, na.rm = TRUE),
  N_NA = rowSums(is.na(lt_pcor_df)),
  Min_pval = apply(lt_pcor_df, 1, function(x) min(x, na.rm = TRUE)),
  N_target = filter(tfs_curated, TF_Symbol %in% rownames(lt_pcor_df))$N
)


cor(summ_lt$N_sig, summ_lt$N_target, method = "spearman")


ggplot(summ_lt, aes(x = log10(N_target+1), y = N_sig)) +
  geom_jitter(size = 4, shape = 21, colour = "black", fill = "slategrey") +
  xlab("Log10 count of curated targets") +
  ylab("Count of significant tests (RNA correlation lists)") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))


# inspect extremes


tf <- "MYC"
tf_cor <- cor_df(cor_list, tf)
lt_tf <- filter(lt, TF_Symbol == tf)
tf_cor$Curated <- tf_cor$Symbol %in% lt_tf$Target_Symbol

l <- lapply(names(mat_list), function(x) {
  
  if (! tf %in% rownames(mat_list[[x]])) {
    return(NA)
  }
  
  df <- tf_cor[, c("Symbol", "Curated", x)] %>%
    dplyr::rename(Cor = all_of(x))

  ggplot(df, aes(x = Cor, fill = Curated)) +
    geom_density(alpha = 0.4) +
    ggtitle(x) +
    scale_fill_manual(values = c("grey", "goldenrod")) +
    theme_classic() +
    theme(legend.position = "none")
})
l <- l[!is.na(l)]

plot_grid(plotlist = l, ncol = 3)


# Coexpr ~ binding


tfs_bind <- names(diff_bind_hg$Group_mean)
tf <- "ASCL1"
tf_cor <- cor_df(cor_list, tf)
fit_df <- diff_bind_hg$Fit[[tf]] %>% rownames_to_column(var = "Symbol")

bind_df <- diff_bind_hg$Group_mean[[tf]] %>% 
  dplyr::rename(Mean_bind = Mean) %>% 
  left_join(tf_cor, by = "Symbol") %>% 
  left_join(fit_df[, c("Symbol", "adj.P.Val")], by = "Symbol") %>% 
  mutate(Sig = adj.P.Val < 0.05)

cor(dplyr::select_if(bind_df, is.numeric), use = "pairwise.complete.obs")
plot(bind_df$Mean_bind, bind_df$HPA)
boxplot(bind_df$HPA ~ bind_df$Sig)

bind_df %>% 
  filter(!is.na(HPA) & !is.na(Sig)) %>% 
  ggplot(aes(x = HPA, colour = Sig)) +
  geom_density()



# regularized regression: must first scale matrix






gene <- "DLL3"
y <- mat_hpa[, gene]
x <- mat_hpa[, colnames(mat_hpa) != gene]


ridge_model <- cv.glmnet(x, y, alpha = 0, trace.it = 1, standardize = FALSE, standardize.response = FALSE)
lasso_model <- cv.glmnet(x, y, alpha = 1, trace.it = 1, standardize = FALSE, standardize.response = FALSE)


# plot(ridge_model)
# plot(lasso_model)


ridge_coef <- coef(ridge_model, s = "lambda.min")
lasso_coef <- coef(lasso_model, s = "lambda.min")

ridge_df <- data.frame(Symbol = ridge_coef@Dimnames[[1]],
                       Ridge_coef = ridge_coef@x) %>%
  slice(-1) %>%
  mutate(Ridge_rank = rank(-Ridge_coef))


lasso_df <- data.frame(Symbol = lasso_coef@Dimnames[[1]][lasso_coef@i+1],
                       Lasso_coef = lasso_coef@x) %>% 
  slice(-1) %>% 
  mutate(Lasso_rank = rank(-Lasso_coef))


coef_df <- left_join(ridge_df, lasso_df, by = "Symbol")


test_mat <- mat_gtex[, colnames(x)]

pred <- predict(lasso_model, 
                newx = test_mat, 
                s = "lambda.min")

stopifnot(identical(rownames(pred), rownames(mat_gtex)))
plot(mat_gtex[, gene], pred)
cor(mat_gtex[, gene], pred)


doMC::registerDoMC(cores = 8)



model_predict <- function(gene, train_mat, test_mat) {
  
  #
  y <- train_mat[, gene]
  x <- train_mat[, colnames(train_mat) != gene]
  
  #
  ridge_model <-
    cv.glmnet(
      x,
      y,
      alpha = 0,
      trace.it = 1,
      standardize = FALSE,
      standardize.response = FALSE,
      parallel = TRUE
    )
  
  lasso_model <-
    cv.glmnet(
      x,
      y,
      alpha = 1,
      trace.it = 1,
      standardize = FALSE,
      standardize.response = FALSE,
      parallel = TRUE
    )
  
  #
  ridge_coef <- coef(ridge_model, s = "lambda.min")
  lasso_coef <- coef(lasso_model, s = "lambda.min")
  
  #
  ridge_df <- data.frame(Symbol = ridge_coef@Dimnames[[1]],
                         Ridge_coef = ridge_coef@x) %>%
    slice(-1) %>%
    mutate(Ridge_rank = rank(-Ridge_coef))
  
  
  lasso_df <- data.frame(Symbol = lasso_coef@Dimnames[[1]][lasso_coef@i+1],
                         Lasso_coef = lasso_coef@x) %>% 
    slice(-1) %>% 
    mutate(Lasso_rank = rank(-Lasso_coef))
  
  
  coef_df <- left_join(ridge_df, lasso_df, by = "Symbol")
  
  
  
  #
  test_gene <- test_mat[, gene]
  test_mat <- test_mat[, colnames(x)]
  
  
  pred_ridge <- predict(ridge_model, 
                        newx = test_mat, 
                        s = "lambda.min")
  
  pred_lasso <- predict(lasso_model, 
                        newx = test_mat, 
                        s = "lambda.min")
  
  
  stopifnot(identical(rownames(pred_ridge), rownames(test_mat)))
  cor_ridge <- cor(test_gene, pred_ridge)
  cor_lasso <- cor(test_gene, pred_lasso)
  
  list(
    Coef_df = coef_df,
    Pred_ridge = pred_ridge,
    Pred_lasso = pred_lasso,
    Cor_ridge = cor_ridge,
    Cor_lasso = cor_lasso
  )
  
}


train_mat <- scale(t(mat_list$HPA), center = TRUE, scale = TRUE)
test_mat <- scale(t(mat_list$GTEX), center = TRUE, scale = TRUE)
train_mat <- train_mat[, intersect(colnames(train_mat), colnames(test_mat))]
test_mat <- test_mat[, colnames(train_mat)]



# genomics
tf <- "ASCL1"
genes <- allrank$Human[[tf]]$Symbol[1:topn]

pred_list <- lapply(genes, function(x) {
  try(model_predict(x, train_mat, test_mat))
})
names(pred_list) <- genes



pred_df <- data.frame(
  Symbol = genes,
  Ridge_rank = unlist(lapply(pred_list, function(x) filter(x$Coef_df, Symbol == tf)$Ridge_rank)),
  Lasso_rank = unlist(lapply(pred_list, function(x) filter(x$Coef_df, Symbol == tf)$Lasso_rank)),
  Ridge_cor = unlist(lapply(pred_list, function(x) x$Cor_ridge)),
  Lasso_cor = unlist(lapply(pred_list, function(x) x$Cor_lasso)),
  TFcor_HPA = cor_list$HPA[tf, genes],
  TFcor_GTEX = cor_list$GTEX[tf, genes]
)


data.frame(Observed = test_mat[, "DLL3"],
           Predicted = as.numeric(pred_list$DLL3$Pred_lasso)) %>%
  ggplot(aes(x = Predicted, y = Observed)) +
  geom_point(size = 3, shape = 21, fill = "#a6cee3") +
  geom_smooth(method = "lm", col = "red") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))



# low-throughput
tf <- "ASCL1"
genes <- unique(filter(lt, TF_Symbol == tf)$Target_Symbol)

pred_list <- lapply(genes, function(x) {
  try(model_predict(x, train_mat, test_mat))
})
names(pred_list) <- genes



pred_df <- data.frame(
  Symbol = genes,
  Ridge_rank = unlist(lapply(pred_list, function(x) filter(x$Coef_df, Symbol == tf)$Ridge_rank)),
  Lasso_rank = unlist(lapply(pred_list, function(x) filter(x$Coef_df, Symbol == tf)$Lasso_rank)),
  Ridge_cor = unlist(lapply(pred_list, function(x) x$Cor_ridge)),
  Lasso_cor = unlist(lapply(pred_list, function(x) x$Cor_lasso)),
  TFcor_HPA = cor_list$HPA[tf, genes],
  TFcor_GTEX = cor_list$GTEX[tf, genes]
)


data.frame(Observed = test_mat[, "DLL3"],
           Predicted = as.numeric(pred_list$DLL3$Pred_lasso)) %>%
  ggplot(aes(x = Predicted, y = Observed)) +
  geom_point(size = 3, shape = 21, fill = "#a6cee3") +
  geom_smooth(method = "lm", col = "red") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))
