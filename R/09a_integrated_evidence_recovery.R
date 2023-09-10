## TODO: point to missing genes
## TODO: missing genes excluded from rank df makes a difference?
## TODO: common process for curated and top k evidence for get_rank_df
## TODO: remove tf?
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

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# For each dataset, load the subset gene x TF aggregation matrix 
agg_tf_hg <- load_or_generate_agg(path = agg_tf_hg_path, ids = ids_hg, genes = pc_hg$Symbol, sub_genes = tfs_hg$Symbol)
agg_tf_mm <- load_or_generate_agg(path = agg_tf_mm_path, ids = ids_mm, genes = pc_mm$Symbol, sub_genes = tfs_mm$Symbol)

# Saved list RDS of the summaries
rank_tf_hg <- readRDS(rank_tf_hg_path)
rank_tf_mm <- readRDS(rank_tf_mm_path)

# TODO: Unibind and low throughput
evidence_l <- readRDS(evidence_path)



# 
# ------------------------------------------------------------------------------




get_colwise_performance <- 





auprc_hg_curated <- get_tf_performance(agg_l = agg_tf_hg,
                                       msr_mat = msr_hg,
                                       evidence_df = evidence_l$Human[[tf_hg]],
                                       evidence_col = "Curated_target",
                                       tf = tf_hg,
                                       measure = "AUPRC")



pr_df_hg_curated <- get_all_perf_df(agg_l = agg_tf_hg,
                                    msr_mat = msr_hg,
                                    evidence_df = evidence_l$Human[[tf_hg]],
                                    evidence_col = "Curated_target",
                                    tf = tf_hg,
                                    measure = "PR")





# Demo a single TF
# ------------------------------------------------------------------------------


tf_mm <- "Mef2c"
tf_hg <- "MEF2C"
evidence_col <- "Rank_integrated"
k <- 500


# agg_df <- arrange(rank_tf_mm[[tf]], Rank_RSR)
# agg_vec <- agg_df$Avg_RSR
# names(agg_vec) <- agg_df$Symbol


# Variable nrow of performance df due to ties, not sure how ROCR handles

# rank_df <- get_rank_df(gene_vec = agg_tf_hg$Velmeshev[, tf_hg],
#                        tf = tf_hg,
#                        evidence_df = evidence_l$Human[[tf_hg]],
#                        evidence_col = "Curated_target")
# 
# 
# perf_df <- get_perf_df(rank_df = rank_df,
#                        label_col = "Label",
#                        score_col = "Score",
#                        measure = "PR")



# AUPRCs for genomic evidence

auprc_hg_genomic <- get_tf_performance(agg_l = agg_tf_hg,
                                       msr_mat = msr_hg,
                                       evidence_df = evidence_l$Human[[tf_hg]],
                                       evidence_col = evidence_col,
                                       tf = tf_hg,
                                       measure = "AUPRC",
                                       k = k)


auprc_mm_genomic <- get_tf_performance(agg_l = agg_tf_mm,
                                       msr_mat = msr_mm,
                                       evidence_df = evidence_l$Mouse[[tf_mm]],
                                       evidence_col = evidence_col,
                                       tf = tf_mm,
                                       measure = "AUPRC",
                                       k = k)


pr_df_hg_genomic <- get_all_perf_df(agg_l = agg_tf_hg,
                                    msr_mat = msr_hg,
                                    evidence_df = evidence_l$Human[[tf_hg]],
                                    evidence_col = evidence_col,
                                    tf = tf_hg,
                                    measure = "PR",
                                    k = k)


pr_df_mm_genomic <- get_all_perf_df(agg_l = agg_tf_mm,
                                    msr_mat = msr_mm,
                                    evidence_df = evidence_l$Mouse[[tf_mm]],
                                    evidence_col = evidence_col,
                                    tf = tf_mm,
                                    measure = "PR",
                                    k = k)



# AUROCs for genomic evidence

auroc_hg_genomic <- get_tf_performance(agg_l = agg_tf_hg,
                                       msr_mat = msr_hg,
                                       evidence_df = evidence_l$Human[[tf_hg]],
                                       evidence_col = evidence_col,
                                       tf = tf_hg,
                                       measure = "AUROC",
                                       k = k)


auroc_mm_genomic <- get_tf_performance(agg_l = agg_tf_mm,
                                       msr_mat = msr_mm,
                                       evidence_df = evidence_l$Mouse[[tf_mm]],
                                       evidence_col = evidence_col,
                                       tf = tf_mm,
                                       measure = "AUROC",
                                       k = k)


roc_df_hg_genomic <- get_all_perf_df(agg_l = agg_tf_hg,
                                     msr_mat = msr_hg,
                                     evidence_df = evidence_l$Human[[tf_hg]],
                                     evidence_col = evidence_col,
                                     tf = tf_hg,
                                     measure = "ROC",
                                     k = k)


roc_df_mm_genomic <- get_all_perf_df(agg_l = agg_tf_mm,
                                     msr_mat = msr_mm,
                                     evidence_df = evidence_l$Mouse[[tf_mm]],
                                     evidence_col = evidence_col,
                                     tf = tf_mm,
                                     measure = "ROC",
                                     k = k)




# AUPRCs for curated evidence

auprc_hg_curated <- get_tf_performance(agg_l = agg_tf_hg,
                                       msr_mat = msr_hg,
                                       evidence_df = evidence_l$Human[[tf_hg]],
                                       evidence_col = "Curated_target",
                                       tf = tf_hg,
                                       measure = "AUPRC")


auprc_mm_curated <- get_tf_performance(agg_l = agg_tf_mm,
                                       msr_mat = msr_mm,
                                       evidence_df = evidence_l$Mouse[[tf_mm]],
                                       evidence_col = "Curated_target",
                                       tf = tf_mm,
                                       measure = "AUPRC")


pr_df_hg_curated <- get_all_perf_df(agg_l = agg_tf_hg,
                                    msr_mat = msr_hg,
                                    evidence_df = evidence_l$Human[[tf_hg]],
                                    evidence_col = "Curated_target",
                                    tf = tf_hg,
                                    measure = "PR")


pr_df_mm_curated <- get_all_perf_df(agg_l = agg_tf_mm,
                                    msr_mat = msr_mm,
                                    evidence_df = evidence_l$Mouse[[tf_mm]],
                                    evidence_col = "Curated_target",
                                    tf = tf_mm,
                                    measure = "PR")



# AUROCs for curated evidence

auroc_hg_curated <- get_tf_performance(agg_l = agg_tf_hg,
                                       msr_mat = msr_hg,
                                       evidence_df = evidence_l$Human[[tf_hg]],
                                       evidence_col = "Curated_target",
                                       tf = tf_hg,
                                       measure = "AUROC")


auroc_mm_curated <- get_tf_performance(agg_l = agg_tf_mm,
                                       msr_mat = msr_mm,
                                       evidence_df = evidence_l$Mouse[[tf_mm]],
                                       evidence_col = "Curated_target",
                                       tf = tf_mm,
                                       measure = "AUROC")


roc_df_hg_curated <- get_all_perf_df(agg_l = agg_tf_hg,
                             msr_mat = msr_hg,
                             evidence_df = evidence_l$Human[[tf_hg]],
                             evidence_col = "Curated_target",
                             tf = tf_hg,
                             measure = "ROC")


roc_df_mm_curated <- get_all_perf_df(agg_l = agg_tf_mm,
                             msr_mat = msr_mm,
                             evidence_df = evidence_l$Mouse[[tf_mm]],
                             evidence_col = "Curated_target",
                             tf = tf_mm,
                             measure = "ROC")



# Plotting

# Human genomic

p_auroc_hg_genomic <- plot_perf(df = roc_df_hg_genomic,
                                auc_l = roc_hg_genomic,
                                measure = "ROC",
                                # cols,
                                title = tf_hg,
                                ncol_legend = 1)


p_auprc_hg_genomic <- plot_perf(df = pr_df_hg_genomic,
                                auc_l = auprc_hg_genomic,
                                measure = "PR",
                                # cols,
                                title = tf_hg,
                                ncol_legend = 1)


p_hg_genomic <- plot_grid(p_auroc_hg_genomic, p_auprc_hg_genomic)



# Mouse genomic

p_auroc_mm_genomic <- plot_perf(df = roc_df_mm_genomic,
                                auc_l = roc_mm_genomic,
                                measure = "ROC",
                                # cols,
                                title = tf_mm,
                                ncol_legend = 1)


p_auprc_mm_genomic <- plot_perf(df = pr_df_mm_genomic,
                                auc_l = auprc_mm_genomic,
                                measure = "PR",
                                # cols,
                                title = tf_mm,
                                ncol_legend = 1)


p_mm_genomic <- plot_grid(p_auroc_mm_genomic, p_auprc_mm_genomic)



# Human curated

p_auroc_hg_curated <- plot_perf(df = roc_df_hg_curated,
                                auc_l = roc_hg_curated,
                                measure = "ROC",
                                # cols,
                                title = tf_hg,
                                ncol_legend = 1)


p_auprc_hg_curated <- plot_perf(df = pr_df_hg_curated,
                                auc_l = auprc_hg_curated,
                                measure = "PR",
                                # cols,
                                title = tf_hg,
                                ncol_legend = 1)


p_hg_curated <- plot_grid(p_auroc_hg_curated, p_auprc_hg_curated)



# Mouse curated

p_auroc_mm_curated <- plot_perf(df = roc_df_mm_curated,
                                auc_l = roc_mm_curated,
                                measure = "ROC",
                                # cols,
                                title = tf_mm,
                                ncol_legend = 1)


p_auprc_mm_curated <- plot_perf(df = pr_df_mm_curated,
                                auc_l = auprc_mm_curated,
                                measure = "PR",
                                # cols,
                                title = tf_mm,
                                ncol_legend = 1)


p_mm_curated <- plot_grid(p_auroc_mm_curated, p_auprc_mm_curated)



# Look at ability of genomic + coexpr to recover curated

species <- "Human"
measure <- "ROC"
tf <- tf_hg

integrated <- get_perf_df(rank_df = arrange(evidence_l[[species]][[tf]], Rank_integrated),
                          "Curated_target",
                          # score_col = "Rank_integrated",
                          measure = measure)

integrated$ID <- "Integrated"


get_au_perf(rank_df = arrange(evidence_l[[species]][[tf]], Rank_integrated),
            "Curated_target",
            # score_col = "Rank_integrated",
            measure = "AUROC")


perturb <- get_perf_df(rank_df = arrange(evidence_l[[species]][[tf]], Rank_perturbation),
                       "Curated_target",
                       # score_col = "Rank_perturbation",
                       measure = measure)

perturb$ID <- "Perturbation"


get_au_perf(rank_df = arrange(evidence_l[[species]][[tf]], Rank_perturbation),
            "Curated_target",
            # score_col = "Rank_integrated",
            measure = "AUROC")


binding <- get_perf_df(rank_df = arrange(evidence_l[[species]][[tf]], Rank_binding),
                       "Curated_target",
                       # score_col = "Rank_binding",
                       measure = measure)

binding$ID <- "Binding"

get_au_perf(rank_df = arrange(evidence_l[[species]][[tf]], Rank_binding),
            "Curated_target",
            # score_col = "Rank_integrated",
            measure = "AUROC")




sort(auroc_mm_curated)
median(auroc_mm_curated[names(auroc_mm_curated) != "Aggregate"])





all <- rbind(roc_df_mm_curated, 
            do.call(rbind, list(integrated, perturb, binding)))



all$ID <- ifelse(
  all$ID %in% c("Aggregate", "Integrated", "Binding", "Perturbation"),
  all$ID,
  "Single_coexpression"
)


cols <- c("Aggregate" = "black",
          "Binding" = "darkblue",
          "Perturbation" = "firebrick",
          "Integrated" = "darkgreen",
          "Single_coexpression" = "lightgrey")


all$ID <- factor(all$ID, levels = rev(unique(names(cols))))


plot_perf2 <- function(df, 
                      measure = NULL, 
                      cols,
                      title,
                      ncol_legend = 1) {
  
  stopifnot(measure %in% c("ROC", "PR"), "ID" %in% colnames(df))
  
  if (measure == "ROC") {
    p <- ggplot(df, aes(x = FPR, y = TPR, col = ID))
  } else {
    p <- ggplot(df, aes(x = Recall, y = Precision, col = ID))
  }
  
  p <- p +
    geom_path() +
    ggtitle(title) +
    scale_color_manual(values = cols) +
    guides(colour = guide_legend(ncol = ncol_legend)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 30),
          plot.title = element_text(size = 30),
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.position = c(0.8, 0.25))
  
  return(p)
}


p_all <- plot_perf2(df = all,
                    measure = measure,
                    cols = cols,
                    title = tf,
                    ncol_legend = 1)




# Checking examples where individual experiments start with high precision

check_p1_1 <- split(perf_df1, perf_df1$ID)
check_p1_2 <- lapply(check_p1_1, function(x) x$Precision[2] == 1)
check_p1_3 <- names(check_p1_1[unlist(check_p1)])

check_p1_4 <- lapply(check_p1_3, function(x) {
  topk_genes <- topk_sort(agg_tf_hg[[x]][, tf], k = 6)
  topk_genes <- topk_genes[topk_genes != tf]
  # topk_genes %in% evidence_df$Symbol[1:1000]
})



# Individual versus aggregate recovery
# ------------------------------------------------------------------------------


# Human

# tf_performance_hg <- lapply(tfs_hg, function(tf) {
#   
#   get_tf_performance(agg_l = agg_hg, 
#                      evidence_df = evidence_l$Human[[tf]], 
#                      evidence_col = "Rank_integrated", 
#                      tf = tf, 
#                      k = 1000)
# })
# names(tf_performance_hg) <- tfs_hg



# Mouse

# tf_performance_mm <- lapply(tfs_mm, function(tf) {
#   
#   get_tf_performance(agg_l = agg_mm, 
#                      evidence_df = evidence_l$Mouse[[tf]], 
#                      evidence_col = "Rank_integrated", 
#                      tf = tf, 
#                      k = 1000)
# })
# names(tf_performance_mm) <- tfs_mm





# Every gene recovery of evidence (slow!)
# ------------------------------------------------------------------------------


# outfile_hg1 <- "/space/scratch/amorin/R_objects/05-06-2023_coexpr_agg_recovery_integrated_evidence_human.RDS"
# outfile_hg2 <- "/space/scratch/amorin/R_objects/05-06-2023_coexpr_agg_recovery_curated_evidence_human.RDS"
# outfile_mm1 <- "/space/scratch/amorin/R_objects/05-06-2023_coexpr_agg_recovery_integrated_evidence_mouse.RDS"
# outfile_mm2 <- "/space/scratch/amorin/R_objects/05-06-2023_coexpr_agg_recovery_curated_evidence_mouse.RDS"
# 

# Human

# if (!file.exists(outfile_hg1)) {
#   
#   hg1 <- lapply(tfs_hg, function(tf) {
#     
#     get_all_performance(agg_l = agg_hg,
#                         evidence_df = evidence_l$Human[[tf]],
#                         evidence_col = "Rank_integrated",
#                         genes = genes_hg,
#                         k = 1000,
#                         ncores = ncore)
#   })
#   names(hg1) <- tfs_hg
#   
#   saveRDS(hg1, outfile_hg1)
# 
# } else {
#   
#   hg1 <- readRDS(outfile_hg1)
#   
# }


# 
# if (!file.exists(outfile_hg2)) {
#   
#   hg2 <- lapply(tfs_hg, function(tf) {
#     
#     get_all_performance(agg_l = agg_hg,
#                         evidence_df = evidence_l$Human[[tf]],
#                         evidence_col = "Curated_target",
#                         genes = genes_hg,
#                         k = 1000,
#                         ncores = ncore)
#   })
#   names(hg2) <- tfs_hg
#   
#   saveRDS(hg2, outfile_hg2)
#   
# } else {
#   
#   hg2 <- readRDS(outfile_hg2)
#   
# }



# Mouse

# 
# if (!file.exists(outfile_mm1)) {
#   
#   mm1 <- lapply(tfs_mm, function(tf) {
#     
#     get_all_performance(agg_l = agg_mm,
#                         evidence_df = evidence_l$Mouse[[tf]],
#                         evidence_col = "Rank_integrated",
#                         genes = genes_mm,
#                         k = 1000,
#                         ncores = ncore)
#   })
#   names(mm1) <- tfs_mm
#   
#   saveRDS(mm1, outfile_mm1)
#   
# } else {
#   
#   mm1 <- readRDS(outfile_mm1)
#   
# }



# if (!file.exists(outfile_mm2)) {
#   
#   mm2 <- lapply(tfs_mm, function(tf) {
#     
#     get_all_performance(agg_l = agg_mm,
#                         evidence_df = evidence_l$Mouse[[tf]],
#                         evidence_col = "Curated_target",
#                         genes = genes_mm,
#                         k = 1000,
#                         ncores = ncore)
#   })
#   names(mm2) <- tfs_mm
#   
#   saveRDS(mm2, outfile_mm2)
#   
# } else {
#   
#   mm2 <- readRDS(outfile_mm2)
#   
# }



# TODO:
# ------------------------------------------------------------------------------


# which_hg1 <- which_rank_mat(hg1)
# which_hg2 <- which_rank_mat(hg2)
# which_mm1 <- which_rank_mat(mm1)
# which_mm2 <- which_rank_mat(mm2)


# tf <- "ASCL1"
# 
# hist(hg1[[tf]]$Tabula_Sapiens$Topk, breaks = 100)
# abline(v = filter(hg1[[tf]]$Tabula_Sapiens, Symbol == tf)$Topk, col = "red")
# head(hg1[[tf]]$Tabula_Sapiens, 10)
# 
# hist(hg2[[tf]]$Tabula_Sapiens$Topk, breaks = 100)
# abline(v = filter(hg2[[tf]]$Tabula_Sapiens, Symbol == tf)$Topk, col = "red")
# head(hg2[[tf]]$Tabula_Sapiens, 10)
# 
# 
# plot(hg1$ASCL1$Tabula_Sapiens$Topk, hg1$ASCL1$Tabula_Sapiens$AUPRC)
# cor(hg1$ASCL1$Tabula_Sapiens$Topk, hg1$ASCL1$Tabula_Sapiens$AUPRC, method = "spearman")
# 
# 
# tf <- "Pax6"
# 
# hist(mm2[[tf]]$HypoMap$Topk, breaks = 100)
# abline(v = filter(mm2[[tf]]$HypoMap, Symbol == tf)$Topk, col = "red")
# head(mm2[[tf]]$HypoMap, 10)
# 
# 
# # Tally the genes whose aggregate expression vector was among the top k at 
# # recovering evidence (also top k)
# 
# tally_topk <- lapply(hg1$HES1, function(x) x$Symbol[1:k]) %>%
#   unlist() %>%
#   table() %>%
#   sort(decreasing = TRUE)




# Notice that while RUNX1 coexpression itself is not highly ranked, the genes
# that are highly ranked (by their coexpression recovry of RUNX1 targets) are
# also highly ranked RUNX1 targets...

# filter(evidence, Symbol %in% names(tally_topk)[1:100]

       