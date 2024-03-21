## Organize and describe the ability of the aggregate coexpression and binding 
## rankings to recover literature curated targets
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# The minimum count of curated targets for a TF for reporting
min_targets <- 5

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Measurement matrices used for filtering when a gene was never expressed
msr_hg <- readRDS(msr_mat_hg_path)
msr_mm <- readRDS(msr_mat_mm_path)

# Saved list RDS of the integrated rankings
rank_tf_hg <- readRDS(rank_int_hg_path)
rank_tf_mm <- readRDS(rank_int_mm_path)

# Curated low throughput targets
curated <- read.delim(curated_all_path, stringsAsFactors = FALSE)

# Aggregate coexpr profiles versus individual dataset profiles
avi_auc_hg <- readRDS(avg_vs_ind_auc_hg_path)
avi_auc_mm <- readRDS(avg_vs_ind_auc_mm_path)

# Aggregate coexpr profiles versus null 
coexpr_auc_hg <- readRDS(coexpr_auc_hg_path)
coexpr_auc_mm <- readRDS(coexpr_auc_mm_path)

# Aggregate reversed coexpr (prioritize negative cor) profiles versus null 
rev_coexpr_auc_hg <- readRDS(rev_coexpr_auc_hg_path)
rev_coexpr_auc_mm <- readRDS(rev_coexpr_auc_mm_path)

# Aggregate binding profiles versus null
bind_auc_hg <- readRDS(unibind_auc_hg_path)
bind_auc_mm <- readRDS(unibind_auc_mm_path)


# Functions
# ------------------------------------------------------------------------------


# Give the proportion of data within the supplied numeric columns of df that are
# in the extreme ranges

extreme_proportions <- function(df, stat_cols, upper = 0.9, lower = 0.1) {
  
  res <- lapply(df[, stat_cols], function(col) {
    
    col <- col[!is.na(col)]
    n <- length(col)
    
    stats <- c(
      N = n,
      Eq1 = sum(col == 1) / n,
      Gt_upper = sum(col > upper) / n,
      Lt_lower = sum(col < lower) / n
    )
    
    stats <- round(stats, 3)
  })
  
  df <- do.call(rbind, res)
  row.names(df) <- stat_cols
  colnames(df) <- c("N", "Eq1", paste0("Gt", upper), paste0("Lt", lower))
  
  return(df)
}



# Join the summary dataframes from the list of coexpression and bind AUCs

join_auc_df <- function(coexpr_l, bind_l, min_targets) {
  
  coexpr_df <- do.call(rbind, lapply(coexpr_l, `[[`, "Perf_df"))
  bind_df <- do.call(rbind, lapply(bind_l, `[[`, "Perf_df"))
  
  left_join(coexpr_df, bind_df,
            by = c("Symbol", "N_targets"),
            suffix = c("_coexpr", "_bind")) %>% 
    filter(N_targets >= min_targets)
}



# The AUC quantile of the averaged/final network relative to the distribution
# of the individual TFs
# ------------------------------------------------------------------------------


avi_df_hg <-  lapply(avi_auc_hg, `[[`, "Summary_df") %>% 
  do.call(rbind, .) %>%
  rownames_to_column(var = "Symbol") %>% 
  filter(N_targets >= min_targets)


avi_df_mm <-  lapply(avi_auc_mm, `[[`, "Summary_df") %>% 
  do.call(rbind, .) %>%
  rownames_to_column(var = "Symbol") %>% 
  filter(N_targets >= min_targets)


# Summaries

avi_cols <- c("AUROC_quantile", "AUPRC_quantile")

avi_summ_hg <- summary(Filter(is.numeric, avi_df_hg))
avi_summ_mm <- summary(Filter(is.numeric, avi_df_mm))

ext_avi_hg <- extreme_proportions(avi_df_hg, avi_cols)
ext_avi_mm <- extreme_proportions(avi_df_mm, avi_cols)


# Inspect TFs where the average was outperformed by most individual

avi_underperf_hg <- filter(avi_df_hg, AUROC_quantile < 0.1 | AUPRC_quantile < 0.1)
avi_underperf_mm <- filter(avi_df_mm, AUROC_quantile < 0.1 | AUPRC_quantile < 0.1)

avi_underperf_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% avi_underperf_hg$Symbol & Symbol_mm %in% avi_underperf_mm$Symbol) %>% 
  pull(Symbol_hg)



# Join and summarize the aggregated coexpression and binding AUCs
# ------------------------------------------------------------------------------


agg_df_hg <- join_auc_df(coexpr_l = coexpr_auc_hg, 
                         bind_l = bind_auc_hg, 
                         min_targets = min_targets)


agg_df_mm <- join_auc_df(coexpr_l = coexpr_auc_mm, 
                         bind_l = bind_auc_mm,
                         min_targets = min_targets)


agg_df_hg$N_coexpr_msr <- rowSums(msr_hg[agg_df_hg$Symbol, ])
agg_df_mm$N_coexpr_msr <- rowSums(msr_mm[agg_df_mm$Symbol, ])


# Subset of TFs with data for both binding and coexpr along with min targets

agg_df_common_hg <- filter(agg_df_hg, !is.na(AUPRC_coexpr) & !is.na(AUPRC_bind))
agg_df_common_mm <- filter(agg_df_mm, !is.na(AUPRC_coexpr) & !is.na(AUPRC_bind))


# Summaries

agg_summ_hg <- summary(Filter(is.numeric, agg_df_hg))
agg_summ_mm <- summary(Filter(is.numeric, agg_df_mm))

agg_summ_common_hg <- summary(Filter(is.numeric, agg_df_common_hg))
agg_summ_common_mm <- summary(Filter(is.numeric, agg_df_common_mm))


agg_cols <- c(
  "AUPRC_quantile_coexpr",
  "AUPRC_quantile_bind",
  "AUROC_quantile_coexpr",
  "AUROC_quantile_bind"
)


ext_agg_hg <- extreme_proportions(agg_df_hg, agg_cols)
ext_agg_mm <- extreme_proportions(agg_df_mm, agg_cols)

ext_agg_common_hg <- extreme_proportions(agg_df_common_hg, agg_cols)
ext_agg_common_mm <- extreme_proportions(agg_df_common_mm, agg_cols)



# Difference between aggregate coexpression AUCs to binding AUCs
# ------------------------------------------------------------------------------


# More positive values mean that coexpression AUC was more performant
delta_agg_auc <- list(
  AUPRC_human = agg_df_common_hg$AUPRC_coexpr - agg_df_common_hg$AUPRC_bind,
  AUPRC_mouse = agg_df_common_mm$AUPRC_coexpr - agg_df_common_mm$AUPRC_bind,
  AUROC_human = agg_df_common_hg$AUROC_coexpr - agg_df_common_hg$AUROC_bind,
  AUROC_mouse = agg_df_common_mm$AUROC_coexpr - agg_df_common_mm$AUROC_bind
)



# More positive values mean that coexpression better separated from its null than binding
diff_delta_agg_auc <- list(
  AUPRC_human = agg_df_common_hg$AUPRC_diff_coexpr - agg_df_common_hg$AUPRC_diff_bind,
  AUPRC_mouse = agg_df_common_mm$AUPRC_diff_coexpr - agg_df_common_mm$AUPRC_diff_bind,
  AUROC_human = agg_df_common_hg$AUROC_diff_coexpr - agg_df_common_hg$AUROC_diff_bind,
  AUROC_mouse = agg_df_common_mm$AUROC_diff_coexpr - agg_df_common_mm$AUROC_diff_bind
)


delta_summ <- lapply(delta_agg_auc, summary)
diff_delta_auc_summ <- lapply(diff_delta_agg_auc, summary)



# Isolating TFs that performed well for for both aggregations
# ------------------------------------------------------------------------------



filter_both_performant <- function(df, type1, type2, upper = 0.9) {
  
  type1_auprc <- !!sym(paste0("AUPRC_quantile_", type1))
  type2_auprc <- !!sym(paste0("AUPRC_quantile_", type2))
  type1_auroc <- !!sym(paste0("AUROC_quantile_", type1))
  type2_auroc <- !!sym(paste0("AUROC_quantile_", type2))
  
  filter(df, 
        (type1_auprc > upper & type2_auprc > upper) |
        (type1_auroc > upper & type2_auroc > upper))
  
}


both_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_quantile_coexpr > 0.9 & AUPRC_quantile_bind > 0.9) |
    (AUROC_quantile_coexpr > 0.9 & AUROC_quantile_bind > 0.9)
)


both_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_quantile_coexpr > 0.9 & AUPRC_quantile_bind > 0.9) |
    (AUROC_quantile_coexpr > 0.9 & AUROC_quantile_bind > 0.9)
)


both_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% both_hg$Symbol & Symbol_mm %in% both_mm$Symbol) %>% 
  pull(Symbol_hg)


# Proportion of TFs who had good performance in both data type aggregations
both_prop_hg <- length(both_ortho) / nrow(both_hg)
both_prop_mm <- length(both_ortho) / nrow(both_mm)



# Looking at TFs that performed well in one aggregation but not the other
# ------------------------------------------------------------------------------


# Coexpression

only_coexpr_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_quantile_coexpr > 0.9 & AUPRC_quantile_bind < 0.5) |
  (AUROC_quantile_coexpr > 0.9 & AUROC_quantile_bind < 0.5)
)


only_coexpr_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_quantile_coexpr > 0.9 & AUPRC_quantile_bind < 0.5) |
    (AUROC_quantile_coexpr > 0.9 & AUROC_quantile_bind < 0.5)
)


only_coexpr_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% only_coexpr_hg$Symbol & Symbol_mm %in% only_coexpr_mm$Symbol) %>% 
  pull(Symbol_hg)



# Binding

only_binding_hg <- filter(
  agg_df_common_hg, 
  (AUPRC_quantile_bind > 0.9 & AUPRC_quantile_coexpr < 0.5) |
    (AUROC_quantile_bind > 0.9 & AUROC_quantile_coexpr < 0.5)
)


only_binding_mm <- filter(
  agg_df_common_mm, 
  (AUPRC_quantile_bind > 0.9 & AUPRC_quantile_coexpr < 0.5) |
    (AUROC_quantile_bind > 0.9 & AUROC_quantile_coexpr < 0.5)
)


only_binding_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% only_binding_hg$Symbol & Symbol_mm %in% only_binding_mm$Symbol) %>% 
  pull(Symbol_hg)



# Inspecting behavior of null AUCs
# ------------------------------------------------------------------------------


# Instances where a TFs aggregation was not performant relative to 
# the null, and yet its AUROC is greater than 0.5. This would suggest that its
# aggregate ranking is somewhat "generic" in its ability to recover targets.


check_coexpr_auroc <- agg_df_hg %>% 
  filter(AUROC_quantile_coexpr < 0.5 & AUROC_coexpr > 0.5) %>% 
  slice_max(AUROC_coexpr)


check_binding_auroc <- agg_df_hg %>% 
  filter(AUROC_quantile_bind < 0.5 & AUROC_bind > 0.5) %>% 
  slice_max(AUROC_bind)


# Summarize and order null AUROCs

summarize_null_auc <- function(auc_l, agg_df) {
  
  lapply(auc_l, function(x) summary(unlist(lapply(x$Null, `[[`, "AUROC")))) %>% 
    do.call(rbind, .) %>%
    as.data.frame() %>% 
    rownames_to_column(var = "Symbol") %>% 
    filter(Symbol %in% agg_df$Symbol) %>% 
    arrange(desc(Median))
}



null_coexpr_auroc_hg <- summarize_null_auc(coexpr_auc_hg, agg_df_hg)
null_coexpr_auroc_mm <- summarize_null_auc(coexpr_auc_mm, agg_df_mm)

null_bind_auroc_hg <- summarize_null_auc(bind_auc_hg, agg_df_hg)
null_bind_auroc_mm <- summarize_null_auc(bind_auc_mm, agg_df_mm)


# Note that the TFs with the highest null AUROCs still exceeded the null
# expectation in their respective aggregations.

top_nulls <- filter(agg_df_hg, Symbol %in% c("GATA6", "FOXC1"))



# Are there TFs where reversing ranks (neg coexpr) increased performance?
# ------------------------------------------------------------------------------


make_rev_coexpr_df <- function(rev_l, agg_df, min_targets) {
  
  lapply(rev_l, `[[`, "Perf_df") %>%
    do.call(rbind, .) %>% 
    as.data.frame() %>% 
    rename_with(~paste0(., "_reverse_coexpr"), -c("Symbol", "N_targets")) %>% 
    left_join(agg_df, by = c("Symbol", "N_targets")) %>% 
    filter(N_targets >= min_targets)
  
}


rev_coexpr_hg <- make_rev_coexpr_df(rev_coexpr_auc_hg, agg_df_hg, min_targets)
rev_coexpr_mm <- make_rev_coexpr_df(rev_coexpr_auc_mm, agg_df_mm, min_targets)


only_rev_hg <- filter(
  rev_coexpr_hg, 
  (AUPRC_quantile_reverse_coexpr > 0.9 & AUPRC_quantile_coexpr < 0.5) |
    (AUROC_quantile_reverse_coexpr > 0.9 & AUROC_quantile_coexpr < 0.5)
)


only_rev_mm <- filter(
  rev_coexpr_mm, 
  (AUPRC_quantile_reverse_coexpr > 0.9 & AUPRC_quantile_coexpr < 0.5) |
    (AUROC_quantile_reverse_coexpr > 0.9 & AUROC_quantile_coexpr < 0.5)
)


only_rev_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% only_rev_hg$Symbol & Symbol_mm %in% only_rev_mm$Symbol) %>% 
  pull(Symbol_hg)



# Elevated performance in reverse coexpression and binding: repression?

rev_binding_hg <- filter(only_rev_hg, Symbol %in% only_binding_hg$Symbol)

rev_binding_mm <- filter(only_rev_mm, Symbol %in% only_binding_mm$Symbol)

rev_binding_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% rev_binding_hg$Symbol & 
         Symbol_mm %in% rev_binding_mm$Symbol) %>% 
  pull(Symbol_hg)



# Performant in coexpression and reverse coexpression... suggests targets are
# at extremes of rankings


both_coexpr_hg <- filter(
  rev_coexpr_hg, 
  (AUPRC_quantile_reverse_coexpr > 0.9 & AUPRC_quantile_coexpr > 0.9) |
    (AUROC_quantile_reverse_coexpr > 0.9 & AUROC_quantile_coexpr > 0.9)
)


both_coexpr_mm <- filter(
  rev_coexpr_mm, 
  (AUPRC_quantile_reverse_coexpr > 0.9 & AUPRC_quantile_coexpr > 0.9) |
    (AUROC_quantile_reverse_coexpr > 0.9 & AUROC_quantile_coexpr > 0.9)
)



# Checking the ranks of curated targets for each aggregate for a given TF
# ------------------------------------------------------------------------------

# LEF1 used in paper: performant in coexpr but not ChIP-seq in both species
# filter(agg_df_common_hg, Symbol == "LEF1")
# filter(agg_df_common_mm, Symbol == "Lef1")


check_tf <- "LEF1"


targets <- get_curated_labels(tf = check_tf,
                              curated_df = curated,
                              ortho_df = pc_ortho,
                              pc_df = pc_hg,
                              species = "Human",
                              remove_self = TRUE)

rank_df <- filter(rank_tf_hg[[check_tf]], Symbol %in% targets)



# N targets at k=500 cutoff
n_coexpr <- sum(rank_df$Rank_aggr_coexpr <= 500)
n_bind <- sum(rank_df$Rank_bind <= 500)



# Plotting
# ------------------------------------------------------------------------------


# Histograms of Average versus individual AUCs

p1a <- plot_hist(avi_df_hg, stat_col = paste0("AUPRC", "_quantile"), xlab = "AUPRC quantile", title = "Human")
p1b <- plot_hist(avi_df_hg, stat_col = paste0("AUROC", "_quantile"), xlab = "AUROC quantile", title = "Human")
p1c <- plot_hist(avi_df_mm, stat_col = paste0("AUPRC", "_quantile"), xlab = "AUPRC quantile", title = "Mouse")
p1d <- plot_hist(avi_df_mm, stat_col = paste0("AUROC", "_quantile"), xlab = "AUROC quantile", title = "Mouse")

p1 <- plot_grid(p1a, p1b, p1c, p1d, nrow = 2)

ggsave(p1, height = 9, width = 12, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_average_vs_individual_hists.png")))



# Density plot of a TF's individual AUCs overlaid with the average AUC

plot_auc_density <- function(plot_df, vline, stat, tf) {
  
  xlim <- c(min(plot_df[[stat]], vline), max(plot_df[[stat]], vline))
  
  ggplot(plot_df, aes(x = !!sym(stat))) +
    geom_density(fill = "lightgrey") +
    geom_vline(xintercept = vline, col = "firebrick", linewidth = 1.4) +
    ylab("Density") +
    ggtitle(tf) +
    xlim(xlim) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          plot.margin = margin(c(10, 20, 10, 10)))
}



p2_tf <- "ASCL1"
p2_stat <- "AUROC"
pdf2 <- avi_auc_hg[[p2_tf]]$AUC_df
vline2 <- filter(pdf2, ID == "Average")[[p2_stat]]
pdf2 <- filter(pdf2, ID != "Average")

p2 <- plot_auc_density(pdf2, vline = vline2, stat = p2_stat, tf = p2_tf)

ggsave(p2, height = 4, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, p2_tf, "_coexpr_average_vs_individual_density.png")))



# Histogram of quantiles of aggregate AUCs versus null


# Using all TFs available for coexpression
plot_df3_hg <- agg_df_hg
plot_df3_mm <- agg_df_mm

# Using only TFs common to coexpression and bind
# plot_df3_hg <- agg_df_common_hg
# plot_df3_mm <- agg_df_common_mm


p3a <- plot_hist(plot_df3_hg, stat_col = paste0("AUPRC_quantile_coexpr"), xlab = "Coexpression AUPRC quantile", title = "Human")
p3b <- plot_hist(plot_df3_mm, stat_col = paste0("AUPRC_quantile_coexpr"), xlab = "Coexpression AUPRC quantile", title = "Mouse")

p3c <- plot_hist(plot_df3_hg, stat_col = paste0("AUROC_quantile_coexpr"), xlab = "Coexpression AUROC quantile", title = "Human")
p3d <- plot_hist(plot_df3_mm, stat_col = paste0("AUROC_quantile_coexpr"), xlab = "Coexpression AUROC quantile", title = "Mouse")

p3e <- plot_hist(plot_df3_hg, stat_col = paste0("AUPRC_quantile_bind"), xlab = "Binding AUPRC quantile", title = "Human")
p3f <- plot_hist(plot_df3_mm, stat_col = paste0("AUPRC_quantile_bind"), xlab = "Binding AUPRC quantile", title = "Mouse")

p3g <- plot_hist(plot_df3_hg, stat_col = paste0("AUROC_quantile_bind"), xlab = "Binding AUROC quantile", title = "Human")
p3h <- plot_hist(plot_df3_mm, stat_col = paste0("AUROC_quantile_bind"), xlab = "Binding AUROC quantile", title = "Mouse")


# TODO: fussing with fit in figure. when finalized clean up save 
# # Focus on AUPRC, showing both species and both aggregations
# p3 <- plot_grid(p3a, p3b, p3e, p3f, ncol = 1)


ggsave(p3a, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_quant_auprc_human.png")))

ggsave(p3b, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_quant_auprc_mouse.png")))

ggsave(p3c, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_quant_auroc_human.png")))

ggsave(p3d, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_quant_auroc_mouse.png")))

ggsave(p3e, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "binding_quant_auprc_human.png")))

ggsave(p3f, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "binding_quant_auprc_mouse.png")))

ggsave(p3g, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "binding_quant_auroc_human.png")))

ggsave(p3h, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "binding_quant_auroc_mouse.png")))




# Density plot of a TF's null AUCs overlaid with the aggregate AUC

p4_tf <- "ASCL1"
p4_stat <- "AUROC"
p4_list <- coexpr_auc_hg
pdf4 <- data.frame(AUC = unlist(lapply(p4_list[[p4_tf]]$Null, `[[`, p4_stat)))
colnames(pdf4) <- p4_stat
vline4 <- p4_list[[p4_tf]]$Perf_df[[p4_stat]]

p4 <- plot_auc_density(pdf4, vline = vline4, stat = p4_stat, tf = p4_tf)

ggsave(p4, height = 6, width = 9, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, p4_tf, "_coexpr_aggregate_vs_null_density.png")))


# Scatterplot of coexpr vs bind quantiles, with boxes indicating if one, both,
# or no method was performant

p5_stat <- "AUROC"

p5 <- 
  ggplot(
  agg_df_common_hg,
  aes(x = !!sym(paste0(p5_stat, "_quantile_coexpr")),
      y = !!sym(paste0(p5_stat, "_quantile_bind")))) +
  
  # Squares with no fill
  # geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = 0.9, ymax = 1.05), fill = NA, colour = "forestgreen") +
  # geom_rect(aes(xmin = -0.05, xmax = 0.1, ymin = 0.9, ymax = 1.05), fill = NA, colour = "grey") +
  # geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = -0.05, ymax = 0.1), fill = NA, colour = "grey") +
  # geom_rect(aes(xmin = -0.05, xmax = 0.5, ymin = -0.05, ymax = 0.5), fill = NA, colour = "firebrick") +

  # Squares with fill
  geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = 0.9, ymax = 1.05), fill = "forestgreen", colour = NA, alpha = 0.002) +
  geom_rect(aes(xmin = -0.05, xmax = 0.1, ymin = 0.9, ymax = 1.05), fill = "grey", colour = NA, alpha = 0.006) +
  geom_rect(aes(xmin = 0.9, xmax = 1.05, ymin = -0.05, ymax = 0.1), fill = "grey", colour = NA, alpha = 0.006) +
  geom_rect(aes(xmin = -0.05, xmax = 0.5, ymin = -0.05, ymax = 0.5), fill = "firebrick", colour = NA, alpha = 0.002) +

  geom_point(shape = 21, size = 3.4) +
  xlab(paste("Coexpression", p5_stat, "quantile")) +
  ylab(paste("Binding", p5_stat, "quantile")) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 10, 10, 10)))


ggsave(p5, height = 9, width = 9, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_vs_binding_", p5_stat, "_perc_scatter_human.png")))



# Relationship between AUC and quantile AUC
p6a <- qplot(agg_df_hg, xvar = "AUPRC_coexpr", yvar = "AUPRC_quantile_coexpr")
p6b <- qplot(agg_df_hg, xvar = "AUROC_coexpr", yvar = "AUROC_quantile_coexpr")

# Relationship between coexpression and binding AUC
p7a <- qplot(agg_df_common_hg, xvar = "AUPRC_coexpr", yvar = "AUPRC_bind")
p7b <- qplot(agg_df_common_mm, xvar = "AUPRC_coexpr", yvar = "AUPRC_bind")
p7c <- qplot(agg_df_common_hg, xvar = "AUROC_coexpr", yvar = "AUROC_bind")
p7d <- qplot(agg_df_common_mm, xvar = "AUROC_coexpr", yvar = "AUROC_bind")



# Vbplot of AUC quantiles for coexpression and binding aggregations 

auc_vbplot <- function(agg_df, title) {
  
  pdf <- data.frame(
    Quantile = c(
      agg_df$AUPRC_quantile_coexpr,
      agg_df$AUPRC_quantile_bind,
      agg_df$AUPRC_quantile_reverse_coexpr,
      agg_df$AUROC_quantile_coexpr,
      agg_df$AUROC_quantile_bind,
      agg_df$AUROC_quantile_reverse_coexpr
    ),
    Group1 = c(
      rep("(+) Coexpression", length(agg_df$AUPRC_quantile_coexpr)),
      rep("Binding", length(agg_df$AUPRC_quantile_bind)),
      rep("(-) Coexpression", length(agg_df$AUPRC_quantile_bind)),
      rep("(+) Coexpression", length(agg_df$AUROC_quantile_coexpr)),
      rep("Binding", length(agg_df$AUROC_quantile_bind)),
      rep("(-) Coexpression", length(agg_df$AUROC_ratio_reverse_coexpr))
    )
  )
  
  pdf$Group1 <- factor(pdf$Group1, levels = unique(pdf$Group1))
  pdf$Group2 <- c(rep("AUPRC", nrow(pdf) / 2), rep("AUROC", nrow(pdf) / 2))
  
  
  ggplot(pdf, aes(x = Group1, y = Quantile)) +
    facet_wrap(~Group2, scales = "free") +
    geom_violin(fill = "slategrey") +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    ggtitle(title) +
    ylab("Quantile") +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 25),
          strip.text = element_text(size = 25),
          plot.title = element_text(size = 25))
  
  
}


p8a <- auc_vbplot(rev_coexpr_hg, "Human")
p8b <- auc_vbplot(rev_coexpr_hg, "Mouse")

ggsave(p8a, height = 9, width = 18, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_vs_binding_quantile_vbplot_human.png")))

ggsave(p8b, height = 9, width = 18, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "coexpr_vs_binding_quantile_vbplot_mouse.png")))



# Plotting distribution of the median null versus observed AUROCs
# Note that this currently fixed to AUROC, which is shown in the paper.

auc_density <- function(agg_df, null_df, title, colour) {
  
  pdf <- data.frame(
    AUROC = c(agg_df$AUROC_coexpr, null_df$Median),
    Group = c(rep("Aggregate", nrow(agg_df)), 
              rep("Median null", nrow(agg_df))))
 
    ggplot(pdf, aes(x = AUROC, fill = Group)) +
    geom_density() +
    ylab("Density") +
    ggtitle(title) +
    scale_fill_manual(values = c(colour, "lightgrey")) +
    scale_x_continuous(breaks = seq(0, 1, 0.25)) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          legend.position = c(0.8, 0.8),
          legend.title = element_blank(),
          legend.text = element_text(size = 25),
          plot.margin = margin(c(10, 10, 10, 10)))
  
  
}


p9a <- auc_density(agg_df_hg, null_coexpr_auroc_hg, "Human coexpression", "royalblue")
p9b <- auc_density(agg_df_mm, null_coexpr_auroc_mm, "Mouse coexpression", "goldenrod")

p9 <- plot_grid(p9a, p9b, ncol = 1)

ggsave(p9, height = 12, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, "coexpr_obs_vs_median_null_density.png"))


# Scatter showing the range of null AUCs

p10a <-
  left_join(agg_df_hg, null_coexpr_auroc_hg, by = "Symbol") %>%
  qplot(., xvar = "Median", yvar = "AUROC_coexpr") +
  geom_hline(yintercept = 0.5, col = "red") +
  geom_vline(xintercept = 0.5, col = "red") +
  ylab("Aggregate coexpression AUROC") +
  xlab("Median null coexpression AUROC")


p10b <-
  left_join(agg_df_hg, null_bind_auroc_hg, by = "Symbol") %>%
  qplot(., xvar = "Median", yvar = "AUROC_bind") +
  geom_hline(yintercept = 0.5, col = "red") +
  geom_vline(xintercept = 0.5, col = "red") +
  ylab("Aggregate binding AUROC") +
  xlab("Median null binding AUROC")



# Histograms of the difference in raw AUCs (used in paper) between coexpression
# and binding aggregates for available TFs

p11a_l <- list(
  plot_hist(data.frame(AUPRC_human = delta_agg_auc$AUPRC_human), stat_col = "AUPRC_human", xlab = "Coexpression - Binding AUPRC", title = "Human"),
  plot_hist(data.frame(AUROC_human = delta_agg_auc$AUROC_human), stat_col = "AUROC_human", xlab = "Coexpression - Binding AUROC", title = "Human"),
  plot_hist(data.frame(AUPRC_mouse = delta_agg_auc$AUPRC_mouse), stat_col = "AUPRC_mouse", xlab = "Coexpression - Binding AUPRC", title = "Mouse"),
  plot_hist(data.frame(AUROC_mouse = delta_agg_auc$AUROC_mouse), stat_col = "AUROC_mouse", xlab = "Coexpression - Binding AUROC", title = "Mouse")
)


p11a <- plot_grid(plotlist = p11a_l, nrow = 2)

ggsave(p11a, height = 9, width = 12, device = "png", dpi = 300,
       filename = file.path(paste0(plot_dir, "delta_coexpr_bind_histograms.png")))
