## Inspect the results of GO enrichment and export barplots of sig. terms
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/utils/plot_functions.R")


# The score/significance and GO process column names of the results
score_col <- "CorrectedMFPvalue"
term_col <- "Name"

# Loading ortho gene table and GO result lists
pc_ortho <- read.delim(pc_ortho_path)
go_coexpr_hg <- readRDS(erminer_coexpr_hg_path)
go_coexpr_mm <- readRDS(erminer_coexpr_mm_path)
go_coexpr_ortho <- readRDS(erminer_coexpr_ortho_path)



# Only consider results reaching sig. cut-off

filter_results <- function(go_l, cutoff_col, cutoff = 0.05) {
  lapply(go_l, function(x) filter(x, !!sym(cutoff_col) < cutoff))
}


go_coexpr_filt_hg <- filter_results(go_coexpr_hg, score_col)
go_coexpr_filt_mm <- filter_results(go_coexpr_mm, score_col)
go_coexpr_filt_ortho <- filter_results(go_coexpr_ortho, score_col)


# Count of significant terms per TF

n_terms_per_tf <- function(go_l) {
  n_l <- lapply(go_l, nrow)
  n_df <- data.frame(Symbol = names(go_l), N = unlist(n_l))
  n_df <- arrange(n_df, desc(N))
  return(n_df)
}

tf_n_hg <- n_terms_per_tf(go_coexpr_filt_hg)
tf_n_mm <- n_terms_per_tf(go_coexpr_filt_mm)
tf_n_ortho <- n_terms_per_tf(go_coexpr_filt_ortho)


# TFs with no associated terms

zero_hg <- tf_n_hg[tf_n_hg$N == 0, "Symbol"]
zero_mm <- tf_n_hg[tf_n_mm$N == 0, "Symbol"]
zero_ortho <- tf_n_hg[tf_n_ortho$N == 0, "Symbol"]


# Are ortho TFs correlated in their count of terms?

tf_n_all <- pc_ortho %>% 
  left_join(tf_n_hg, by = c("Symbol_hg" = "Symbol")) %>% 
  dplyr::rename(N_hg = N) %>% 
  left_join(tf_n_mm, by = c("Symbol_mm" = "Symbol")) %>%
  dplyr::rename(N_mm = N) %>% 
  filter(!is.na(N_hg), !is.na(N_mm)) %>% 
  left_join(tf_n_ortho, by = c("Symbol_hg" = "Symbol")) %>% 
  dplyr::rename(N_ortho = N) 


human0 <- tf_n_all$N_hg == 0
mouse0 <- tf_n_all$N_mm == 0

cor.test(tf_n_all$N_hg, tf_n_all$N_mm)
fisher.test(table(human0, mouse0))


# How well overlapping are the terms for ortho TFs?

tf_n_all$Jaccard <- unlist(lapply(1:nrow(tf_n_all), function(x) {
  tf_hg <- tf_n_all$Symbol_hg[x]
  tf_mm <- tf_n_all$Symbol_mm[x]
  terms_hg <- go_coexpr_filt_hg[[tf_hg]]$Name
  terms_mm <- go_coexpr_filt_mm[[tf_mm]]$Name
  length(intersect(terms_hg, terms_mm)) / length(union(terms_hg, terms_mm))
}))


summary(filter(tf_n_all, N_hg > 1 | N_mm > 1)$Jaccard)
summary(filter(tf_n_all, N_hg > 1 & N_mm > 1)$Jaccard)



# What are the most frequent/rare terms? Include the associated TFs

n_tfs_per_term <- function(go_l, term_col) {
  
  tf_terms <- lapply(names(go_l), function(tf) {
    term <- go_l[[tf]][[term_col]]
    if (length(term) == 0) return(NA)
    data.frame(Term = term, TF = tf, stringsAsFactors = FALSE)
  })
  
  tf_terms <- tf_terms[!is.na(tf_terms)]
  
  all_terms_df <- do.call(rbind, tf_terms)
  
  # Group by Term and tally occurrences
  result <- all_terms_df %>%
    group_by(Term) %>%
    summarise(
      N = n(),  # Count occurrences of each term
      TF = paste(unique(TF), collapse = ", ")  # Concat unique TFs
    ) %>% 
    arrange(desc(N))
  
  return(result)
}



term_n_hg <- n_tfs_per_term(go_coexpr_filt_hg, term_col)
term_n_mm <- n_tfs_per_term(go_coexpr_filt_mm, term_col)
term_n_ortho <- n_tfs_per_term(go_coexpr_filt_ortho, term_col)



# Are there terms that show species enrichment?

term_n_all <- data.frame(Term = union(term_n_hg$Term, term_n_mm$Term)) %>% 
  left_join(term_n_hg, by = "Term") %>% 
  dplyr::rename(N_hg = N, TF_hg = TF) %>% 
  left_join(term_n_mm, by = "Term") %>% 
  dplyr::rename(N_mm = N, TF_mm = TF) %>% 
  left_join(term_n_ortho, by = "Term") %>% 
  dplyr::rename(N_ortho = N, TF_ortho = TF) %>% 
  replace(is.na(.), 0)


cor.test(term_n_all$N_hg, term_n_all$N_mm)



# Plots
# ------------------------------------------------------------------------------


plot_bar <- function(df, count_col, group_col, xlab, ylab, title) {
  
  ggplot(df, aes(x = !!sym(count_col), 
                 y = reorder(!!sym(group_col), !!sym(count_col)))) +
    geom_bar(stat = "identity", colour = "slategrey", fill = "slategrey") +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
  
}



# Count of terms per TR

p1a <- plot_bar(tf_n_hg, 
                count_col = "N", 
                group_col = "Symbol", 
                xlab = "Count of GO terms (biological process)", 
                ylab = "TR", 
                title = "Human")


p1b <- plot_bar(tf_n_mm, 
                count_col = "N", 
                group_col = "Symbol", 
                xlab = "Count of GO terms (biological process)", 
                ylab = "TR", 
                title = "Mouse")


# Count of TRs per term

p2a <- plot_bar(term_n_hg, 
                count_col = "N", 
                group_col = "Term", 
                xlab = "Count of TRs belonging to term", 
                ylab = "GO terms (biological process)",
                title = "Human")


p2b <- plot_bar(term_n_mm, 
                count_col = "N", 
                group_col = "Term", 
                xlab = "Count of TRs belonging to term",
                ylab = "GO terms (biological process)", 
                title = "Mouse")



# Scatter of count of terms between ortho TRs

p3 <- 
  ggplot(tf_n_all, aes(x = N_hg, y = N_mm)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = "lm") +
  xlab("Count of significant terms (human)") +
  ylab("Count of significant terms (mouse)") +
  ggtitle("N=1,241 orthologous TRs") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))


# Scatter of count of TRs per term

p4 <- 
  ggplot(term_n_all, aes(x = N_hg, y = N_mm)) +
  geom_point(shape = 21, size = 2.4) +
  geom_smooth(method = "lm") +
  xlab("Count of associated TRs (human)") +
  ylab("Count of associated TRs (mouse)") +
  ggtitle("N=5,692 unique terms between species") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        plot.title = element_text(size = 20),
        plot.margin = margin(c(10, 20, 10, 10)))
  



# Simple bar chart of GO terms and -log10(adj.pval)

go_barplot <- function(df, topn = 15, title) {
  
  df %>% 
    mutate(Score = -log10(CorrectedMFPvalue)) %>% 
    slice_max(Score, n = topn) %>% 
    arrange(Score) %>%
    mutate(Name = factor(Name, levels = unique(Name))) %>% 
    ggplot(., aes(x = Score, y = Name)) +
    geom_bar(stat = "identity", col = "black", fill = "slategrey") +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
    xlab(expr("-log"[!!10] ~ "adjusted p-value")) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text.y = element_text(size = 15),
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
  
}



generate_and_save_plot <- function(params) {
  
  plot <- go_barplot(params$data, title = params$title)
  
  ggsave(
    plot, height = 7, width = 14, device = "png", dpi = 300,
    filename = file.path(plot_dir, params$filename)
  )
}


# Examples called in paper
plot_params <- list(
  list(data = go_coexpr_filt_hg$ASCL1, title = "ASCL1 human ranking", filename = "GO_ASCL1_human_coexpr.png"),
  list(data = go_coexpr_filt_mm$Ascl1, title = "ASCL1 mouse ranking", filename = "GO_ASCL1_mouse_coexpr.png"),
  list(data = go_coexpr_filt_ortho$ASCL1, title = "ASCL1 orthologous ranking", filename = "GO_ASCL1_ortho_coexpr.png"),
  list(data = go_coexpr_filt_hg$OLIG1, title = "Human OLIG1", filename = "GO_OLIG1_human_coexpr.png"),
  list(data = go_coexpr_filt_mm$Irf8, title = "Mouse Irf8", filename = "GO_Irf8_mouse_coexpr.png"),
  list(data = go_coexpr_filt_hg$NEUROD6, title = "Human NEUROD6", filename = "GO_NEUROD6_human_coexpr.png"),
  list(data = go_coexpr_filt_mm$Gata1, title = "Mouse Gata1", filename = "GO_Gata1_mouse_coexpr.png"),
  list(data = go_coexpr_filt_mm$Pax6, title = "Mouse Pax6", filename = "GO_Pax6_mouse_coexpr.png"),
  list(data = go_coexpr_filt_hg$SOX4, title = "Human SOX4", filename = "GO_SOX4_human_coexpr.png"),
  list(data = go_coexpr_filt_hg$E2F8, title = "Human E2F8", filename = "GO_E2F8_human_coexpr.png")
)


# Saving out

invisible(lapply(plot_params, generate_and_save_plot))
