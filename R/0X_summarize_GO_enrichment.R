## TODO
## -----------------------------------------------------------------------------

library(clusterProfiler)
library(parallel)
library(tidyverse)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
source("R/00_config.R")

pc_ortho <- read.delim(pc_ortho_path)
go_coexpr_hg <- readRDS(go_coexpr_hg_path)
go_coexpr_mm <- readRDS(go_coexpr_mm_path)


# Example of one table - ASCL1 used in paper
example_hg <- "ASCL1"
example_mm <- "Ascl1"
go_coexpr_hg[[example_hg]]
go_coexpr_mm[[example_mm]]
go_coexpr_ortho[[example_hg]]



# Only consider results reaching sig. cut-off

filter_results <- function(go_l, cutoff = 0.05) {
  lapply(go_l, function(x) filter(x, p.adjust < cutoff))
}


go_coexpr_hg <- filter_results(go_coexpr_hg)
go_coexpr_mm <- filter_results(go_coexpr_mm)
go_coexpr_ortho <- filter_results(go_coexpr_ortho)


# Count of significant terms

n_terms <- function(go_l) {
  n_l <- lapply(go_l, nrow)
  n_df <- data.frame(Symbol = names(go_l), N = unlist(n_l))
  n_df <- arrange(n_df, desc(N))
  return(n_df)
}

n_hg <- n_terms(go_coexpr_hg)
n_mm <- n_terms(go_coexpr_mm)
n_ortho <- n_terms(go_coexpr_ortho)


# TFs with no associated terms

zero_hg <- n_hg[n_hg$N == 0, "Symbol"]
zero_mm <- n_hg[n_mm$N == 0, "Symbol"]
zero_ortho <- n_hg[n_ortho$N == 0, "Symbol"]


# Are ortho TFs correlated in their count of terms?

n_all <- pc_ortho %>% 
  left_join(n_hg, by = c("Symbol_hg" = "Symbol")) %>% 
  dplyr::rename(N_hg = N) %>% 
  left_join(n_mm, by = c("Symbol_mm" = "Symbol")) %>%
  dplyr::rename(N_mm = N) %>% 
  filter(!is.na(N_hg), !is.na(N_mm))


cor.test(n_all$N_hg, n_all$N_mm)









n_term_hist <- function(n_df, title) {
  
  ggplot(n_df, aes(x = N)) + 
    geom_histogram(bins = 100, colour = "slategrey", fill = "slategrey") +
    ylab("Count of TRs") +
    xlab("Count of significant GO terms (biological process)") +
    ggtitle(title) +
    theme_classic() +
    theme(axis.text = element_text(size = 20),
          axis.title = element_text(size = 20),
          plot.title = element_text(size = 20),
          plot.margin = margin(c(10, 20, 10, 10)))
  
}


p1a <- n_term_hist(n_hg, title = "Human")
p1b <- n_term_hist(n_mm, title = "Mouse")
p1c <- n_term_hist(n_ortho, title = "Ortho")



ggplot(n_all, aes(x = N_hg, y = N_mm)) +
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
