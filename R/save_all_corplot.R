## Save the correlation of the two input genes for the given species across
## all cell types and datasets.
## Usage: Rscript save_all_corplot.R human ASCL1 DLL1
## -----------------------------------------------------------------------------

library(tidyverse, quietly = TRUE)
library(WGCNA, quietly = TRUE)
library(parallel, quietly = TRUE)
source("R/utils/functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

out_dir <- file.path(plot_dir, "All_corplots")

args <- commandArgs(trailingOnly = TRUE)
species <- args[1] 
gene1 <- args[2] 
gene2 <- args[3] 


if (length(args) != 3) {
  stop("Incorrect number of arguments. Usage: Rscript save_all_corplot.R <species> <gene1> <gene2>")
}


species <- str_to_lower(species)


if (!(species %in% c("human", "hg", "mouse", "mm"))) {
  stop("Species not recognized")
}


# Return a list of the necessary data.
# species: one of human, hg, mouse, mm
# gene1 and gene2: a valid ensembl protein coding symbol

load_data <- function(species, gene1, gene2) {
  
  if (species %in% c("human", "hg")) {
    
    pc <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
    sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)
    ids <- filter(sc_meta, Species == "Human")$ID
    
    gene1 <- str_to_upper(gene1)
    gene2 <- str_to_upper(gene2)
    stopifnot(all(c(gene1, gene2) %in% pc$Symbol))
    
    cor_path <- file.path(out_dir, paste0("Human_", gene1, "_", gene2, ".RDS"))
    plot_path <- file.path(out_dir, paste0("Human_", gene1, "_", gene2, ".png"))
    
  } else if (species %in% c("mouse", "mm")) {
    
    pc <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
    sc_meta <- read.delim(sc_meta_path, stringsAsFactors = FALSE)
    ids <- filter(sc_meta, Species == "Mouse")$ID
    
    gene1 <- str_to_title(gene1)
    gene2 <- str_to_title(gene2)
    stopifnot(all(c(gene1, gene2) %in% pc$Symbol))
    
    cor_path <- file.path(out_dir, paste0("Mouse_", gene1, "_", gene2, ".RDS"))
    plot_path <- file.path(out_dir, paste0("Mouse_", gene1, "_", gene2, ".png"))
    
  } 
  
  return(list(
    ids = ids,
    gene1 = gene1,
    gene2 = gene2,
    cor_path = cor_path,
    plot_path = plot_path
  ))
}



dat <- load_data(species, gene1, gene2)



# Generate a list of gene-gene cor for every cell type in each dataset


if (!file.exists(dat$cor_path)) {
  cor_l <- get_all_cor_l(ids = dat$ids, gene1 = dat$gene1, dat$gene2)
  saveRDS(cor_l, dat$cor_path)
} else {
  cor_l <- readRDS(dat$cor_path)
}



# Plot


# Akin to a forest plot, show the spread of cell-type correlations for a given
# gene pair across all experiments
# cor_l: a named list of named numeric vectors, where list names are the dataset
# IDs, vector names are the cell types, and the values are cell type correlations


all_corplot <- function(cor_l) {
  
  cor_df <- do.call(
    rbind, 
    lapply(names(cor_l), function(x) data.frame(Cor = cor_l[[x]], ID = x))
  )
  
  ggplot(cor_df, aes(x = Cor, y = reorder(ID, Cor, FUN = median))) +
    geom_point(alpha = 0.4, shape = 21, size = 2.4) +
    geom_boxplot(outlier.shape = NA, coef = 0, fill = "slategrey") +
    geom_vline(xintercept = 0, colour = "black") +
    xlab("Pearson's correlation") +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 30),
          axis.text.x = element_text(size = 25),
          axis.text.y = element_text(size = 10),
          plot.margin = margin(c(10, 20, 10, 10)))
  
}



p <- all_corplot(cor_l)


ggsave(p, height = 14, width = 9, device = "png", dpi = 300,
       filename = dat$plot_path)
