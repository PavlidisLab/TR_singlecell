## Load a google sheets metadata table tracking scRNA-seq experiments, then
## format and save a local copy for available experiments
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")


meta <- read_sheet(gsheets_id, sheet = "Main", trim_ws = TRUE)


loaded <- lapply(meta$ID, function(x) {
  file.exists(file.path(amat_dir, x, paste0(x, "_RSR_allrank.tsv")))
})


meta_loaded <- meta[unlist(loaded), ]


# For inspecting those that are not loaded/failed
missing <- setdiff(meta$ID, meta_loaded$ID)


stopifnot(all(meta_loaded$Species %in% c("Human", "Mouse")))


# Inspect assay/technology for consistency
n_platform <- sort(table(meta$Platform))



# Load a dataset's count matrix and meta and return a list of two data.frames, 
# the first which counts the cells and cell types, and the second which tallies
# the count of cells for each unique cell type

# NOTE: as.character() call in ct_count is done to prevent counting cell types 
# with 0 counts (cell type is a factor), consistent with how NA counts are 
# tracked in the main RSR_allrank() function. Cell types with counts above 0 but
# still under the min filter (default 20) are still tallied, however.

add_meta <- function(id) {
  
  dat <- load_dat_list(id)[[1]]
  
  ct_count <- data.frame(table(as.character(dat$Meta$Cell_type)))
  colnames(ct_count) <- c("Cell_type", "N_cells")
  
  dat_count <- data.frame(
    ID = id,
    N_cells = ncol(dat$Mat),
    N_celltypes = nrow(ct_count)
  )
  
  return(list(Dat_count = dat_count, Ct_count = ct_count))
  
}



meta_counts <- lapply(meta_loaded$ID, add_meta)



# List of cell counts per cell type for each dataset for later inspection
counts_l <- lapply(meta_counts, `[`, "Ct_count")
names(counts_l) <- meta_loaded$ID


# Bind cell counts and enter into meta
meta_cols <- do.call(rbind, lapply(meta_counts, `[[`, "Dat_count"))
stopifnot(identical(meta_cols$ID, meta_loaded$ID))
meta_loaded$N_cells <- meta_cols$N_cells
meta_loaded$N_celltypes <- meta_cols$N_celltypes


# Summarize counts by species
summ_cells <- lapply(split(meta_loaded, meta_loaded$Species), function(x) summary(x$N_cells))
summ_cts <- lapply(split(meta_loaded, meta_loaded$Species), function(x) summary(x$N_celltypes))


# Save meta and list of cell types

write.table(meta_loaded,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE,
            file = sc_meta_path)


saveRDS(counts_l, file = celltype_list_path)



# Plot 
# ------------------------------------------------------------------------------


# Barchart of data set counts by species

p1 <- count(meta_loaded, Species) %>% 
  ggplot(., aes(x = Species, y = n)) +
  geom_bar(stat = "identity", fill = "#1c9099", col = "black", width = 0.6) +
  ylab("Count of data sets") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.title.x = element_blank(),
        plot.margin = margin(10, 5, 5, 5))



ggsave(p1, height = 6, width = 4.5, device = "png", dpi = 300,
       filename = file.path(plot_dir, "dataset_counts_by_species.png"))


# Histogram of log10 cell counts by species

p2 <- ggplot(meta_loaded, aes(x = log10(N_cells), fill = Species)) +
  geom_histogram(bins = 30) +
  ylab("Count of data sets") +
  xlab("Log10 count of cells") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = c(0.85, 0.75))


ggsave(p2, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(plot_dir, "cell_counts_by_species.png"))


# Histogram of cell type counts by species

p3 <- ggplot(meta_loaded, aes(x = N_celltypes, fill = Species)) +
  geom_histogram(bins = 50) +
  ylab("Count of data sets") +
  xlab("Count of cell types") +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = c(0.85, 0.75))


ggsave(p3, height = 6, width = 6, device = "png", dpi = 300,
       filename = file.path(plot_dir, "celltype_counts_by_species.png"))
