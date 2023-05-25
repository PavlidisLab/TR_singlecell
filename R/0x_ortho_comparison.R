## Examining aggregate coexpression between mouse and human
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/00_config.R")
source("R/utils/functions.R")
source("R/utils/plot_functions.R")

pc_ortho <- read.delim(pc_ortho_path, stringsAsFactors = FALSE)
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# Load mouse/human data
# Summarized ranks across data sets
# Access to organized species datasets will require cleaned metadata


dat_hg_agg <- matrix(0)
dat_mm_agg <- matrix(0)
dat_hg_all <- list()
dat_mm_all <- list()

common <- filter(pc_ortho, Symbol_hg %in% pc_hg$Symbol, Symbol_mm %in% pc_mm$Symbol)


# Most analysis would focus on the 1:1 DIOPT orthologs.
# But it would be interesting to look for gene family members presence

# For each TR, get the ability of one species to recover the other species's TR.
# AUPRC, ROC. Have to decide on a cut-off. Could show plot of recovery at different k
# Null would be... sample of matched genes?

