## Export regions bound by a TF around a given TF
# TODO: range table functions...
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(igvR)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")
source("/home/amorin/Projects/TR_aggregation/R/utils/range_table_functions.R")

# Protein coding genes 
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# Load GR objects of bound regions 
bind_gr_hg <- readRDS("/space/scratch/amorin/R_objects/unibind_grlist_perm_human.RDS")
bind_gr_mm <- readRDS("/space/scratch/amorin/R_objects/unibind_grlist_perm_mouse.RDS")
bind_meta <- readRDS("/space/scratch/amorin/R_objects/Unibind_metadata.RDS")

# Convert protein coding to GR
pc_gr_hg <- pc_to_gr(filter(pc_hg, !is.na(Start)))
pc_gr_mm <- pc_to_gr(filter(pc_mm, !is.na(Start)))


# Functions
# ------------------------------------------------------------------------------ 


# Subset tf_gr to peaks within dist of gene_gr

filter_nearest_peaks <- function(tf_gr, gene_gr, dist = 500e3) {
  
  tf_gr <- 
    tf_gr[tf_gr@seqnames == intersect(seqnames(tf_gr), seqnames(gene_gr))]
  
  tf_gr$Distance <- 
    GenomicRanges::distance(tf_gr, gene_gr, ignore.strand = TRUE)
  
  tf_gr <- tf_gr[!is.na(tf_gr$Distance) & tf_gr$Distance < dist]
  tf_gr <- tf_gr[order(tf_gr$Distance)]
  
  return(tf_gr)
}



# Given duplicated ids, reduce the corresponding GRs into a single object.

reduce_duplicates <- function(ids, meta, gr_l, pad) {
  
  dedup <- lapply(ids, function(x) {
  file_id <- filter(meta, ID == x)$File
  stopifnot(length(file_id) == 2)
  reduce(c(gr_l[[file_id[1]]] + pad, gr_l[[file_id[2]]] + pad))
  })
  names(dedup) <- ids
  
  return(dedup)
}



# Organizing peak overlap for given TF/gene
# ------------------------------------------------------------------------------


# TF for ChIP-seq peaks, gene for neighbourhood of gene TSS to inspect
tf <- "ASCL1"
gene <- "DLL3"


# How much to add to both sides of 1bp peak centre for overlaps 
pad <- 150


# Subset protein coding to gene of interest
gene_gr_hg <- pc_gr_hg[pc_gr_hg$Symbol == gene]
gene_gr_mm <- pc_gr_mm[pc_gr_mm$Symbol == str_to_title(gene)]


# Experiment IDs of TF
meta_hg <- filter(bind_meta$Permissive_hg, Symbol == tf)
meta_mm <- filter(bind_meta$Permissive_mm, Symbol == tf)


# Subset GR list to TF of interest
tf_hg <- bind_gr_hg[meta_hg$File]
tf_mm <- bind_gr_mm[meta_mm$File]


# Reduce (union) peaks of duplicated experiments that have different scored motifs
dup_ids_hg <- unique(filter(meta_hg, Duplicate)$ID)
nodup_ids_hg <- unique(filter(meta_hg, !Duplicate)$File)
dedup_hg <- reduce_duplicates(dup_ids_hg, meta_hg, tf_hg, pad)

dup_ids_mm <- unique(filter(meta_mm, Duplicate)$ID)
nodup_ids_mm <- unique(filter(meta_mm, !Duplicate)$File)
dedup_mm <- reduce_duplicates(dup_ids_mm, meta_mm, tf_mm, pad)


# Combine non-duplicated GR list with deduplicated
dedup_tf_hg <- c(tf_hg[nodup_ids_hg], dedup_hg)
dedup_tf_mm <- c(tf_mm[nodup_ids_mm], dedup_mm)


# Pad non-duplicated GRs to same size as those that were duplicated/reduced
dedup_tf_hg[nodup_ids_hg] <- lapply(dedup_tf_hg[nodup_ids_hg], `+`, pad)
dedup_tf_mm[nodup_ids_mm] <- lapply(dedup_tf_mm[nodup_ids_mm], `+`, pad)


# Another reduction in case padding caused overlap (don't want to double count)
dedup_tf_hg <- lapply(dedup_tf_hg, reduce)
dedup_tf_mm <- lapply(dedup_tf_mm, reduce)


# Identify the peaks within 500kb of the gene TSS
peaks_hg <- lapply(dedup_tf_hg, filter_nearest_peaks, gene_gr = gene_gr_hg)
peaks_mm <- lapply(dedup_tf_mm, filter_nearest_peaks, gene_gr = gene_gr_mm)


# Create a single GR of all bound regions
all_hg <- reduce(unlist(reduce(GRangesList(peaks_hg))))
all_mm <- reduce(unlist(reduce(GRangesList(peaks_mm))))


# Get distance of the all set to gene TSS
all_hg$Distance <- GenomicRanges::distance(all_hg, gene_gr_hg, ignore.strand = TRUE)
all_hg <- all_hg[order(all_hg$Distance)]

all_mm$Distance <- GenomicRanges::distance(all_mm, gene_gr_mm, ignore.strand = TRUE)
all_mm <- all_mm[order(all_mm$Distance)]


# Count how many dataset overlap with the all set
all_hg$Count <- GenomicRanges::countOverlaps(all_hg, GRangesList(peaks_hg), ignore.strand = TRUE)
all_hg <- all_hg[order(-all_hg$Count)]

all_mm$Count <- GenomicRanges::countOverlaps(all_mm, GRangesList(peaks_mm), ignore.strand = TRUE)
all_mm <- all_mm[order(-all_mm$Count)]


# Convert to table and format for igvR and export
out_hg <- GenomicRanges::as.data.frame(all_hg) %>% 
  dplyr::select(seqnames, start, end, Count, Distance) %>% 
  mutate(seqnames = paste0("chr", seqnames))


out_mm <- GenomicRanges::as.data.frame(all_mm) %>% 
  dplyr::select(seqnames, start, end, Count, Distance) %>% 
  mutate(seqnames = paste0("chr", seqnames))


write.table(out_hg, sep = "\t", quote = FALSE, row.names = FALSE,
            file = paste0("/space/scratch/amorin/R_objects/Human_", tf, "_", gene, "_count_regions.tsv"))
                         

write.table(out_mm, sep = "\t", quote = FALSE, row.names = FALSE,
            file = paste0("/space/scratch/amorin/R_objects/Mouse_", tf, "_", gene, "_count_regions.tsv"))


# Inspecting the overlap of notable peaks from the reduced all set to original
# names(peaks_hg)[findOverlaps(all_hg[1], GRangesList(peaks_hg))@to]
# names(peaks_mm)[findOverlaps(all_mm[1], GRangesList(peaks_mm))@to]
