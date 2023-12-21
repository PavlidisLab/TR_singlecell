library(tidyverse)
library(GenomicRanges)
library(igvR)
source("R/utils/functions.R")
source("R/utils/vector_comparison_functions.R")
source("R/utils/plot_functions.R")
source("R/00_config.R")

# Protein coding genes 
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

# 
collection <- "Permissive"
stopifnot(collection %in% c("Robust", "Permissive"))
bind_summary_path <- paste0("/space/scratch/amorin/R_objects/unibind_", collection, "_bindscore_summary.RDS")
bind_summary <- readRDS(bind_summary_path)
bind_dat <- readRDS("/space/scratch/amorin/R_objects/processed_unibind_data.RDS")
bind_gr_hg <- readRDS("/space/scratch/amorin/R_objects/unibind_grlist_perm_human.RDS")
bind_gr_mm <- readRDS("/space/scratch/amorin/R_objects/unibind_grlist_perm_mouse.RDS")
bind_meta <- readRDS("/space/scratch/amorin/R_objects/Unibind_metadata.RDS")

source("~/TRagg/R/utils/range_table_functions.R")
pc_gr_hg <- pc_to_gr(filter(pc_hg, !is.na(Start)))
pc_gr_mm <- pc_to_gr(filter(pc_mm, !is.na(Start)))


# Functions
# ------------------------------------------------------------------------------ 


filter_nearest_peaks <- function(tf_gr, gene_gr, dist = 500e3) {
  
  tf_gr <- 
    tf_gr[tf_gr@seqnames == intersect(seqnames(tf_gr), seqnames(gene_gr))]
  
  tf_gr$Distance <- 
    GenomicRanges::distance(tf_gr, gene_gr, ignore.strand = TRUE)
  
  tf_gr <- tf_gr[!is.na(tf_gr$Distance) & tf_gr$Distance < dist]
  tf_gr <- tf_gr[order(tf_gr$Distance)]
  
  return(tf_gr)
}




#
tf <- "ASCL1"
gene <- "DLL3"


#
gene_gr_hg <- pc_gr_hg[pc_gr_hg$Symbol == gene]
gene_gr_mm <- pc_gr_mm[pc_gr_mm$Symbol == str_to_title(gene)]

#
meta_hg <- filter(bind_meta$Permissive_hg, Symbol == tf)
meta_mm <- filter(bind_meta$Permissive_mm, Symbol == tf)

#
l_hg <- bind_gr_hg[meta_hg$File]
l_mm <- bind_gr_mm[meta_mm$File]


# Reduce duplicated experiments that have different scored motifs
dup_ids_hg <- unique(filter(meta_hg, Duplicate)$ID)
nodup_ids_hg <- unique(filter(meta_hg, !Duplicate)$File)

dup_ids_mm <- unique(filter(meta_mm, Duplicate)$ID)
nodup_ids_mm <- unique(filter(meta_mm, !Duplicate)$File)



dedup_hg <- lapply(dup_ids_hg, function(x) {
  file_id <- filter(meta_hg, ID == x)$File
  reduce(c(l_hg[[file_id[1]]] + 150, l_hg[[file_id[2]]] + 150))
})
names(dedup_hg) <- dup_ids_hg



dedup_mm <- lapply(dup_ids_mm, function(x) {
  file_id <- filter(meta_mm, ID == x)$File
  reduce(c(l_mm[[file_id[1]]] + 150, l_mm[[file_id[2]]] + 150))
})
names(dedup_mm) <- dup_ids_mm


#
dedup_l_hg <- c(l_hg[nodup_ids_hg], dedup_hg)
dedup_l_mm <- c(l_mm[nodup_ids_mm], dedup_mm)

#
peaks_hg <- lapply(dedup_l_hg, filter_nearest_peaks, gene_gr = gene_gr_hg)
peaks_mm <- lapply(dedup_l_mm, filter_nearest_peaks, gene_gr = gene_gr_mm)

#
peaks_add_hg <- lapply(peaks_hg, `+`, 150)
peaks_add_mm <- lapply(peaks_mm, `+`, 150)

#
peaks_red_hg <- lapply(peaks_add_hg, reduce)
peaks_red_mm <- lapply(peaks_add_mm, reduce)

#
n_peaks_hg <- data.frame(
  Original = unlist(lapply(peaks_hg, length)),
  Reduce = unlist(lapply(peaks_red_hg, length)))


n_peaks_mm <- data.frame(
  Original = unlist(lapply(peaks_mm, length)),
  Reduce = unlist(lapply(peaks_red_mm, length)))


#
all_red_hg <- reduce(unlist(reduce(GRangesList(peaks_red_hg))))
all_red_mm <- reduce(unlist(reduce(GRangesList(peaks_red_mm))))

#
all_red_hg$Distance <- GenomicRanges::distance(all_red_hg, gene_gr_hg, ignore.strand = TRUE)
all_red_hg <- all_red_hg[order(all_red_hg$Distance)]

#
all_red_mm$Distance <- GenomicRanges::distance(all_red_mm, gene_gr_mm, ignore.strand = TRUE)
all_red_mm <- all_red_mm[order(all_red_mm$Distance)]

#
all_ol_hg <- findOverlaps(all_red_hg, GRangesList(peaks_red_hg), ignore.strand = TRUE)
all_ol_mm <- findOverlaps(all_red_mm, GRangesList(peaks_red_mm), ignore.strand = TRUE)

#
all_red_hg$Count <- GenomicRanges::countOverlaps(all_red_hg, GRangesList(peaks_red_hg), ignore.strand = TRUE)
all_red_mm$Count <- GenomicRanges::countOverlaps(all_red_mm, GRangesList(peaks_red_mm), ignore.strand = TRUE)

head(all_red_hg[order(-all_red_hg$Count)])
head(all_red_mm[order(-all_red_mm$Count)])



out_hg <- GenomicRanges::as.data.frame(all_red_hg[order(-all_red_hg$Count)])
out_mm <- GenomicRanges::as.data.frame(all_red_mm[order(-all_red_mm$Count)])


write.table(out_hg, sep = "\t", quote = FALSE, row.names = FALSE,
            file = paste0("~/Robj/Human_", tf, "_", gene, "_count_regions.tsv"))
                         

write.table(out_mm, sep = "\t", quote = FALSE, row.names = FALSE,
            file = paste0("~/Robj/Mouse_", tf, "_", gene, "_count_regions.tsv"))


# findOverlaps(all_red_hg[18], GRangesList(peaks_red_hg))
# findOverlaps(all_red_mm[35], GRangesList(peaks_red_mm))
