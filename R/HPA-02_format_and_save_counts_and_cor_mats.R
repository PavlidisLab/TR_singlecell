## Load HPA expression datasets, convert from long/skinny data frames into count
## matrices, generate gene-gene cor matrices for each dataset, and export.
## -----------------------------------------------------------------------------

library(tidyverse)
library(WGCNA)
library(parallel)
library(reshape2)

# Output paths
expr_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_expression_mat_list.RDS"
cor_mat_l_path <- "/home/amorin/scratch/R_objects/HPA_cor_mat_list.RDS"

# Loading the HPA expression data sets
hpa_dir <- "/space/grp/amorin/Expression_files/HPA/"
hpa <- read.delim(file.path(hpa_dir, "rna_tissue_hpa.tsv"), stringsAsFactors = FALSE)
gtex <- read.delim(file.path(hpa_dir, "rna_tissue_gtex.tsv"), stringsAsFactors = FALSE)
fantom <- read.delim(file.path(hpa_dir, "rna_tissue_fantom.tsv"), stringsAsFactors = FALSE)
allen <- read.delim(file.path(hpa_dir, "rna_mouse_brain_allen.tsv"), stringsAsFactors = FALSE)
sc <- read.delim(file.path(hpa_dir, "rna_single_cell_type.tsv"), stringsAsFactors = FALSE)
sc_tissue <- read.delim(file.path(hpa_dir, "rna_single_cell_type_tissue.tsv"), stringsAsFactors = FALSE)
pigbrain <- read.delim(file.path(hpa_dir, "rna_pig_brain_hpa.tsv"), stringsAsFactors = FALSE)
mousebrain <- read.delim(file.path(hpa_dir, "rna_mouse_brain_hpa.tsv"), stringsAsFactors = TRUE)
blood_hpa <- read.delim(file.path(hpa_dir, "rna_blood_cell.tsv"), stringsAsFactors = TRUE)
blood_monaco <- read.delim(file.path(hpa_dir, "rna_blood_cell_monaco.tsv"), stringsAsFactors = TRUE)
blood_schmiedel <- read.delim(file.path(hpa_dir, "rna_blood_cell_schmiedel.tsv"), stringsAsFactors = TRUE)
cline <- read.delim(file.path(hpa_dir, "rna_celline.tsv"), stringsAsFactors = TRUE)

# Ensembl protein coding genes
pc_hg <- read.delim("/space/grp/amorin/Metadata/ensembl_human_protein_coding_105.tsv", stringsAsFactors = FALSE)


# Collapsing single cell tissue columns to get unique ID
sc_tissue <- sc_tissue %>%
  mutate(ID = str_replace_all(paste(Tissue, Cluster, Cell.type, sep = "_"), " ", "_"))

sc$Cell.type <- str_replace_all(sc$Cell.type, " ", "_")


# TCGA cancer processing - collapse samples grouped by cancer type by mean expression.
# Note that FPKM - others were some form of TPM. Also ensmebl IDs so must get symbol
cancer <- read.delim(file.path(hpa_dir, "rna_cancer_sample.tsv"), stringsAsFactors = FALSE)
cancer <- aggregate(cancer$FPKM, list(cancer$Gene, cancer$Cancer), FUN = mean) 
colnames(cancer) <- c("Gene.name", "Cancer", "FPKM")
cancer$Gene.name <- pc_hg$Symbol[match(cancer$Gene.name, pc_hg$Gene_ID)]
cancer <- filter(cancer, Gene.name != "")


# Process df before casting to matrix: Only keep genes present in all tissues,
# remove genes with all 0 counts, and log2 transform + 1
# https://stackoverflow.com/questions/9617348/
# ------------------------------------------------------------------------------


process_mat <- function(expr_df, 
                        var = "nTPM", 
                        group = "Tissue",
                        gene_name = "Gene.name") {
  
  dup <- split(expr_df, expr_df$Gene.name)
  dup <- unlist(lapply(dup, nrow))
  
  # hacky means of getting the most represented count
  dom_count <- as.integer(names(sort(table(dup), decreasing = TRUE)[1]))
  
  # Remove duplicated gene symbols
  dup <- dup[dup != dom_count]
  expr_df <- filter(expr_df, !(Gene.name %in% names(dup)))

  # To matrix and remove all 0s
  mat <- acast(expr_df, 
               as.formula(paste0(gene_name, " ~", group)),
               value.var = var, 
               drop = FALSE)
  
  all0 <- which(rowSums(mat) == 0)
  mat <- mat[-all0, ]
  mat <- log2(mat + 1)
  
  return(mat)
}



# All processed matrics in a list
# ------------------------------------------------------------------------------


mat_list <- list(
  HPA = process_mat(hpa),
  GTEX = process_mat(gtex),
  FANTOM = process_mat(fantom, var = "Normalized.tags.per.million"),
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



# Gene by gene cor mats for each data set


cor_list <- mclapply(mat_list, function(x) {
  WGCNA::cor(t(x), use = "pairwise.complete.obs")
}, mc.cores = 8)


# Save out

saveRDS(mat_list, expr_mat_l_path)
saveRDS(cor_list, cor_mat_l_path)
