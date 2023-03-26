## Download organized data provided by Human Protein Atlas
## https://www.proteinatlas.org/about/download
## -----------------------------------------------------------------------------

library(tidyverse)
library(assertthat)


# Download, unzip, and ensure the TSV is readable

dl_check <- function(url, out_dir = "/space/grp/amorin/Expression_files/HPA/") {
  
  path_zip <-
    paste0(out_dir,
           str_replace(url, "https://www.proteinatlas.org/download/", ""))
  
  path_tsv <- str_replace(path_zip, ".zip", "")
  
  if (!file.exists(path_tsv)) {
    download.file(url, path_zip)
    unzip(path_zip, exdir = out_dir)
    stopifnot(is.readable(path_tsv))
    file.remove(path_zip)
  }
}


# 1: Normal tissue data
# Expression profiles for proteins in human tissues based on immunohistochem
# using tissue micro arrays. The tab-separated file includes Ensembl gene 
# identifier ("Gene"), tissue name ("Tissue"), annotated cell type ("Cell type"), 
# expression value ("Level"), and the gene reliability of the expression value 
# ("Reliability"). The data is based on The Human Protein Atlas version 21.1 and 
# Ensembl version 103.38.

url_1 <- "https://www.proteinatlas.org/download/normal_tissue.tsv.zip"
dl_check(url_1)


# 2: Subcellular location data
# Subcellular location of proteins based on immunofluorescently stained cells. The
# tab-separated file includes the following columns: Ensembl gene identifier 
# ("Gene"), name of gene ("Gene name"), gene reliability score ("Reliability"), 
# enhanced locations ("Enhanced"), supported locations ("Supported"), Approved 
# locations ("Approved"), uncertain locations ("Uncertain"), locations with 
# single-cell variation in intensity ("Single-cell variation intensity"), 
# locations with spatial single-cell variation ("Single-cell variation spatial"), 
# locations with observed cell cycle dependency (type can be one or more of 
# biological definition, custom data or correlation) ("Cell cycle dependency"), 
# Gene Ontology Cellular Component term identifier ("GO id") The data is based on
# The Human Protein Atlas version 21.1 and Ensembl version 103.38.

url_2 <- "https://www.proteinatlas.org/download/subcellular_location.tsv.zip"
dl_check(url_2)


# 3: RNA consensus tissue gene data
# Consensus transcript expression levels summarized per gene in 54 tissues based 
# on transcriptomics data from HPA and GTEx. The consensus normalized expression 
# ("nTPM") value is calculated as the maximum nTPM value for each gene in the two
# data sources. For tissues with multiple sub-tissues (brain regions, lymphoid 
# tissues and intestine) the maximum of all sub-tissues is used for the tissue 
# type. The tab-separated file includes Ensembl gene identifier ("Gene"), analysed
# sample ("Tissue") and normalized expression ("nTPM"). The data is based on The 
# Human Protein Atlas version 21.1 and Ensembl version 103.38.

url_3 <- "https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"
dl_check(url_3)


# 4: RNA HPA tissue gene data
# Transcript expression levels summarized per gene in 256 tissues based on RNA-seq. 
# The tab-separated file includes Ensembl gene identifier ("Gene"), analysed sample 
# ("Tissue"), transcripts per million ("TPM"), protein-transcripts per million 
# ("pTPM") and normalized expression ("nTPM"). The data is based on The Human 
# Protein Atlas version 21.1 and Ensembl version 103.38.

url_4 <- "https://www.proteinatlas.org/download/rna_tissue_hpa.tsv.zip"
dl_check(url_4)


# 5: RNA GTEx tissue gene data
# Transcript expression levels summarized per gene in 37 tissues based on RNA-seq.
# The tab-separated file includes Ensembl gene identifier ("Gene"), analysed 
# sample ("Tissue"), transcripts per million ("TPM"), protein-transcripts per 
# million ("pTPM") and normalized expression ("nTPM"). The data was obtained from 
# GTEx and is based on The Human Protein Atlas version 21.1 and Ensembl version 103.38.

url_5 <- "https://www.proteinatlas.org/download/rna_tissue_gtex.tsv.zip"
dl_check(url_5)


# 6: RNA FANTOM tissue gene data
# Transcript expression levels summarized per gene in 60 tissues based on CAGE 
# data. The tab-separated file includes Ensembl gene identifier ("Gene"), analysed
# sample ("Tissue"), tags per million ("Tags per million"), scaled-tags per million
# ("Scaled tags per million") and normalized expression ("nTPM"). The data was 
# obtained from FANTOM5 and is based on The Human Protein Atlas version 21.1 and 
# Ensembl version 103.38.

url_6 <- "https://www.proteinatlas.org/download/rna_tissue_fantom.tsv.zip"
dl_check(url_6)


# 7: RNA single cell type data
# Transcript expression levels summarized per gene in 76 cell types from 26 
# datasets. The tab-separated file includes Ensembl gene identifier ("Gene"), gene
# name ("Gene name"), analysed sample ("Cell type") and normalized expresion ("nTPM").

url_7 <- "https://www.proteinatlas.org/download/rna_single_cell_type.tsv.zip"
dl_check(url_7)


# 8: RNA single cell type tissue cluster data
# Transcript expression levels summarized per gene and cluster in 26 datasets. The
# tab-separated file includes Ensembl gene identifier ("Gene"), gene name 
# ("Gene name"), tissue ("Tissue"), analysed sample ("Cell type"), cluster
# ("Cluster"), read count ("Read count") and protein-transcripts per million ("pTPM")

url_8 <- "https://www.proteinatlas.org/download/rna_single_cell_type_tissue.tsv.zip"
dl_check(url_8)


# 9: RNA single cell read count data
# Read count per gene and cell in 26 datasets. The tab-separated file is in matrix
# format with Ensembl gene identifiers as columns and single cells as rows. 
# Columns included are tissue ("Tissue"), cell ("Cell") and cluster ("Cluster"). 

# NOTE: times out. download manually
#url_9 <- "https://www.proteinatlas.org/download/rna_single_cell_read_count.tsv.zip"
#dl_check(url_9)


# 10: RNA GTEx brain region gene data
# Transcript expression levels summarized per gene in 10 brain regions based on
# RNA-seq. The tab-separated file includes Ensembl gene identifier ("Gene"), 
# analysed sample ("Brain region"), transcripts per million ("TPM"), 
# protein-transcripts per million ("pTPM") and normalized expression ("nTPM"). 
# The data was obtained from GTEx and is based on The Human Protein Atlas version 
# 21.1 and Ensembl version 103.38.

url_10 <- "https://www.proteinatlas.org/download/rna_brain_gtex.tsv.zip"
dl_check(url_10)


# 11: RNA FANTOM brain region gene data
# Transcript expression levels summarized per gene in 14 brain regions based on 
# CAGE. The tab-separated file includes Ensembl gene identifier ("Gene"), analysed
# sample ("Brain region"), tags per million ("Tags per million"), scaled-tags per 
# million ("Scaled tags per million") and normalized expression ("nTPM"). The data
# was obtained from FANTOM5 and is based on The Human Protein Atlas version 21.1 
# and Ensembl version 103.38.

url_11 <- "https://www.proteinatlas.org/download/rna_brain_fantom.tsv.zip"
dl_check(url_11)


# 12: RNA pig brain region gene data
# Transcript expression levels summarized per gene in 15 brain regions based on 
# RNA-seq. The tab-separated file includes Ensembl gene identifier ("Gene"), 
# analysed sample ("Brain region") and transcripts per million ("TPM") and 
# protein-coding transcripts per million ("pTPM") and normalized expression 
# ("nTPM"). The data is based on The Human Protein Atlas version 21.1 and Ensembl 
# version 103.38

url_12 <- "https://www.proteinatlas.org/download/rna_pig_brain_hpa.tsv.zip"
dl_check(url_12)


# 13: RNA pig brain subregion sample gene data
# Transcript expression levels summarized per gene in 32 brain subregions per 
# sample based on RNA-seq. The tab-separated file includes Ensembl gene identifier
# for pig gene ("Gene"), main region ("Main region"), sub region ("Sub region"), 
# animal ("Animal"), transcripts per million ("TPM") and protein-coding 
# transcripts per million ("pTPM"). The data is based on The Human Protein Atlas 
# version 21.1 and Ensembl version 103.38.

url_13 <- "https://www.proteinatlas.org/download/rna_pig_brain_sample_hpa.tsv.zip"
dl_check(url_13)


# 14: RNA mouse brain region gene data
# Transcript expression levels summarized per gene in 13 brain regions based on 
# RNA-seq. The tab-separated file includes Ensembl gene identifier ("Gene"), 
# analysed sample ("Brain region") and transcripts per million ("TPM") and 
# protein-coding transcripts per million ("pTPM") and normalized expression 
# ("nTPM"). The data is based on The Human Protein Atlas version 21.1 and Ensembl 
# version 103.38.

url_14 <- "https://www.proteinatlas.org/download/rna_mouse_brain_hpa.tsv.zip"
dl_check(url_14)


# 15: RNA mouse brain subregion sample gene data
# Transcript expression levels summarized per gene in 19 brain subregions per 
# sample based on RNA-seq. The tab-separated file includes Ensembl gene identifier
# for mouse gene ("Gene"), main brain region ("Main region"), subregion 
# ("Subregion"), animal ("Animal"), transcripts per million ("TPM") and 
# protein-coding transcripts per million ("pTPM"). The data is based on The Human 
# Protein Atlas version 21.1 and Ensembl version 103.38.

url_15 <- "https://www.proteinatlas.org/download/rna_mouse_brain_sample_hpa.tsv.zip"
dl_check(url_15)


# 16: RNA Allen mouse brain region gene data
# Transcript expression levels summarized per gene in 11 brain regions based on 
# ISH. The tab-separated file includes Ensembl gene identifier ("Gene"), 
# analysed sample ("Tissue") and expression energy ("Expression energy"). The 
# data was obtained from Allen brain atlas and is based on The Human Protein 
# Atlas version 21.1 and Ensembl version 103.38.

url_16 <- "https://www.proteinatlas.org/download/rna_mouse_brain_allen.tsv.zip"
dl_check(url_16)


# 17: RNA HPA blood cell gene data
# Transcript expression levels summarized per gene in 18 cell types and total PBMC.
# The tab-separated file includes Ensembl gene identifier ("Gene"), analysed 
# sample ("Blood cell"), transcripts per million ("TPM"), protein-coding 
# transcripts per million ("pTPM") and normalized expression ("nTPM"). The data is
# based on The Human Protein Atlas version 21.1 and Ensembl version 103.38.

url_17 <- "https://www.proteinatlas.org/download/rna_blood_cell.tsv.zip"
dl_check(url_17)


# 18: RNA HPA blood cell sample gene data
# Transcript expression levels summarized per gene in 109 blood cell samples. 
# "rna_blood_cell_sample.tsv.zip" includes Ensembl gene identifier ("Gene"), 
# analysed sample ("Blood cell"), donor ("Donor"), transcripts per million 
# ("TPM"), protein-coding transcripts per million ("pTPM"). The data is based on 
# The Human Protein Atlas version 21.1 and Ensembl version 103.38.

url_18 <- "https://www.proteinatlas.org/download/rna_blood_cell_sample.tsv.zip"
dl_check(url_18)


# 19: RNA Monaco blood cell gene data
# Transcript expression levels summarized per gene in 30 blood cell types. The 
# tab-separated file includes Ensembl gene identifier ("Gene"), analysed sample 
# ("Blood cell"), transcripts per million ("TPM") and protein-coding transcripts 
# per million ("pTPM"). The data was obtained from Monaco publication and is based
# on The Human Protein Atlas version 21.1 and Ensembl version 103.38.

url_19 <- "https://www.proteinatlas.org/download/rna_blood_cell_monaco.tsv.zip"
dl_check(url_19)


# 20: RNA Schmiedel blood cell gene data
# Transcript expression levels summarized per gene in 15 blood cell types. The 
# tab-separated file includes Ensembl gene identifier ("Gene"), analysed sample 
# ("Blood cell") and transcripts per million ("TPM"). The data was obtained from
# Schmiedel publication and is based on The Human Protein Atlas version 21.1 and 
# Ensembl version 103.38.

url_20 <- "https://www.proteinatlas.org/download/rna_blood_cell_schmiedel.tsv.zip"
dl_check(url_20)


# 21: RNA HPA cell line gene data
# Transcript expression levels summarized per gene in 69 cell lines. The 
# tab-separated file includes Ensembl gene identifier ("Gene"), analysed sample 
# ("Cell line"), transcripts per million ("TPM"), protein-coding transcripts per 
# million ("pTPM") and normalized expression ("nTPM"). The data is based on The 
# Human Protein Atlas version 21.1 and Ensembl version 103.38.

url_21 <- "https://www.proteinatlas.org/download/rna_celline.tsv.zip"
dl_check(url_21)


# 22: RNA TCGA cancer sample gene data
# Transcript expression levels summarized per gene in 7932 samples from 17 
# different cancer types. The tab-separated file includes Ensembl gene identifier
# ("Gene"), analysed sample ("Sample"), cancer type ("Cancer") and fragments per 
# kilobase million ("FPKM"). The data is based on The Human Protein Atlas version 
# 21.1 and Ensembl version 103.38.

# NOTE: times out. download manually
# url_22 <- "https://www.proteinatlas.org/download/rna_cancer_sample.tsv.zip"
# dl_check(url_22)
