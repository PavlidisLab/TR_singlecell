## Download the necessary tables for downstream analysis
## 1) Protein coding genes (both ensembl and refseq)
## 2) List of TFs from Animal TF DB
## 3) DIOPT 1:1 orthologous genes between mouse and human [Pavlab server only]
## 4) L/S cytosolic ribosomal genes from HUGO
## -----------------------------------------------------------------------------

source("R/00_config.R")
library(biomaRt)
library(tidyverse)


# 1) Protein coding tables
# Ensembl has every TSS for a gene. Refseq select is only using one TSS per gene.
# All analysis was done on hg38/mm10. Fixing ensembl at V98.
# Using ensembl range formatting: no 'chr' prefix and strand as 1/-1
# https://www.ncbi.nlm.nih.gov/refseq/refseq_select/
# ------------------------------------------------------------------------------


download_refseq <- function(outfile, 
                            species) {  # Human|Mouse
  
  if (file.exists(outfile)) return(message(outfile, " already exists!"))
  
  if (species == "Mouse") {
    link <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ncbiRefSeqSelect.txt.gz"
    chr <- c(1:19, "MT", "X", "Y")
  } else if (species == "Human") {
    link <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz"
    chr <- c(1:22, "MT", "X", "Y")
  }
  
  download.file(link, outfile)
  
  refseq <- read.delim(outfile, stringsAsFactors = FALSE, header = FALSE)
  
  refseq <- dplyr::select(refseq, c(V3, V5, V6, V4, V2, V13))
  
  colnames(refseq) <- c("Chromosome",
                        "Start",
                        "End",
                        "Strand",
                        "Refseq_ID",
                        "Symbol")
  
  refseq <- refseq %>% 
    mutate(
      Strand = ifelse(Strand == "+", 1, -1),
      Transcription_start_site = ifelse(Strand == 1, Start, End),
      Chromosome = str_replace(Chromosome, "chr", "")) %>% 
    filter(Chromosome %in% chr) %>% 
    arrange(match(Chromosome, chr), Transcription_start_site) %>% 
    dplyr::relocate(Transcription_start_site, .after = Chromosome)
  
  
  write.table(refseq,
              quote = FALSE,
              row.names = FALSE,
              sep = "\t",
              file = outfile)
  
}


download_ensembl_pcoding <- function(outfile,
                                     species,  # Human|Mouse
                                     version = "105") {
  
  if (file.exists(outfile)) return(message(outfile, " already exists!"))
  
  if (species == "Mouse") {
    symbol = "mgi_symbol"
    species_data = "mmusculus_gene_ensembl"
    chr_filter <- c(1:19, "MT", "X", "Y")
    
  } else if (species == "Human") {
    symbol = "hgnc_symbol"
    species_data = "hsapiens_gene_ensembl"
    chr_filter <- c(1:22, "MT", "X", "Y")
  }
  
  attributes <- c(
    "chromosome_name",
    "transcription_start_site",
    "transcript_start",
    "transcript_end",
    "strand",
    "ensembl_gene_id",
    symbol,
    "ucsc",
    "gene_biotype"
  )
  
  ens_mart <- useEnsembl(biomart = "ensembl",
                         dataset = species_data,
                         version = version)
  
  anno_table <- getBM(
    attributes = attributes,
    filters = "chromosome_name",
    values = chr_filter,
    mart = ens_mart,
    useCache = FALSE
  )
  
  # only protein coding gene type and order the table by chromosome then by TSS
  anno_table <- anno_table %>% 
    filter(gene_biotype == "protein_coding") %>%
    arrange(match(chromosome_name, chr_filter), transcription_start_site)
  
  anno_table$gene_biotype <- NULL
  
  colnames(anno_table) <- c(
    "Chromosome",
    "Transcription_start_site",
    "Start",
    "End",
    "Strand",
    "Gene_ID",
    "Symbol",
    "Transcript_ID"
  )
  
  write.table(anno_table,
              quote = FALSE,
              row.names = FALSE,
              sep = "\t",
              file = outfile)
  
}



download_refseq(outfile = ref_hg_path, species = "Human")
download_refseq(outfile = ref_mm_path, species = "Mouse")


download_ensembl_pcoding(outfile = ens_hg_path, species = "Human")
download_ensembl_pcoding(outfile = ens_mm_path, species = "Mouse")


# Ensuring Refseq and Ensembl share the same symbols. Refseq used for most 
# analyses, Ensembl to map Ensembl IDs to symbols for datasets with IDs only


ref_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
ref_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)
ens_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
ens_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)


# Keep only X copies of duplicated genes in human

dupl_genes <- ref_hg %>% group_by(Symbol) %>% filter(n() > 1)

ref_hg <- filter(ref_hg, !(Symbol %in% dupl_genes$Symbol & Chromosome == "Y"))

  
# Reduce Ensembl to Refseq select, keeping only a single entry per symbol, and
# only the X copy of duplicated counts. As using gene-summarized counts, the 
# exact transcript is unimportant


ens_hg <- ens_hg %>%
  filter(Symbol %in% ref_hg$Symbol & Symbol != "") %>%
  right_join(ref_hg[, c("Symbol", "Refseq_ID")], by = "Symbol") %>%
  filter(!(Symbol %in% dupl_genes$Symbol & Chromosome == "Y")) %>%
  distinct(Symbol, .keep_all = TRUE) %>%
  arrange(match(Symbol, ref_hg$Symbol))


ens_mm <- ens_mm %>% 
  filter(Symbol %in% ref_mm$Symbol & Symbol != "") %>% 
  right_join(ref_mm[, c("Symbol", "Refseq_ID")], by = "Symbol") %>% 
  distinct(Symbol, .keep_all = TRUE) %>%
  arrange(match(Symbol, ref_mm$Symbol))


stopifnot(identical(ens_hg$Symbol, ref_hg$Symbol))
stopifnot(identical(ens_mm$Symbol, ref_mm$Symbol))


write.table(ref_hg, sep = "\t", quote = FALSE, file = ref_hg_path)
write.table(ref_mm, sep = "\t", quote = FALSE, file = ref_mm_path)
write.table(ens_hg, sep = "\t", quote = FALSE, file = ens_hg_path)
write.table(ens_mm, sep = "\t", quote = FALSE, file = ens_mm_path)



# 2) List of TFs http://bioinfo.life.hust.edu.cn/AnimalTFDB4/#/
# NOTE: Analysis was carried out using V3; V4 has since been released.
# https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Homo_sapiens_TF
# https://guolab.wchscu.cn/AnimalTFDB4_static/download/TF_list_final/Mus_musculus_TF
# ------------------------------------------------------------------------------


tfs_hg_url <- "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF"
tfs_mm_url <- "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Mus_musculus_TF"


download_tfs <- function(tfs_url,
                         tfs_path,
                         pc_path) {
  
  if (!file.exists(tfs_path)) {
    
    download.file(url = tfs_url, destfile = tfs_path)
    
    pc <- read.delim(pc_path, stringsAsFactors = FALSE)
    tfs <- read.delim(tfs_path, stringsAsFactors = FALSE)
    
    tfs <- tfs %>% 
      filter(Symbol %in% pc$Symbol) %>%
      distinct(Symbol, .keep_all = TRUE)
    
    write.table(tfs, sep = "\t", quote = FALSE, row.names = FALSE, file = tfs_path)
  
  }
}


download_tfs(tfs_url = tfs_hg_url, tfs_path = tfs_hg_path, pc_path = ens_hg_path)
download_tfs(tfs_url = tfs_mm_url, tfs_path = tfs_mm_path, pc_path = ens_mm_path)



# 3) Filter/organize high-confidence 1:1 orthologs between mouse and human.
# https://www.flyrnai.org/diopt
# Note that this was a provided data dump that lives on Pavlab servers.
# Recommended heuristic (passed on by Sanja): keep scores >= 5, require mutual 
# bestscore, and remove symbols with more than one match
# ------------------------------------------------------------------------------


diopt <- read.delim(diopt_path, stringsAsFactors = FALSE)

# Protein coding genes
pc_hg <- read.delim(ref_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm_path, stringsAsFactors = FALSE)

symbol_hg <- unique(pc_hg$Symbol)
symbol_mm <- unique(pc_mm$Symbol)

diopt_filt <- filter(diopt,
                     score >= 5 &
                     human_symbol %in% pc_hg$Symbol &
                     symbol2 %in% pc_mm$Symbol &
                     ID_type2 == "MGI" &
                     best_score == "Yes" &
                     best_score_rev == "Yes")

split_hg <- split(diopt_filt, diopt_filt$human_symbol)
split_mm <- split(diopt_filt, diopt_filt$symbol2)

which_gt1_hg <- which(unlist(lapply(split_hg, nrow)) > 1)
which_gt1_mm <- which(unlist(lapply(split_mm, nrow)) > 1)

diopt_filt <- filter(diopt_filt, 
                        !(human_symbol) %in% names(which_gt1_hg) &
                          !(symbol2) %in% names(which_gt1_mm))

stopifnot(all(diopt_filt$human_symbol %in% pc_hg$Symbol))
stopifnot(all(diopt_filt$symbol2 %in% pc_mm$Symbol))

stopifnot(identical(n_distinct(diopt_filt$symbol2), 
                    n_distinct(diopt_filt$human_symbol)))

# create a df with an ID/key column as not all symbols have an exact 1:1 naming match

symbols <- data.frame(
  Symbol_hg = diopt_filt$human_symbol,
  Symbol_mm = diopt_filt$symbol2,
  ID = paste(diopt_filt$human_symbol, diopt_filt$symbol2, sep = "_")
)

stopifnot(identical(n_distinct(symbols$ID), nrow(symbols)))


if (!file.exists(pc_ortho_path)) {
  write.table(symbols,
              sep = "\t",
              quote = FALSE,
              file = pc_ortho_path)
}


# 4) L/S cytosolic ribosomal genes from HUGO
# https://www.genenames.org/data/genegroup/#!/group/728
# Group 728 == (S)mall, Group 729 = (L)arge
# ------------------------------------------------------------------------------


sribo_url <- "https://www.genenames.org/cgi-bin/genegroup/download?id=728&type=node"
lribo_url <- "https://www.genenames.org/cgi-bin/genegroup/download?id=729&type=node"


if (!file.exists(ribo_path)) {
  
  pc_ortho <- read.delim(pc_ortho_path, stringsAsFactors = FALSE)
  sribo <- read.delim(url(sribo_url), sep = "\t", stringsAsFactors = FALSE)
  lribo <- read.delim(url(lribo_url), sep = "\t", stringsAsFactors = FALSE)
  
  ribo <- filter(pc_ortho, Symbol_hg %in% c(sribo$Approved.symbol, lribo$Approved.symbol))
  
  write.table(ribo,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              file = ribo_path)

}
