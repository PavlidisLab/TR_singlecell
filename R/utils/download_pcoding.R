## Ensembl has every TSS for a gene. Refseq select is only using one TSS per gene.
## -----------------------------------------------------------------------------

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

library(biomaRt)
library(tidyverse)


# Usage: download_ensembl_pcoding("ensembl_pcoding_human_V109.tsv", "Human")
# This will return the most recent ensembl version: at time of writing V109
# (Human GRCh38.p13, Mouse GRCm39)  https://useast.ensembl.org/index.html
# ------------------------------------------------------------------------------


download_ensembl_pcoding <- function(outfile, 
                                     species) {  # Human|Mouse 
  
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
  
  ens_mart <- useEnsembl(biomart = "ensembl", dataset = species_data)
  
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




# Usage: download_refseq("refseq_pcoding_human.tsv", "Human")
# Using ensembl range formatting: no 'chr' prefix and strand as 1/-1
# Overview of selection: https://www.ncbi.nlm.nih.gov/refseq/refseq_select/
# Table schema (mouse) http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1827432270_fzD31oCmH8TcANzQ82Yzz7h0RKlz&hgta_doSchemaDb=mm10&hgta_doSchemaTable=ncbiRefSeqSelect
# ---
# bin	
# V2 == name == Refseq_ID	
# V3 == chrom == Chromsome	
# V4 == strand == Strand	
# V5 == txStart == Start	
# V6 == txEnd == End	
# cdsStart	
# cdsEnd	
# exonCount	
# exonStarts	
# exonEnds	
# score	
# V13 == name2 == Symbol	
# cdsStartStat	
# cdsEndStat	
# exonFrames
#
# ------------------------------------------------------------------------------


download_refseq <- function(outfile, 
                            species) {  # Human|Mouse
  
  if (file.exists(outfile)) return(message(outfile, " already exists!"))
  
  if (species == "Mouse") {
    link <- "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ncbiRefSeqSelect.txt.gz"
    chr <- c(1:19, "MT", "X", "Y")
  } else if (species == "Human") {
    link <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz"
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
