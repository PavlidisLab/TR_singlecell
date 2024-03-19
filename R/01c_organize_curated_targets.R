## Organize curated low-throughput TF-target interactions from Chu 2021, which
## also aggregates from other low throughput databases, and from on going Pavlab
## curation on google sheets.
## -----------------------------------------------------------------------------

library(tidyverse)
library(googlesheets4)
source("R/00_config.R")
source("R/utils/functions.R")

# Protein coding genes 
pc_hg <- read.delim(ens_hg_path, stringsAsFactors = FALSE)
pc_mm <- read.delim(ens_mm_path, stringsAsFactors = FALSE)
pc_ortho <- read.delim(pc_ortho_path)

# Transcription Factors
tfs_hg <- filter(read.delim(tfs_hg_path, stringsAsFactors = FALSE), Symbol %in% pc_hg$Symbol)
tfs_mm <- filter(read.delim(tfs_mm_path, stringsAsFactors = FALSE), Symbol %in% pc_mm$Symbol)
tfs_mm <- distinct(tfs_mm, Symbol, .keep_all = TRUE)

# Tables downloaded from supplement of Chu 2021 (all includes external dbs)
lt_chu2021 <- read.delim(chu2021_records_path, stringsAsFactors = FALSE, skip = 1)
lt_all <- read.delim(chu2021_all_path, stringsAsFactors = FALSE, skip = 1)

# On-going curation on googlesheets

if (!file.exists(pavlab_curation_path)) {
  pavlab <- read_sheet(ss = gsheets_curated, sheet = "Master_Curation", trim_ws = TRUE, col_types = "c", range = "A:M")
  write.table(pavlab, quote = FALSE, row.names = FALSE, sep = "\t", file = pavlab_curation_path)
} else {
  pavlab <- read.delim(pavlab_curation_path, stringsAsFactors = FALSE)
}


# Format and join targets curated in Chu2021 and aggregated from external dbs. 
# Note that external interactions are missing details like pubmed ID or species.
# ------------------------------------------------------------------------------


keep_cols <- c("DTRI_ID", 
               "TF_Symbol", 
               "Target_Symbol", 
               "TF_Species", 
               "Target_Species",
               "Cell_Type",
               "Experiment_Type",
               "PubMed_ID",
               "Databases")

# Chu2021: keep only columns of interest and reduce column names

lt_chu2021 <- lt_chu2021 %>% 
  dplyr::rename(
    TF_Symbol = TF_Symbol_Human,
    Target_Symbol = Target_Symbol_Human) %>% 
  dplyr::select(any_of(keep_cols)) %>% 
  mutate(
    TF_Species = str_to_title(TF_Species),
    Target_Species = str_to_title(Target_Species),
    Databases = "Chu2021"
  )


# External dbs aggregated in Chu2021: db name from Current -> Chu2021

lt_all <- lt_all %>%
  dplyr::rename(
    TF_Symbol = TF_Symbol_Human,
    Target_Symbol = Target_Symbol_Human) %>%
  dplyr::select(any_of(keep_cols)) %>% 
  mutate(Databases = str_replace(Databases, "Current", "Chu2021"))


# Join Chu2021 curated targets and those aggregated in Chu2021. Keep only genes
# found in human or mouse Refseq select protein coding; coerce to mouse symbol
# casing if recorded as mouse experiment

lt_all <-
  left_join(lt_all, lt_chu2021, by = "DTRI_ID", suffix = c("", ".y")) %>%
  select(-c(ends_with(".y"), "DTRI_ID")) %>%
  filter(TF_Symbol %in% tfs_hg$Symbol | TF_Symbol %in% tfs_mm$Symbol) %>%
  filter(Target_Symbol %in% pc_hg$Symbol | Target_Symbol %in% pc_mm$Symbol) %>%
  mutate(
    TF_Symbol = ifelse(
      TF_Species == "Mouse" & !is.na(TF_Species),
      str_to_title(TF_Symbol),
      TF_Symbol
    ),
    Target_Symbol = ifelse(
      Target_Species == "Mouse" & !is.na(Target_Species),
      str_to_title(Target_Symbol),
      Target_Symbol
    )
  )


# For pavlab/updated resource, only keep/add relevant cols, and only keep
# protein coding genes. Add species info, which has been encoded by the casing
# of the symbol name. Db name -> "Pavlab" (distinguished from prior Chu2021)
# ------------------------------------------------------------------------------


pavlab <- pavlab %>% 
  dplyr::rename(
    TF_Symbol = TF_Gene_Name, 
    Target_Symbol = Target_Gene_Name
    ) %>% 
  mutate(PubMed_ID = as.integer(PubMed_ID)) %>% 
  dplyr::select(any_of(keep_cols)) %>% 
  filter(TF_Symbol %in% tfs_hg$Symbol | TF_Symbol %in% tfs_mm$Symbol) %>%
  filter(Target_Symbol %in% pc_hg$Symbol | Target_Symbol %in% pc_mm$Symbol) %>% 
  mutate(
    Databases = "Pavlab",
    TF_Species = ifelse(TF_Symbol %in% pc_hg$Symbol, "Human", "Mouse"),
    Target_Species = ifelse(Target_Symbol %in% pc_hg$Symbol, "Human", "Mouse"))



# Combine Chu2021 (+ aggregated external) and updated curation
lt_final <- rbind(lt_all, pavlab)



# Certain low-throughput only:
# TRRUST, Pavlab, Chu2021, TFe, InnateDB
# HTRIdb_LC is only the literature curated part of HTRIdb, same for ORegAnno_LC

#  Uncertain if contains genomic evidence: 
# CytReg (note 2020 update using eY1H, not included here!)


# Remove TCF4 from non Pavlab sources as heavily contaminated with TCF7L2.
# Keep TCF4 itself, autoreg shown in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5243153/
# Also remove TCF4 targets from Pavlab that are actually TCF7L2.
# Also format the perturbation experiment type (absent from external)
# ------------------------------------------------------------------------------


rm_target <- c("lef1", "cd36", "sox9", "glce", "vegfa")


lt_final <- lt_final %>% 
  mutate(Rm = !str_detect(Databases, "Chu2021") & TF_Symbol == "TCF4") %>% 
  mutate(
    Rm = ifelse(str_to_lower(TF_Symbol) == "tcf4" & str_to_lower(Target_Symbol) == "tcf4", FALSE, Rm),
    Experiment_Type = case_when(
      str_detect(str_to_lower(Experiment_Type), ".*perturbation,*") ~ "Perturbation",
      str_detect(str_to_lower(Experiment_Type), ".*binding*") ~ "Binding",
      str_detect(str_to_lower(Experiment_Type), ".*reporter*") ~ "Reporter",
      TRUE ~ "NA"
    )
  ) %>% 
  filter(!Rm) %>%
  filter(!(str_to_lower(TF_Symbol) == "tcf4" & str_to_lower(Target_Symbol) %in% rm_target)) %>% 
  select(-Rm)


# Get an estimate of count of targets per TF, collapsing species/casing

n_target <- lt_final %>% 
  mutate(
    TF_Symbol = str_to_upper(TF_Symbol),
    Target_Symbol = str_to_upper(Target_Symbol)
  ) %>% 
  group_by(TF_Symbol) %>% 
  summarise(N_target = n_distinct(Target_Symbol)) %>% 
  arrange(N_target) %>% 
  mutate(TF_Symbol = factor(TF_Symbol, levels = unique(TF_Symbol)))
  

# Summarize: minimum of 5 targets used for downstream benchmarking

n_target_min <- filter(n_target, N_target >= 5)

n_genes <- lt_final %>% 
  filter(TF_Symbol %in% n_target_min$TF_Symbol) %>% 
  summarise(Distinct_targets = n_distinct(Target_Symbol),
            Distinct_TFs = n_distinct(TF_Symbol))

summ_n <- summary(n_target_min$N_target)


# Inspect genes that are not in the 1:1 ortho set. Note that downstream
# benchmarking uses genes with 1:1 orthologs for either species. Workflow grabs 
# genes regardless of casing AND that have non-equivalent naming. Eg, curated
# table has Tp53, TP53, and Trp53. All affiliated targets will be
# considered for TP53/Trp53

not_ortho <- lt_final %>% 
  filter(Target_Symbol %!in% pc_ortho$Symbol_hg & 
        Target_Symbol %!in% pc_ortho$Symbol_mm) %>% 
  pull(Target_Symbol) %>% 
  unique()


# Save cleaned curation table

write.table(lt_final, 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = curated_all_path)


# Plots
# ------------------------------------------------------------------------------


# Barchart of the count of targets per TF (collapsing species/casing), adding a
# vertical line showing the cut-off of minimum 5 targets

p1 <- ggplot(n_target, aes(x = TF_Symbol, y = N_target)) +
  geom_bar(stat = "identity", fill = "slategrey") +
  geom_vline(xintercept = which(n_target$N_target == 5)[1], col = "firebrick") +
  ylab("Count of distinct curated targets") +
  xlab("TF") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())


ggsave(p1, height = 6, width = 9, device = "png", dpi = 300,
       filename = file.path(plot_dir, "count_of_distinct_curated_targets.png"))
