## Organize curated low-throughput TF-target interactions from Chu 2021, which
## also aggregates from other low throughput databases, and from on going Pavlab
## curation on google sheets.
## -----------------------------------------------------------------------------

library(tidyverse)
library(googlesheets4)
source("R/00_config.R")

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



lt_all <- lt_all %>%
  dplyr::rename(
    TF_Symbol = TF_Symbol_Human,
    Target_Symbol = Target_Symbol_Human) %>%
  dplyr::select(any_of(keep_cols)) %>% 
  mutate(Databases = str_replace(Databases, "Current", "Chu2021"))



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
# of the symbol name
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



lt_final <- rbind(lt_all, pavlab)



# Keep:
# TRRUST, Pavlab, Chu2021, TFe, InnateDB
# TFactS may have some non-low throughput
# HTRIdb_LC is only the literature curated part of HTRIdb
sort(table(lt_final$Databases), decreasing = TRUE)
unique(unlist(str_split(lt_final$Databases, ", ")))



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



n_target <- lt_final %>% 
  mutate(
    TF_Symbol = str_to_upper(TF_Symbol),
    Target_Symbol = str_to_upper(Target_Symbol)
  ) %>% 
  group_by(TF_Symbol) %>% 
  summarise(N_target = n_distinct(Target_Symbol)) %>% 
  arrange(N_target) %>% 
  mutate(TF_Symbol = factor(TF_Symbol, levels = unique(TF_Symbol)))
  


#
# ------------------------------------------------------------------------------


ggplot(n_target, aes(x = TF_Symbol, y = N_target)) +
  geom_bar(stat = "identity") +
  ylab("Count of distinct curated targets") +
  xlab("TF") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 20),
        axis.ticks.x = element_blank())


n_distinct(str_to_upper(lt_final$Target_Symbol))
n_distinct(str_to_upper(lt_final$TF_Symbol))


# table(n_target$TF_Symbol %in% c(tfs_hg$Symbol, str_to_upper(tfs_mm$Symbol)))
# non_tfs <- setdiff(n_target$TF_Symbol, c(tfs_hg$Symbol, str_to_upper(tfs_mm$Symbol)))
# filter(n_target, TF_Symbol %in% non_tfs) %>% view
# sort(table(filter(lt_final, str_to_upper(TF_Symbol) == "CTNNB1")$Databases))

sort(table(filter(lt_final, str_to_upper(TF_Symbol) == "SP1")$Databases))


write.table(lt_final, 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = curated_all_path)
