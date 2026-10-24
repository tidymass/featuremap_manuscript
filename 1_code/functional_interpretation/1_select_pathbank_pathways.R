library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(tidygraph)

data_dir <- "2_data/functional_interpretation"
out_dir <- "3_data_analysis/functional_interpretation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# PathBank downloads: https://pathbank.org/downloads (all pathways / all metabolites)
clean_names <- function(df) rename_with(df, ~ str_replace_all(.x, "[ -]", "_"))
pathways <- read_csv(file.path(data_dir, "pathbank_pathways.csv")) %>% clean_names()
metabolites <- read_csv(file.path(data_dir, "pathbank_metabolites.csv")) %>%
  clean_names() %>%
  filter(str_detect(Species, regex("Homo sapiens", ignore_case = TRUE)))

load("2_data/metabolic_reaction_network.rda")
mrn_hmdb <- metabolic_reaction_network %>%
  activate(nodes) %>%
  pull(HMDB_ID) %>%
  na.omit() %>%
  unique()

# Human metabolic pathways with > 50% of HMDB metabolites mapped to the MRN
pathway_reference <- pathways %>%
  filter(SMPDB_ID %in% metabolites$PathBank_ID, Subject == "Metabolic", nchar(Description) > 10) %>%
  inner_join(metabolites %>%
               filter(!is.na(HMDB_ID)) %>%
               distinct(SMPDB_ID = PathBank_ID, HMDB_ID),
             by = "SMPDB_ID") %>%
  group_by(SMPDB_ID, PW_ID, Name, Subject, Description) %>%
  summarise(pathway_hmdb_n = n_distinct(HMDB_ID),
            mrn_overlap_n = n_distinct(HMDB_ID[HMDB_ID %in% mrn_hmdb]),
            mrn_overlap_hmdb_ids = paste(sort(unique(HMDB_ID[HMDB_ID %in% mrn_hmdb])), collapse = "; "),
            .groups = "drop") %>%
  mutate(coverage = mrn_overlap_n / pathway_hmdb_n) %>%
  filter(coverage > 0.5) %>%
  # Merge lipid-species variants of the same pathway, keep the best-covered one
  mutate(Name = str_trim(str_remove(Name, "\\s+(TG|DG|MG|PC|PE|PS|PI|PG|CL|SM|Cer|LPC|LPE|LPI|LPS|LPG)\\(.*\\)$"))) %>%
  group_by(Name) %>%
  arrange(desc(coverage), desc(mrn_overlap_n), desc(pathway_hmdb_n), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

write_csv(pathway_reference, file.path(out_dir, "pathway_reference.csv"))
