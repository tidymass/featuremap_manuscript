library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/phenomap_analysis/phenomap_tools.R")

disease_similarity <- read_similarity(file.path(out_dir, "disease_similarity_matrix.csv"))

# CTD disease-pathway associations: http://ctdbase.org/reports/CTD_diseases_pathways.tsv.gz
ctd_pathway <- read_tsv(file.path(data_dir, "CTD_diseases_pathways.tsv.gz"), comment = "#",
                        col_names = c("DiseaseName", "DiseaseID", "PathwayName", "PathwayID", "InferenceGeneSymbol")) %>%
  distinct(DiseaseID, PathwayID)
ctd_kegg <- ctd_pathway %>%
  filter(str_detect(PathwayID, regex("^KEGG:", ignore_case = TRUE))) %>%
  split(.$DiseaseID) %>%
  map("PathwayID")

# Disease phenotypes mapped to CTD MEDIC IDs (name / ontology ID matching, manually reviewed)
ctd_mapping <- read_csv(file.path(data_dir, "disease_ctd_mapping.csv")) %>%
  filter(Use_For_CTD_Validation, !is.na(CTD_Disease_ID_Final),
         Phenotype_ID %in% rownames(disease_similarity),
         CTD_Disease_ID_Final %in% ctd_pathway$DiseaseID)
ctd_id <- set_names(ctd_mapping$CTD_Disease_ID_Final, ctd_mapping$Phenotype_ID)

# Shared mechanism: number of shared CTD KEGG pathways above the mean across pairs
combn(ctd_mapping$Phenotype_ID, 2) %>%
  t() %>%
  as_tibble(.name_repair = ~ c("Phenotype_ID_1", "Phenotype_ID_2")) %>%
  filter(ctd_id[Phenotype_ID_1] != ctd_id[Phenotype_ID_2]) %>%
  mutate(BMA_Similarity = disease_similarity[cbind(Phenotype_ID_1, Phenotype_ID_2)],
         Shared_KEGG_N = map2_int(ctd_id[Phenotype_ID_1], ctd_id[Phenotype_ID_2],
                                  ~ length(intersect(ctd_kegg[[.x]], ctd_kegg[[.y]]))),
         Shared_Mechanism = Shared_KEGG_N > mean(Shared_KEGG_N)) %>%
  write_csv(file.path(out_dir, "disease_ctd_pairs.csv"))
