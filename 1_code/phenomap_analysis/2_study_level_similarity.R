library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/phenomap_analysis/phenomap_tools.R")

# Unit = one study-phenotype pair, represented by its module set;
# phenotypes investigated in >= 2 independent studies
units <- read_csv(file.path(data_dir, "selected_module_phenotype_map.csv")) %>%
  group_by(Phenotype_ID, Phenotype, Study_ID) %>%
  summarise(modules = list(sort(unique(Module_ID))), .groups = "drop") %>%
  group_by(Phenotype_ID) %>%
  filter(n_distinct(Study_ID) >= 2) %>%
  ungroup() %>%
  arrange(Phenotype_ID, Study_ID)

S <- module_similarity_matrix(read_csv(file.path(data_dir, "selected_module_module_similarity_edges.csv")),
                              unique(unlist(units$modules)))

# Same-phenotype pairs come from different studies; different-phenotype pairs
# from the same study are excluded
combn(nrow(units), 2) %>%
  t() %>%
  as_tibble(.name_repair = ~ c("i", "j")) %>%
  transmute(Phenotype_ID_1 = units$Phenotype_ID[i], Study_ID_1 = units$Study_ID[i],
            Phenotype_ID_2 = units$Phenotype_ID[j], Study_ID_2 = units$Study_ID[j],
            Same_Phenotype = Phenotype_ID_1 == Phenotype_ID_2,
            Study_Level_BMA_Similarity = map2_dbl(i, j, ~ bma(units$modules[[.x]], units$modules[[.y]], S))) %>%
  filter(Same_Phenotype | Study_ID_1 != Study_ID_2) %>%
  write_csv(file.path(out_dir, "study_level_similarity_pairs.csv"))
