library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/phenomap_analysis/phenomap_tools.R")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Inputs exported from the PhenoMAP knowledge graph: high-confidence Host-Endogenous
# modules mapped to phenotypes, and module-module similarity edges
phenotype_similarity <- function(module_map, edges) {
  sets <- split(module_map$Module_ID, module_map$Phenotype_ID) %>%
    map(unique) %>%
    keep(~ length(.x) >= 2)
  S <- module_similarity_matrix(edges, unique(unlist(sets)))
  ids <- sort(names(sets))
  M <- outer(ids, ids, Vectorize(function(i, j) bma(sets[[i]], sets[[j]], S)))
  dimnames(M) <- list(ids, ids)
  M
}

write_matrix <- function(M, file) {
  as.data.frame(M) %>%
    rownames_to_column("Phenotype_ID") %>%
    write_csv(file)
}

module_map <- read_csv(file.path(data_dir, "selected_module_phenotype_map.csv"))

module_map %>%
  group_by(Phenotype_ID, Phenotype, Phenotype_Type) %>%
  summarise(Module_N = n_distinct(Module_ID), .groups = "drop") %>%
  write_csv(file.path(out_dir, "phenotype_info.csv"))

# All phenotypes with >= 2 modules
phenotype_similarity(module_map,
                     read_csv(file.path(data_dir, "selected_module_module_similarity_edges.csv"))) %>%
  write_matrix(file.path(out_dir, "phenotype_similarity_matrix.csv"))

# Disease phenotypes with >= 2 modules (for the CTD comparison)
phenotype_similarity(read_csv(file.path(data_dir, "selected_disease_module_phenotype_map.csv")),
                     read_csv(file.path(data_dir, "selected_disease_module_module_similarity_edges.csv"))) %>%
  write_matrix(file.path(out_dir, "disease_similarity_matrix.csv"))
