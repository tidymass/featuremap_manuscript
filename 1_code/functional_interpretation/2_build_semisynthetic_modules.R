library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(tidygraph)
library(igraph)

out_dir <- "3_data_analysis/functional_interpretation"

load("2_data/metabolic_reaction_network.rda")
graph <- as.igraph(metabolic_reaction_network)
nodes <- metabolic_reaction_network %>%
  activate(nodes) %>%
  as_tibble() %>%
  mutate(vid = row_number(),
         HMDB_ID = str_extract(HMDB_ID, "HMDB[0-9]{5,}"),
         KEGG_ID = str_extract(KEGG_ID, "C[0-9]{5}"))

# HMDB IDs within 1-5 reaction steps of the given vertices
within_steps <- function(vids) {
  d <- distances(graph, v = unique(vids))
  nodes$HMDB_ID[apply(d, 2, min) %in% 1:5]
}

sample_ids <- function(x, n) {
  x <- unique(na.omit(x))
  if (length(x) <= n) x else sample(x, n)
}

# One module = pathway metabolites + equal numbers of neighbor and noise metabolites
build_module <- function(pathway, difficulty, index) {
  pathway_ids <- intersect(str_split(pathway$mrn_overlap_hmdb_ids, "; ")[[1]], nodes$HMDB_ID)
  n <- length(pathway_ids)
  supporting_ids <- sample_ids(setdiff(within_steps(nodes$vid[nodes$HMDB_ID %in% pathway_ids]),
                                       pathway_ids), n)
  noise_ids <- sample_ids(setdiff(within_steps(nodes$vid[nodes$HMDB_ID %in% supporting_ids]),
                                  c(pathway_ids, supporting_ids)), n)

  # Easy: noise weighted low; hard: noise weighted like pathway metabolites
  noise_range <- if (difficulty == "easy") c(0.4, 0.6) else c(0.8, 1)
  weight <- function(role) {
    case_when(role == "pathway" ~ runif(length(role), 0.8, 1),
              role == "supporting" ~ runif(length(role), 0.6, 0.8),
              role == "noise" ~ runif(length(role), noise_range[1], noise_range[2]),
              TRUE ~ runif(length(role), 0.2, 1))
  }
  module_id <- sprintf("HMDB_SIM_%03d_%s", index, toupper(substr(difficulty, 1, 1)))

  module_table <- nodes %>%
    filter(HMDB_ID %in% c(pathway_ids, supporting_ids, noise_ids)) %>%
    distinct(HMDB_ID, .keep_all = TRUE) %>%
    mutate(role = case_when(HMDB_ID %in% pathway_ids ~ "pathway",
                            HMDB_ID %in% supporting_ids ~ "supporting",
                            TRUE ~ "noise"),
           Importance = weight(role),
           Node_freq = weight(role)) %>%
    arrange(factor(role, levels = c("pathway", "supporting", "noise")),
            desc(0.7 * Importance + 0.3 * Node_freq)) %>%
    transmute(Module = module_id, INCHI_ID, HMDB_ID, KEGG_ID, Compound_name,
              Node_freq = round(Node_freq, 4), Type = "annotated",
              Importance = round(Importance, 4), Use_for_annotation = TRUE)

  list(module_table = module_table,
       module_score = tibble(Module = module_id, Module_size = nrow(module_table),
                             Annotated_nodes = nrow(module_table),
                             Module_score = round(sum(module_table$Importance), 4),
                             p = 0.001, fdr = 0.01),
       gold = tibble(Module = module_id, pathway_id = pathway$SMPDB_ID, pathway_name = pathway$Name,
                     difficulty = difficulty, n_pathway = n,
                     n_supporting = length(supporting_ids), n_noise = length(noise_ids)))
}

pathways <- read_csv(file.path(out_dir, "pathway_reference.csv")) %>%
  filter(!is.na(mrn_overlap_hmdb_ids)) %>%
  arrange(SMPDB_ID)

set.seed(123)
modules <- crossing(i = seq_len(nrow(pathways)), difficulty = c("easy", "hard")) %>%
  mutate(index = row_number()) %>%
  pmap(function(i, difficulty, index) build_module(pathways[i, ], difficulty, index))

semisynthetic_res <- list(object_type = "featuremap_result",
                          stage = "module_identified",
                          module_table = map_dfr(modules, "module_table"),
                          module_score = map_dfr(modules, "module_score"))
module_gold <- map_dfr(modules, "gold")

save(semisynthetic_res, module_gold, file = file.path(out_dir, "semisynthetic_modules.rda"))
write_csv(module_gold, file.path(out_dir, "module_gold.csv"))
