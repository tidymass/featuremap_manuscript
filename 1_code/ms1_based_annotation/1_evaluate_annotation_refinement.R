library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(tidygraph)

case_ids <- c("MTBLS13267", "MTBLS267", "MTBLS13855", "MTBLS13095", "MTBLS12329",
              "MTBLS13578", "MTBLS11987", "MTBLS417", "MTBLS5028", "MTBLS12854")

data_dir <- "2_data/ms1_based_annotation"
out_dir <- "3_data_analysis/ms1_based_annotation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load("2_data/metabolic_reaction_network.rda")
mrnref <- metabolic_reaction_network %>%
  activate(nodes) %>%
  as_tibble() %>%
  dplyr::select(HMDB_ID, KEGG_ID)

# Tie-adjusted hit: 1/n if the MS/MS reference is among the n tied top-k candidates
hit_score <- function(id, pred, ref, score_col, top) {
  cand <- pred %>%
    filter(variable_id == id) %>%
    arrange(desc(.data[[score_col]])) %>%
    distinct(HMDB.ID, KEGG.ID, .keep_all = TRUE)
  if (nrow(cand) == 0) return(0)
  top_scores <- head(sort(unique(cand[[score_col]]), decreasing = TRUE), top)
  cand <- cand %>% filter(.data[[score_col]] %in% top_scores)
  truth <- ref %>%
    filter(variable_id == id) %>%
    arrange(desc(Total.score)) %>%
    slice(1)
  hit <- any(na.omit(cand$HMDB.ID) %in% na.omit(truth$HMDB.ID)) |
    any(na.omit(cand$KEGG.ID) %in% na.omit(truth$KEGG.ID))
  ifelse(hit, 1 / nrow(cand), 0)
}

tie_adjusted_accuracy <- function(pred, ref, score_col, top) {
  ids <- unique(ref$variable_id)
  mean(map_dbl(ids, hit_score, pred = pred, ref = ref, score_col = score_col, top = top))
}

evaluate_case <- function(case_id) {
  message(case_id)
  res <- load_rda_first(file.path(data_dir, case_id, "res.rda"))
  anno1 <- res$anno_feature_table1
  anno2 <- res$anno_feature_table2

  ref <- load_rda_first(file.path(data_dir, case_id, "annotated_table_ms2.rda")) %>%
    prepare_ms2_reference(normalize_id = FALSE) %>%
    filter(variable_id %in% anno1$variable_id,
           variable_id %in% anno2$variable_id,
           HMDB.ID %in% mrnref$HMDB_ID | KEGG.ID %in% mrnref$KEGG_ID,
           HMDB.ID %in% anno1$HMDB.ID | KEGG.ID %in% anno1$KEGG.ID)

  accuracy <- map_dfr(1:3, function(top) {
    tibble(
      case_id = case_id,
      top = top,
      ms1_based = tie_adjusted_accuracy(anno1, ref, "score_norm", top),
      network_refined = tie_adjusted_accuracy(anno2, ref, "score_norm", top),
      network_weighted = tie_adjusted_accuracy(anno2, ref, "weighted_anno_score", top)
    )
  })

  # Positive pair: candidate matching the top-scoring MS/MS reference of a feature
  positive_key <- ref %>%
    group_by(variable_id) %>%
    slice_max(Total.score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(key = paste(variable_id, coalesce(HMDB.ID, KEGG.ID), sep = "_")) %>%
    pull(key)

  anno2_ref <- anno2 %>%
    filter(variable_id %in% ref$variable_id) %>%
    mutate(key = paste(variable_id, coalesce(HMDB.ID, KEGG.ID), sep = "_"),
           type = ifelse(key %in% positive_key, "positive", "negative"))

  edge_frequency <- anno2_ref %>%
    group_by(type) %>%
    summarise(n = n(),
              mean_freq = mean(edge_freq, na.rm = TRUE),
              se_freq = sd(edge_freq, na.rm = TRUE) / sqrt(n()),
              .groups = "drop") %>%
    mutate(case_id = case_id, .before = 1)

  positive_pair_rank <- anno2_ref %>%
    filter(type == "positive") %>%
    transmute(case_id = case_id, variable_id, key, rank = confidence_rank)

  redundancy <- tibble(
    case_id = case_id,
    stage = c("ms1_based", "network_refined"),
    redundancy = c(nrow(anno1) / n_distinct(anno1$variable_id),
                   nrow(anno2) / n_distinct(anno2$variable_id))
  )

  list(annotation_accuracy = accuracy,
       edge_frequency = edge_frequency,
       positive_pair_rank = positive_pair_rank,
       annotation_redundancy = redundancy)
}

result <- map(case_ids, evaluate_case)

for (name in names(result[[1]])) {
  map_dfr(result, name) %>%
    write_csv(file.path(out_dir, paste0(name, ".csv")))
}
