library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(tidygraph)

case_ids <- c("MTBLS13267", "MTBLS267", "MTBLS13855", "MTBLS13095", "MTBLS12329",
              "MTBLS13578", "MTBLS11987", "MTBLS417", "MTBLS5028", "MTBLS12854")

data_dir <- "2_data/ms1_based_annotation"
out_dir <- "3_data_analysis/ms1_based_annotation"

load("2_data/metabolic_reaction_network.rda")
mrnref <- metabolic_reaction_network %>%
  activate(nodes) %>%
  as_tibble() %>%
  transmute(HMDB.ID = normalize_hmdb(HMDB_ID), KEGG.ID = normalize_kegg(KEGG_ID))

# Top-k candidates (ties kept together) are compared with the top MS/MS reference
evaluate_tool <- function(pred, ref, score_col, top) {
  id_cols <- intersect(c("HMDB.ID", "KEGG.ID"), colnames(pred))
  hits <- map_dfr(unique(ref$variable_id), function(id) {
    cand <- pred %>%
      filter(variable_id == id) %>%
      arrange(desc(.data[[score_col]])) %>%
      distinct(across(all_of(id_cols)), .keep_all = TRUE)
    if (nrow(cand) == 0) return(NULL)
    top_scores <- head(sort(unique(cand[[score_col]]), decreasing = TRUE), top)
    cand <- cand %>% filter(.data[[score_col]] %in% top_scores)
    truth <- ref %>%
      filter(variable_id == id) %>%
      arrange(desc(Total.score)) %>%
      slice(1)
    hit <- any(map_lgl(id_cols, ~ any(na.omit(cand[[.x]]) %in% na.omit(truth[[.x]]))))
    tibble(hit = hit, n = nrow(cand))
  })
  tibble(
    n_feature = n_distinct(ref$variable_id),
    efficient_score = sum(hits$hit / hits$n) / n_feature,
    precision = sum(hits$hit) / sum(hits$n),
    F1 = 2 * sum(hits$hit) / (2 * sum(hits$hit) + sum(hits$n - hits$hit) + sum(!hits$hit))
  )
}

benchmark_case <- function(case_id) {
  message(case_id)
  res <- load_rda_first(file.path(data_dir, case_id, "res.rda"))
  anno1 <- normalize_id_cols(res$anno_feature_table1)
  featuremap <- normalize_id_cols(res$anno_feature_table2)

  tool_dir <- file.path(out_dir, "tool_annotation", case_id)
  mwise <- read_csv(file.path(tool_dir, "mWISE_anno.csv")) %>%
    transmute(variable_id, KEGG.ID = normalize_kegg(KEGG_ID), score = -Ranking) %>%
    filter(!is.na(KEGG.ID))
  xms <- read_csv(file.path(tool_dir, "xms_anno.csv")) %>%
    transmute(variable_id, HMDB.ID = normalize_hmdb(chemical_ID), score = Confidence) %>%
    filter(!is.na(HMDB.ID))

  ms2 <- load_rda_first(file.path(data_dir, case_id, "annotated_table_ms2.rda")) %>%
    prepare_ms2_reference()

  ref_featuremap <- ms2 %>%
    filter(HMDB.ID %in% mrnref$HMDB.ID | KEGG.ID %in% mrnref$KEGG.ID,
           HMDB.ID %in% anno1$HMDB.ID | KEGG.ID %in% anno1$KEGG.ID,
           variable_id %in% featuremap$variable_id)
  ref_mwise <- ms2 %>% filter(!is.na(KEGG.ID), variable_id %in% mwise$variable_id)
  ref_xms <- ms2 %>% filter(!is.na(HMDB.ID), variable_id %in% xms$variable_id)

  # Common set: MS/MS-supported features evaluable by all three tools
  common <- Reduce(intersect, list(ref_featuremap$variable_id,
                                   ref_mwise$variable_id,
                                   ref_xms$variable_id))

  metrics <- map_dfr(1:3, function(top) {
    bind_rows(
      evaluate_tool(featuremap, filter(ref_featuremap, variable_id %in% common),
                    "weighted_anno_score", top) %>% mutate(tool = "FeatureMAP"),
      evaluate_tool(mwise, filter(ref_mwise, variable_id %in% common),
                    "score", top) %>% mutate(tool = "mWISE"),
      evaluate_tool(xms, filter(ref_xms, variable_id %in% common),
                    "score", top) %>% mutate(tool = "xMSannotator")
    ) %>%
      mutate(case_id = case_id, top = top, .before = 1)
  })

  # Redundancy: mean candidates per feature in each tool's full output
  candidates <- list(FeatureMAP = featuremap, mWISE = mwise, xMSannotator = xms) %>%
    imap_dfr(function(pred, tool) {
      id_cols <- intersect(c("HMDB.ID", "KEGG.ID"), colnames(pred))
      pred %>%
        filter(if_any(all_of(id_cols), ~ !is.na(.x))) %>%
        distinct(across(c(variable_id, all_of(id_cols)))) %>%
        count(variable_id) %>%
        summarise(case_id = case_id, tool = tool, mean_candidates_per_feature = mean(n))
    })

  list(tool_benchmark_metrics = metrics,
       tool_candidates_per_feature = candidates)
}

result <- map(case_ids, benchmark_case)

for (name in names(result[[1]])) {
  map_dfr(result, name) %>%
    write_csv(file.path(out_dir, paste0(name, ".csv")))
}
