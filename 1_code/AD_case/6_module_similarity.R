library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/functional_interpretation/llm_tools.R")

out_dir <- "3_data_analysis/AD_case"

split_terms <- function(x) {
  if (is.na(x) || x == "") return(character())
  terms <- str_squish(str_split(x, ";")[[1]])
  unique(terms[terms != "" & !toupper(terms) %in% c("NA", "NULL", "NAN", "<NA>")])
}
normalize_token <- function(x) str_squish(str_replace_all(str_to_lower(x), "[[:punct:]]+", " "))
jaccard <- function(a, b) if (length(union(a, b)) == 0) 0 else length(intersect(a, b)) / length(union(a, b))

# Accepted, high-confidence Host-Endogenous modules with nominal P < 0.05
modules <- crossing(sample = c("csf", "ipsc"), column = c("hilic", "rp")) %>%
  pmap_dfr(function(sample, column) {
    read_csv(file.path(out_dir, paste0("module_function_", sample, "_", column, ".csv"))) %>%
      mutate(side = ifelse(sample == "csf", "CSF", "iPSC"), method = column,
             Module_Key = paste(ifelse(sample == "csf", "CSF_AD_biotree", "IPSC_AD_NNI"),
                                column, module_id, hypothesis_id, sep = "_"))
  }) %>%
  filter(annotation_mode == "advanced", hypothesis_type == "Host-Endogenous",
         confidence_label == "high", verdict_status == "accepted", p < 0.05) %>%
  mutate(anchor_set = map(anchor_metabolite_ids, split_terms),
         top_set = map(top_metabolite_ids, split_terms),
         class_set = map(class_names, ~ unique(normalize_token(split_terms(.x)))),
         pathway_set = map2(top_kegg_pathways, top_hmdb_pathways,
                            ~ unique(normalize_token(str_remove(c(split_terms(.x), split_terms(.y)), "\\s*\\([^()]*\\)\\s*$")))))

write_csv(modules %>% dplyr::select(-ends_with("_set")), file.path(out_dir, "selected_modules.csv"))

# Same module-module similarity as in PhenoMAP: functional text (title, summary) + composition
title_emb <- openai_embedding(coalesce(modules$title, " "), "text-embedding-3-small")
summary_emb <- openai_embedding(coalesce(modules$functional_summary, " "), "text-embedding-3-small")

csf <- which(modules$side == "CSF")
ipsc <- which(modules$side == "iPSC")

crossing(i = ipsc, j = csf) %>%
  mutate(
    IPSC_Module_Key = modules$Module_Key[i],
    CSF_Module_Key = modules$Module_Key[j],
    title_similarity = map2_dbl(i, j, ~ cosine_sim(title_emb[[.x]], title_emb[[.y]])),
    summary_similarity = map2_dbl(i, j, ~ cosine_sim(summary_emb[[.x]], summary_emb[[.y]])),
    type_score = as.numeric(tolower(modules$hypothesis_type[i]) == tolower(modules$hypothesis_type[j])),
    metabolite_score = 0.6 * map2_dbl(i, j, ~ jaccard(modules$anchor_set[[.x]], modules$anchor_set[[.y]])) +
      0.4 * map2_dbl(i, j, ~ jaccard(modules$top_set[[.x]], modules$top_set[[.y]])),
    Compositional_Score = 0.40 * metabolite_score +
      0.25 * map2_dbl(i, j, ~ jaccard(modules$pathway_set[[.x]], modules$pathway_set[[.y]])) +
      0.25 * map2_dbl(i, j, ~ jaccard(modules$class_set[[.x]], modules$class_set[[.y]])) +
      0.10 * type_score,
    Biological_Score = 0.50 * summary_similarity + 0.40 * title_similarity + 0.10 * type_score,
    Similarity_Score = 0.75 * pmax(Biological_Score, Compositional_Score) +
      0.25 * (Biological_Score + Compositional_Score) / 2
  ) %>%
  dplyr::select(-i, -j) %>%
  write_csv(file.path(out_dir, "module_similarity_pairs.csv"))
