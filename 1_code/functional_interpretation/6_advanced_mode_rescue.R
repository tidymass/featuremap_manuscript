library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/functional_interpretation/llm_tools.R")
library(featuremap)

out_dir <- "3_data_analysis/functional_interpretation"
load(file.path(out_dir, "semisynthetic_modules.rda"))

# Weak fast-mode titles: no word overlap, or partial overlap with below-median ROUGE-L recall
candidates <- read_csv(file.path(out_dir, "title_evaluation.csv")) %>%
  group_by(difficulty) %>%
  mutate(partial_median = median(rougeL_recall[overlap_class == "Partial word overlap"])) %>%
  ungroup() %>%
  mutate(selection_reason = case_when(
    overlap_class == "No word overlap" ~ "no_word_overlap",
    overlap_class == "Partial word overlap" & rougeL_recall < partial_median ~ "partial_word_overlap_low_rougeL"
  )) %>%
  filter(!is.na(selection_reason))

run_advanced <- function(module_ids) {
  res <- semisynthetic_res
  res$module_table <- filter(res$module_table, Module %in% module_ids)
  res$module_score <- filter(res$module_score, Module %in% module_ids)
  annotate_module(res = res, p = "fdr", pcutoff = 0.05, mode = "advanced", species = "human",
                  biospecimen = NULL, phenotype = NULL, pubmed = FALSE, provider = "openai",
                  llm_model = "gpt-4.1-mini", llm_api_key = Sys.getenv("OPENAI_API_KEY"),
                  synthesize = FALSE, debate_html = FALSE)
}

res_advanced <- map(c(easy = "easy", hard = "hard"),
                    ~ run_advanced(candidates$module_id[candidates$difficulty == .x]))
save(res_advanced, file = file.path(out_dir, "advanced_mode_result.rda"))

# Rescued: no overlap -> partial/exact overlap, or higher ROUGE-L recall for partial overlap
advanced_comparison <- candidates %>%
  left_join(map_dfr(res_advanced, ~ .x$module_function %>%
                      transmute(module_id, advanced_title = title, advanced_summary = functional_summary)),
            by = "module_id") %>%
  mutate(advanced_title_norm = str_to_lower(advanced_title),
         advanced_overlap_class = map2_chr(gold_title, advanced_title_norm, word_overlap),
         advanced_rougeL_recall = map2_dbl(gold_title, advanced_title_norm, rouge_l_recall),
         title_rescued = (overlap_class == "No word overlap" & advanced_overlap_class != "No word overlap") |
           (overlap_class == "Partial word overlap" & advanced_rougeL_recall > rougeL_recall))
write_csv(advanced_comparison, file.path(out_dir, "advanced_comparison.csv"))

# LLM judge of fast- and advanced-mode summaries for rescued modules
rescued <- advanced_comparison %>%
  filter(title_rescued) %>%
  left_join(read_csv(file.path(out_dir, "pathway_reference.csv")) %>%
              dplyr::select(gold_pathway_id = SMPDB_ID, Description),
            by = "gold_pathway_id")

rescued %>%
  dplyr::select(module_id, difficulty) %>%
  left_join(read_csv(file.path(out_dir, "summary_evaluation.csv")) %>%
              dplyr::select(module_id, fast_judge_overall = judge_overall),
            by = "module_id") %>%
  mutate(advanced_judge_overall = map2_dfr(rescued$Description, rescued$advanced_summary,
                                           judge_summary)$judge_overall) %>%
  write_csv(file.path(out_dir, "rescued_summary_judge.csv"))
