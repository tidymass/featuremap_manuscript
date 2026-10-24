library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/functional_interpretation/llm_tools.R")

out_dir <- "3_data_analysis/functional_interpretation"

fast <- read_csv(file.path(out_dir, "fast_annotation.csv")) %>%
  left_join(read_csv(file.path(out_dir, "pathway_reference.csv")) %>%
              dplyr::select(pathway_id = SMPDB_ID, Description),
            by = "pathway_id")

# Embedding similarity between generated summary and PathBank description
embedding <- openai_embedding(c(fast$Description, fast$ai_summary), "text-embedding-3-large")
n <- nrow(fast)

summary_evaluation <- fast %>%
  transmute(module_id = Module, difficulty,
            score_openai_embedding = map_dbl(seq_len(n), ~ cosine_sim(embedding[[.x]], embedding[[n + .x]]))) %>%
  bind_cols(map2_dfr(fast$Description, fast$ai_summary, judge_summary))

write_csv(summary_evaluation, file.path(out_dir, "summary_evaluation.csv"))
