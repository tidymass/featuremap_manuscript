library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/functional_interpretation/llm_tools.R")

out_dir <- "3_data_analysis/functional_interpretation"

titles <- read_csv(file.path(out_dir, "fast_annotation.csv")) %>%
  transmute(module_id = Module, difficulty, gold_pathway_id = pathway_id,
            gold_title = normalize_title(pathway_name), ai_title = normalize_title(ai_title)) %>%
  filter(!is.na(ai_title), ai_title != "")

all_titles <- unique(c(titles$gold_title, titles$ai_title))
embedding <- set_names(openai_embedding(all_titles, "text-embedding-3-small"), all_titles)

# Random control: AI titles shuffled within each setting so that no title keeps its own module
derangement <- function(n, seed) {
  set.seed(seed)
  repeat {
    idx <- sample.int(n)
    if (all(idx != seq_len(n))) return(idx)
  }
}

title_evaluation <- titles %>%
  group_by(difficulty) %>%
  mutate(shuffled_ai_title = ai_title[derangement(n(), 123 + match(first(difficulty), c("easy", "medium", "hard")))]) %>%
  ungroup() %>%
  mutate(title_similarity = map2_dbl(gold_title, ai_title, ~ cosine_sim(embedding[[.x]], embedding[[.y]])),
         shuffled_title_similarity = map2_dbl(gold_title, shuffled_ai_title,
                                              ~ cosine_sim(embedding[[.x]], embedding[[.y]])),
         overlap_class = map2_chr(gold_title, ai_title, word_overlap),
         rougeL_recall = map2_dbl(gold_title, ai_title, rouge_l_recall))

write_csv(title_evaluation, file.path(out_dir, "title_evaluation.csv"))
