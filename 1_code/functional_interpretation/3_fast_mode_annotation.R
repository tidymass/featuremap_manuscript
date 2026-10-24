library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(featuremap)

out_dir <- "3_data_analysis/functional_interpretation"
load(file.path(out_dir, "semisynthetic_modules.rda"))

# Fast mode without phenotype or literature input
res_fast <- annotate_module(
  res = semisynthetic_res,
  p = "fdr",
  pcutoff = 0.05,
  mode = "fast",
  species = "human",
  biospecimen = NULL,
  phenotype = NULL,
  pubmed = FALSE,
  provider = "openai",
  llm_model = "gpt-4.1-mini",
  llm_api_key = Sys.getenv("OPENAI_API_KEY")
)
save(res_fast, file = file.path(out_dir, "fast_mode_result.rda"))

module_gold %>%
  left_join(res_fast$module_function %>%
              transmute(Module = module_id, ai_type = hypothesis_type,
                        ai_title = title, ai_summary = functional_summary),
            by = "Module") %>%
  write_csv(file.path(out_dir, "fast_annotation.csv"))
