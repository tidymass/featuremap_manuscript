library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(featuremap)

out_dir <- "3_data_analysis/AD_case"

datasets <- crossing(sample = c("csf", "ipsc"), column = c("hilic", "rp")) %>%
  mutate(biospecimen = ifelse(sample == "csf", "Cerebrospinal fluid", "induced pluripotent stem cell"))

# HILIC and RP analysed separately; modules with nominal P < 0.05 interpreted in advanced mode
pwalk(datasets, function(sample, column, biospecimen) {
  name <- paste(sample, column, sep = "_")
  message(name)
  feature_table <- read_csv(file.path(out_dir, paste0("feature_table_", name, ".csv")))

  res <- annotate_feature(feature_table = feature_table, ms2_table = NULL, column = column,
                          ms1_match_ppm = 15, mfc_rt_tol = 5, isotope_number = 3, workers = 6)
  res <- extract_subnetwork(res = res, feature = "all", workers = 6, times = 100, b = 20, w = 25, mu = 1e-4)
  res <- identify_module(res = res, method = "louvain", louvain_resolution = 1, min_size = 3,
                         n_perm = 100, seed = 123, workers = 4, node_freq = 0.05)

  res_adv <- annotate_module(
    res = res,
    p = "p",
    pcutoff = 0.05,
    mode = "advanced",
    species = "human",
    biospecimen = biospecimen,
    phenotype = "Alzheimer's disease",
    pubmed = TRUE,
    NCBI_api_key = Sys.getenv("NCBI_API_KEY"),
    debate_html = TRUE,
    html_dir = file.path(out_dir, "advanced_html", name),
    synthesize = TRUE,
    synthesis_model = "gpt-4.1",
    provider = "openai",
    llm_model = "gpt-4.1",
    llm_api_key = Sys.getenv("OPENAI_API_KEY")
  )
  save(res, res_adv, file = file.path(out_dir, paste0("featuremap_result_", name, ".rda")))

  res_adv$module_function %>%
    as.data.frame() %>%
    left_join(res_adv$module_score, by = c("module_id" = "Module")) %>%
    write_csv(file.path(out_dir, paste0("module_function_", name, ".csv")))
})
