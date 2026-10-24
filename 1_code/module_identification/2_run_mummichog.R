library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

out_dir <- "3_data_analysis/module_identification/arginine_biosynthesis"

# mummichog 2.7.0, run separately for positive and negative features
run_mummichog <- function(input_file, run_out) {
  feature <- read_tsv(input_file) %>%
    filter(!is.na(`p-value`), !is.na(`t-score`))
  for (mode in c("pos", "neg")) {
    mode_dir <- file.path(run_out, mode)
    dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)
    feature %>%
      filter(grepl("_NEG", toupper(custom_id)) == (mode == "neg")) %>%
      write_tsv(file.path(mode_dir, "input.txt"))
    withr::with_dir(mode_dir, system2("mummichog", c(
      "-f", "input.txt", "-o", "result", "-c", "0.05", "-u", "5", "-z", "0", "-p", "100",
      "-m", ifelse(mode == "pos", "positive", "negative")
    )))
    result_file <- list.files(mode_dir, pattern = "modularanalysis.*\\.tsv$",
                              recursive = TRUE, full.names = TRUE)[1]
    file.copy(result_file, file.path(mode_dir, "module_analysis.tsv"), overwrite = TRUE)
  }
}

walk(1:100, function(i) {
  run_mummichog(
    file.path(out_dir, "mummichog_input", paste0("feature_mummichog_", i, ".txt")),
    file.path(out_dir, "mummichog_result", paste0("feature_mummichog_", i))
  )
})
