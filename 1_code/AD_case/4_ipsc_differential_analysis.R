library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

out_dir <- "3_data_analysis/AD_case"

# iPSC samples: two AD and two control donor lines, three culture replicates each
sample_info <- readxl::read_xlsx("2_data/AD_case/ipsc_sample_info.xlsx") %>%
  set_names(c("sampleID", "cell_type", "cell_count", "date", "ID", "group")) %>%
  mutate(ID = as.character(ID))

# Wilcoxon rank-sum test on log2 intensities, AD vs control
for (column in c("hilic", "rp")) {
  column_dir <- file.path(out_dir, "ipsc_processing", column)
  long <- read_csv(file.path(column_dir, "feature_expression.csv")) %>%
    dplyr::select(variable_id, all_of(sample_info$ID)) %>%
    pivot_longer(-variable_id, names_to = "ID", values_to = "intensity") %>%
    left_join(sample_info %>% dplyr::select(ID, group), by = "ID") %>%
    mutate(log2_intensity = log2(intensity + 1))

  diff <- long %>%
    group_by(variable_id) %>%
    summarise(mean_intensity = mean(intensity, na.rm = TRUE),
              condition = mean(log2_intensity[group == "case"], na.rm = TRUE) -
                mean(log2_intensity[group == "control"], na.rm = TRUE),
              pvalue = tryCatch(wilcox.test(log2_intensity[group == "case"], log2_intensity[group == "control"],
                                            exact = FALSE)$p.value,
                                error = function(e) NA_real_),
              .groups = "drop") %>%
    mutate(padj = p.adjust(pvalue, method = "BH"))

  read_csv(file.path(column_dir, "feature_info.csv")) %>%
    mutate(polarity = ifelse(str_detect(variable_id, "_POS$"), "positive", "negative")) %>%
    left_join(diff, by = "variable_id") %>%
    write_csv(file.path(out_dir, paste0("feature_table_ipsc_", column, ".csv")))
}
