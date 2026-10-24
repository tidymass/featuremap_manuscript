library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

data_dir <- "2_data/AD_case"
out_dir <- "3_data_analysis/AD_case"

matched <- read_csv(file.path(out_dir, "csf_matched_samples.csv")) %>%
  mutate(ID = as.character(ID))

# CSF feature table from LC-MS processing (HILIC = Amide, RP = C18), with MS/MS annotations (BiotreeDB)
features <- readxl::read_xlsx(file.path(data_dir, "csf_feature_table.xlsx")) %>%
  mutate(variable_id = paste0("M", gsub("\\.", "", sprintf("%.4f", mz)),
                              "T", gsub("\\.", "", sprintf("%.1f", rt)), "_", toupper(type)),
         polarity = ifelse(toupper(type) == "POS", "positive", "negative")) %>%
  distinct(variable_id, .keep_all = TRUE)

# Paired Wilcoxon signed-rank test on log2 intensities of matched AD-control pairs
paired_test <- function(df) {
  long <- df %>%
    dplyr::select(variable_id, all_of(matched$ID)) %>%
    pivot_longer(-variable_id, names_to = "ID", values_to = "intensity") %>%
    left_join(matched %>% dplyr::select(ID, group, match_set), by = "ID") %>%
    mutate(log2_intensity = log2(intensity + 1))

  diff <- long %>%
    dplyr::select(variable_id, match_set, group, log2_intensity) %>%
    pivot_wider(names_from = group, values_from = log2_intensity) %>%
    group_by(variable_id) %>%
    summarise(condition = mean(case - control, na.rm = TRUE),
              pvalue = tryCatch(wilcox.test(case, control, paired = TRUE, exact = FALSE)$p.value,
                                error = function(e) NA_real_),
              .groups = "drop") %>%
    mutate(padj = p.adjust(pvalue, method = "BH"))

  df %>%
    transmute(variable_id, mz, rt, polarity) %>%
    left_join(long %>% group_by(variable_id) %>% summarise(mean_intensity = mean(intensity, na.rm = TRUE)),
              by = "variable_id") %>%
    left_join(diff, by = "variable_id")
}

for (column in c("hilic", "rp")) {
  column_features <- features %>% filter(type2 == ifelse(column == "hilic", "Amide", "RP"))
  write_csv(paired_test(column_features), file.path(out_dir, paste0("feature_table_csf_", column, ".csv")))

  # Level 1/2 MS/MS annotations
  column_features %>%
    filter(level %in% c("level1", "level2")) %>%
    transmute(variable_id, HMDB.ID = HMDB, KEGG.ID = `KEGG COMPOUND ID`, INCHIKEY, level) %>%
    write_csv(file.path(out_dir, paste0("ms2_table_csf_", column, ".csv")))
}
