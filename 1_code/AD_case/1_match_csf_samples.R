library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

data_dir <- "2_data/AD_case"
out_dir <- "3_data_analysis/AD_case"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# CSF cohort: ID, group (case = AD, control = non-AD), sex, age
sample_info <- readxl::read_xlsx(file.path(data_dir, "csf_sample_info.xlsx")) %>%
  mutate(ID = as.character(ID), age = as.numeric(age))

# 1:1 matching within sex, minimizing the total absolute age difference
# (dynamic programming over age-sorted controls and cases, without replacement)
match_by_age <- function(sex_level) {
  controls <- sample_info %>% filter(group == "control", sex == sex_level) %>% arrange(age, ID)
  cases <- sample_info %>% filter(group == "case", sex == sex_level) %>% arrange(age, ID)
  m <- nrow(controls)
  n <- nrow(cases)
  cost <- matrix(Inf, m + 1, n + 1)
  take <- matrix(FALSE, m + 1, n + 1)
  cost[1, ] <- 0
  for (i in seq_len(m)) for (j in seq_len(n)) {
    match_cost <- cost[i, j] + abs(controls$age[i] - cases$age[j])
    take[i + 1, j + 1] <- match_cost <= cost[i + 1, j]
    cost[i + 1, j + 1] <- min(match_cost, cost[i + 1, j])
  }
  pairs <- list()
  i <- m
  j <- which.min(cost[m + 1, ])
  while (i > 0 && j > 1) {
    if (take[i + 1, j]) {
      pairs[[i]] <- tibble(match_set = paste0(sex_level, "_", controls$ID[i]),
                           control_ID = controls$ID[i], case_ID = cases$ID[j - 1])
      i <- i - 1
    }
    j <- j - 1
  }
  bind_rows(pairs)
}

pairs <- map_dfr(c("female", "male"), match_by_age)

sample_info %>%
  inner_join(bind_rows(pairs %>% transmute(ID = control_ID, match_set),
                       pairs %>% transmute(ID = case_ID, match_set)),
             by = "ID") %>%
  arrange(match_set, group) %>%
  write_csv(file.path(out_dir, "csf_matched_samples.csv"))
