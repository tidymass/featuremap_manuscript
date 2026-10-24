library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/module_identification/module_tools.R")

out_dir <- "3_data_analysis/module_identification/two_pathways"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bg <- load_background()
targets <- list(hsa00220 = pathway_targets(bg, "hsa00220"),
                hsa01040 = pathway_targets(bg, "hsa01040"))
seeds <- list(hsa00220 = c(fc = 40000, spike = 30000),
              hsa01040 = c(fc = 50000, spike = 60000))
iwalk(targets, ~ write_csv(.x$target, file.path(out_dir, paste0(.y, "_spiked_features.csv"))))

# Same pseudo-case / pseudo-control split as the single-pathway simulation
sample_info <- read_csv("3_data_analysis/module_identification/arginine_biosynthesis/sample_split.csv") %>%
  mutate(ID = as.character(ID))
case_idx <- rownames(bg$intensity) %in% sample_info$ID[sample_info$group == "case"]

# Arginine biosynthesis and unsaturated fatty acid biosynthesis spiked together
result <- map(1:100, function(i) {
  message("Simulation ", i)
  mat <- as.matrix(bg$intensity)
  for (id in names(targets)) {
    mat <- spike_features(mat, targets[[id]]$features, case_idx,
                          fc_seed = seeds[[id]][["fc"]] + i,
                          spike_seed = seeds[[id]][["spike"]] + i)
  }
  res <- run_featuremap(bg$base_res, diff_test(bg, mat, sample_info$group))
  list(module_table = res$module_table,
       module_score = res$module_score %>% mutate(iter_id = i, .before = 1),
       recovery = imap_dfr(targets, ~ count_recovery(res$module_table, .x) %>%
                             mutate(pathway_id = .y)) %>%
         mutate(iter_id = i, .before = 1))
})

write_csv(map_dfr(result, "module_score"), file.path(out_dir, "fmap_module_score.csv"))
write_csv(map_dfr(result, "recovery"), file.path(out_dir, "fmap_module_recovery.csv"))
saveRDS(map(result, "module_table"), file.path(out_dir, "fmap_module_tables.rds"))
