library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/module_identification/module_tools.R")

out_dir <- "3_data_analysis/module_identification/arginine_biosynthesis"
dir.create(file.path(out_dir, "mummichog_input"), recursive = TRUE, showWarnings = FALSE)

bg <- load_background()
arginine <- pathway_targets(bg, "hsa00220")
write_csv(arginine$target, file.path(out_dir, "arginine_spiked_features.csv"))

# Random split of 40 controls into pseudo-case and pseudo-control groups
set.seed(20003)
sample_info <- bg$sample_info %>%
  mutate(group = sample(rep(c("case", "control"), each = n() / 2)))
write_csv(sample_info, file.path(out_dir, "sample_split.csv"))
case_idx <- rownames(bg$intensity) %in% sample_info$ID[sample_info$group == "case"]

# 100 simulations: spike arginine-biosynthesis features 2-10 fold in pseudo-cases
result <- map(1:100, function(i) {
  message("Simulation ", i)
  mat <- spike_features(as.matrix(bg$intensity), arginine$features, case_idx,
                        fc_seed = 40000 + i, spike_seed = 30000 + i)
  feature_table <- diff_test(bg, mat, sample_info$group)

  feature_table %>%
    transmute(mz, rtime = rt, `p-value` = pvalue, `t-score` = condition, custom_id = variable_id) %>%
    write_tsv(file.path(out_dir, "mummichog_input", paste0("feature_mummichog_", i, ".txt")))

  res <- run_featuremap(bg$base_res, feature_table)
  list(module_table = res$module_table,
       module_score = res$module_score %>%
         left_join(count_recovery(res$module_table, arginine), by = "Module") %>%
         mutate(iter_id = i, .before = 1))
})

module_score <- map_dfr(result, "module_score") %>%
  mutate(is_sig = !is.na(fdr) & fdr < 0.05)

# Main module per run: the significant module recovering the most spiked metabolites
run_summary <- module_score %>%
  group_by(iter_id) %>%
  mutate(n_sig_fdr = sum(is_sig)) %>%
  arrange(desc(is_sig), desc(n_spiked), desc(n_pathway), desc(Module_score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(iter_id, main_module = Module, main_module_p = p, main_module_fdr = fdr,
            main_module_score = Module_score, n_sig_fdr,
            spiked_rec = n_spiked, hidden_rec = n_hidden)

write_csv(module_score, file.path(out_dir, "fmap_module_score.csv"))
write_csv(run_summary, file.path(out_dir, "fmap_run_summary.csv"))
saveRDS(map(result, "module_table"), file.path(out_dir, "fmap_module_tables.rds"))
