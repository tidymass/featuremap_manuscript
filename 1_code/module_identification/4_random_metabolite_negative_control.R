library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/module_identification/module_tools.R")

out_dir <- "3_data_analysis/module_identification/negative_control"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

bg <- load_background()
pathway_kegg <- read_csv(file.path(data_dir, "hsa00220_ipmap.csv"))$KEGG

# Candidate pool: MS/MS-supported, MRN-mapped metabolites outside arginine biosynthesis
mrn_nodes <- bg$mrn %>%
  activate(nodes) %>%
  as_tibble() %>%
  filter(!is.na(INCHI_ID), !is.na(KEGG_ID)) %>%
  distinct(KEGG_ID, .keep_all = TRUE) %>%
  dplyr::select(INCHI_ID, KEGG_ID, Compound_name)

candidate <- bg$ms2 %>%
  filter(!is.na(KEGG.ID), !KEGG.ID %in% pathway_kegg, variable_id %in% colnames(bg$intensity)) %>%
  distinct(variable_id, KEGG.ID) %>%
  inner_join(mrn_nodes, by = c("KEGG.ID" = "KEGG_ID"))

mrn_graph <- as.igraph(bg$mrn)
dist_mat <- distances(mrn_graph,
                      v = V(mrn_graph)[name %in% candidate$INCHI_ID],
                      to = V(mrn_graph)[name %in% candidate$INCHI_ID])

# Seven metabolites separated by at least 6 reaction steps
pick_distant <- function(seed, n = 7, min_dist = 6) {
  set.seed(seed)
  for (attempt in 1:2000) {
    picked <- character()
    for (id in sample(rownames(dist_mat))) {
      if (all(dist_mat[id, picked] >= min_dist)) picked <- c(picked, id)
      if (length(picked) == n) return(picked)
    }
  }
  stop("Not enough distant metabolites")
}

result <- map(1:100, function(i) {
  message("Negative control ", i)
  set.seed(190000 + i)
  sample_info <- bg$sample_info %>%
    mutate(group = sample(rep(c("case", "control"), each = n() / 2)))
  case_idx <- rownames(bg$intensity) %in% sample_info$ID[sample_info$group == "case"]

  picked <- pick_distant(seed = 160000 + i)
  features <- candidate %>%
    filter(INCHI_ID %in% picked) %>%
    arrange(match(INCHI_ID, picked), variable_id) %>%
    pull(variable_id) %>%
    unique()

  mat <- spike_features(as.matrix(bg$intensity), features, case_idx,
                        fc_seed = 180000 + i, spike_seed = 170000 + i)
  res <- run_featuremap(bg$base_res, diff_test(bg, mat, sample_info$group))

  list(metabolites = tibble(iter_id = i, INCHI_ID = picked),
       module_score = res$module_score %>% mutate(iter_id = i, .before = 1),
       sig_module_table = res$module_table %>%
         filter(Module %in% res$module_score$Module[res$module_score$fdr < 0.05]) %>%
         mutate(iter_id = i, .before = 1))
})

module_score <- map_dfr(result, "module_score")

module_score %>%
  group_by(iter_id) %>%
  summarise(n_sig_modules = sum(fdr < 0.05, na.rm = TRUE)) %>%
  write_csv(file.path(out_dir, "negative_control_runs.csv"))

module_score %>%
  filter(fdr < 0.05) %>%
  transmute(iter_id, Module, module_p = p, module_fdr = fdr, module_score = Module_score) %>%
  write_csv(file.path(out_dir, "negative_control_sig_modules.csv"))

write_csv(map_dfr(result, "metabolites"), file.path(out_dir, "negative_control_metabolites.csv"))
saveRDS(map_dfr(result, "sig_module_table"), file.path(out_dir, "negative_control_sig_module_tables.rds"))
