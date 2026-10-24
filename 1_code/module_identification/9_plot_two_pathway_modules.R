library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/module_identification/module_tools.R")

out_dir <- "3_data_analysis/module_identification/two_pathways"
load("2_data/metabolic_reaction_network.rda")

# Module per pathway: the one recovering the most pathway metabolites
selected <- read_csv(file.path(out_dir, "fmap_module_recovery.csv")) %>%
  left_join(read_csv(file.path(out_dir, "fmap_module_score.csv")) %>%
              dplyr::select(iter_id, Module, fdr, Module_score),
            by = c("iter_id", "Module")) %>%
  group_by(iter_id, pathway_id) %>%
  arrange(desc(n_pathway), desc(n_spiked), desc(n_hidden), desc(Module_score), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

# Success: both pathways recovered as two distinct significant modules
success <- selected %>%
  group_by(iter_id) %>%
  filter(all(!is.na(fdr) & fdr < 0.05), n_distinct(Module) == 2) %>%
  ungroup()

module_tables <- readRDS(file.path(out_dir, "fmap_module_tables.rds"))

# Fig. 3g: merged modules of the two pathways across successful runs
for (pathway_id in c("hsa00220", "hsa01040")) {
  node_file <- file.path(out_dir, paste0("two_pathway_module_nodes_", pathway_id, ".csv"))
  module_node_table(
    module_tables,
    success %>% filter(pathway_id == !!pathway_id) %>% dplyr::select(iter_id, Module),
    read_targets(file.path(out_dir, paste0(pathway_id, "_spiked_features.csv")), pathway_id,
                 metabolic_reaction_network),
    pathway_id,
    metabolic_reaction_network
  ) %>%
    write_csv(node_file)

  ggsave(file.path(out_dir, paste0("fig3g_module_network_", pathway_id, ".pdf")),
         plot_module_network(read_csv(node_file), metabolic_reaction_network),
         width = 9, height = 7)
}
