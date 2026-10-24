library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

out_dir <- "3_data_analysis/module_identification/arginine_biosynthesis"

target <- read_csv(file.path(out_dir, "arginine_spiked_features.csv"))
idmap <- read_csv("2_data/module_identification/hsa00220_ipmap.csv") %>% filter(!is.na(KEGG))
spiked_kegg <- unique(na.omit(target$KEGG.ID))
hidden_kegg <- setdiff(idmap$KEGG, spiked_kegg)

fmap_run_summary <- read_csv(file.path(out_dir, "fmap_run_summary.csv"))
fmap_module_score <- read_csv(file.path(out_dir, "fmap_module_score.csv"))
fmap_module_tables <- readRDS(file.path(out_dir, "fmap_module_tables.rds"))

# mummichog: primary module = smallest p across positive and negative modes
read_mummichog <- function(i) {
  map_dfr(c("pos", "neg"), function(mode) {
    f <- file.path(out_dir, "mummichog_result", paste0("feature_mummichog_", i), mode, "module_analysis.tsv")
    if (!file.exists(f)) return(NULL)
    read_tsv(f, col_types = cols(.default = "c"))
  }) %>%
    mutate(p = as.numeric(`p-value`)) %>%
    filter(!is.na(p)) %>%
    arrange(p)
}

mummichog_members <- map(1:100, function(i) {
  modules <- read_mummichog(i)
  list(p = modules$p,
       members = if (nrow(modules) > 0) unique(str_trim(str_split(modules$`members (id)`[1], ",")[[1]])) else character())
})

mummichog_run_summary <- imap_dfr(mummichog_members, function(x, i) {
  tibble(iter_id = i,
         main_module_p = x$p[1],
         n_sig = sum(x$p < 0.05),
         spiked_rec = sum(x$members %in% spiked_kegg),
         hidden_rec = sum(x$members %in% hidden_kegg))
})
write_csv(mummichog_run_summary, file.path(out_dir, "mummichog_run_summary.csv"))

# Success: exactly one significant module containing >= k spiked metabolites
success_rate <- bind_rows(
  fmap_run_summary %>% transmute(tool = "FeatureMAP", n_sig = n_sig_fdr, spiked_rec),
  mummichog_run_summary %>% transmute(tool = "Mummichog", n_sig, spiked_rec)
) %>%
  crossing(k = 1:7) %>%
  group_by(tool, k) %>%
  summarise(success_rate = mean(n_sig == 1 & spiked_rec >= k), .groups = "drop")
write_csv(success_rate, file.path(out_dir, "success_rate.csv"))

# Recovery frequency of each pathway metabolite in the primary module
fmap_main <- fmap_module_score %>%
  filter(!is.na(p)) %>%
  group_by(iter_id) %>%
  slice_min(p, n = 1, with_ties = FALSE) %>%
  ungroup()

fmap_hits <- map_dfr(seq_len(nrow(fmap_main)), function(j) {
  module <- fmap_module_tables[[fmap_main$iter_id[j]]] %>%
    filter(Module == fmap_main$Module[j])
  tibble(KEGG.ID = unique(c(target$KEGG.ID[target$INCHI_ID %in% module$INCHI_ID],
                            intersect(module$KEGG_ID, hidden_kegg))))
})

mummichog_hits <- map_dfr(mummichog_members, ~ tibble(KEGG.ID = intersect(.x$members, idmap$KEGG)))

node_recovery <- full_join(count(fmap_hits, KEGG.ID, name = "n_fmap"),
                           count(mummichog_hits, KEGG.ID, name = "n_mummi"),
                           by = "KEGG.ID") %>%
  left_join(idmap %>% distinct(KEGG.ID = KEGG, label = Match), by = "KEGG.ID") %>%
  transmute(KEGG.ID, label,
            fmap_freq = coalesce(n_fmap, 0L) / 100,
            mummi_freq = coalesce(n_mummi, 0L) / 100,
            node_type = ifelse(KEGG.ID %in% spiked_kegg, "Spiked", "Unspiked"))
write_csv(node_recovery, file.path(out_dir, "node_recovery_frequency.csv"))
