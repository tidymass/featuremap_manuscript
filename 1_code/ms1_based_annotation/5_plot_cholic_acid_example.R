library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(patchwork)

out_dir <- "3_data_analysis/ms1_based_annotation"

# Fig. 2e: feature M407T446_NEG in MTBLS11987, MS/MS reference cholic acid (HMDB0000619)
res <- load_rda_first("2_data/ms1_based_annotation/MTBLS11987/res.rda")

ms1_candidate <- res$anno_feature_table1 %>%
  normalize_id_cols() %>%
  filter(variable_id == "M407T446_NEG") %>%
  arrange(desc(score_norm)) %>%
  distinct(HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  dplyr::select(Compound.name, HMDB.ID, KEGG.ID, score_norm)

network_candidate <- res$anno_feature_table2 %>%
  normalize_id_cols() %>%
  filter(variable_id == "M407T446_NEG") %>%
  arrange(desc(weighted_anno_score)) %>%
  distinct(HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
  dplyr::select(HMDB.ID, KEGG.ID, weighted_anno_score)

ms1_candidate %>%
  left_join(network_candidate, by = c("HMDB.ID", "KEGG.ID")) %>%
  write_csv(file.path(out_dir, "cholic_acid_candidates.csv"))

cholic_acid <- read_csv(file.path(out_dir, "cholic_acid_candidates.csv")) %>%
  mutate(order = row_number(),
         is_reference = HMDB.ID %in% "HMDB0000619")

p_ms1 <- ggplot(cholic_acid, aes(order, score_norm, fill = is_reference)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("#7F95A0", "#D7301F")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Annotation candidates", y = "MS1-based annotation score") +
  theme_base +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_network <- cholic_acid %>%
  filter(!is.na(weighted_anno_score)) %>%
  mutate(Compound.name = fct_reorder(Compound.name, weighted_anno_score)) %>%
  ggplot(aes(weighted_anno_score, Compound.name, fill = is_reference)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = c("#EFC7A7", "#D7301F")) +
  scale_x_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  labs(x = "Network-weighted\nannotation score", y = NULL) +
  theme_base

ggsave(file.path(out_dir, "fig2e_cholic_acid.pdf"),
       p_ms1 + p_network + plot_layout(widths = c(1.4, 1)),
       width = 11, height = 3.2)
