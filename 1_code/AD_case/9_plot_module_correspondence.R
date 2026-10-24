library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(ggforce)

out_dir <- "3_data_analysis/AD_case"

theme_colors <- c("Membrane lipid remodeling & sphingolipid homeostasis" = "#3C5488",
                  "Arachidonic acid–eicosanoid inflammatory signaling" = "#E64B35",
                  "Mitochondrial bioenergetics & redox metabolism" = "#F28E2B",
                  "Aromatic amino-acid & neuroactive metabolite metabolism" = "#8A4198",
                  "Sugar-acid & glycoconjugate metabolism" = "#7E6148",
                  "Purine/xanthine–methylxanthine metabolism" = "#C2185B")

# Modules manually grouped into metabolic themes from their titles and summaries
themes <- read_csv(file.path(out_dir, "module_theme_classification.csv")) %>%
  inner_join(read_csv(file.path(out_dir, "selected_modules.csv")) %>% dplyr::select(module_key = Module_Key, fdr),
             by = "module_key")

# Fig. 6c: CSF-iPSC module pairs with similarity >= 0.70
pairs <- read_csv(file.path(out_dir, "module_similarity_pairs.csv")) %>%
  filter(Similarity_Score >= 0.7,
         CSF_Module_Key %in% themes$module_key[themes$side == "CSF"],
         IPSC_Module_Key %in% themes$module_key[themes$side == "iPSC"])

modules <- themes %>%
  filter(module_key %in% c(pairs$CSF_Module_Key, pairs$IPSC_Module_Key)) %>%
  mutate(manual_theme = factor(manual_theme, levels = names(theme_colors))) %>%
  group_by(side) %>%
  arrange(manual_theme, module_label, .by_group = TRUE) %>%
  mutate(x = ifelse(side == "CSF", 1, 2),
         y = -(row_number() - 1) / (n() - 1)) %>%
  ungroup()

theme_nodes <- modules %>%
  group_by(side, manual_theme) %>%
  summarise(y = mean(y), .groups = "drop") %>%
  mutate(x = ifelse(side == "CSF", 0, 3))

theme_links <- modules %>%
  left_join(theme_nodes %>% dplyr::select(side, manual_theme, x_theme = x, y_theme = y),
            by = c("side", "manual_theme"))

module_links <- pairs %>%
  left_join(modules %>% dplyr::select(CSF_Module_Key = module_key, y_csf = y), by = "CSF_Module_Key") %>%
  left_join(modules %>% dplyr::select(IPSC_Module_Key = module_key, y_ipsc = y), by = "IPSC_Module_Key")

p_correspondence <- ggplot() +
  geom_diagonal(data = theme_links, aes(x = x_theme, y = y_theme, xend = x, yend = y, color = manual_theme),
                linewidth = 0.5) +
  geom_diagonal(data = module_links, aes(x = 1, y = y_csf, xend = 2, yend = y_ipsc, linewidth = Similarity_Score),
                color = "grey60", alpha = 0.7) +
  geom_point(data = modules, aes(x, y, size = -log10(fdr), fill = side), shape = 22) +
  geom_point(data = theme_nodes, aes(x, y, fill = manual_theme), shape = 23, size = 5) +
  geom_text(data = modules, aes(x + ifelse(side == "CSF", -0.05, 0.05), y, label = module_label,
                                hjust = ifelse(side == "CSF", 1, 0)), size = 2.5, vjust = -0.6) +
  scale_color_manual(values = theme_colors, guide = "none") +
  scale_fill_manual(values = c(theme_colors, CSF = "#3A9BC9", iPSC = "#5FA83A"), name = NULL) +
  scale_linewidth(range = c(0.5, 3), name = "Module similarity") +
  scale_size(range = c(2, 5), name = expression(-Log[10](FDR))) +
  scale_x_continuous(breaks = 0:3, labels = c("Metabolic themes", "CSF modules", "iPSC modules", "Metabolic themes")) +
  theme_void() +
  theme(axis.text.x = element_text(), legend.position = "bottom")

ggsave(file.path(out_dir, "fig6c_module_correspondence.pdf"), p_correspondence, width = 10, height = 12)
