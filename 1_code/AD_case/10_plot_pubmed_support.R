library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

out_dir <- "3_data_analysis/AD_case"

# Fig. 6d: AD-related PubMed records per module in each metabolic theme
pubmed <- read_csv(file.path(out_dir, "pubmed_module_counts.csv")) %>%
  mutate(manual_theme = fct_reorder(manual_theme, pubmed_n, .fun = mean))

p_pubmed <- ggplot(pubmed, aes(pubmed_n, manual_theme)) +
  stat_summary(aes(fill = manual_theme), fun = mean, geom = "col", width = 0.6, color = "black") +
  stat_summary(fun = mean, fun.min = function(x) mean(x) - sd(x), fun.max = function(x) mean(x) + sd(x),
               geom = "errorbar", width = 0.2) +
  geom_point(shape = 21, fill = "white", size = 2.5) +
  scale_fill_manual(values = c("Membrane lipid remodeling & sphingolipid homeostasis" = "#3C5488",
                               "Arachidonic acid–eicosanoid inflammatory signaling" = "#E64B35",
                               "Mitochondrial bioenergetics & redox metabolism" = "#F28E2B",
                               "Aromatic amino-acid & neuroactive metabolite metabolism" = "#8A4198",
                               "Sugar-acid & glycoconjugate metabolism" = "#7E6148",
                               "Purine/xanthine–methylxanthine metabolism" = "#C2185B")) +
  scale_y_discrete(labels = ~ str_wrap(.x, 30)) +
  labs(x = "PubMed records per module (mean ± SD)", y = NULL) +
  theme_base

ggsave(file.path(out_dir, "fig6d_pubmed_support.pdf"), p_pubmed, width = 6, height = 4)
