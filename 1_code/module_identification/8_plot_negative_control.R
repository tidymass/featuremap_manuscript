library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(ggpubr)
library(patchwork)

out_dir <- "3_data_analysis/module_identification/negative_control"

# Fig. 3f (left): runs with at least one significant module
p_donut <- read_csv(file.path(out_dir, "negative_control_runs.csv")) %>%
  count(group = ifelse(n_sig_modules > 0, "With significant module", "No significant module")) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(x = 2, y = prop, fill = group)) +
  geom_col(color = "white", width = 1) +
  geom_text(aes(label = scales::percent(prop)), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  scale_fill_manual(values = c("#E7E1F2", "#6A51A3"), name = NULL) +
  theme_void()

# Fig. 3f (right): significant modules in negative runs vs successful positive runs
positive <- read_csv("3_data_analysis/module_identification/arginine_biosynthesis/fmap_module_score.csv") %>%
  filter(!is.na(p)) %>%
  group_by(iter_id) %>%
  slice_min(p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(n_spiked >= 5) %>%
  transmute(source = "Positive", module_fdr = fdr, module_score = Module_score)

strength <- read_csv(file.path(out_dir, "negative_control_sig_modules.csv")) %>%
  transmute(source = "Negative", module_fdr, module_score) %>%
  bind_rows(positive) %>%
  mutate(neglog10_fdr = -log10(module_fdr))

plot_strength <- function(y, y_lab) {
  ggplot(strength, aes(source, .data[[y]], fill = source)) +
    geom_jitter(aes(color = source), width = 0.1, size = 1.4, alpha = 0.7) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8) +
    stat_compare_means(comparisons = list(c("Negative", "Positive")),
                       method = "wilcox.test", label = "p.signif") +
    scale_fill_manual(values = c(Negative = "#6A51A3", Positive = "#5B8E62")) +
    scale_color_manual(values = c(Negative = "#6A51A3", Positive = "#5B8E62")) +
    labs(x = NULL, y = y_lab) +
    theme_base
}

ggsave(file.path(out_dir, "fig3f_negative_control.pdf"),
       p_donut + plot_strength("neglog10_fdr", expression(-Log[10](FDR))) +
         plot_strength("module_score", "Module score") +
         plot_layout(widths = c(1.6, 1, 1)),
       width = 10, height = 3.5)
