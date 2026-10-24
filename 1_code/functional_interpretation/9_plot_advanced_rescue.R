library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(ggpubr)
library(ggalluvial)
library(patchwork)

out_dir <- "3_data_analysis/functional_interpretation"

state_levels <- c("Exact match", "Partial word overlap", "No word overlap", "Rescued", "Not rescued")
state_colors <- c("#4E79A7", "#F2C14E", "#A7A9AC", "#D95F02", "#B8B8D9")
names(state_colors) <- state_levels

# Fig. 4e,f (left): title overlap of weak modules in fast and advanced mode
sankey <- read_csv(file.path(out_dir, "advanced_comparison.csv")) %>%
  count(difficulty,
        Fast = factor(overlap_class, levels = state_levels),
        Advanced = factor(advanced_overlap_class, levels = state_levels),
        Outcome = factor(ifelse(title_rescued, "Rescued", "Not rescued"), levels = state_levels))

plot_sankey <- function(setting) {
  sankey %>%
    filter(difficulty == setting) %>%
    ggplot(aes(axis1 = Fast, axis2 = Advanced, axis3 = Outcome, y = n)) +
    geom_alluvium(aes(fill = Fast), width = 0.2, alpha = 0.6) +
    geom_stratum(aes(fill = after_stat(stratum)), width = 0.2, color = "white") +
    geom_text(stat = "stratum", aes(label = paste(after_stat(stratum), after_stat(count))),
              angle = 90, size = 3) +
    scale_x_discrete(limits = c("Fast", "Advanced", "Outcome")) +
    scale_fill_manual(values = state_colors) +
    theme_void() +
    theme(legend.position = "none", axis.text.x = element_text())
}

# Fig. 4e,f (right): LLM judge score of rescued modules, fast vs advanced
rescued_judge <- read_csv(file.path(out_dir, "rescued_summary_judge.csv")) %>%
  transmute(module_id, difficulty,
            Fast = fast_judge_overall / 4,
            Advanced = advanced_judge_overall / 4) %>%
  pivot_longer(c(Fast, Advanced), names_to = "mode", values_to = "score") %>%
  mutate(mode = factor(mode, levels = c("Fast", "Advanced"))) %>%
  arrange(module_id)

plot_judge <- function(setting) {
  rescued_judge %>%
    filter(difficulty == setting) %>%
    ggplot(aes(mode, score)) +
    geom_line(aes(group = module_id), color = "grey70") +
    geom_boxplot(aes(fill = mode), width = 0.5, outlier.shape = NA) +
    geom_point(shape = 21, fill = "white") +
    stat_compare_means(comparisons = list(c("Fast", "Advanced")), method = "wilcox.test",
                       paired = TRUE, label = "p.signif") +
    scale_fill_manual(values = c(Fast = "grey70", Advanced = "#D95F02")) +
    scale_y_continuous(limits = c(0, 1.1)) +
    labs(x = NULL, y = "LLM judge score") +
    theme_base
}

for (setting in c("easy", "hard")) {
  ggsave(file.path(out_dir, paste0("fig4_advanced_rescue_", setting, ".pdf")),
         plot_sankey(setting) + plot_judge(setting) + plot_layout(widths = c(2.5, 1)),
         width = 9, height = 4.5)
}
