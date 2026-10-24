library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(ggpubr)
library(geomtextpath)
library(patchwork)

out_dir <- "3_data_analysis/ms1_based_annotation"

case_ids <- c("MTBLS13267", "MTBLS267", "MTBLS13855", "MTBLS13095", "MTBLS12329",
              "MTBLS13578", "MTBLS11987", "MTBLS417", "MTBLS5028", "MTBLS12854")
tool_levels <- c("FeatureMAP", "mWISE", "xMSannotator")
tool_colors <- c(FeatureMAP = "#CB4B2A", mWISE = "#3D87BF", xMSannotator = "#2E7A5C")

# Fig. 2f: top-1 accuracy and F1 on the common feature set, one sector per dataset
metrics <- read_csv(file.path(out_dir, "tool_benchmark_metrics.csv")) %>%
  filter(top == 1)

plot_circular <- function(metric, title) {
  gap <- 2
  n_slot <- length(case_ids) * (length(tool_levels) + gap)
  df <- metrics %>%
    transmute(case_id = factor(case_id, levels = case_ids),
              tool = factor(tool, levels = tool_levels),
              value = .data[[metric]]) %>%
    mutate(pos = (as.integer(case_id) - 1) * (length(tool_levels) + gap) + as.integer(tool) + 1)
  y_max <- max(df$value)
  y_refs <- round(c(0, y_max / 2, y_max), 2)
  ring <- df %>%
    group_by(case_id) %>%
    summarise(xmin = min(pos) - 1.4, xmax = max(pos) + 1.4)
  label_path <- ring %>%
    group_by(case_id) %>%
    reframe(x = seq(xmin + 0.05, xmax - 0.05, length.out = 100))
  separator <- tibble(x = seq_len(length(case_ids) - 1) * (length(tool_levels) + gap) + 0.5)

  ggplot(df, aes(pos, value, fill = tool)) +
    geom_hline(yintercept = y_refs, color = "grey82", linewidth = 0.25) +
    geom_col(width = 0.82) +
    geom_segment(data = separator, aes(x = x, xend = x, y = 0, yend = y_max),
                 inherit.aes = FALSE, linetype = "dashed", color = "grey55", linewidth = 0.3) +
    geom_rect(data = ring, aes(xmin = xmin, xmax = xmax, ymin = 1.1 * y_max, ymax = 1.3 * y_max),
              inherit.aes = FALSE, fill = NA, color = "grey45") +
    geom_textpath(data = label_path,
                  aes(x = x, y = 1.2 * y_max, label = as.character(case_id), group = case_id),
                  inherit.aes = FALSE, text_only = TRUE, upright = TRUE, size = 4) +
    annotate("text", x = n_slot + 0.3, y = y_refs, label = y_refs, color = "grey45", vjust = -0.25) +
    coord_polar() +
    scale_x_continuous(limits = c(0.5, n_slot + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.5 * y_max, 1.37 * y_max), expand = c(0, 0)) +
    scale_fill_manual(values = tool_colors, name = NULL) +
    labs(title = title) +
    theme_void(base_size = 14) +
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5))
}

ggsave(file.path(out_dir, "fig2f_tool_benchmark.pdf"),
       plot_circular("efficient_score", "Tie-adjusted annotation accuracy") +
         plot_circular("F1", "F1 score"),
       width = 12, height = 6.5)

# Fig. 2g: mean candidates per feature in each tool's full output
candidates <- read_csv(file.path(out_dir, "tool_candidates_per_feature.csv")) %>%
  mutate(tool = factor(tool, levels = tool_levels)) %>%
  arrange(case_id)

p_candidates <- ggplot(candidates, aes(tool, mean_candidates_per_feature)) +
  stat_summary(aes(fill = tool), fun = mean, geom = "bar", width = 0.55) +
  stat_summary(fun = mean,
               fun.min = function(x) mean(x) - sd(x),
               fun.max = function(x) mean(x) + sd(x),
               geom = "errorbar", width = 0.18) +
  geom_jitter(width = 0.1, shape = 21, fill = "white", size = 2) +
  stat_compare_means(
    comparisons = list(tool_levels[1:2], tool_levels[2:3], tool_levels[c(1, 3)]),
    method = "wilcox.test", paired = TRUE, label = "p.signif"
  ) +
  scale_fill_manual(values = tool_colors) +
  labs(x = NULL, y = "Mean candidates per feature") +
  theme_base

ggsave(file.path(out_dir, "fig2g_candidates_per_feature.pdf"), p_candidates, width = 3.6, height = 4.2)
