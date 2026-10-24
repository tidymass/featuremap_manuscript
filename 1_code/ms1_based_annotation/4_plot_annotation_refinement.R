library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(ggpubr)
library(patchwork)

out_dir <- "3_data_analysis/ms1_based_annotation"

# Fig. 2b: top-1 tie-adjusted accuracy across annotation stages
stage_levels <- c("ms1_based", "network_refined", "network_weighted")

accuracy <- read_csv(file.path(out_dir, "annotation_accuracy.csv")) %>%
  filter(top == 1) %>%
  pivot_longer(all_of(stage_levels), names_to = "stage", values_to = "accuracy") %>%
  mutate(stage = factor(stage, levels = stage_levels)) %>%
  arrange(case_id)

p_accuracy <- ggplot(accuracy, aes(stage, accuracy)) +
  geom_line(aes(group = case_id), color = "grey60", linewidth = 0.3) +
  geom_boxplot(aes(fill = stage), width = 0.4, outlier.shape = NA) +
  geom_point(color = "grey30", size = 1.8) +
  stat_compare_means(
    comparisons = list(stage_levels[1:2], stage_levels[2:3], stage_levels[c(1, 3)]),
    method = "wilcox.test", paired = TRUE, label = "p.signif"
  ) +
  scale_x_discrete(labels = c("MS1-based\nonly", "Network\nrefined", "Network\nweighted")) +
  scale_fill_manual(values = c("#CCEBC5", "#FDB462", "#7FB1D3")) +
  labs(x = NULL, y = "Tie-adjusted annotation accuracy") +
  theme_base

ggsave(file.path(out_dir, "fig2b_annotation_accuracy.pdf"), p_accuracy, width = 3.6, height = 4.2)

# Fig. 2c: edge frequency of positive and negative feature-metabolite pairs
edge_frequency <- read_csv(file.path(out_dir, "edge_frequency.csv")) %>%
  mutate(type = factor(str_to_title(type), levels = c("Negative", "Positive"))) %>%
  arrange(case_id)

p_edge <- ggplot(edge_frequency, aes(type, mean_freq)) +
  geom_line(aes(group = case_id), color = "grey60", linewidth = 0.3) +
  geom_boxplot(aes(fill = type), width = 0.45, outlier.shape = NA, alpha = 0.7) +
  geom_point(color = "grey20", size = 2) +
  stat_compare_means(comparisons = list(c("Negative", "Positive")),
                     method = "wilcox.test", paired = TRUE, label = "p.signif") +
  scale_fill_manual(values = c("#8C8C8C", "#B65A5A")) +
  labs(x = NULL, y = "Mean edge frequency") +
  theme_base

ggsave(file.path(out_dir, "fig2c_edge_frequency.pdf"), p_edge, width = 2.8, height = 4)

# Fig. 2d: annotation redundancy (left) and rank of positive pairs (right)
rank_dist <- read_csv(file.path(out_dir, "positive_pair_rank.csv")) %>%
  mutate(rank_bin = case_when(rank == 1 ~ "Rank 1",
                              rank == 2 ~ "Rank 2",
                              rank == 3 ~ "Rank 3",
                              TRUE ~ "Rank >3"),
         rank_bin = factor(rank_bin, levels = c("Rank 1", "Rank 2", "Rank 3", "Rank >3"))) %>%
  count(case_id, rank_bin, .drop = FALSE) %>%
  group_by(case_id) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

case_order <- rank_dist %>%
  filter(rank_bin == "Rank 1") %>%
  arrange(prop) %>%
  pull(case_id)

p_rank <- ggplot(rank_dist, aes(factor(case_id, levels = case_order), prop, fill = rank_bin)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7, color = "white") +
  geom_text(aes(label = ifelse(prop >= 0.06, scales::percent(prop, accuracy = 1), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            color = "white", size = 3.8) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  scale_fill_manual(values = c("#C94040", "#E8854A", "#F2C35A", "#7BAEC8"), name = NULL) +
  labs(x = NULL, y = "Proportion of positive pairs") +
  theme_base +
  theme(legend.position = "right", axis.text.y = element_blank(), axis.ticks.y = element_blank())

redundancy <- read_csv(file.path(out_dir, "annotation_redundancy.csv")) %>%
  pivot_wider(names_from = stage, values_from = redundancy) %>%
  mutate(case_id = factor(case_id, levels = case_order))

p_redundancy <- ggplot(redundancy, aes(y = case_id)) +
  geom_vline(xintercept = mean(redundancy$ms1_based), linetype = "dashed", color = "#8FBBDA") +
  geom_vline(xintercept = mean(redundancy$network_refined), linetype = "dashed", color = "#E8854A") +
  geom_segment(aes(x = ms1_based, xend = network_refined, yend = case_id),
               color = "grey65", linewidth = 1) +
  geom_point(aes(x = ms1_based), shape = 21, size = 3.4, fill = "#8FBBDA") +
  geom_point(aes(x = network_refined), shape = 21, size = 3.4, fill = "#E8854A") +
  scale_x_reverse() +
  labs(x = "Annotation redundancy", y = NULL) +
  theme_base

ggsave(file.path(out_dir, "fig2d_redundancy_rank.pdf"),
       p_redundancy + p_rank + plot_layout(widths = c(0.72, 1.45)),
       width = 11, height = 4)
