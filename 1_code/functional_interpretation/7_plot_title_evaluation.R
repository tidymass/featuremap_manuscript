library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(ggpubr)
library(patchwork)

out_dir <- "3_data_analysis/functional_interpretation"

title_evaluation <- read_csv(file.path(out_dir, "title_evaluation.csv")) %>%
  mutate(difficulty = str_to_title(difficulty))

# Fig. 4b: title embedding similarity to the reference vs randomly paired titles
title_long <- title_evaluation %>%
  pivot_longer(c(shuffled_title_similarity, title_similarity),
               names_to = "comparison", values_to = "similarity") %>%
  mutate(comparison = recode(comparison, shuffled_title_similarity = "Random", title_similarity = "True"))

title_stat <- title_long %>%
  group_by(difficulty) %>%
  summarise(p = wilcox.test(similarity[comparison == "True"], similarity[comparison == "Random"],
                            paired = TRUE)$p.value) %>%
  mutate(xmin = match(difficulty, c("Easy", "Hard")) - 0.2, xmax = xmin + 0.4, y.position = 1.08) %>%
  bind_rows(tibble(p = with(title_evaluation, wilcox.test(title_similarity[difficulty == "Easy"],
                                                          title_similarity[difficulty == "Hard"])$p.value),
                   xmin = 1.2, xmax = 2.2, y.position = 1.18)) %>%
  mutate(p.signif = sig_label(p), group1 = NA, group2 = NA)

p_title <- ggplot(title_long, aes(difficulty, similarity, fill = comparison)) +
  geom_violin(position = position_dodge(0.8), width = 0.7, alpha = 0.6) +
  geom_boxplot(position = position_dodge(0.8), width = 0.15, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 0.8, alpha = 0.6) +
  stat_pvalue_manual(title_stat, label = "p.signif", tip.length = 0) +
  scale_fill_manual(values = c(Random = "#9E9E9E", True = "#4DB6AC"), name = NULL) +
  labs(x = NULL, y = "Title embedding similarity") +
  theme_base +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig4b_title_similarity.pdf"), p_title, width = 4.5, height = 4)

# Fig. 4c: word overlap between generated and reference titles
overlap_colors <- c("Exact match" = "#4E79A7", "Partial word overlap" = "#F2C14E", "No word overlap" = "#A7A9AC")

p_overlap <- title_evaluation %>%
  count(difficulty, overlap_class) %>%
  group_by(difficulty) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(overlap_class = factor(overlap_class, levels = names(overlap_colors))) %>%
  ggplot(aes(x = 2, y = prop, fill = overlap_class)) +
  geom_col(color = "white", width = 1) +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)), position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_wrap(~ difficulty) +
  scale_fill_manual(values = overlap_colors, name = NULL) +
  theme_void() +
  theme(legend.position = "bottom")

# ROUGE-L recall of partially matched titles; lines link the same pathway
rouge <- title_evaluation %>%
  filter(overlap_class == "Partial word overlap")

p_rouge <- ggplot(rouge, aes(difficulty, rougeL_recall)) +
  geom_line(data = ~ group_by(.x, gold_pathway_id) %>% filter(n() == 2),
            aes(group = gold_pathway_id), color = "grey70") +
  geom_boxplot(aes(fill = difficulty), width = 0.5, outlier.shape = NA) +
  geom_point(shape = 21, fill = "white") +
  stat_compare_means(comparisons = list(c("Easy", "Hard")), method = "wilcox.test", label = "p.signif") +
  scale_fill_manual(values = c(Easy = "#4DB6AC", Hard = "#F28E2B")) +
  labs(x = NULL, y = "ROUGE-L recall") +
  theme_base

ggsave(file.path(out_dir, "fig4c_title_overlap.pdf"),
       p_overlap + p_rouge + plot_layout(widths = c(2, 1)),
       width = 8, height = 4)
