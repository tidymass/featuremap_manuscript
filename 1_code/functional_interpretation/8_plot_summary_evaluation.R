library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(ggpubr)

out_dir <- "3_data_analysis/functional_interpretation"

# Fig. 4d: summary agreement with PathBank descriptions (embedding and LLM judge)
summary_long <- read_csv(file.path(out_dir, "summary_evaluation.csv")) %>%
  transmute(difficulty = str_to_title(difficulty),
            `OpenAI embedding` = score_openai_embedding,
            `LLM judge` = judge_score) %>%
  pivot_longer(-difficulty, names_to = "metric", values_to = "score") %>%
  mutate(metric = factor(metric, levels = c("OpenAI embedding", "LLM judge")))

summary_stat <- summary_long %>%
  group_by(metric) %>%
  summarise(p = wilcox.test(score[difficulty == "Easy"], score[difficulty == "Hard"])$p.value) %>%
  mutate(p.signif = sig_label(p), xmin = as.numeric(metric) - 0.2, xmax = xmin + 0.4,
         y.position = 1.08, group1 = NA, group2 = NA)

p_summary <- ggplot(summary_long, aes(metric, score, fill = difficulty)) +
  geom_violin(position = position_dodge(0.8), width = 0.7, alpha = 0.6) +
  geom_boxplot(position = position_dodge(0.8), width = 0.15, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size = 0.8, alpha = 0.6) +
  stat_pvalue_manual(summary_stat, label = "p.signif", tip.length = 0) +
  scale_fill_manual(values = c(Easy = "#4DB6AC", Hard = "#F28E2B"), name = NULL) +
  labs(x = NULL, y = "Summary similarity") +
  theme_base +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig4d_summary_similarity.pdf"), p_summary, width = 4.5, height = 4)
