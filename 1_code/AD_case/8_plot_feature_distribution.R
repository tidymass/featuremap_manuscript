library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")

out_dir <- "3_data_analysis/AD_case"

# Fig. 6b: log2 fold change of all features; rugs mark nominally differential features
datasets <- c(csf_hilic = "CSF - HILIC", csf_rp = "CSF - RP", ipsc_hilic = "iPSC - HILIC", ipsc_rp = "iPSC - RP")

features <- imap_dfr(datasets, ~ read_csv(file.path(out_dir, paste0("feature_table_", .y, ".csv"))) %>%
                       transmute(dataset = .x, condition, pvalue)) %>%
  filter(is.finite(condition), is.finite(pvalue)) %>%
  mutate(dataset = factor(dataset, levels = datasets),
         direction = case_when(pvalue < 0.05 & condition < 0 ~ "Down (P < 0.05)",
                               pvalue < 0.05 & condition > 0 ~ "Up (P < 0.05)"))

counts <- features %>%
  group_by(dataset) %>%
  summarise(label = sprintf("All %s | P < 0.05: %s\nDown %s | Up %s",
                            scales::comma(n()), scales::comma(sum(pvalue < 0.05)),
                            sum(direction == "Down (P < 0.05)", na.rm = TRUE),
                            sum(direction == "Up (P < 0.05)", na.rm = TRUE)))

p_features <- ggplot(features, aes(condition)) +
  geom_density(aes(y = after_stat(count)), fill = "#DCE2E8", color = "#8995A3") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  geom_rug(data = ~ filter(.x, !is.na(direction)), aes(color = direction), sides = "b", length = unit(0.075, "npc")) +
  geom_text(data = counts, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.2, size = 3) +
  facet_grid(dataset ~ ., scales = "free_y") +
  scale_color_manual(values = c(`Down (P < 0.05)` = "#59678C", `Up (P < 0.05)` = "#B98552"), name = NULL) +
  labs(x = expression(log[2]~"(fold change)"), y = "Feature count") +
  theme_base +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "fig6b_feature_distribution.pdf"), p_features, width = 6, height = 7)
