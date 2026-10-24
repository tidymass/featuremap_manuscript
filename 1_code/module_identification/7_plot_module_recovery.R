library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(patchwork)

out_dir <- "3_data_analysis/module_identification/arginine_biosynthesis"

run_summary <- read_csv(file.path(out_dir, "fmap_run_summary.csv")) %>%
  mutate(success = n_sig_fdr == 1 & spiked_rec >= 5)
module_score <- read_csv(file.path(out_dir, "fmap_module_score.csv")) %>%
  mutate(is_sig = !is.na(fdr) & fdr < 0.05, neglog10_fdr = -log10(fdr))
mummichog_run_summary <- read_csv(file.path(out_dir, "mummichog_run_summary.csv"))

# Fig. 3b: recovered metabolites, success and module significance per run
run_x <- scale_x_continuous(limits = c(0.5, 100.5), expand = c(0, 0))

p_recovered <- run_summary %>%
  pivot_longer(c(spiked_rec, hidden_rec), names_to = "type", values_to = "n") %>%
  mutate(type = factor(type, levels = c("hidden_rec", "spiked_rec"), labels = c("Unspiked", "Spiked"))) %>%
  ggplot(aes(iter_id, n, fill = type)) +
  geom_col(width = 0.75) +
  run_x +
  scale_fill_manual(values = c(Unspiked = "#2171B5", Spiked = "#D94F00"), name = NULL) +
  labs(x = NULL, y = "Recovered\ncompounds") +
  theme_base +
  theme(legend.position = "top", axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_success <- ggplot(run_summary, aes(iter_id, 1, fill = success)) +
  geom_tile() +
  run_x +
  scale_fill_manual(values = c(`FALSE` = "#E0E0E0", `TRUE` = "#238B45")) +
  theme_void() +
  theme(legend.position = "none")

p_fdr <- ggplot(module_score, aes(iter_id, neglog10_fdr)) +
  geom_jitter(data = ~ filter(.x, !is_sig), color = "grey70", size = 0.8, width = 0.3, height = 0) +
  geom_point(data = ~ filter(.x, is_sig), color = "#990000", shape = 17, size = 1.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#B5563A") +
  run_x +
  labs(x = "Simulation runs", y = expression(-Log[10](FDR))) +
  theme_base

ggsave(file.path(out_dir, "fig3b_simulation_runs.pdf"),
       p_recovered / p_success / p_fdr + plot_layout(heights = c(0.4, 0.05, 0.6)),
       width = 11, height = 3.5)

# Fig. 3c: spiked and unspiked metabolites recovered per run
recovery <- bind_rows(
  run_summary %>% transmute(iter_id, tool = "FeatureMAP", Spiked = spiked_rec, Unspiked = hidden_rec),
  mummichog_run_summary %>% transmute(iter_id, tool = "Mummichog", Spiked = spiked_rec, Unspiked = hidden_rec)
) %>%
  pivot_longer(c(Spiked, Unspiked), names_to = "type", values_to = "n")

spiked_p <- recovery %>%
  filter(type == "Spiked") %>%
  pivot_wider(names_from = tool, values_from = n) %>%
  with(wilcox.test(FeatureMAP, Mummichog, paired = TRUE, exact = FALSE)$p.value)

recovery_mean <- recovery %>%
  group_by(tool, type) %>%
  summarise(mean = mean(n), se = sd(n) / sqrt(n()), .groups = "drop")
y_bracket <- max(recovery_mean$mean + recovery_mean$se) + 0.85

p_recovery <- ggplot(recovery_mean, aes(tool, mean, fill = type)) +
  geom_col(position = position_dodge(0.56), width = 0.46, color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(0.56), width = 0.1) +
  geom_text(aes(y = mean + se + 0.15, label = sprintf("%.2f", mean)), position = position_dodge(0.56)) +
  annotate("segment", x = 0.86, xend = 1.86, y = y_bracket, yend = y_bracket) +
  annotate("text", x = 1.36, y = y_bracket + 0.15, label = ifelse(spiked_p < 0.001, "***", signif(spiked_p, 2))) +
  scale_fill_manual(values = c(Spiked = "#D94F00", Unspiked = "#2171B5"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Recovered metabolites per run") +
  theme_base +
  theme(legend.position = c(0.8, 0.9))

ggsave(file.path(out_dir, "fig3c_recovered_metabolites.pdf"), p_recovery, width = 4, height = 3.5)

# Fig. 3d: recovery frequency of each pathway metabolite; overlapping nodes placed on a ring
node_recovery <- read_csv(file.path(out_dir, "node_recovery_frequency.csv")) %>%
  arrange(node_type, desc(fmap_freq), desc(mummi_freq), label) %>%
  group_by(node_type) %>%
  mutate(node_id = paste0(ifelse(node_type == "Spiked", "S", "U"), row_number())) %>%
  group_by(fmap_freq, mummi_freq) %>%
  mutate(n_overlap = n(),
         angle = 2 * pi * (row_number() - 1) / n_overlap,
         x = fmap_freq + ifelse(n_overlap > 1, 0.034 * cos(angle), 0),
         y = mummi_freq + ifelse(n_overlap > 1, 0.034 * sin(angle), 0)) %>%
  ungroup()

p_node <- ggplot(node_recovery) +
  geom_abline(linetype = "dashed", color = "grey70") +
  geom_point(data = ~ distinct(filter(.x, n_overlap > 1), fmap_freq, mummi_freq),
             aes(fmap_freq, mummi_freq), shape = 21, size = 16, fill = NA) +
  geom_point(aes(x, y, color = node_type), size = 6) +
  geom_text(aes(x, y, label = node_id), color = "white", size = 2.6) +
  scale_color_manual(values = c(Spiked = "#D94F00", Unspiked = "#2171B5"), name = NULL) +
  scale_x_continuous(labels = scales::percent, limits = c(-0.08, 1.08)) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.08, 1.08)) +
  coord_fixed() +
  labs(x = "FeatureMAP recovery frequency", y = "Mummichog recovery frequency") +
  theme_base +
  theme(legend.position = c(0.15, 0.9))

ggsave(file.path(out_dir, "fig3d_node_recovery.pdf"), p_node, width = 5, height = 5)

# Fig. 3e: success rate when >= k spiked metabolites are required
p_success_rate <- read_csv(file.path(out_dir, "success_rate.csv")) %>%
  ggplot(aes(k, success_rate, color = tool)) +
  geom_vline(xintercept = 5, linetype = "dashed") +
  geom_line() +
  geom_point() +
  geom_text(data = ~ filter(.x, k == 5), aes(label = scales::percent(success_rate, accuracy = 1)),
            hjust = -0.3, color = "black", size = 3) +
  scale_color_manual(values = c(FeatureMAP = "#2F7D4F", Mummichog = "#7A5A2A"), name = NULL) +
  scale_x_continuous(breaks = 1:7) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Required recovered spike-in nodes (k)", y = "Success rate") +
  theme_base +
  theme(legend.position = c(0.3, 0.5))

ggsave(file.path(out_dir, "fig3e_success_rate.pdf"), p_success_rate, width = 4, height = 3.5)
