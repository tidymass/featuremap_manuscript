library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/phenomap_analysis/phenomap_tools.R")

plot_violin <- function(df, group, value, y_lab, colors) {
  p <- wilcox.test(df[[value]][df[[group]] == names(colors)[1]],
                   df[[value]][df[[group]] == names(colors)[2]],
                   alternative = "greater", exact = FALSE)$p.value
  ggplot(df, aes(.data[[group]], .data[[value]], fill = .data[[group]])) +
    geom_violin(scale = "width", width = 0.8, color = NA) +
    geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
    annotate("segment", x = 1, xend = 2, y = 1.02, yend = 1.02) +
    annotate("text", x = 1.5, y = 1.07, label = p_label(p)) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(limits = c(0, 1.1)) +
    labs(x = NULL, y = y_lab) +
    theme_base
}

# Fig. 5e: module-level similarity between independent studies of the same vs different phenotypes
study_pairs <- read_csv(file.path(out_dir, "study_level_similarity_pairs.csv")) %>%
  mutate(group = factor(ifelse(Same_Phenotype, "Within phenotype", "Between phenotypes"),
                        levels = c("Within phenotype", "Between phenotypes")))

ggsave(file.path(out_dir, "fig5e_study_level_similarity.pdf"),
       plot_violin(study_pairs, "group", "Study_Level_BMA_Similarity", "Study-level BMA similarity",
                   c(`Within phenotype` = "#A07BBF", `Between phenotypes` = "#DDD5E8")),
       width = 3.6, height = 4.2)

# Fig. 5f: disease pairs with vs without a shared CTD KEGG mechanism
ctd_pairs <- read_csv(file.path(out_dir, "disease_ctd_pairs.csv")) %>%
  mutate(group = factor(ifelse(Shared_Mechanism, "Shared mechanism", "Non-shared mechanism"),
                        levels = c("Shared mechanism", "Non-shared mechanism")))

ggsave(file.path(out_dir, "fig5f_ctd_shared_mechanism.pdf"),
       plot_violin(ctd_pairs, "group", "BMA_Similarity", "Disease BMA similarity",
                   c(`Shared mechanism` = "#5BB5A2", `Non-shared mechanism` = "#D5EBE5")),
       width = 3.6, height = 4.2)
