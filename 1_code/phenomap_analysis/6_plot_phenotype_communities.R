library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/phenomap_analysis/phenomap_tools.R")
library(igraph)
library(scatterpie)
library(ggrepel)
library(patchwork)

similarity <- read_similarity(file.path(out_dir, "phenotype_similarity_matrix.csv"))
members <- read_csv(file.path(out_dir, "community_members.csv")) %>%
  filter(!is.na(Community))
communities <- sort(unique(members$Community))

# Fig. 5g (left): communities as pies of phenotype types; links = phenotype pairs above the network cutoff
cutoff <- quantile(similarity[upper.tri(similarity)], 0.95, na.rm = TRUE)
community_edges <- combn(communities, 2) %>%
  t() %>%
  as_tibble(.name_repair = ~ c("from", "to")) %>%
  mutate(block = map2(from, to, ~ similarity[members$Phenotype_ID[members$Community == .x],
                                             members$Phenotype_ID[members$Community == .y]]),
         mean_bma = map_dbl(block, mean, na.rm = TRUE),
         edge_n = map_int(block, ~ sum(.x >= cutoff, na.rm = TRUE)))

g <- graph_from_data_frame(community_edges %>% transmute(from = as.character(from), to = as.character(to), weight = mean_bma),
                           directed = FALSE, vertices = tibble(name = as.character(communities)))
set.seed(20260921)
xy <- apply(layout_with_fr(g, weights = E(g)$weight, niter = 3000), 2,
            function(v) (v - min(v)) / (max(v) - min(v)) * 2 - 1)

types <- intersect(names(type_fill), unique(members$Phenotype_Type))
nodes <- members %>%
  count(Community, Phenotype_Type) %>%
  pivot_wider(names_from = Phenotype_Type, values_from = n, values_fill = 0) %>%
  mutate(size = rowSums(across(all_of(types))),
         r = 0.16 * sqrt(size / max(size)),
         x = xy[match(as.character(Community), V(g)$name), 1],
         y = xy[match(as.character(Community), V(g)$name), 2])

edge_xy <- community_edges %>%
  filter(edge_n >= 1) %>%
  left_join(nodes %>% dplyr::select(from = Community, x, y), by = "from") %>%
  left_join(nodes %>% dplyr::select(to = Community, xend = x, yend = y), by = "to")

p_communities <- ggplot() +
  geom_segment(data = edge_xy, aes(x, y, xend = xend, yend = yend, linewidth = edge_n), color = "grey60") +
  geom_scatterpie(data = nodes, aes(x, y, r = r), cols = types, color = "white") +
  geom_text(data = nodes, aes(x, y - r - 0.06, label = paste0("C", Community)), size = 4.4) +
  scale_fill_manual(values = type_fill, name = NULL) +
  scale_linewidth(range = c(0.4, 1.8), name = "Cross-cluster links") +
  coord_fixed() +
  theme_void()

# Fig. 5g (right): bootstrap stability of each community
p_stability <- read_csv(file.path(out_dir, "community_stability.csv")) %>%
  mutate(Community = fct_reorder(paste0("C", Community), Bootstrap_Jaccard)) %>%
  ggplot(aes(Bootstrap_Jaccard, Community, fill = Bootstrap_Jaccard >= 0.8)) +
  geom_col(width = 0.68, color = "black") +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "#D64541") +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "grey70")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Bootstrap Jaccard", y = NULL) +
  theme_base

ggsave(file.path(out_dir, "fig5g_phenotype_communities.pdf"),
       p_communities + p_stability + plot_layout(widths = c(2, 1)),
       width = 11, height = 4.6)

# Fig. 5h: communities C7, C9 and C5; edges curated against the literature (Supplementary Data 6-8)
case_edges <- read_csv(file.path(out_dir, "case_cluster_edges.csv"))

layout_cluster <- function(cluster_id, offset) {
  cluster_nodes <- members %>% filter(Community == cluster_id)
  cluster_graph <- graph_from_data_frame(case_edges %>% filter(cluster == cluster_id) %>%
                                           transmute(from, to, weight = bma),
                                         directed = FALSE, vertices = cluster_nodes %>% dplyr::select(Phenotype_ID))
  set.seed(20260919)
  xy <- apply(layout_with_fr(cluster_graph, weights = E(cluster_graph)$weight, niter = 3000), 2,
              function(v) (v - min(v)) / (max(v) - min(v)) * 2 - 1)
  cluster_nodes %>%
    mutate(x = xy[match(Phenotype_ID, V(cluster_graph)$name), 1] + offset,
           y = xy[match(Phenotype_ID, V(cluster_graph)$name), 2])
}

node_xy <- map2_dfr(c(7, 9, 5), c(0, 3, 6), layout_cluster)
edge_xy <- case_edges %>%
  left_join(node_xy %>% dplyr::select(from = Phenotype_ID, x, y), by = "from") %>%
  left_join(node_xy %>% dplyr::select(to = Phenotype_ID, xend = x, yend = y), by = "to")

p_cases <- ggplot() +
  geom_segment(data = edge_xy, aes(x, y, xend = xend, yend = yend, linewidth = bma, linetype = class),
               color = "grey55") +
  geom_point(data = node_xy, aes(x, y, fill = Phenotype_Type), shape = 21, size = 7) +
  geom_text_repel(data = node_xy, aes(x, y, label = Phenotype), size = 3.9, point.size = 7,
                  min.segment.length = Inf, seed = 20260919) +
  scale_fill_manual(values = type_fill, name = "Phenotype type") +
  scale_linewidth(range = c(0.6, 1.4), name = "BMA similarity") +
  scale_linetype_manual(values = c(Direct = "solid", Indirect = "22"), name = "Literature support") +
  coord_fixed(clip = "off") +
  theme_void() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig5h_case_communities.pdf"), p_cases, width = 13, height = 5.6)
