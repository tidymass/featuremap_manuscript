library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/module_identification/module_tools.R")
library(KEGGREST)
library(xml2)

out_dir <- "3_data_analysis/module_identification/arginine_biosynthesis"
load("2_data/metabolic_reaction_network.rda")

idmap <- read_csv(file.path(data_dir, "hsa00220_ipmap.csv")) %>%
  filter(!is.na(KEGG)) %>%
  distinct(KEGG, .keep_all = TRUE)
arginine <- read_targets(file.path(out_dir, "arginine_spiked_features.csv"), "hsa00220",
                         metabolic_reaction_network)

# Fig. 3a (left): reactions among arginine-biosynthesis metabolites from KEGG
kgml <- read_xml(keggGet("hsa00220", option = "kgml")[[1]])
reaction_ids <- c(xml_attr(xml_find_all(kgml, ".//reaction"), "name"),
                  xml_attr(xml_find_all(kgml, ".//entry[@reaction]"), "reaction")) %>%
  str_extract_all("R\\d{5}") %>%
  unlist() %>%
  unique()

map_dfr(reaction_ids, function(rid) {
  equation <- tryCatch(keggGet(paste0("rn:", rid))[[1]]$EQUATION, error = function(e) NULL)
  if (is.null(equation)) return(NULL)
  side <- str_split(equation, "<=>|=>|=", n = 2)[[1]]
  crossing(from = intersect(str_extract_all(side[1], "C\\d{5}")[[1]], idmap$KEGG),
           to = intersect(str_extract_all(side[2], "C\\d{5}")[[1]], idmap$KEGG))
}) %>%
  filter(from != to) %>%
  transmute(a = pmin(from, to), b = pmax(from, to)) %>%
  distinct(a, b) %>%
  rename(from = a, to = b) %>%
  write_csv(file.path(out_dir, "kegg_hsa00220_edges.csv"))

pathway_nodes <- idmap %>%
  transmute(name = KEGG, is_spiked = KEGG %in% arginine$target$KEGG.ID) %>%
  left_join(arginine$target %>% distinct(KEGG.ID, .keep_all = TRUE) %>%
              dplyr::select(name = KEGG.ID, label = Compound.name),
            by = "name")

set.seed(20260608)
p_pathway <- graph_from_data_frame(read_csv(file.path(out_dir, "kegg_hsa00220_edges.csv")),
                                   directed = FALSE, vertices = pathway_nodes) %>%
  ggraph(layout = "kk") +
  geom_edge_link(color = "grey60", linewidth = 0.4) +
  geom_node_point(aes(fill = is_spiked), shape = 21, size = 4.5) +
  geom_node_text(aes(label = label), repel = TRUE, size = 4) +
  scale_fill_manual(values = c(`FALSE` = "grey86", `TRUE` = "#D55E00"), guide = "none") +
  theme_void()

ggsave(file.path(out_dir, "fig3a_arginine_pathway.pdf"), p_pathway, width = 2.8, height = 3.5)

# Fig. 3a (right): merged module across successful simulations (>= 5 spiked metabolites)
success_modules <- read_csv(file.path(out_dir, "fmap_run_summary.csv")) %>%
  filter(spiked_rec >= 5) %>%
  dplyr::select(iter_id, Module = main_module)

module_node_table(readRDS(file.path(out_dir, "fmap_module_tables.rds")),
                  success_modules, arginine, "hsa00220", metabolic_reaction_network) %>%
  write_csv(file.path(out_dir, "arginine_module_nodes.csv"))

p_module <- plot_module_network(read_csv(file.path(out_dir, "arginine_module_nodes.csv")),
                                metabolic_reaction_network)

ggsave(file.path(out_dir, "fig3a_arginine_module_network.pdf"), p_module, width = 9, height = 7)
