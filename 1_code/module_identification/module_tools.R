library(featuremap)
library(tidygraph)
library(igraph)
library(ggraph)

data_dir <- "2_data/module_identification"

# Background: 40 control samples, features passing missing/zero/intensity filters
load_background <- function() {
  object_pos <- load_rda_first(file.path(data_dir, "annotated_object_pos.rda"))
  object_neg <- load_rda_first(file.path(data_dir, "annotated_object_neg.rda"))
  sample_info <- readxl::read_xlsx(file.path(data_dir, "sampleID.xlsx")) %>%
    set_names(c("ID", "name", "group_name")) %>%
    filter(group_name == "Control_E")

  trim_names <- function(df) {
    colnames(df) <- substr(colnames(df), 10, 25)
    df
  }
  intensity <- bind_rows(trim_names(object_pos@expression_data),
                         trim_names(object_neg@expression_data)) %>%
    t() %>%
    as.data.frame()
  intensity <- intensity[sample_info$ID, ]
  keep <- colMeans(is.na(intensity)) < 0.2 &
    colMeans(intensity == 0 | is.na(intensity)) < 0.3 &
    apply(intensity, 2, median, na.rm = TRUE) > 1e3

  load("2_data/metabolic_reaction_network.rda")

  list(
    sample_info = sample_info,
    intensity = intensity[, keep],
    feature_meta = bind_rows(
      object_pos@variable_info %>% mutate(polarity = "positive"),
      object_neg@variable_info %>% mutate(polarity = "negative")
    ) %>%
      dplyr::select(variable_id, mz, rt, polarity),
    ms2 = load_rda_first(file.path(data_dir, "annotated_table_ms2.rda")) %>%
      group_by(variable_id) %>%
      slice_max(SS, n = 3, with_ties = FALSE) %>%
      ungroup(),
    mrn = metabolic_reaction_network,
    base_res = load_rda_first(file.path(data_dir, "featuremap_result.rda"))
  )
}

# MS/MS-supported features of a KEGG pathway, used as spike-in targets
pathway_targets <- function(bg, pathway_id) {
  pathway_kegg <- read_csv(file.path(data_dir, paste0(pathway_id, "_ipmap.csv"))) %>%
    pull(KEGG) %>%
    na.omit() %>%
    unique()
  pathway_node <- bg$mrn %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(KEGG_ID %in% pathway_kegg, !is.na(INCHI_ID)) %>%
    distinct(KEGG_ID, .keep_all = TRUE) %>%
    dplyr::select(INCHI_ID, Compound_name, KEGG_ID)
  target <- bg$ms2 %>%
    filter(KEGG.ID %in% pathway_kegg, variable_id %in% colnames(bg$intensity)) %>%
    distinct(variable_id, KEGG.ID, .keep_all = TRUE) %>%
    dplyr::select(variable_id, KEGG.ID, HMDB.ID, Compound.name, Total.score) %>%
    left_join(pathway_node, by = c("KEGG.ID" = "KEGG_ID")) %>%
    arrange(KEGG.ID, variable_id)
  spiked_inchi <- unique(na.omit(target$INCHI_ID))
  list(
    pathway_id = pathway_id,
    target = target,
    features = sort(unique(target$variable_id)),
    spiked_inchi = spiked_inchi,
    hidden_inchi = setdiff(pathway_node$INCHI_ID, spiked_inchi)
  )
}

# Spike-in targets saved by the simulation scripts
read_targets <- function(target_file, pathway_id, mrn) {
  target <- read_csv(target_file)
  pathway_kegg <- read_csv(file.path(data_dir, paste0(pathway_id, "_ipmap.csv")))$KEGG
  pathway_inchi <- mrn %>%
    activate(nodes) %>%
    as_tibble() %>%
    filter(KEGG_ID %in% pathway_kegg, !is.na(INCHI_ID)) %>%
    pull(INCHI_ID) %>%
    unique()
  spiked_inchi <- unique(na.omit(target$INCHI_ID))
  list(target = target, spiked_inchi = spiked_inchi, hidden_inchi = setdiff(pathway_inchi, spiked_inchi))
}

# Multiply case intensities by a feature-specific fold change (2-10) with log-normal noise
spike_features <- function(mat, features, case_idx, fc_seed, spike_seed) {
  set.seed(fc_seed)
  fc <- runif(length(features), 2, 10)
  set.seed(spike_seed)
  for (j in seq_along(features)) {
    mat[case_idx, features[j]] <- mat[case_idx, features[j]] *
      rlnorm(sum(case_idx), meanlog = log(fc[j]), sdlog = 0.15)
  }
  mat
}

# Welch t-test and log2 fold change for every feature
diff_test <- function(bg, mat, group) {
  stat <- map_dfr(colnames(mat), function(f) {
    v <- mat[, f]
    ok <- is.finite(v)
    p <- tryCatch(t.test(v[ok] ~ group[ok])$p.value, error = function(e) NA_real_)
    tibble(variable_id = f,
           condition = log2(mean(v[ok & group == "case"]) / mean(v[ok & group == "control"])),
           pvalue = p)
  })
  bg$feature_meta %>% inner_join(stat, by = "variable_id")
}

run_featuremap <- function(base_res, feature_table) {
  res <- base_res$clone(deep = TRUE)
  res$feature_table <- feature_table
  res <- extract_subnetwork(res, feature = "sig", workers = 4, times = 100,
                            b = 20, mu = 1e-4, w = 25, r = 0.3)
  identify_module(res, method = "louvain", min_size = 3, max_size = 1e5,
                  n_perm = 100, seed = 13579, workers = 4, node_freq = 0.05)
}

# Spiked / unspiked pathway metabolites recovered in each module
count_recovery <- function(module_table, target) {
  module_table %>%
    group_by(Module) %>%
    summarise(
      n_spiked = n_distinct(INCHI_ID[INCHI_ID %in% target$spiked_inchi]),
      n_hidden = n_distinct(INCHI_ID[INCHI_ID %in% target$hidden_inchi]),
      .groups = "drop"
    ) %>%
    mutate(n_pathway = n_spiked + n_hidden)
}

# Node frequency across selected modules (spiked, unspiked and their MRN neighbors)
module_node_table <- function(module_tables, modules, target, pathway_id, mrn) {
  kegg_name <- read_tsv(file.path(data_dir, pathway_id), col_names = c("KEGG_ID", "name")) %>%
    mutate(KEGG_ID = str_trim(KEGG_ID), name = str_trim(name))
  spiked <- target$target %>% filter(!is.na(INCHI_ID)) %>% distinct(INCHI_ID, .keep_all = TRUE)
  mrn_nodes <- mrn %>% activate(nodes) %>% as_tibble()
  target_vid <- which(mrn_nodes$INCHI_ID %in% c(target$spiked_inchi, target$hidden_inchi))
  neighbor_inchi <- mrn_nodes$INCHI_ID[unique(unlist(neighborhood(as.igraph(mrn), 1, target_vid)))]

  map_dfr(seq_len(nrow(modules)), function(j) {
    module_tables[[modules$iter_id[j]]] %>%
      filter(Module == modules$Module[j], !is.na(INCHI_ID)) %>%
      distinct(INCHI_ID, .keep_all = TRUE) %>%
      mutate(iter_id = modules$iter_id[j])
  }) %>%
    group_by(INCHI_ID) %>%
    summarise(n_success_runs = n_distinct(iter_id),
              KEGG_ID = first(na.omit(KEGG_ID)),
              Compound_name = first(na.omit(Compound_name)),
              .groups = "drop") %>%
    mutate(node_class = case_when(INCHI_ID %in% target$spiked_inchi ~ "Spiked",
                                  INCHI_ID %in% target$hidden_inchi ~ "Unspiked",
                                  TRUE ~ "Other")) %>%
    filter(node_class != "Other" | INCHI_ID %in% neighbor_inchi) %>%
    left_join(spiked %>% dplyr::select(INCHI_ID, spiked_kegg = KEGG.ID), by = "INCHI_ID") %>%
    mutate(KEGG_ID = coalesce(spiked_kegg, KEGG_ID),
           plot_id = ifelse(is.na(KEGG_ID), paste0("INCHI:", INCHI_ID), paste0("KEGG:", KEGG_ID))) %>%
    left_join(kegg_name, by = "KEGG_ID") %>%
    mutate(plot_label = ifelse(node_class == "Other", NA, coalesce(name, Compound_name))) %>%
    group_by(plot_id) %>%
    summarise(node_class = ifelse(any(node_class == "Spiked"), "Spiked",
                                  ifelse(any(node_class == "Unspiked"), "Unspiked", "Other")),
              KEGG_ID = first(KEGG_ID),
              Compound_name = first(na.omit(Compound_name)),
              plot_label = first(na.omit(plot_label)),
              n_success_runs = max(n_success_runs),
              .groups = "drop") %>%
    mutate(freq = n_success_runs / nrow(modules))
}

# Module network: node size = recovery frequency, edges from the MRN
plot_module_network <- function(nodes, mrn) {
  mrn_nodes <- mrn %>%
    activate(nodes) %>%
    as_tibble() %>%
    mutate(plot_id = case_when(paste0("KEGG:", KEGG_ID) %in% nodes$plot_id ~ paste0("KEGG:", KEGG_ID),
                               paste0("INCHI:", INCHI_ID) %in% nodes$plot_id ~ paste0("INCHI:", INCHI_ID)))
  edges <- mrn %>%
    activate(edges) %>%
    as_tibble() %>%
    transmute(a = mrn_nodes$plot_id[from], b = mrn_nodes$plot_id[to]) %>%
    filter(!is.na(a), !is.na(b), a != b) %>%
    transmute(from = pmin(a, b), to = pmax(a, b)) %>%
    distinct()
  graph <- graph_from_data_frame(edges, directed = FALSE, vertices = relocate(nodes, plot_id))

  set.seed(20240608)
  ggraph(graph, layout = "stress") +
    geom_edge_link(color = "grey78", linewidth = 0.4, alpha = 0.6) +
    geom_node_point(aes(size = freq, fill = node_class), shape = 21, color = "black") +
    geom_node_text(aes(label = plot_label), repel = TRUE, size = 5) +
    scale_fill_manual(values = c(Spiked = "#D94F00", Unspiked = "#2171B5", Other = "grey82"),
                      name = "Node type") +
    scale_size_continuous(range = c(3, 13), limits = c(0, 1), breaks = c(0.25, 0.5, 0.75, 1),
                          labels = scales::percent, name = "Frequency") +
    theme_void(base_size = 14)
}
