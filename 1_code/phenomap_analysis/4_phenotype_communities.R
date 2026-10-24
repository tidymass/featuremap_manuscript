library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
source("1_code/phenomap_analysis/phenomap_tools.R")
library(igraph)

similarity <- read_similarity(file.path(out_dir, "phenotype_similarity_matrix.csv"))
diag(similarity) <- NA
ids <- rownames(similarity)
phenotype_info <- read_csv(file.path(out_dir, "phenotype_info.csv"))

# Phenotype network: keep pairs above the q-th quantile of BMA similarity
build_graph <- function(sim, q) {
  adj <- sim
  adj[is.na(adj) | adj < quantile(sim[upper.tri(sim)], q, na.rm = TRUE)] <- 0
  graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
}

# Consensus Louvain: co-assignment across reshuffled node orders, then Louvain on the consensus graph
consensus_communities <- function(g, resolution, runs) {
  nodes <- V(g)$name
  co <- Reduce(`+`, map(seq_len(runs), function(r) {
    m <- membership(cluster_louvain(permute(g, order(sample.int(length(nodes)))), resolution = resolution))[nodes]
    outer(m, m, "==")
  })) / runs
  co[co < 0.5] <- 0
  diag(co) <- 0
  if (sum(co) == 0) return(set_names(seq_along(nodes), nodes))
  m <- membership(cluster_louvain(graph_from_adjacency_matrix(co, mode = "undirected", weighted = TRUE, diag = FALSE),
                                  resolution = resolution))
  set_names(as.integer(m), nodes)
}

# Communities smaller than 3 phenotypes are left unassigned (0)
size_floor <- function(m) {
  small <- names(which(table(m) < 3))
  m[as.character(m) %in% small] <- 0L
  m
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# Bootstrap stability: best Jaccard match of each community in 85% subsamples
bootstrap_stability <- function(q, resolution, reference) {
  groups <- split(names(reference)[reference > 0], reference[reference > 0])
  scores <- map(groups, ~ numeric())
  for (b in 1:50) {
    picked <- sort(sample(seq_along(ids), floor(0.85 * length(ids))))
    sub <- size_floor(consensus_communities(build_graph(similarity[picked, picked], q), resolution, 20))
    sub_groups <- split(names(sub)[sub > 0], sub[sub > 0])
    if (length(sub_groups) == 0) next
    for (k in names(groups)) {
      present <- intersect(groups[[k]], ids[picked])
      if (length(present) > 0) scores[[k]] <- c(scores[[k]], max(map_dbl(sub_groups, ~ jaccard(present, .x))))
    }
  }
  map_dbl(scores, ~ if (length(.x) == 0) NA_real_ else mean(.x))
}

# Sweep edge density and Louvain resolution
set.seed(20260917)
sweep <- crossing(q = c(0.80, 0.85, 0.90, 0.95), resolution = c(0.6, 0.8, 1.0, 1.2, 1.5, 2.0))
result <- pmap(sweep, function(q, resolution) {
  message("q = ", q, ", resolution = ", resolution)
  g <- build_graph(similarity, q)
  m <- size_floor(consensus_communities(g, resolution, 50))
  stability <- bootstrap_stability(q, resolution, m)
  list(membership = m, stability = stability,
       summary = tibble(q, resolution,
                        n_communities = n_distinct(m[m > 0]),
                        coverage = mean(m > 0),
                        modularity = modularity(g, ifelse(m == 0, max(m) + seq_along(m), m), weights = E(g)$weight),
                        mean_stability = mean(stability, na.rm = TRUE)))
})
sweep <- map_dfr(result, "summary")
write_csv(sweep, file.path(out_dir, "community_sweep.csv"))

# Rule fixed in advance: coverage >= 70%, 3-20 communities, highest mean stability
chosen <- sweep %>%
  mutate(index = row_number()) %>%
  filter(coverage >= 0.7, between(n_communities, 3, 20)) %>%
  arrange(desc(mean_stability), desc(modularity)) %>%
  slice(1) %>%
  pull(index)
membership <- result[[chosen]]$membership
stability <- result[[chosen]]$stability

community_stability <- tibble(Community = as.integer(names(stability)),
                              Size = as.integer(table(membership)[names(stability)]),
                              Bootstrap_Jaccard = stability) %>%
  mutate(Reportable = Bootstrap_Jaccard >= 0.6) %>%
  arrange(desc(Bootstrap_Jaccard))
write_csv(community_stability, file.path(out_dir, "community_stability.csv"))

tibble(Phenotype_ID = names(membership), Community = na_if(membership, 0L)) %>%
  left_join(phenotype_info, by = "Phenotype_ID") %>%
  left_join(community_stability %>% dplyr::select(Community, Community_Stability = Bootstrap_Jaccard),
            by = "Community") %>%
  arrange(Community, desc(Module_N)) %>%
  write_csv(file.path(out_dir, "community_members.csv"))
