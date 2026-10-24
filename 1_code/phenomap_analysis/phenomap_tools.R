data_dir <- "2_data/phenomap_analysis"
out_dir <- "3_data_analysis/phenomap_analysis"

# Module-module similarity from PhenoMAP "similar to" edges; unlinked pairs = 0
module_similarity_matrix <- function(edges, modules) {
  S <- matrix(0, length(modules), length(modules), dimnames = list(modules, modules))
  diag(S) <- 1
  edges <- edges %>% filter(`:START_ID` %in% modules, `:END_ID` %in% modules)
  S[cbind(edges$`:START_ID`, edges$`:END_ID`)] <- edges$Similarity_Score
  S[cbind(edges$`:END_ID`, edges$`:START_ID`)] <- edges$Similarity_Score
  S
}

# Best-match average between two module sets
bma <- function(a, b, S) {
  m <- S[unique(a), unique(b), drop = FALSE]
  mean(c(mean(apply(m, 1, max)), mean(apply(m, 2, max))))
}

read_similarity <- function(file) {
  as.matrix(read.csv(file, row.names = 1, check.names = FALSE))
}

p_label <- function(p) ifelse(p < 1e-4, "p < 0.0001", paste("p =", signif(p, 2)))

type_fill <- c(disease = "#A07BBF", exposure = "#5BB5A2", chemical = "#E8896B",
               intervention = "#D9B26F", trait = "#E7A9C4", biological_process = "#6FA3D2",
               cell_type_or_line = "#9C8A7A", genetic_perturbation = "#4E6E9E", symptom = "#B7B7B7")
