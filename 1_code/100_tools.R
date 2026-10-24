library(tidyverse)
library(ggplot2)

load_rda_first <- function(path) {
  env <- new.env()
  name <- load(path, envir = env)
  env[[name[1]]]
}

normalize_hmdb <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  idx <- which(grepl("^HMDB\\d+", x))
  out[idx] <- sprintf("HMDB%07d", as.integer(gsub("\\D", "", x[idx])))
  out
}

normalize_kegg <- function(x) {
  x <- as.character(x)
  x[x %in% c("", "<NA>")] <- NA_character_
  x
}

normalize_id_cols <- function(df) {
  df %>%
    mutate(HMDB.ID = normalize_hmdb(HMDB.ID),
           KEGG.ID = normalize_kegg(KEGG.ID))
}

# MS/MS annotations used as reference: one row per feature-metabolite pair
prepare_ms2_reference <- function(ms2_table, normalize_id = TRUE) {
  ms2_table <- ms2_table %>%
    filter(!is.na(HMDB.ID) | !is.na(KEGG.ID)) %>%
    separate_rows(HMDB.ID, sep = "\\{\\}") %>%
    separate_rows(KEGG.ID, sep = "\\{\\}")
  if (normalize_id) ms2_table <- normalize_id_cols(ms2_table)
  ms2_table %>%
    filter(!is.na(HMDB.ID) | !is.na(KEGG.ID)) %>%
    distinct(variable_id, HMDB.ID, KEGG.ID, .keep_all = TRUE) %>%
    dplyr::select(variable_id, HMDB.ID, KEGG.ID, Total.score, Level)
}

theme_base <- theme_classic(base_size = 15) +
  theme(
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.65, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black")
  )

sig_label <- function(p) {
  case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "ns")
}
