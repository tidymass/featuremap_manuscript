library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(mWISE)
library(mWiseData)
library(xMSannotator)

case_ids <- c("MTBLS13267", "MTBLS267", "MTBLS13855", "MTBLS13095", "MTBLS12329",
              "MTBLS13578", "MTBLS11987", "MTBLS417", "MTBLS5028", "MTBLS12854")

data_dir <- "2_data/ms1_based_annotation"
out_dir <- "3_data_analysis/ms1_based_annotation/tool_annotation"

# Peak table (variable_id, mz, rt, sample intensities) for one ionization mode
get_peak_table <- function(case_id, mode) {
  object_file <- file.path(data_dir, case_id, paste0("object_", mode, "2.rda"))
  if (!file.exists(object_file)) return(NULL)
  object <- load_rda_first(object_file)
  load_rda_first(file.path(data_dir, case_id, "feature_df.rda")) %>%
    filter(tolower(substr(polarity, 1, 3)) == mode) %>%
    dplyr::select(variable_id, mz, rt) %>%
    left_join(rownames_to_column(as.data.frame(object@expression_data), "variable_id"),
              by = "variable_id")
}

run_mwise <- function(peak_table, mode) {
  mWISE.annotation(
    Peak.List = dplyr::rename(peak_table, Peak.Id = variable_id),
    polarity = ifelse(mode == "pos", "positive", "negative"),
    diffusion.input.type = "binary",
    score = "z",
    Cpd.Add = Cpd.Add,
    graph.name = "RClass",
    Unique.Annotation = TRUE,
    Intensity.idx = 4:ncol(peak_table),
    nClust = 2,
    do.Par = TRUE
  )$Ranked.Tab %>%
    transmute(variable_id = Peak.Id, KEGG_ID = Compound, Ranking)
}

run_xms <- function(peak_table, mode, outloc) {
  dir.create(outloc, recursive = TRUE, showWarnings = FALSE)
  xms_input <- peak_table %>%
    column_to_rownames("variable_id") %>%
    dplyr::rename(time = rt) %>%
    as.matrix()
  num_nodes <<- 6
  multilevelannotation(
    xms_input,
    max.mz.diff = 10, max.rt.diff = 10, cormethod = "pearson", num_nodes = 6,
    queryadductlist = "all", mode = mode, db_name = "HMDB", adduct_weights = NA,
    num_sets = 3000, allsteps = TRUE, corthresh = 0.7, NOPS_check = TRUE,
    customIDs = NA, missing.value = NA, deepsplit = 2, networktype = "unsigned",
    outloc = outloc, minclustsize = 10, module.merge.dissimilarity = 0.2,
    filter.by = ifelse(mode == "pos", "M+H", "M-H"), redundancy_check = TRUE,
    min_ions_perchem = 1, biofluid.location = NA, origin = NA, status = NA,
    boostIDs = NA, max_isp = 5, MplusH.abundance.ratio.check = FALSE,
    customDB = NA, HMDBselect = "union", mass_defect_window = 0.01,
    mass_defect_mode = mode, pathwaycheckmode = "pm"
  )
  # Map xMSannotator output back to feature IDs by m/z and RT
  read_csv(file.path(outloc, "Stage5.csv")) %>%
    mutate(key = paste(round(mz, 5), round(time, 1), sep = "_")) %>%
    inner_join(peak_table %>%
                 transmute(variable_id, key = paste(round(mz, 5), round(rt, 1), sep = "_")),
               by = "key") %>%
    dplyr::select(variable_id, chemical_ID, Confidence, score)
}

Cpd.Add <- CpdaddPreparation(KeggDB = NULL, do.Par = TRUE, nClust = 2)

for (case_id in case_ids) {
  message(case_id)
  case_out <- file.path(out_dir, case_id)
  dir.create(case_out, recursive = TRUE, showWarnings = FALSE)

  peak_tables <- map(c(pos = "pos", neg = "neg"), ~ get_peak_table(case_id, .x)) %>%
    compact()

  imap_dfr(peak_tables, run_mwise) %>%
    write_csv(file.path(case_out, "mWISE_anno.csv"))

  imap_dfr(peak_tables, ~ run_xms(.x, .y, file.path(case_out, "xMSannotator", .y))) %>%
    write_csv(file.path(case_out, "xms_anno.csv"))
}
