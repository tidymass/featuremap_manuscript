library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(tidymass)

# Thermo RAW files: <ipsc_raw>/<Amide|C18>/; samples 1-6 = control, 7-12 = AD, QC*, BLK*
raw_dir <- normalizePath("2_data/AD_case/ipsc_raw")
work_dir <- "3_data_analysis/AD_case/ipsc_processing"

docker <- function(...) {
  status <- system2("docker", c("run", "--rm", "--platform", "linux/amd64", "-v", paste0(raw_dir, ":/data"), ...))
  if (status != 0) stop("docker command failed")
}

sample_group <- function(file) {
  stem <- toupper(tools::file_path_sans_ext(basename(file)))
  case_when(grepl("^BLK[0-9]+$", stem) ~ "Blank",
            grepl("^QC[0-9]+$", stem) ~ "QC",
            stem %in% as.character(1:6) ~ "Control",
            stem %in% as.character(7:12) ~ "Case")
}

for (column in c("hilic", "rp")) {
  column_dir <- file.path(work_dir, column)
  raw_column <- ifelse(column == "hilic", "Amide", "C18")

  # RAW -> mzML (ThermoRawFileParser 1.4.5), then split by polarity (OpenMS 3.5.0 FileFilter)
  dir.create(file.path(raw_dir, "mzML", column), recursive = TRUE, showWarnings = FALSE)
  for (raw_file in list.files(file.path(raw_dir, raw_column), pattern = "\\.raw$", ignore.case = TRUE)) {
    docker("-e", "MONO_ENV_OPTIONS=--interpreter", "quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1",
           "ThermoRawFileParser.sh", paste0("-i=/data/", raw_column, "/", raw_file),
           paste0("-o=/data/mzML/", column), "-f=1", "-m=2", "-L=1")
  }
  for (mzml in list.files(file.path(raw_dir, "mzML", column), pattern = "\\.mzML$")) {
    group <- sample_group(mzml)
    if (is.na(group)) next
    for (mode in c("POS", "NEG")) {
      dir.create(file.path(raw_dir, "MS1", column, mode, group), recursive = TRUE, showWarnings = FALSE)
      docker("ghcr.io/openms/openms-executables:3.5.0", "FileFilter",
             "-in", paste0("/data/mzML/", column, "/", mzml),
             "-out", paste0("/data/MS1/", column, "/", mode, "/", group, "/", mzml),
             "-spectra:select_polarity", ifelse(mode == "POS", "positive", "negative"))
    }
  }

  # Peak detection and alignment (massprocesser), then cleaning (tidymass)
  objects <- map(c(POS = "positive", NEG = "negative"), function(polarity) {
    mode_dir <- file.path(raw_dir, "MS1", column, ifelse(polarity == "positive", "POS", "NEG"))
    process_data(path = mode_dir, polarity = polarity, ppm = 10, peakwidth = c(10, 60), threads = 6,
                 output_tic = FALSE, output_bpc = FALSE, output_rt_correction_plot = FALSE,
                 min_fraction = 0.5, group_for_figure = "QC")
    load(file.path(mode_dir, "Result/object"))
    object %>%
      activate_mass_dataset(what = "sample_info") %>%
      filter(class != "Blank") %>%
      impute_mv(method = "knn") %>%
      normalize_data(method = "median") %>%
      integrate_data(method = "subject_median")
  })

  dir.create(column_dir, recursive = TRUE, showWarnings = FALSE)
  save(objects, file = file.path(column_dir, "mass_dataset.rda"))
  map_dfr(objects, ~ .x@variable_info %>% dplyr::select(variable_id, mz, rt)) %>%
    write_csv(file.path(column_dir, "feature_info.csv"))
  map_dfr(objects, ~ .x@expression_data %>% rownames_to_column("variable_id")) %>%
    dplyr::select(-starts_with("QC")) %>%
    write_csv(file.path(column_dir, "feature_expression.csv"))
}
