library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source("1_code/100_tools.R")
library(jsonlite)

out_dir <- "3_data_analysis/AD_case"

ad_query <- paste0('("Alzheimer Disease"[Mesh] OR "Alzheimer disease"[Title/Abstract] ',
                   'OR "Alzheimer\'s disease"[Title/Abstract] OR "Alzheimer dementia"[Title/Abstract])')

# Search terms from a module title; "X and Y metabolism" also searched as "X metabolism" and "Y metabolism"
module_terms <- function(title) {
  if (str_to_lower(title) == "ubiquinone and other terpenoid-quinone biosynthesis") {
    return(c(title, "Ubiquinone biosynthesis", "terpenoid-quinone biosynthesis"))
  }
  if (str_to_lower(title) == "pentose and glucuronate interconversions" ||
      !str_detect(str_to_lower(title), "metabolism$")) return(title)
  parts <- str_squish(str_split(str_remove(title, "\\s+metabolism$"), "\\s+(?:and|&)\\s+")[[1]])
  if (length(parts) <= 1) return(title)
  unique(c(title, ifelse(str_detect(str_to_lower(parts), "metabolism$"), parts, paste(parts, "metabolism"))))
}

pubmed_count <- function(title) {
  terms <- map_chr(module_terms(title), ~ sprintf('("%s"[Title/Abstract] OR "%s"[MeSH Terms])', .x, .x))
  query <- paste0(ad_query, " AND (", paste(terms, collapse = " OR "), ")")
  url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=0",
                "&api_key=", Sys.getenv("NCBI_API_KEY"), "&term=", URLencode(query, reserved = TRUE))
  Sys.sleep(0.15)
  as.integer(fromJSON(url)$esearchresult$count)
}

# Modules in CSF-iPSC correspondences, manually grouped into metabolic themes
read_csv(file.path(out_dir, "module_theme_classification.csv")) %>%
  mutate(pubmed_n = map_int(title, pubmed_count)) %>%
  write_csv(file.path(out_dir, "pubmed_module_counts.csv"))
