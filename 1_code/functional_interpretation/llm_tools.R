library(httr2)
library(jsonlite)

openai_request <- function(endpoint, body) {
  request(paste0("https://api.openai.com/v1/", endpoint)) %>%
    req_auth_bearer_token(Sys.getenv("OPENAI_API_KEY")) %>%
    req_body_json(body) %>%
    req_retry(max_tries = 5) %>%
    req_timeout(120) %>%
    req_perform() %>%
    resp_body_json()
}

openai_embedding <- function(texts, model) {
  texts <- unname(as.character(texts))
  split(texts, ceiling(seq_along(texts) / 100)) %>%
    map(~ map(openai_request("embeddings", list(model = model, input = .x))$data,
              ~ unlist(.x$embedding))) %>%
    do.call(c, .)
}

cosine_sim <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))

normalize_title <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("\\([^)]*\\)", " ") %>%
    str_replace_all(",", " ") %>%
    str_squish()
}

tokenize_title <- function(x) {
  words <- str_split(str_squish(str_replace_all(str_to_lower(x), "[^a-z0-9]+", " ")), " ")[[1]]
  words[words != ""]
}

word_overlap <- function(gold, ai) {
  case_when(gold == ai ~ "Exact match",
            length(intersect(tokenize_title(gold), tokenize_title(ai))) > 0 ~ "Partial word overlap",
            TRUE ~ "No word overlap")
}

# ROUGE-L recall: longest common subsequence / reference length
rouge_l_recall <- function(reference, prediction) {
  a <- tokenize_title(reference)
  b <- tokenize_title(prediction)
  if (length(a) == 0 || length(b) == 0) return(0)
  dp <- matrix(0, length(a) + 1, length(b) + 1)
  for (i in seq_along(a)) for (j in seq_along(b)) {
    dp[i + 1, j + 1] <- if (a[i] == b[j]) dp[i, j] + 1 else max(dp[i, j + 1], dp[i + 1, j])
  }
  dp[length(a) + 1, length(b) + 1] / length(a)
}

# LLM judge (gpt-4o, two samples averaged); scores 0-4 on four axes
JUDGE_PROMPT <- paste(readLines("1_code/functional_interpretation/judge_prompt.txt"), collapse = "\n")

judge_summary <- function(reference, summary) {
  res <- openai_request("chat/completions", list(
    model = "gpt-4o",
    temperature = 0.2,
    n = 2,
    response_format = list(type = "json_object"),
    messages = list(
      list(role = "system", content = JUDGE_PROMPT),
      list(role = "user", content = sprintf("Text A (reference):\n%s\n\nText B (to evaluate):\n%s",
                                            reference, summary))
    )
  ))
  scores <- map(res$choices, ~ fromJSON(.x$message$content))
  axes <- c("topic_alignment", "entity_support", "mechanism_consistency", "unsupported_claims")
  axis_mean <- map_dbl(axes, function(a) mean(map_dbl(scores, a)))
  tibble(!!!set_names(as.list(axis_mean), axes),
         judge_overall = mean(axis_mean),
         judge_score = mean(axis_mean) / 4,
         judge_rationale = scores[[1]]$rationale)
}
