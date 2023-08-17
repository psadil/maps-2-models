
r2_score <- function(y_true, y_pred){
  numerator <- sum((y_true - y_pred)^2) 
  denominator <-  sum((y_true - mean(y_true))^2 )
  1 - numerator / denominator
}

make_data_model_study_to_study <- function(dataset){
  arrow::open_dataset(dataset) |>
    dplyr::select(sub, fold, n_sub, task, y_hat) |>
    dplyr::collect() |>
    dplyr::group_nest(n_sub, task) |>
    dplyr::mutate(
      rex = purrr::map(
        data,
        ~ReX::lme_ICC_2wayR(data=.x$y_hat, subID = .x$sub, session = .x$fold) |>
          tibble::as_tibble())) |>
    dplyr::select(-data) |>
    tidyr::unnest(rex) |>
    dplyr::select(n_sub, task, tidyselect::starts_with("ICC")) |>
    tidyr::pivot_longer(tidyselect::starts_with("ICC"), names_to = "Model", values_to = "ICC") |>
    dplyr::mutate(
      Measure = dplyr::case_match(
        Model,
        "ICC.a" ~ "Agreement",
        "ICC.c" ~ "Consistency",
        "ICCk.a" ~ "Agreement",
        "ICCk.c" ~ "Consistency"
      ),
      Model = dplyr::case_match(
        Model,
        "ICC.a" ~ "ICC(2,1)",
        "ICC.c" ~ "ICC(2,1)",
        "ICCk.a" ~ "ICC(2,k)",
        "ICCk.c" ~ "ICC(2,k)"
      ))
}

make_data_model_study_to_study2 <- function(dataset){
  arrow::open_dataset(dataset) |>
    dplyr::select(sub, fold, n_sub, task, y_hat) |>
    dplyr::collect() |>
    summarise(s = sd(y_hat), .by = c(sub, n_sub, task)) |>
    dplyr::mutate(n_sub = factor(n_sub))
}

make_data_model_gold_gold_to_study <- function(
    dataset_all, 
    dataset,
    n_sub_gold_thresh = 200,
    method="spearman"){
  
  total <- arrow::open_dataset(dataset_all) |> 
    dplyr::collect() |>
    dplyr::mutate(n_sub = dplyr::n_distinct(sub), .by = task) |>
    dplyr::summarize(
      R2 = cor(g, y_hat, method = method),
      .by = c(fold, task, n_sub)
    ) |>
    dplyr::summarise(
      avg = mean(R2),
      sd = sd(R2),
      N = dplyr::n(),
      .by = c(task, n_sub)
    ) 
  
  arrow::open_dataset(dataset) |>
    dplyr::select(g, sub, fold, n_sub, task, y_hat) |>
    dplyr::collect() |>
    dplyr::summarize(
      R2 = cor(g, y_hat, method=method),
      .by = c(fold, task, n_sub)
    ) |>
    dplyr::summarise(
      avg = mean(R2),
      sd = sd(R2),
      N = dplyr::n(),
      .by = c(task, n_sub)
    ) |>
    dplyr::bind_rows(total) |>
    dplyr::mutate(
      type = dplyr::if_else(n_sub < n_sub_gold_thresh, "simulation", "gold"))
}

make_data_model_sub_to_sub <- function(parquet_files){
  difumos <- purrr::map_dfr(parquet_files, arrow::read_parquet, .id = "task") |>
    dplyr::mutate(task = stringr::str_extract(task, "(?<=task=)[[:alpha:]]+")) |>
    dplyr::filter(stringr::str_detect(task, "EMOTION", TRUE)) |>
    tidyr::pivot_longer(c(-sub, -task)) |> 
    tidyr::pivot_wider(names_from = sub) |>
    dplyr::select(-name) |>
    dplyr::group_nest(task) |>
    dplyr::mutate(
      rho = purrr::map(
        data,
        ~ .x |>
          corrr::correlate(method = "spearman") |> 
          corrr::shave() |> 
          corrr::stretch(na.rm = TRUE))
    ) |>
    dplyr::select(task, rho) |>
    tidyr::unnest(rho)
}

