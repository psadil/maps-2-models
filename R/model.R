r2_score <- function(y_true, y_pred) {
  numerator <- sum((y_true - y_pred)^2)
  denominator <- sum((y_true - mean(y_true))^2)
  1 - numerator / denominator
}

make_data_model_study_to_study <- function(dataset) {
  arrow::open_dataset(dataset) |>
    dplyr::filter(dimension == 64, model == "RIDGE_CV") |>
    dplyr::select(sub, n_sub, task, y_hat, measure, study, confounds) |>
    dplyr::collect() |>
    dplyr::group_nest(n_sub, task, measure, confounds) |>
    dplyr::mutate(
      rex = purrr::map(
        data,
        ~ .x |>
          tidyr::pivot_wider(names_from = study, values_from = y_hat) |>
          dplyr::select(-sub) |>
          irr::icc(model = "t", type = "a")
      ),
      icc_agreement = purrr::map_dbl(rex, purrr::pluck, "value"),
      lower_agreement = purrr::map_dbl(rex, purrr::pluck, "lbound"),
      upper_agreement = purrr::map_dbl(rex, purrr::pluck, "ubound"),
      rex = purrr::map(
        data,
        ~ .x |>
          tidyr::pivot_wider(names_from = study, values_from = y_hat) |>
          dplyr::select(-sub) |>
          irr::icc(model = "t", type = "consistency")
      ),
      icc_consistency = purrr::map_dbl(rex, purrr::pluck, "value"),
      lower_consistency = purrr::map_dbl(rex, purrr::pluck, "lbound"),
      upper_consistency = purrr::map_dbl(rex, purrr::pluck, "ubound")
    ) |>
    dplyr::select(-data, -rex) |>
    tidyr::pivot_longer(
      tidyselect::ends_with(
        c("agreement", "consistency")
      ),
      names_to = c("name", "type"),
      names_sep = "_"
    ) |>
    tidyr::pivot_wider()
}


make_data_model_gold_gold_to_study <- function(dataset_gold, dataset) {
  gold <- arrow::open_dataset(dataset_gold) |>
    dplyr::filter(dimension == 64, model == "RIDGE_CV") |>
    dplyr::distinct(statistic_rep, pvalue_rep, r2_rep_p, r2_rep, mae_rep, measure, task, sub, confounds) |>
    dplyr::collect() |>
    dplyr::mutate(n_sub = dplyr::n_distinct(sub), .by = c(task, measure, confounds)) |>
    dplyr::distinct(statistic_rep, pvalue_rep, r2_rep_p, r2_rep, mae_rep, measure, task, n_sub, confounds) |>
    dplyr::mutate(type = "gold")

  arrow::open_dataset(dataset) |>
    dplyr::filter(dimension == 64, model == "RIDGE_CV") |>
    dplyr::distinct(statistic_rep, pvalue_rep, r2_rep_p, r2_rep, mae_rep, measure, task, n_sub, study, confounds) |>
    dplyr::collect() |>
    dplyr::group_nest(n_sub, task, measure, confounds) |>
    dplyr::mutate(
      m = purrr::map2(
        n_sub, data,
        ~ meta::metacor(
          .y$statistic_rep,
          rep(.x, times = nrow(.y)),
          random = TRUE
        )
      ),
      avg = purrr::map_dbl(m, purrr::pluck, "TE.random"),
      lower = purrr::map_dbl(m, purrr::pluck, "lower.random"),
      upper = purrr::map_dbl(m, purrr::pluck, "upper.random"),
    ) |>
    dplyr::select(-data, -m) |>
    dplyr::mutate(type = "simulation") |>
    dplyr::bind_rows(gold)
}

make_data_model_gold_gold_to_study_ukb <- function(dataset_gold, dataset) {
  gold <- arrow::open_dataset(dataset_gold) |>
    dplyr::filter(dimension == 64, model == "RIDGE_CV") |>
    dplyr::select(-dimension, -model) |>
    dplyr::collect() |>
    dplyr::summarise(
      statistic_rep = cor(g, y_hat, method = "spearman"),
      .by = c(measure, task, confounds)
    ) |>
    dplyr::mutate(type = "gold")

  arrow::open_dataset(dataset) |>
    dplyr::filter(dimension == 64, model == "RIDGE_CV") |>
    dplyr::collect() |>
    dplyr::summarise(
      statistic_rep = cor(g, y_hat, method = "spearman"),
      .by = c(measure, task, model, dimension, study, n_sub, confounds)
    ) |>
    dplyr::group_nest(n_sub, task, measure, confounds) |>
    dplyr::mutate(
      m = purrr::map2(
        n_sub, data,
        ~ meta::metacor(
          .y$statistic_rep,
          rep(.x, times = nrow(.y)),
          random = TRUE
        )
      ),
      avg = purrr::map_dbl(m, purrr::pluck, "TE.random"),
      lower = purrr::map_dbl(m, purrr::pluck, "lower.random"),
      upper = purrr::map_dbl(m, purrr::pluck, "upper.random"),
    ) |>
    dplyr::select(-data, -m) |>
    dplyr::mutate(type = "simulation") |>
    dplyr::bind_rows(gold)
}

make_data_model_sub_to_sub <- function(features) {
  arrow::open_dataset(features) |>
    dplyr::filter(dimension == 64) |>
    dplyr::select(-dimension) |>
    dplyr::collect() |>
    tidyr::pivot_longer(c(-sub, -task, -confounds)) |>
    tidyr::pivot_wider(names_from = sub) |>
    dplyr::group_nest(task, confounds) |>
    dplyr::mutate(
      rho = purrr::map(
        data,
        ~ .x |>
          dplyr::select(-name) |>
          corrr::correlate(method = "spearman", quiet = TRUE) |>
          corrr::shave() |>
          corrr::stretch(na.rm = TRUE)
      )
    ) |>
    dplyr::select(task, rho, confounds) |>
    tidyr::unnest(rho)
}
