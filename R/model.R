r2_score <- function(y_true, y_pred) {
  numerator <- sum((y_true - y_pred)^2)
  denominator <- sum((y_true - mean(y_true))^2)
  1 - numerator / denominator
}

make_data_model_study_to_study <- function(dataset, measures) {
  files <- fs::dir_ls(
    dataset,
    recurse = TRUE,
    glob = "*parquet"
  ) |>
    stringr::str_subset("replacement=True") |>
    stringr::str_subset("UKB_SMALL", TRUE) |>
    stringr::str_subset("model=RIDGE_CV")

  ms <- duckplyr::read_parquet_duckdb(files) |>
    dplyr::distinct(measure) |>
    dplyr::filter(measure %in% measures) |>
    dplyr::collect() |>
    purrr::pluck("measure")

  out <- list()
  for (m in ms) {
    out[[m]] <- duckplyr::read_parquet_duckdb(files, prudence = "lavish") |>
      dplyr::filter(measure == .env$m) |>
      dplyr::select(
        sub,
        n_sub,
        task,
        y_hat,
        measure,
        study,
        confounds,
        type,
        model,
        replacement
      ) |>
      dplyr::collect() |>
      na.omit() |>
      dplyr::group_nest(
        n_sub,
        task,
        measure,
        confounds,
        type,
        model,
        replacement
      ) |>
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
        upper_consistency = purrr::map_dbl(rex, purrr::pluck, "ubound"),
        rex = purrr::map(
          data,
          ~ ReX::lme_ICC_2wayM(.x$y_hat, .x$sub, .x$study)
        ),
        sigma2_b = purrr::map_dbl(
          rex,
          ~ .x[1, "sigma2_b"]
        ),
        sigma2_w = purrr::map_dbl(
          rex,
          ~ .x[1, "sigma2_w"]
        ),
        var.data = purrr::map_dbl(
          rex,
          ~ .x[1, "var.data"]
        )
      ) |>
      dplyr::select(-data, -rex) |>
      tidyr::pivot_longer(
        tidyselect::ends_with(
          c("agreement", "consistency")
        ),
        names_to = c("name", "method"),
        names_sep = "_"
      ) |>
      tidyr::pivot_wider() |>
      dplyr::select(-measure)
  }
  dplyr::bind_rows(out, .id = "measure")
}


make_data_model_gold_gold_to_study <- function(
  dataset_gold,
  dataset,
  measures
) {
  gold <- arrow::open_dataset(dataset_gold) |>
    na.omit() |>
    dplyr::distinct(
      statistic_rep,
      measure,
      task,
      sub,
      confounds,
      type,
      model
    ) |>
    dplyr::collect() |>
    dplyr::mutate(
      n_sub = dplyr::n_distinct(sub),
      .by = c(task, measure, confounds, type, model)
    ) |>
    dplyr::distinct(
      avg = statistic_rep,
      measure,
      task,
      n_sub,
      confounds,
      type,
      model
    ) |>
    dplyr::mutate(sim = "gold", replacement = "False")

  arrow::open_dataset(dataset) |>
    dplyr::distinct(
      statistic_rep,
      measure,
      task,
      n_sub,
      study,
      confounds,
      type,
      model,
      replacement
    ) |>
    dplyr::collect() |>
    na.omit() |>
    dplyr::summarise(
      avg = mean(statistic_rep),
      sem = sd(statistic_rep) / sqrt(dplyr::n()),
      lower = quantile(statistic_rep, 0.025),
      upper = quantile(statistic_rep, 0.975),
      .by = c(n_sub, task, measure, confounds, type, model, replacement)
    ) |>
    dplyr::mutate(sim = "simulation") |>
    dplyr::bind_rows(gold) |>
    dplyr::filter(measure %in% measures)
}

make_data_model_gold_gold_to_study_ukb <- function(
  dataset_gold,
  dataset,
  measures
) {
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
    na.omit() |>
    dplyr::collect() |>
    dplyr::summarise(
      statistic_rep = cor(g, y_hat, method = "spearman"),
      .by = c(measure, task, model, dimension, study, n_sub, confounds)
    ) |>
    dplyr::group_nest(n_sub, task, measure, confounds, model) |>
    dplyr::mutate(
      m = purrr::map2(
        n_sub,
        data,
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

make_data_model_gold_gold_to_study2 <- function(
  dataset_gold,
  dataset,
  measures
) {
  # prediction significance
  gold <- arrow::open_dataset(dataset_gold) |>
    dplyr::distinct(pvalue_rep, measure, task, sub, confounds, type, model) |>
    na.omit() |>
    dplyr::collect() |>
    dplyr::mutate(
      n_sub = dplyr::n_distinct(sub),
      .by = c(task, measure, confounds, type, model)
    ) |>
    dplyr::distinct(
      sig = pvalue_rep,
      measure,
      task,
      n_sub,
      confounds,
      type,
      model
    ) |>
    dplyr::mutate(
      sim = "gold",
      avg = ceiling(sig < 0.05),
      replacement = "False"
    ) |>
    dplyr::select(-sig)

  arrow::open_dataset(dataset) |>
    dplyr::distinct(
      pvalue_rep,
      measure,
      task,
      n_sub,
      study,
      confounds,
      type,
      model,
      replacement
    ) |>
    dplyr::mutate(sig = pvalue_rep < 0.05) |>
    dplyr::select(-pvalue_rep) |>
    na.omit() |>
    dplyr::collect() |>
    dplyr::summarise(
      avg = mean(sig),
      lower = qbeta(0.025, 1 / 2 + sum(sig), dplyr::n() - sum(sig) + 1 / 2),
      upper = qbeta(0.975, 1 / 2 + sum(sig), dplyr::n() - sum(sig) + 1 / 2),
      .by = c(n_sub, task, measure, confounds, type, model, replacement)
    ) |>
    dplyr::mutate(sim = "simulation") |>
    dplyr::bind_rows(gold) |>
    dplyr::filter(measure %in% measures)
}

make_data_model_gold_gold_to_study3 <- function(
  dataset_gold,
  dataset,
  measures
) {
  # features
  gold <- duckplyr::read_parquet_duckdb(
    fs::dir_ls(
      dataset_gold,
      recurse = TRUE,
      glob = "*parquet"
    ),
    prudence = "lavish"
  )
  ms <- gold |>
    dplyr::distinct(measure) |>
    dplyr::filter(measure %in% measures) |>
    dplyr::collect() |>
    purrr::pluck("measure")

  mm <- list()
  for (m in ms) {
    mm[[m]] <- duckplyr::read_parquet_duckdb(
      fs::dir_ls(
        dataset,
        recurse = TRUE,
        glob = "*parquet"
      ),
      prudence = "lavish"
    ) |>
      dplyr::filter(measure == m) |>
      dplyr::left_join(
        gold,
        by = dplyr::join_by(index, type, measure, confounds, task, model)
      ) |>
      na.omit() |>
      dplyr::collect() |>
      dplyr::summarise(
        .estimate = cor(coef.x, coef.y, method = "spear"),
        .by = c(
          type,
          measure,
          confounds,
          task,
          study,
          n_sub,
          model,
          replacement
        )
      ) |>
      na.omit() |>
      dplyr::summarise(
        .lower = quantile(.estimate, 0.025),
        .upper = quantile(.estimate, 0.975),
        .estimate = tanh(mean(atanh(.estimate))),
        .by = c(type, measure, confounds, task, n_sub, model, replacement)
      )
  }
  dplyr::bind_rows(mm)
}

do_iccbin <- function(.data) {
  tryCatch(
    expr = {
      aod::iccbin(n = n, y = y, data = .data, method = "B")@rho[1]
    },
    error = function(e) {
      NA_real_
    }
  )
}


make_data_model_study_to_study2 <- function(
  dataset,
  measures,
  n_boot = 100,
  n_workers = 8
) {
  arrow::open_dataset(dataset) |>
    dplyr::distinct(
      pvalue_rep,
      measure,
      task,
      n_sub,
      study,
      confounds,
      type,
      model,
      replacement
    ) |>
    dplyr::filter(measure %in% measures, model == "RIDGE_CV") |>
    na.omit() |>
    dplyr::summarise(
      n = dplyr::n(),
      y = sum(pvalue_rep < 0.05),
      .by = c(n_sub, task, type, measure, confounds, study, model, replacement)
    ) |>
    dplyr::collect() |>
    dplyr::group_nest(
      task,
      measure,
      confounds,
      n_sub,
      type,
      model,
      replacement
    ) |>
    dplyr::mutate(
      .estimate = purrr::map_dbl(data, do_iccbin)
    ) |>
    dplyr::select(-data)
}

do_lme_ICC_2wayM <- function(.data) {
  tmp <- .data |>
    tidyr::pivot_longer(-index)

  out <- suppressMessages(ReX::lme_ICC_2wayM(tmp$value, tmp$index, tmp$name))
  out
}

make_data_model_study_to_study3 <- function(dataset, measures) {
  files <- fs::dir_ls(dataset, recurse = TRUE, glob = "*parquet")

  ms <- duckplyr::read_parquet_duckdb(files) |>
    dplyr::distinct(measure) |>
    dplyr::filter(measure %in% measures) |>
    dplyr::collect() |>
    purrr::pluck("measure")

  mm <- list()
  for (m in ms) {
    mm[[m]] <- duckplyr::read_parquet_duckdb(files, prudence = "lavish") |>
      dplyr::filter(measure == .env$m) |>
      dplyr::collect() |>
      na.omit() |>
      tidyr::pivot_wider(names_from = study, values_from = coef) |>
      dplyr::group_nest(
        type,
        measure,
        confounds,
        model,
        task,
        n_sub,
        replacement
      ) |>
      dplyr::mutate(
        fit = purrr::map(
          data,
          ~ .x |>
            dplyr::select(-index) |>
            irr::icc(model = "t", type = "c")
        ),
        .estimate = purrr::map_dbl(fit, ~ .x$value),
        .lower = purrr::map_dbl(fit, ~ .x$lbound),
        .upper = purrr::map_dbl(fit, ~ .x$ubound)
      ) |>
      dplyr::select(-fit, -data)
  }
  dplyr::bind_rows(mm)
}


make_data_model_sub_to_sub <- function(features) {
  arrow::open_dataset(features) |>
    dplyr::collect() |>
    tidyr::pivot_longer(c(-sub, -task, -confounds, -type)) |>
    tidyr::pivot_wider(names_from = sub) |>
    dplyr::group_nest(task, confounds, type) |>
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


make_data_model_gold_gold_to_study_popsize <- function(
  dataset_gold,
  dataset,
  measures
) {
  gold <- arrow::open_dataset(dataset_gold) |>
    na.omit() |>
    dplyr::distinct(statistic_rep, sub, model, popsize) |>
    dplyr::collect() |>
    dplyr::mutate(n_sub = dplyr::n_distinct(sub), .by = c(model, popsize)) |>
    dplyr::distinct(avg = statistic_rep, n_sub, popsize, model) |>
    dplyr::mutate(sim = "gold", replacement = "False")

  arrow::open_dataset(dataset) |>
    dplyr::distinct(statistic_rep, n_sub, study, popsize, model, replacement) |>
    dplyr::collect() |>
    na.omit() |>
    dplyr::summarise(
      avg = mean(statistic_rep),
      sem = sd(statistic_rep) / sqrt(dplyr::n()),
      lower = quantile(statistic_rep, 0.025),
      upper = quantile(statistic_rep, 0.975),
      .by = c(n_sub, model, replacement, popsize)
    ) |>
    dplyr::mutate(sim = "simulation") |>
    dplyr::bind_rows(gold) |>
    dplyr::filter(measure %in% measures)
}

make_data_model_study_to_study_popsize <- function(dataset, measures) {
  arrow::open_dataset(dataset) |>
    na.omit() |>
    dplyr::select(sub, y_hat, study, model, replacement, popsize) |>
    dplyr::collect() |>
    dplyr::group_nest(model, replacement, popsize) |>
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
      upper_consistency = purrr::map_dbl(rex, purrr::pluck, "ubound"),
      rex = purrr::map(
        data,
        ~ ReX::lme_ICC_2wayM(.x$y_hat, .x$sub, .x$study)
      ),
      sigma2_b = purrr::map_dbl(
        rex,
        ~ .x[1, "sigma2_b"]
      ),
      sigma2_w = purrr::map_dbl(
        rex,
        ~ .x[1, "sigma2_w"]
      ),
      var.data = purrr::map_dbl(
        rex,
        ~ .x[1, "var.data"]
      )
    ) |>
    dplyr::select(-data, -rex) |>
    tidyr::pivot_longer(
      tidyselect::ends_with(
        c("agreement", "consistency")
      ),
      names_to = c("name", "method"),
      names_sep = "_"
    ) |>
    tidyr::pivot_wider() |>
    dplyr::filter(measure %in% measures)
}


make_data_model_gold_gold_to_study_r2 <- function(
  dataset_gold,
  dataset,
  measures
) {
  gold <- arrow::open_dataset(dataset_gold) |>
    na.omit() |>
    dplyr::distinct(
      r2_rep,
      measure,
      task,
      sub,
      confounds,
      type,
      model
    ) |>
    dplyr::collect() |>
    dplyr::mutate(
      n_sub = dplyr::n_distinct(sub),
      .by = c(task, measure, confounds, type, model)
    ) |>
    dplyr::distinct(
      avg = r2_rep,
      measure,
      task,
      n_sub,
      confounds,
      type,
      model
    ) |>
    dplyr::mutate(sim = "gold", replacement = "False")

  arrow::open_dataset(dataset) |>
    dplyr::distinct(
      r2_rep,
      measure,
      task,
      n_sub,
      study,
      confounds,
      type,
      model,
      replacement
    ) |>
    dplyr::collect() |>
    na.omit() |>
    dplyr::summarise(
      avg = mean(r2_rep),
      sem = sd(r2_rep) / sqrt(dplyr::n()),
      lower = quantile(r2_rep, 0.025),
      upper = quantile(r2_rep, 0.975),
      .by = c(n_sub, task, measure, confounds, type, model, replacement)
    ) |>
    dplyr::mutate(sim = "simulation") |>
    dplyr::bind_rows(gold) |>
    dplyr::filter(measure %in% measures)
}

get_measures <- function(
  ukb_f = "data-raw/cognitive.parquet",
  hcp_f = "data-raw/hcp.parquet",
  ukb_n = 15000,
  hcp_n = 150
) {
  ukb <- duckplyr::read_parquet_duckdb(ukb_f) |>
    dplyr::summarise(
      dplyr::across(
        tidyselect::starts_with("f."),
        ~ sum(!is.na(.x))
      )
    ) |>
    dplyr::collect() |>
    tidyr::pivot_longer(tidyselect::everything()) |>
    dplyr::filter(value > ukb_n)

  hcp <- duckplyr::read_parquet_duckdb(hcp_f) |>
    dplyr::select(
      -sub,
      -Release,
      -Acquisition,
      -Gender,
      -Age,
      -QC_Issue,
      -tidyselect::starts_with(c("3T", "7T", "MEG", "fMRI_")),
      -tidyselect::ends_with(c("_Count", "_Compl", "_3T", "_7T", "_Raw")),
      -tidyselect::contains(c("_3T_", "Peak")),
      -tidyselect::matches(c("_Comp[[:digit:]]")),
    ) |>
    dplyr::summarise(
      dplyr::across(
        tidyselect::everything(),
        ~ sum(!is.na(.x))
      )
    ) |>
    dplyr::collect() |>
    tidyr::pivot_longer(tidyselect::everything()) |>
    dplyr::filter(value > hcp_n)

  vctrs::vec_c(hcp$name, ukb$name)
}
