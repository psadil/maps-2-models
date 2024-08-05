get_subs_tasks <- function(dataset, contrasts) {
  arrow::open_dataset(dataset) |>
    dplyr::distinct(sub, task) |> 
    dplyr::filter(
      !(sub %in% not_avail()),
      task %in% unique(contrasts$Task)) |>
    dplyr::collect()
}

do_loo_cor <- function(dataset, subtask, method="spearman") {
  
  m <- MNITemplate::getMNISegPath(res = "2mm") |>
    to_tbl(measure="seg") |>
    dplyr::filter(seg == 2) |>
    dplyr::mutate(i = x - 1, j = y - 1, k = z - 1) |>
    dplyr::select(i, j, k) |>
    arrow::as_arrow_table(
      schema = arrow::schema(
        arrow::field("i", arrow::uint8()),
        arrow::field("j", arrow::uint8()),
        arrow::field("k", arrow::uint8())
      )
    )
  .task <- unique(subtask$task)
  .sub <- unique(subtask$sub)
  
  d <- arrow::open_dataset(dataset) |>
    dplyr::select(i, j, k, cope1, sub, task) |>
    dplyr::filter(task %in% .task) |>
    dplyr::semi_join(m, by = dplyr::join_by(i, j, k))
  
  probe <- arrow::open_dataset(glue::glue("{dataset}/sub={.sub}")) |>
    dplyr::select(i, j, k, src = cope1, task)
  
  d |>
    dplyr::filter(sub != .sub) |>
    dplyr::left_join(probe, dplyr::join_by(i, j, k, task)) |>
    dplyr::collect() |>
    dplyr::summarise(
      rho = cor(cope1, src, method = method),
      .by = c(sub, task)
    ) |>
    dplyr::mutate(probe = .sub)
}


do_loo_cor2 <- function(dataset, subtask, method="spearman") {
  
  m <- MNITemplate::getMNISegPath(res = "2mm") |>
    to_tbl(measure="seg") |>
    dplyr::filter(seg == 2) |>
    dplyr::mutate(i = x - 1, j = y - 1, k = z - 1) |>
    dplyr::select(i, j, k) |>
    polars::as_polars_lf(
      schema = list(
        i = polars::pl$UInt8,
        j = polars::pl$UInt8,
        k = polars::pl$UInt8
      )
    )
  .task <- unique(subtask$task)
  .sub <- unique(subtask$sub)
  .cols <- c("i", "j", "k", "cope1", "sub", "task")
  d <- polars::pl$scan_parquet(fs::path(dataset, "*/*parquet"))$select(
    .cols
  )$join(m, how="semi", on = c("i", "j", "k"))
  
  probe <- polars::pl$scan_parquet(
    glue::glue("{dataset}/sub={.sub}/*parquet")
  )$select(.cols)$drop("sub")$rename(cope1="src")
  
  d$filter(
    polars::pl$col("sub")$neq(.sub)
  )$join(probe$select("task")$unique(), how="semi", on="task"
  )$join(
    probe, how = "left", on = c("i", "j", "k", "task")
  )$collect()$to_data_frame() |>
    dplyr::summarise(
      rho = cor(cope1, src, method = method),
      .by = c(sub, task)
    ) |>
    dplyr::mutate(probe = .sub)
  
  
  # d$filter(polars::pl$col("sub")$neq(.sub))$join(
  #   probe, how = "left", on = c("i", "j", "k", "task")
  # )$group_by(
  #   "sub", "task"
  # )$agg(
  #   rho = polars::pl$corr("cope1", "src", method = method)
  # )$collect()$to_data_frame() |>
  #   tibble::as_tibble() |>
  #   dplyr::mutate(probe = .sub)
}


cor_pairwise_ptfce <- function(tfce, ContrastName, n_sub, storage_dir, method = "spearman") {
  tfce <- tfce |>
    dplyr::filter(
      .data$n_sub == .env$n_sub,
      .data$ContrastName == .env$ContrastName
    )
  
  tmp <- tfce |>
    dplyr::select(Task, CopeNumber, ContrastName, n_sub, ptfce, iter) |>
    dplyr::mutate(
      data2 = purrr::map(
        ptfce,
        ~ to_tbl0(
            qs::qread(fs::path(storage_dir, fs::path_file(.x)))$Z_raw) |> mask_gray()
      )
    ) |>
    tidyr::unnest(data2) |>
    dplyr::select(-ptfce) |>
    dplyr::mutate(value = z_to_g(value, n=n_sub)) |>
    tidyr::pivot_wider(names_from = iter, values_from = value)
  
  tmp |>
    dplyr::distinct(Task, CopeNumber, ContrastName, n_sub) |>
    dplyr::mutate(
      rhos = list(
        tmp |>
          dplyr::select(tidyselect::matches("[[:digit:]]+")) |>
          corrr::correlate(method = .env$method, quiet = TRUE) |>
          corrr::stretch()
      ),
      method = .env$method
    ) |>
    tidyr::unnest(rhos)
}