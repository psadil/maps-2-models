get_subs_tasks <- function(dataset) {
  arrow::open_dataset(dataset) |>
    dplyr::distinct(sub, task) |>
    dplyr::collect()
}

do_loo_cor <- function(dataset, subtask) {
  .subtask <- arrow::as_arrow_table(
    subtask,
    schema = arrow::schema(
      arrow::field("sub", arrow::int32()),
      arrow::field("task", arrow::large_utf8())
    )
  )

  d <- arrow::open_dataset(dataset) |>
    dplyr::filter(
      mask,
      dplyr::between(`aparc+aseg`, 1000, 2999)
    )

  probe <- d |>
    dplyr::semi_join(.subtask, dplyr::join_by(sub, task)) |>
    dplyr::select(i, j, k, src = cope1, task)

  # not all subs have same voxels in GM
  # hence, right_join
  d |>
    dplyr::semi_join(dplyr::distinct(.subtask, task), dplyr::join_by(task)) |>
    dplyr::anti_join(dplyr::distinct(.subtask, sub), dplyr::join_by(sub)) |>
    dplyr::select(i, j, k, cope1, sub, task) |>
    dplyr::right_join(probe, dplyr::join_by(i, j, k, task)) |>
    dplyr::collect() |>
    dplyr::summarise(
      rho = cor(cope1, src),
      .by = c(sub, task)
    ) |>
    dplyr::mutate(probe = unique(dplyr::collect(.subtask)$sub))
}
