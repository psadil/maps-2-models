make_data_peak_sub_to_sub <- function(at, space_sub, gold_peaks) {
  at_ <- at |>
    dplyr::mutate(dplyr::across(c(x, y, z), as.integer))

  gold_peaks_ <- gold_peaks |>
    dplyr::select(Task, m) |>
    tidyr::unnest(m) |>
    dplyr::left_join(at_) |>
    dplyr::group_by(Task, label) |>
    dplyr::slice_max(order_by = Value, n = 1, with_ties = FALSE) |> # grab highest peak from each label
    dplyr::group_by(Task) |>
    dplyr::slice_max(order_by = Value, n = 10, with_ties = FALSE) |> # grab highest 10 peaks (distinct labels)
    dplyr::ungroup() |>
    dplyr::distinct(Task, x, y, z) |>
    dplyr::mutate(peak_i = 1:dplyr::n(), .by = Task)

  space_sub |>
    dplyr::filter(!is.na(study_ind)) |> # can happen in no voxels pass threshold
    dplyr::inner_join(gold_peaks_) |>
    dplyr::select(x = x.study, y = y.study, z = z.study, Task, peak_i) |>
    dplyr::group_nest(Task, peak_i) |>
    dplyr::mutate(
      d = purrr::map(
        data,
        ~ .x |>
          dist() |>
          as.matrix() |>
          corrr::as_cordf() |>
          corrr::shave() |>
          corrr::stretch(na.rm = TRUE) |>
          dplyr::rename(d = r)
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(d)
}


bootstrap_r <- function(.data, times = 2000) {
  boots <- boot::boot(
    data = .data$r,
    statistic = function(x, i) mean(x[i]),
    R = 2000
  )

  ci <- boot::boot.ci(boots, type = "perc")

  tibble::tibble(
    rr = boots$t0,
    lower = ci$percent[4],
    upper = ci$percent[5]
  )
}
