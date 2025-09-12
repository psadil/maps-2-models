get_ukb_copes <- function(src = "data-raw/ukb_copes") {
  tibble::tibble(UKB = readr::read_lines(src)) |>
    dplyr::filter(stringr::str_detect(UKB, "ses-2")) |>
    dplyr::mutate(
      Task = "EMOTION",
      ContrastName = "FACES-SHAPES",
      sub = stringr::str_extract(UKB, "[[:digit:]]{7}")
    )
}

get_ukb_copes <- function(src = "data-raw/ukb_copes") {
  tibble::tibble(UKB = readr::read_lines(src)) |>
    dplyr::mutate(
      Task = "EMOTION",
      ContrastName = "FACES-SHAPES",
      sub = stringr::str_extract(UKB, "[[:digit:]]{7}"),
      CopeNumber = 5
    )
}

sample_ukb <- function(test, n_iter, n_subs) {
  test |>
    dplyr::select(-ContrastName) |>
    tidyr::pivot_longer(UKB, names_to = "type") |>
    tidyr::crossing(
      iter = seq_len(n_iter),
      n_sub = n_subs
    ) |>
    dplyr::group_nest(Task, CopeNumber, type, iter, n_sub) |>
    dplyr::mutate(
      data = purrr::map2(
        data,
        n_sub,
        ~ dplyr::sample_n(.x, size = .y, replace = TRUE)
      )
    ) |>
    tidyr::unnest(data) |>
    dplyr::select(-value)
}


avg_roi_ukb <- function(test, n_parcels = 400) {
  vol <- test |>
    dplyr::mutate(VOL = purrr::map(UKB, roi_from_nifti, n_parcels)) |>
    tidyr::unnest(VOL) |>
    dplyr::select(
      Task,
      CopeNumber,
      ContrastName,
      sub,
      label,
      value
    ) |>
    dplyr::mutate(
      n_parcels = .env$n_parcels,
      type = "UKB"
    )
}


prep_ukb_pop <- function(test_ukb) {
  test_ukb |>
    dplyr::mutate(
      UKB = stringr::str_replace(
        UKB,
        "/fastscratch/myscratch/pssadil/ukb_mni/derivatives",
        Sys.getenv("UKBMNI")
      ) |>
        stringr::str_remove("_task-EMOTION_space-MNI152_contrast-5_cope") |>
        stringr::str_remove("sub-"),
      data = purrr::map(UKB, to_tbl)
    ) |>
    tidyr::unnest(data) |>
    dplyr::select(-UKB, -Task, -ContrastName, -CopeNumber) |>
    mask_gray()
}

get_glm_pop_ukb <- function(src) {
  pop <- arrow::read_parquet(src)

  tibble::tibble(
    Task = "EMOTION",
    CopeNumber = 5,
    type = "UKB",
    iter = 0,
    n_sub = unique(pop$n_sub),
    glm = fs::path(src)
  )
}
