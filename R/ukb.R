
get_ukb_copes <- function(src="data-raw/ukb_copes"){
  tibble::tibble(UKB=readr::read_lines(src)) |>
    dplyr::filter(stringr::str_detect(UKB, "ses-2")) |>
    dplyr::mutate(
      Task="EMOTION", 
      ContrastName="FACES-SHAPES",
      sub = stringr::str_extract(UKB, "[[:digit:]]{7}")
    )
}

get_ukb_copes <- function(src="data-raw/ukb_copes"){
  tibble::tibble(UKB=readr::read_lines(src)) |>
    dplyr::mutate(
      Task="EMOTION", 
      ContrastName="FACES-SHAPES",
      sub = stringr::str_extract(UKB, "[[:digit:]]{7}"),
      CopeNumber = 5 
    )
}

sample_ukb <- function(test, n_iter, n_subs){
  
  test |>
    dplyr::select(-ContrastName) |>
    tidyr::pivot_longer(UKB, names_to = "type") |>
    tidyr::crossing(
      iter = seq_len(n_iter),
      n_sub = n_subs) |>
    dplyr::group_nest(Task, CopeNumber, type, iter, n_sub) |>
    dplyr::mutate(
      data = purrr::map2(
        data, n_sub, 
        ~dplyr::sample_n(.x, size=.y, replace=TRUE))
    ) |>
    tidyr::unnest(data) |>
    dplyr::select(-value)
  
}


avg_roi_ukb <- function(test, n_parcels=400){
  
  vol <- test |>
    dplyr::mutate(VOL = purrr::map(UKB, roi_from_nifti, n_parcels))  |>
    tidyr::unnest(VOL) |>
    dplyr::select(
      Task, CopeNumber, ContrastName, sub, label, value
    ) |>
    dplyr::mutate(
      n_parcels = .env$n_parcels,
      type = "UKB"
    )
  
}
