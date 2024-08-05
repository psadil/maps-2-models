not_avail <- function() {
  c(
    110613, 113417, 113821, 120010, 121719, 130518, 139637, 143830, 146836,
    168139, 175035, 176239, 185038, 189652, 199958, 201515, 202820, 385046,
    401422, 415837, 433839, 462139, 465852, 469961, 644246, 656657, 688569,
    723141, 767464, 872764, 943862, 965367, 969476, 987983, 994273, 433839
  )
}

cor_w_pop_by_region <- function(tfce, tfce_pop, at, storage_dir, method = "spearman") {
  tfce <- tfce |>
    dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
  tfce_pop <- tfce_pop |>
    dplyr::semi_join(tfce, by = c("ContrastName"))
  checkmate::assert_data_frame(tfce, nrows = 1)
  checkmate::assert_data_frame(tfce_pop, nrows = 1)

  study <- get_pairs(fs::path(storage_dir, fs::path_file(tfce$tstat)), tfce$n_sub) |>
    mask() |>
    mask_gray() |>
    dplyr::mutate(study = cope / sigma * correct_d(tfce$n_sub)) |>
    dplyr::select(x, y, z, study)

  test <- get_pairs(fs::path(storage_dir, fs::path_file(tfce_pop$tstat)), tfce_pop$n_sub) |>
    mask() |>
    mask_gray() |>
    dplyr::mutate(test = cope / sigma * correct_d(tfce$n_sub)) |>
    dplyr::select(x, y, z, test)

  dplyr::left_join(study, test, by = c("x", "y", "z")) |>
    dplyr::left_join(at, by = c("x", "y", "z")) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::group_by(label, `Label Name`, `Network Name`, `Full component name`, hemi) |>
    dplyr::summarise(rho = cor(study, test, method = .env$method), .groups = "drop") |>
    dplyr::bind_cols(
      tfce |>
        dplyr::select(Task, CopeNumber, ContrastName, iter, n_sub)
    ) |>
    dplyr::mutate(method = .env$method)
}


get_pairs <- function(tstat, n_sub) {
  cope <- to_tbl(stringr::str_replace(tstat, "_tstat1", "_glm_cope"), measure = "cope")
  varcope <- to_tbl(stringr::str_replace(tstat, "_tstat1", "_glm_varcope"), measure = "varcope") |>
    dplyr::mutate(sigma = sqrt(varcope) * sqrt(n_sub)) |>
    dplyr::select(-varcope)

  dplyr::left_join(cope, varcope, by = c("x", "y", "z")) |>
    mask() 
}


get_pop_d <- function(tfce_pop) {
  tfce_pop |>
    dplyr::select(tstat, Task, CopeNumber, ContrastName, n_sub) |>
    dplyr::mutate(data = purrr::map2(tstat, n_sub, get_pairs)) |>
    dplyr::select(-tstat) |>
    tidyr::unnest(data) |>
    mask() |>
    dplyr::mutate(
      d = cut(
        abs(cope / sigma * correct_d(n_sub)),
        right = FALSE,
        ordered = TRUE,
        breaks = c(0, .2, .5, .8, Inf),
        labels = c("0", "small", "medium", "large")
      )
    )
}
