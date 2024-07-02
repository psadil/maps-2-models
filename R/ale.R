do_z <- function(cope_files) {
  z_stats <- cope_files %>%
    dplyr::group_nest(
      study, n_sub, iter, n_study
    ) |>
    dplyr::mutate(
      z_stat = purrr::map(data, ~ calc_z(.x$copes))
    ) |>
    dplyr::mutate(
      z = here::here(
        "data-raw",
        "niis",
        glue::glue("nstudy-{n_study}_nsub-{n_sub}_study-{study}_iter-{iter}_z.nii.gz")
      )
    )

  purrr::walk2(
    z_stats$z_stat,
    z_stats$z,
    neurobase::writenii
  )

  z_stats |>
    dplyr::select(study, n_sub, iter, n_study, z) |>
    dplyr::mutate(z = fs::path_rel(z, here::here()))
}


do_t <- function(cope_files) {
  stopifnot({
    dplyr::n_distinct(cope_files$n_sub) == 1
    dplyr::n_distinct(cope_files$iter) == 1
    dplyr::n_distinct(cope_files$n_study) == 1
    dplyr::n_distinct(cope_files$study) == 1
  })

  iter <- unique(cope_files$iter)
  n_sub <- unique(cope_files$n_sub)
  n_study <- unique(cope_files$n_study)
  study <- unique(cope_files$study)

  t_vcope <- calc_t_vcope(cope_files$copes)
  varcope_file <- neurobase::writenii(
    t_vcope$varcope,
    here::here(
      "data-raw", "niis", glue::glue("nstudy-{n_study}_nsub-{n_sub}_study-{study}_iter-{iter}_varcope.nii.gz")
    )
  ) %>%
    fs::path_rel(here::here())
  t_file <- neurobase::writenii(
    t_vcope$t,
    here::here(
      "data-raw", "niis", glue::glue("nstudy-{n_study}_nsub-{n_sub}_study-{study}_iter-{iter}_t.nii.gz")
    )
  ) %>%
    fs::path_rel(here::here())

  cope_files %>%
    dplyr::distinct(n_sub, iter, study, n_study) %>%
    dplyr::mutate(
      t = t_file,
      varcope = varcope_file
    )
}

calc_t_vcope <- function(cope_files) {
  copes <- RNifti::readNifti(cope_files) |> simplify2array()
  stopifnot(length(dim(copes)) == 4)
  n <- dim(copes)[4]
  means <- apply(copes, 1:3, mean)
  sds <- apply(copes, 1:3, sd)
  t_stat <- sqrt(n) * means / sds
  t_stat[is.na(t_stat)] <- 0
  sds[is.na(sds)] <- 0
  t_img <- neurobase::niftiarr(neurobase::readnii(cope_files[[1]]), t_stat)
  varcope_img <- neurobase::niftiarr(neurobase::readnii(cope_files[[1]]), (sds / sqrt(n))^2)
  list("t" = t_img, "varcope" = varcope_img)
}


calc_z <- function(cope_files) {
  copes <- RNifti::readNifti(cope_files) |> simplify2array()
  stopifnot(length(dim(copes)) == 4)
  n <- dim(copes)[4]
  means <- apply(copes, 1:3, mean)
  sds <- apply(copes, 1:3, sd)
  t_stat <- sqrt(n) * means / sds
  z_stat <- apply(
    t_stat, 1:3,
    FUN = function(x) {
      pt(x, n - 1, log.p = TRUE, lower.tail = FALSE) |>
        qnorm(log.p = TRUE, lower.tail = FALSE)
    }
  )
  z_stat[is.na(z_stat)] <- 0

  neurobase::niftiarr(neurobase::readnii(cope_files[[1]]), z_stat)
}


calc_dice <- function(nii1, nii2, lower = 0.0001, na.rm = TRUE) {
  abs1 <- abs(nii1) > lower
  abs2 <- abs(nii2) > lower
  2 * sum(abs1 * abs2, na.rm = na.rm) / (sum(abs1, na.rm = na.rm) + sum(abs2, na.rm = na.rm))
}
