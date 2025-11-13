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


cimg_to_tbl <- function(.i) {
  as.raster(.i) |>
    as.matrix() |>
    tibble::as_tibble() |>
    dplyr::mutate(y = dplyr::row_number()) |>
    tidyr::pivot_longer(-y, names_to = "x") |>
    dplyr::mutate(
      x = stringr::str_extract(x, "[[:digit:]]+") |> as.integer(),
      y = max(y) - y + 1,
      rgb = purrr::map(value, col2rgb),
      r = purrr::map_dbl(rgb, purrr::pluck, 1),
      g = purrr::map_dbl(rgb, purrr::pluck, 2),
      b = purrr::map_dbl(rgb, purrr::pluck, 3)
    )
}


plot_msmall <- function(nsub, iter, task = "EMOTION") {
  palmdir <- Sys.getenv("PALMDIR")
  dtseries <- fs::dir_ls(
    palmdir,
    glob = glue::glue(
      "*nsub-{nsub}_iter-{iter}_task-{task}*_type-MSMALL.dtseries.nii"
    )
  )
  l <- fs::dir_ls(
    palmdir,
    glob = glue::glue(
      "*iter-0_task-{task}*_type-MSMALL_L_midthickness.32k_fs_LR.gii"
    )
  )

  dtseries_pop <- fs::dir_ls(
    palmdir,
    glob = glue::glue(
      "*ter-0_task-{task}*_type-MSMALL.dtseries.nii"
    )
  )

  xii <- ciftiTools::read_cifti(cifti_fname = dtseries_pop[1])
  xii_mean <- ciftiTools::apply_xifti(xii, 1, mean)
  xii_sd <- ciftiTools::apply_xifti(xii, 1, sd)
  xii_d <- ciftiTools::transform_xifti(xii_mean, `/`, xii_sd)

  a_max <- max(abs(xii_d))

  xii <- ciftiTools::read_cifti(
    cifti_fname = dtseries[1],
    surfL_fname = l,
    surfR_fname = stringr::str_replace(l, "_L_", "_R_")
  )
  xii_mean <- ciftiTools::apply_xifti(xii, 1, mean)
  xii_sd <- ciftiTools::apply_xifti(xii, 1, sd)
  xii_d <- ciftiTools::transform_xifti(xii_mean, `/`, xii_sd)

  fname <- tempfile(fileext = ".png")
  ciftiTools::view_xifti_surface(
    xii_d,
    zlim = c(-a_max, a_max),
    color_mode = "diverging",
    fname = fname
  )
  out <- imager::load.image(fname)
  unlink(fname)
  cimg_to_tbl(out) |>
    dplyr::mutate(n_sub = nsub)
}

plot_msmall_statmap <- function(task) {
  palmdir <- Sys.getenv("PALMDIR")
  dtseries <- fs::dir_ls(
    palmdir,
    glob = glue::glue(
      "*{task}*type-MSMALL.dtseries.nii"
    )
  )

  nsubs <- stringr::str_extract(dtseries, "(?<=sub-)[[:digit:]]+") |>
    as.integer() |>
    unique() |>
    sort()

  figs <- purrr::map2(
    nsubs,
    c(94, 94, 94, 94, 94, 0),
    ~ plot_msmall(.x, .y, task = task),
  )

  dplyr::bind_rows(figs) |>
    dplyr::mutate(
      n_sub = factor(
        n_sub,
        unique(n_sub),
        labels = c(
          "N Sub: 20",
          "N Sub: 40",
          "N Sub: 60",
          "N Sub: 80",
          "N Sub: 100",
          "N Sub: Gold"
        )
      )
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
    ggplot2::facet_wrap(~n_sub, ncol = 2) +
    terrainr::geom_spatial_rgb(ggplot2::aes(r = r, g = g, b = b)) +
    ggplot2::coord_fixed() +
    ggplot2::theme_void() +
    ggplot2::ggtitle(stringr::str_to_lower(glue::glue("msmall, {task}")))
}
