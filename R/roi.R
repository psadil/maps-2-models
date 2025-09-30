test_roi_wrapper <- function(roi_avg, hcp_samples, .fwer = 0.05) {
  hcp_samples |>
    dplyr::group_nest(tar_group) |>
    dplyr::mutate(rois = purrr::map(data, ~ test_roi(roi_avg, .x))) |>
    dplyr::select(-tar_group, -data) |>
    tidyr::unnest(rois)
}

test_roi <- function(rois, hcp_samples, .fwer = 0.05) {
  rois |>
    dplyr::inner_join(
      hcp_samples,
      by = dplyr::join_by(type, Task, CopeNumber, sub),
      relationship = "many-to-many"
    ) |>
    dplyr::group_nest(n_parcels, iter, n_sub, Task, type, label) |>
    dplyr::mutate(
      fit = purrr::map(
        data,
        ~ stats::t.test(value ~ 1, data = .x) |>
          broom::tidy()
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit) |>
    dplyr::mutate(
      p.adjusted = stats::p.adjust(p.value, "holm"),
      .by = c(n_parcels, iter, n_sub, Task, type)
    ) |>
    dplyr::mutate(active = p.adjusted < .fwer)
}

test_roi_pop <- function(rois, .fwer = 0.05) {
  rois |>
    dplyr::group_nest(label, Task, CopeNumber, n_parcels, type) |>
    dplyr::mutate(
      fit = purrr::map(
        data,
        ~ stats::t.test(value ~ 1, data = .x) |>
          broom::tidy()
      ),
      n_sub = purrr::map_dbl(data, nrow)
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit) |>
    dplyr::mutate(
      p.adjusted = stats::p.adjust(p.value, "holm"),
      .by = c(Task, CopeNumber, n_parcels, type)
    ) |>
    dplyr::mutate(
      active = p.adjusted < .fwer,
      iter = 0
    )
}


roi_from_cifti <- function(file, n_parcels) {
  parc <- load_parc(n_parcels = n_parcels)
  parc <- ciftiTools::parc_add_subcortex(parc)
  xii <- ciftiTools::read_xifti(file)
  ciftiTools::apply_parc(xii, parc, FUN = mean, na.rm = TRUE) |>
    tibble::as_tibble(rownames = "label") |>
    dplyr::rename(value = V1) |>
    dplyr::filter(!label == "???") |>
    na.omit()
}


roi_from_nifti <- function(file, n_parcels) {
  at <- make_atlas_full(n_parcels = n_parcels)

  to_tbl(file) |>
    dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::summarise(value = mean(value), .by = c(label))
}

avg_roi <- function(test, n_parcels) {
  msmall <- test |>
    dplyr::select(-VOL, -SURFACE) |>
    na.omit() |>
    dplyr::mutate(MSMALL = purrr::map(MSMALL, roi_from_cifti, n_parcels)) |>
    tidyr::unnest(MSMALL)

  vol <- test |>
    dplyr::select(-MSMALL, -SURFACE) |>
    na.omit()
  if (nrow(vol) > 0) {
    vol <- vol |>
      dplyr::mutate(VOL = purrr::map(VOL, roi_from_nifti, n_parcels)) |>
      tidyr::unnest(VOL)
  } else {
    vol <- dplyr::select(vol, -VOL)
  }

  surface <- test |>
    dplyr::select(-VOL, -MSMALL) |>
    na.omit() |>
    dplyr::mutate(SURFACE = purrr::map(SURFACE, roi_from_cifti, n_parcels)) |>
    tidyr::unnest(SURFACE)

  dplyr::bind_rows(
    list(MSMALL = msmall, VOL = vol, SURFACE = surface),
    .id = "type"
  ) |>
    dplyr::select(
      type,
      Task,
      CopeNumber,
      ContrastName,
      sub,
      label,
      value
    ) |>
    dplyr::mutate(
      n_parcels = .env$n_parcels
    )
}


make_data_roi_study_to_gold <- function(
  gold_tested,
  rois_tested,
  method = "spearman"
) {
  gold_most <- gold_tested |>
    dplyr::mutate(
      r = dplyr::row_number(dplyr::desc(abs(estimate))),
      .by = c(Task, n_parcels, type)
    ) |>
    dplyr::filter(r < 11) |>
    dplyr::mutate(estimate = statistic / sqrt(parameter)) |>
    dplyr::select(Task, n_parcels, label, type, estimate)

  rois_tested |>
    dplyr::right_join(
      gold_most,
      by = dplyr::join_by(Task, n_parcels, label, type)
    ) |>
    dplyr::mutate(estimate.x = statistic / sqrt(parameter)) |>
    dplyr::group_nest(Task, n_parcels, iter, n_sub, type) |>
    dplyr::mutate(
      fit = purrr::map(
        data,
        ~ cor.test(.x$estimate.x, .x$estimate.y, method = .env$method) |>
          broom::tidy() |>
          dplyr::select(-method, -alternative)
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit) |>
    dplyr::summarise(
      lower = quantile(estimate, 0.025),
      avg = mean(estimate),
      upper = quantile(estimate, 0.975),
      .by = c(Task, n_sub, type, n_parcels)
    )
}

make_data_roi_study_to_gold2 <- function(
  gold_tested,
  rois_tested,
  method = "spearman"
) {
  gold_most <- gold_tested |>
    dplyr::mutate(
      r = dplyr::row_number(dplyr::desc(abs(estimate))),
      .by = c(Task, n_parcels, type)
    ) |>
    dplyr::filter(r < 11) |>
    dplyr::mutate(d = statistic / sqrt(parameter)) |>
    dplyr::select(Task, n_parcels, label, type, d)

  rois_tested |>
    dplyr::semi_join(
      gold_most,
      by = dplyr::join_by(Task, n_parcels, label, type)
    ) |>
    dplyr::summarise(
      avg = mean(active),
      lower = qbeta(
        0.025,
        1 / 2 + sum(active),
        dplyr::n() - sum(active) + 1 / 2
      ),
      upper = qbeta(
        0.975,
        1 / 2 + sum(active),
        dplyr::n() - sum(active) + 1 / 2
      ),
      .by = c(Task, n_sub, type, label, n_parcels)
    ) |>
    dplyr::left_join(
      gold_most,
      by = dplyr::join_by(Task, type, label, n_parcels)
    )
}


get_icc_binary <- function(X, ind) {
  n <- X$n[ind]
  y <- X$y[ind]
  X <- data.frame(n = n, y = y)
  tt <- tryCatch(
    .estimate <- aod::iccbin(n = n, y = y, data = X, method = "B")@rho[1],
    error = function(e) e,
    warning = function(w) w
  )
  if (is(tt, "warning")) {
    rho <- NA
  } else {
    rho <- .estimate
  }
  rho
}

do_iccbins <- function(.data, n_boot = 100, n_workers = 8) {
  if (all(.data$y == 0) | all(.data$y == 1)) {
    .lower <- NA
    .upper <- NA
    .estimate <- NA
  } else {
    .estimate <- aod::iccbin(n = n, y = y, data = .data, method = "B")@rho[1]
    boots <- boot::boot(
      data = .data,
      statistic = get_icc_binary,
      R = n_boot,
      parallel = "multicore",
      ncpus = n_workers
    )
    ci <- boot::boot.ci(boots, type = "perc")
    .lower <- ci$percent[[4]]
    .upper <- ci$percent[[5]]
  }
  tibble::tibble(.estimate = .estimate, .lower = .lower, .upper = .upper)
}


do_iccs <- function(.data) {
  i <- .data |>
    tidyr::pivot_wider(names_from = iter, values_from = statistic) |>
    dplyr::select(-label) |>
    irr::icc(model = "t", type = "c")
  tibble::tibble(.estimate = i$value, .lower = i$lbound, .upper = i$ubound)
}

make_data_roi_study_to_study2 <- function(
  rois_tested,
  n_boot = 100,
  n_workers = 8
) {
  rois_tested |>
    dplyr::summarise(
      n = dplyr::n(),
      y = sum(active),
      .by = c(n_parcels, n_sub, Task, type, label)
    ) |>
    dplyr::group_nest(Task, n_parcels, n_sub, type) |>
    dplyr::mutate(
      fit = purrr::map(
        data,
        do_iccbins,
        n_workers = .env$n_workers,
        n_boot = .env$n_boot
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit)
}


make_data_roi_study_to_study <- function(
  rois_tested,
  n_boot = 100,
  n_workers = 8,
  type = "bca"
) {
  rois_tested |>
    dplyr::select(Task, n_parcels, n_sub, iter, statistic, label, type) |>
    dplyr::group_nest(Task, n_parcels, n_sub, type) |>
    dplyr::mutate(
      fit = purrr::map(
        data,
        do_iccs
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(fit)
}

make_data_roi_sub_to_sub <- function(rois_pop) {
  rois_pop |>
    dplyr::mutate(estimate = estimate / sqrt(parameter)) |>
    dplyr::select(Task, n_parcels, n_sub, type, label, estimate) |>
    dplyr::group_nest(Task, n_parcels, type, n_sub) |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~ .x |>
          dplyr::select(label, Z, sub) |>
          tidyr::pivot_wider(names_from = label, values_from = Z) |>
          dplyr::arrange(sub)
      ),
      rho = purrr::map(
        data,
        ~ dplyr::select(.x, -sub) |>
          as.matrix() |>
          sim.rho() |>
          corrr::as_cordf() |>
          corrr::shave() |>
          corrr::stretch(na.rm = TRUE) |>
          dplyr::rename(rho = r)
      )
    ) |>
    dplyr::mutate(n_parcels = factor(n_parcels)) |>
    dplyr::select(-data) |>
    tidyr::unnest(rho)
}
