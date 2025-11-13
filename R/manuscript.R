make_roi <- function(data_roi_study_to_gold, data_roi_study_to_study) {
  a <- data_roi_study_to_gold |>
    dplyr::filter(n_parcels == 400) |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task),
      type = stringr::str_to_lower(type)
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = avg, color = type)
    ) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.3
    ) +
    ggplot2::facet_wrap(~Task, nrow = 2, scales = "free_x") +
    ggplot2::scale_y_continuous(
      "Rank Correlation with\nGold Standard\n(Most Active ROI)",
      breaks = c(0, 0.5, 1),
      labels = c(0, 0.5, 1)
    ) +
    ggplot2::scale_x_continuous(
      "N Sub",
      transform = "log10"
    ) +
    ggplot2::scale_color_viridis_d(name = NULL, option = "turbo") +
    ggplot2::guides(colour = ggplot2::guide_legend(position = "inside"))

  b <- data_roi_study_to_study |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task),
      type = stringr::str_to_lower(type)
    ) |>
    dplyr::filter(n_parcels == 400) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = .estimate, color = type)) +
    ggplot2::facet_wrap(~Task, scales = "free_x", nrow = 2) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .lower, ymax = .upper),
      alpha = 0.3
    ) +
    ggplot2::scale_y_continuous(
      "ICC(C,1) Across\nBootstrap Samples",
      breaks = c(0, 0.5, 1),
      labels = c(0, 0.5, 1)
    ) +
    ggplot2::scale_x_continuous(
      "N Sub",
      transform = "log10"
    ) +
    ggplot2::scale_color_viridis_d(name = NULL, option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
    )

  a +
    b +
    patchwork::plot_layout(ncol = 1) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        legend.margin = ggplot2::margin(0, 0, 0, 0), # turned off for alignment
        legend.justification.top = "left",
        legend.justification.left = "bottom",
        legend.justification.bottom = "right",
        legend.justification.inside = c(1, 0),
        legend.location = "plot"
      )
}

make_roi2 <- function(data_roi_study_to_gold, data_roi_study_to_study) {
  a <- data_roi_study_to_gold |>
    dplyr::filter(n_parcels == 400) |>
    dplyr::mutate(
      d = abs(d),
      Task = stringr::str_to_lower(Task),
      type = stringr::str_to_lower(type)
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = avg, color = d, group = label)
    ) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper)) +
    ggplot2::facet_wrap(~ Task + type, scales = "free_x", nrow = 3) +
    ggplot2::scale_y_continuous(
      "Proportion Bootstrap Samples Active (Most Active ROI)",
      breaks = c(0, 0.5, 1),
      labels = c(0, 0.5, 1)
    ) +
    ggplot2::scale_x_continuous(
      "N Sub",
      transform = "log10"
    ) +
    ggplot2::scale_color_viridis_c(name = "Cohen's d", option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_colorbar(position = "inside"),
    )

  b <- data_roi_study_to_study |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task),
      type = stringr::str_to_lower(type)
    ) |>
    dplyr::filter(n_parcels == 400) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = .estimate, color = type)) +
    ggplot2::facet_wrap(~Task, scales = "free_x", nrow = 2) +
    ggplot2::geom_point(alpha = 0.3) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .lower, ymax = .upper),
      alpha = 0.3
    ) +
    ggplot2::scale_y_continuous(
      "ICC(1) Across Bootstrap Samples",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c(0, 0.5, 1)
    ) +
    ggplot2::scale_x_continuous(
      "N Sub",
      transform = "log10"
    ) +
    ggplot2::scale_color_viridis_d(name = NULL, option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
    )

  a +
    b +
    patchwork::plot_layout(ncol = 1, heights = c(2, 1)) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        legend.margin = ggplot2::margin(0, 0, 0, 0), # turned off for alignment
        legend.justification.top = "left",
        legend.justification.left = "bottom",
        legend.justification.bottom = "right",
        legend.justification.inside = c(1, 0),
        legend.location = "plot"
      )
}

get_max <- function(q) {
  x <- qs::qread(q)
  to_tbl0(x$Z, measure = "Z") |> mask()
}


make_prop_active_most_active_roi_ptfce_null <- function(
  at_list,
  active_null,
  iter,
  gold_tested,
  active_threshold = 0.02
) {
  n_sims <- dplyr::n_distinct(iter)

  # regions with at least one voxel active
  out_null <- active_null |>
    dplyr::collect() |>
    dplyr::left_join(
      at_list,
      by = dplyr::join_by(x, y, z),
      relationship = "many-to-many"
    ) |>
    dplyr::distinct(
      n_sub,
      iter,
      label,
      `Label Name`,
      `Full component name`,
      n_parcels
    ) |>
    dplyr::count(
      n_sub,
      label,
      `Label Name`,
      `Full component name`,
      n_parcels
    ) |>
    dplyr::mutate(prop = n / n_sims)

  gold_null <- gold_tested |>
    dplyr::filter(Task == "WM") |>
    dplyr::mutate(
      d = statistic / sqrt(n_sub),
      active = abs(d) > active_threshold
    ) |>
    dplyr::mutate(
      r = dplyr::row_number(dplyr::desc(abs(estimate))),
      .by = c(Task, n_parcels)
    ) |>
    dplyr::filter(r < 11) |>
    dplyr::select(Task, n_parcels, label, d) |>
    dplyr::mutate(
      l = label |>
        factor() |>
        as.numeric() |>
        factor()
    )

  out_null |>
    dplyr::semi_join(dplyr::distinct(gold_null, l, label, n_parcels)) |>
    dplyr::right_join(
      dplyr::distinct(gold_null, label, n_parcels, l) |>
        tidyr::crossing(
          dplyr::distinct(out_null, n_sub)
        )
    ) |>
    dplyr::filter(n_parcels %in% c(200, 400, 600, 800, 1000)) |>
    dplyr::mutate(
      n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
      n_parcels = forcats::fct_relabel(
        n_parcels,
        .fun = ~ glue::glue("N Parcels: {.x}")
      ),
      propr = dplyr::if_else(is.na(prop), 0, prop)
    ) |>
    dplyr::mutate(prop = dplyr::if_else(is.na(prop), 0, prop)) |>
    ggplot2::ggplot(aes(x = n_sub, y = prop, group = l)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = qbeta(0.05 / 2, 5, 100 - 5 + 1),
        ymax = qbeta(1 - 0.05 / 2, 5, 100 - 5)
      ),
      alpha = 0.01,
      fill = "lightblue",
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::geom_point(show.legend = FALSE, alpha = 0.2) +
    ggplot2::geom_line(show.legend = FALSE, alpha = 0.2) +
    ggplot2::facet_wrap(~n_parcels) +
    ggplot2::scale_y_continuous(
      "Proportion Simulations w/\nActivity in Most Active ROI",
      limits = c(0, 0.12)
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::theme_minimal(base_size = 12)
}

.make_prop_active_most_active_roi <- function(
  d,
  breaks = c(40, 80),
  transform = "identity"
) {
  d |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = n_sub,
        y = avg,
        group = label,
        color = abs(d)
      )
    ) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::geom_line(alpha = 0.2) +
    ggplot2::facet_grid(n_parcels ~ Task) +
    ggplot2::scale_y_continuous(
      "Proportion Simulations w/\nActivity in Most Active ROI",
      limits = c(0, 1),
      labels = c(0, 0.5, 1),
      breaks = c(0, 0.5, 1)
    ) +
    ggplot2::scale_x_continuous(
      "N Sub",
      breaks = breaks,
      labels = breaks,
      transform = transform
    ) +
    ggplot2::scale_color_viridis_c(
      "Abs. Effect Size",
      option = "turbo",
      limits = c(0, 3),
      n.breaks = 3
    )
}

make_prop_active_most_active_roi <- function(data_roi_study_to_gold2) {
  d <- data_roi_study_to_gold2 |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task),
    )
  vol <- dplyr::filter(d, type == "VOL") |>
    .make_prop_active_most_active_roi() +
    ggplot2::ggtitle("VOL")

  ukb <- dplyr::filter(d, type == "UKB") |>
    .make_prop_active_most_active_roi(
      transform = "log10",
      breaks = c(40, 80, 100, 1000, 10000)
    ) +
    ggplot2::ggtitle("UKB")

  msmall <- dplyr::filter(d, type == "MSMALL") |>
    .make_prop_active_most_active_roi() +
    ggplot2::ggtitle("MSMALL")

  surface <- dplyr::filter(d, type == "SURFACE") |>
    .make_prop_active_most_active_roi() +
    ggplot2::ggtitle("SURFACE")

  vol +
    surface +
    msmall +
    ukb +
    patchwork::plot_layout(ncol = 1, guides = "collect") &
    ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(legend.position = "bottom")
}

.make_1_peaks <- function(data_peak_study_to_gold, type) {
  data_peak_study_to_gold |>
    dplyr::filter(type == .env$type) |>
    dplyr::mutate(Task = stringr::str_to_lower(Task)) |>
    ggplot2::ggplot(
      ggplot2::aes(
        x = within,
        group = label,
        y = n_simulations,
        color = d,
      )
    ) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::geom_line(alpha = 0.2) +
    ggplot2::facet_grid(Task ~ n_sub) +
    ggplot2::scale_y_continuous(
      "Proportion Simulations w/ Peak in Radius",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c(0, 0.5, 1)
    ) +
    ggplot2::scale_x_continuous(
      "Radius (mm)",
      limits = c(0, 20),
      breaks = c(0, 10, 20),
      labels = c(0, 10, 20)
    ) +
    ggplot2::scale_color_viridis_c(
      "Cohen's d",
      option = "turbo",
      limits = c(0, 4)
    ) +
    ggplot2::ggtitle(type)
}

make_peaks_validity <- function(
  data_peak_study_to_gold,
  threshold = "unthresholded"
) {
  if (threshold == "unthresholded") {
    thresholded <- dplyr::filter(data_peak_study_to_gold, threshold == 0)
  } else {
    thresholded <- dplyr::filter(data_peak_study_to_gold, threshold > 0)
  }

  vol <- .make_1_peaks(thresholded, "VOL")
  ukb <- .make_1_peaks(thresholded, "UKB")
  msm <- .make_1_peaks(thresholded, "MSMALL")
  surf <- .make_1_peaks(thresholded, "SURFACE")

  vol +
    ukb +
    msm +
    surf +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 8)
}

.make_1_peak_reliability <- function(.d, type, nrow) {
  if (type %in% c("VOL", "UKB")) {
    upper <- 100
    labels <- breaks <- c(0, 50)
  } else {
    upper <- 300
    labels <- breaks <- c(0, 50, 100, 150, 200, 250, 300)
    labels <- c(0, 50, "", "", "", "", "")
  }
  .d |>
    dplyr::filter(type == .env$type) |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task)
    ) |>
    ggplot2::ggplot(ggplot2::aes(y = n_sub, x = d)) +
    ggplot2::facet_wrap(~rank, nrow = nrow, labeller = ggplot2::label_both) +
    ggdist::stat_dots(
      ggplot2::aes(color = Task, fill = Task),
      quantiles = 100,
      position = "dodge",
      show.legend = TRUE
    ) +
    ggplot2::scale_fill_viridis_d(
      option = "turbo",
      drop = FALSE,
    ) +
    ggplot2::scale_color_viridis_d(
      option = "turbo",
      drop = FALSE,
    ) +
    ggplot2::ylab("N Sub") +
    ggplot2::scale_x_continuous(
      "Distance Between\nAssociated Peaks\n(Study-Study)",
      limits = c(0, upper),
      breaks = breaks,
      labels = labels,
      transform = "pseudo_log"
    ) +
    ggplot2::ggtitle(type)
}

make_peaks_reliability <- function(
  data_peak_study_to_study,
  threshold = "unthresholded",
  nrow_subfig = 2,
  base_size = 8
) {
  data_peak_study_to_study <- data_peak_study_to_study |>
    dplyr::mutate(Task = factor(Task))
  if (threshold == "unthresholded") {
    thresholded <- dplyr::filter(data_peak_study_to_study, threshold == 0)
  } else {
    thresholded <- dplyr::filter(data_peak_study_to_study, threshold > 0)
  }

  .make_1_peak_reliability(thresholded, "VOL", nrow = nrow_subfig) +
    .make_1_peak_reliability(thresholded, "UKB", nrow = nrow_subfig) +
    .make_1_peak_reliability(thresholded, "MSMALL", nrow = nrow_subfig) +
    .make_1_peak_reliability(thresholded, "SURFACE", nrow = nrow_subfig) +
    patchwork::plot_layout(nrow = 2, guides = "collect") +
    patchwork::plot_annotation(
      tag_levels = "a",
      tag_suffix = ")"
    ) &
    ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(legend.position = "bottom")
}


make_peak_bysize <- function(study_to_gold_distances, glm_pop2) {
  peaks <- study_to_gold_distances |>
    dplyr::filter(type == "VOL") |>
    dplyr::summarise(
      avg_d = mean(d),
      .by = c(type, Task, n_sub, x, y, z, threshold)
    )

  mapping <- to_tbl(MNITemplate::getMNIPath("Brain", res = "2mm")) |>
    mask() |>
    mask_gray() |>
    mask_atlas() |>
    dplyr::mutate(index = 1:dplyr::n()) |>
    dplyr::select(-value)

  d <- glm_pop2 |>
    dplyr::filter(type == "VOL") |>
    dplyr::mutate(
      glm = stringr::str_remove(glm, "/dcl01/smart/data/psadil/meta/"),
      data = purrr::map(
        glm,
        ~ duckplyr::read_parquet_duckdb(.x) |> dplyr::collect()
      )
    ) |>
    dplyr::select(-glm) |>
    tidyr::unnest(data) |>
    dplyr::mutate(g = abs(pe) / sigma * correct_d(n_sub)) |>
    dplyr::select(-z, -pe, -sigma, -n_sub) |>
    dplyr::left_join(mapping, by = dplyr::join_by(index)) |>
    dplyr::full_join(peaks, by = dplyr::join_by(type, Task, x, y, z)) |>
    na.omit() |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task),
      `N Sub` = glue::glue("N Sub: {n_sub}"),
      `N Sub` = factor(
        `N Sub`,
        levels = c(
          "N Sub: 20",
          "N Sub: 40",
          "N Sub: 60",
          "N Sub: 80",
          "N Sub: 100"
        )
      )
    )

  a <- d |>
    dplyr::filter(threshold == 0) |>
    ggplot2::ggplot(ggplot2::aes(x = g, y = avg_d)) +
    ggplot2::facet_grid(`N Sub` ~ Task) +
    scattermore::geom_scattermore(pointsize = 5, alpha = 0.5) +
    ggplot2::ylab("avg dist(Gold Std., Study) (mm)") +
    ggplot2::scale_x_continuous(
      "Gold Standard Peak Cohen's d",
      breaks = c(0, 1),
      labels = c(0, 1)
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::ggtitle("Unthresholded")

  b <- d |>
    dplyr::filter(threshold > 0) |>
    ggplot2::ggplot(ggplot2::aes(x = g, y = avg_d)) +
    ggplot2::facet_grid(`N Sub` ~ Task) +
    scattermore::geom_scattermore(pointsize = 5, alpha = 0.5) +
    ggplot2::ylab("avg dist(Gold Std., Study) (mm)") +
    ggplot2::scale_x_continuous(
      "Gold Standard Peak Cohen's d",
      breaks = c(0, 1),
      labels = c(0, 1)
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::ggtitle("Thresholded")

  a +
    b +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") +
    patchwork::plot_layout(nrow = 2) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}

.make_fig_bynetwork <- function(.d) {
  .d |>
    ggplot2::ggplot(ggplot2::aes(y = `Network Name`, x = avg_d, color = g)) +
    ggplot2::geom_boxplot(outliers = FALSE) +
    scattermore::geom_scattermore(
      pointsize = 5,
      position = ggplot2::position_jitter(width = 0),
      alpha = 0.5
    ) +
    ggplot2::facet_grid(`N Sub` ~ Task) +
    ggplot2::scale_color_viridis_c(
      option = "turbo",
      guide = ggplot2::guide_colorbar("Cohen's d"),
      limits = c(0, NA)
    ) +
    ggplot2::ylab("Network") +
    ggplot2::scale_x_continuous(
      "avg dist(Gold Standard Peak, Study Peak) (mm)"
    ) +
    ggplot2::theme_minimal(base_size = 8)
}


make_peak_bynetwork <- function(study_to_gold_distances, at, glm_pop2) {
  peaks <- study_to_gold_distances |>
    dplyr::filter(type == "VOL") |>
    dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
    dplyr::mutate(
      `Network Name` = dplyr::if_else(
        is.na(`Network Name`) & !is.na(label),
        "subcortical",
        `Network Name`
      )
    ) |>
    dplyr::filter(!is.na(`Network Name`)) |>
    dplyr::summarise(
      avg_d = mean(d),
      .by = c(type, Task, n_sub, `Network Name`, threshold, iter)
    )

  mapping <- to_tbl(MNITemplate::getMNIPath("Brain", res = "2mm")) |>
    mask() |>
    mask_gray() |>
    mask_atlas() |>
    dplyr::mutate(index = 1:dplyr::n()) |>
    dplyr::select(-value)

  eff_size <- glm_pop2 |>
    dplyr::filter(type == "VOL") |>
    dplyr::mutate(
      glm = stringr::str_remove(glm, "/dcl01/smart/data/psadil/meta/"),
      data = purrr::map(
        glm,
        ~ duckplyr::read_parquet_duckdb(.x) |> dplyr::collect()
      )
    ) |>
    dplyr::select(-glm) |>
    tidyr::unnest(data) |>
    dplyr::mutate(g = abs(pe) / sigma * correct_d(n_sub)) |>
    dplyr::select(-z, -pe, -sigma, -n_sub) |>
    dplyr::left_join(mapping, by = dplyr::join_by(index)) |>
    dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
    dplyr::mutate(
      `Network Name` = dplyr::if_else(
        is.na(`Network Name`) & !is.na(label),
        "subcortical",
        `Network Name`
      )
    ) |>
    dplyr::filter(!is.na(`Network Name`)) |>
    dplyr::summarise(
      g = mean(g),
      .by = c(type, Task, `Network Name`)
    )

  d <- peaks |>
    dplyr::left_join(
      eff_size,
      by = dplyr::join_by(type, Task, `Network Name`)
    ) |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task),
      `N Sub` = glue::glue("N Sub: {n_sub}"),
      `N Sub` = factor(
        `N Sub`,
        levels = c(
          "N Sub: 20",
          "N Sub: 40",
          "N Sub: 60",
          "N Sub: 80",
          "N Sub: 100"
        )
      )
    )

  a <- d |>
    dplyr::filter(threshold == 0) |>
    .make_fig_bynetwork() +
    ggplot2::ggtitle("Unthresholded")

  b <- d |>
    dplyr::filter(threshold > 0) |>
    .make_fig_bynetwork() +
    ggplot2::ggtitle("Thresholded")

  a +
    b +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") +
    patchwork::plot_layout(guides = "collect", nrow = 2) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}

make_topo <- function(data_topo_gold_to_study, data_topo_study_to_study) {
  a <- data_topo_gold_to_study |>
    dplyr::mutate(
      `N Sub` = factor(n_sub),
      Task = stringr::str_to_lower(Task),
      type = stringr::str_to_lower(type)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = rho, y = `N Sub`, color = type)) +
    ggplot2::facet_wrap(~Task, scales = "free_y", nrow = 2) +
    ggplot2::geom_boxplot(outliers = FALSE) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::xlab("Rank Correlation (Gold to Study)") +
    ggplot2::guides(colour = ggplot2::guide_legend(position = "inside"))

  b <- data_topo_study_to_study |>
    dplyr::filter(method == "consistency") |>
    dplyr::mutate(
      Task = stringr::str_to_lower(Task),
      type = stringr::str_to_lower(type)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = estimate, color = type)) +
    ggplot2::facet_wrap(~Task, scales = "free_x", nrow = 2) +
    ggplot2::geom_line(alpha = 0.3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.3
    ) +
    ggplot2::scale_x_continuous("N Sub", transform = "log10") +
    ggplot2::scale_y_continuous(
      "ICC(C,1) Across Bootstrap Samples",
      limits = c(0, 1)
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::guides(colour = ggplot2::guide_legend(position = "inside"))

  a +
    b +
    patchwork::plot_layout(nrow = 2) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        legend.margin = ggplot2::margin(0, 0, 0, 0), # turned off for alignment
        legend.justification.top = "left",
        legend.justification.left = "bottom",
        legend.justification.bottom = "right",
        legend.justification.inside = c(1, 0),
        legend.location = "plot"
      )
}

make_prop_effect_size <- function(glm_pop2) {
  glm_pop2 |>
    dplyr::mutate(
      data = purrr::map(
        glm,
        ~ duckplyr::read_parquet_duckdb(.x) |>
          dplyr::collect()
      )
    ) |>
    dplyr::select(-glm, -n_sub) |>
    tidyr::unnest(data) |>
    dplyr::select(-n_sub) |>
    dplyr::mutate(
      d = cut(abs(pe / sigma), breaks = c(0, 0.2, 0.5, 0.8, Inf))
    ) |>
    dplyr::filter(!is.na(d)) |>
    dplyr::count(Task, d, type, name = "N") |>
    dplyr::group_by(Task, type) |>
    dplyr::mutate(
      Proportion = N / sum(N),
      Task = stringr::str_to_lower(Task)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = Task, fill = d, y = Proportion)) +
    ggplot2::facet_wrap(~type, scales = "free_x") +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::guides(fill = ggplot2::guide_legend("Cohen's d")) +
    ggplot2::xlab(NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
}

make_topo_bynetwork <- function(
  data_topo_gold_to_study_bynetwork,
  glm_pop2,
  at
) {
  mapping <- to_tbl(MNITemplate::getMNIPath("Brain", res = "2mm")) |>
    mask() |>
    mask_gray() |>
    mask_atlas() |>
    dplyr::mutate(index = 1:dplyr::n()) |>
    dplyr::select(-value)

  eff_size <- glm_pop2 |>
    dplyr::filter(type == "VOL") |>
    dplyr::mutate(
      glm = stringr::str_remove(glm, "/dcl01/smart/data/psadil/meta/"),
      data = purrr::map(
        glm,
        ~ duckplyr::read_parquet_duckdb(.x) |> dplyr::collect()
      )
    ) |>
    dplyr::select(-glm) |>
    tidyr::unnest(data) |>
    dplyr::mutate(g = abs(pe) / sigma * correct_d(n_sub)) |>
    dplyr::select(-z, -pe, -sigma, -n_sub) |>
    dplyr::left_join(mapping, by = dplyr::join_by(index)) |>
    dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
    dplyr::mutate(
      `Network Name` = dplyr::if_else(
        is.na(`Network Name`) & !is.na(label),
        "subcortical",
        `Network Name`
      ),
      Task = stringr::str_to_lower(Task)
    ) |>
    dplyr::filter(!is.na(`Network Name`)) |>
    dplyr::summarise(
      g = mean(g),
      .by = c(type, Task, `Network Name`)
    )

  data_topo_gold_to_study_bynetwork |>
    dplyr::left_join(eff_size) |>
    ggplot2::ggplot(ggplot2::aes(
      x = rho,
      y = `Network Name`,
      color = g
    )) +
    ggplot2::facet_grid(`N Sub` ~ Task) +
    ggplot2::geom_boxplot(outlier.alpha = 0.25) +
    ggplot2::scale_color_viridis_c(
      option = "turbo",
      guide = ggplot2::guide_colorbar("Cohen's d"),
      limits = c(0, NA),
    ) +
    ggplot2::scale_x_continuous(
      "Rank Correlation with Reference",
      labels = c(0, 0.5, 1),
      breaks = c(0, 0.5, 1)
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(legend.position = "bottom")
}

make_model <- function(
  data_model_gold_gold_to_study,
  data_model_study_to_study
) {
  prep <- data_model_gold_gold_to_study |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(confounds == "True" | stringr::str_detect(type, "UKB")) |>
    dplyr::filter(
      stringr::str_detect(replacement, "False", TRUE) |
        stringr::str_detect(sim, "gold")
    ) |>
    dplyr::filter(stringr::str_detect(type, "UKB_SMALL", TRUE)) |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0")
    ) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type)
    )

  gold_a <- prep |> dplyr::filter(sim == "gold")

  a <- prep |>
    dplyr::filter(sim == "simulation") |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = avg, color = type)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      linewidth = 0.5,
      width = 0,
      alpha = 0.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = avg - 2 * sem, ymax = avg + 2 * sem),
      linewidth = 3,
      width = 0
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = n_sub, y = avg, fill = type),
      data = gold_a,
      pch = 21,
      color = "gold"
    ) +
    ggplot2::scale_x_log10("N Sub") +
    ggplot2::ylab("Average Rank Correlation\nPrediction-Truth (gF)") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
    )

  b <- data_model_study_to_study |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(stringr::str_detect(replacement, "False", TRUE)) |>
    dplyr::filter(stringr::str_detect(type, "UKB_SMALL", TRUE)) |>
    dplyr::filter(confounds == "False", method == "consistency") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0")
    ) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type)
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = icc, color = type),
      alpha = 0.5
    ) +
    ggplot2::facet_wrap(~task, nrow = 2, scales = "free_x") +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper)) +
    ggplot2::scale_x_log10("N Sub") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::ylab("ICC(C,1) of Predictions") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
    )

  a +
    b +
    patchwork::plot_layout(ncol = 1) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        legend.margin = margin(0, 0, 0, 0), # turned off for alignment
        legend.justification.top = "left",
        legend.justification.left = "bottom",
        legend.justification.bottom = "right",
        legend.justification.inside = c(1, 0),
        legend.location = "plot"
      )
}


make_model2 <- function(data_model_gold_gold_to_study2) {
  prep <- data_model_gold_gold_to_study2 |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(confounds == "True" | stringr::str_detect(type, "UKB")) |>
    dplyr::filter(
      stringr::str_detect(replacement, "False", TRUE) |
        stringr::str_detect(sim, "gold")
    ) |>
    dplyr::filter(stringr::str_detect(type, "UKB_SMALL", TRUE)) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type),
      n_sub = dplyr::case_when(
        sim == "gold" & stringr::str_detect(measure, "^f.") ~ "ukb",
        sim == "gold" & stringr::str_detect(measure, "^f.", TRUE) ~ "hcp",
        TRUE ~ as.character(n_sub)
      ) |>
        factor(
          levels = c(
            "20",
            "40",
            "60",
            "80",
            "100",
            "hcp",
            "1000",
            "10000",
            "ukb"
          ),
          labels = c("20", "40", "60", "80", "100", "hcp", "1k", "10k", "ukb"),
          ordered = TRUE
        ),
      v = (a * b) / ((a + b)^2 * (a + b + 1))
    )

  avgs <- prep |>
    dplyr::summarise(
      avg = mean(avg),
      v = mean(v),
      .by = c(task, measure, n_sub)
    ) |>
    dplyr::mutate(
      dataset = dplyr::if_else(
        stringr::str_detect(measure, "^f."),
        "ukb",
        "hcpya"
      )
    )

  ss <- avgs |>
    dplyr::filter(n_sub %in% c("ukb", "hcp")) |>
    dplyr::filter(avg > 0) |>
    dplyr::distinct(measure, task)

  a <- prep |>
    dplyr::semi_join(ss) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = avg)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = type),
      outliers = FALSE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(group = measure, color = dataset),
      data = dplyr::semi_join(avgs, ss),
      alpha = 0.2
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::ylab(
      "Rate of Significant Rank Correlation"
    ) +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside")
    )

  b <- prep |>
    dplyr::filter(!n_sub %in% c("hcp", "ukb")) |>
    dplyr::semi_join(ss) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = v)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = type),
      outliers = FALSE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(group = measure, color = dataset),
      data = dplyr::semi_join(avgs, ss) |>
        dplyr::filter(!n_sub %in% c("hcp", "ukb")),
      alpha = 0.2
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::ylab(
      "Var of Significant Rank Correlation"
    ) +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside")
    )

  a +
    b +
    patchwork::plot_layout(nrow = 2) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        legend.margin = ggplot2::margin(0, 0, 0, 0), # turned off for alignment
        legend.justification.top = "left",
        legend.justification.left = "bottom",
        legend.justification.bottom = "right",
        legend.justification.inside = c(1, 0),
        legend.box = "horizontal"
      )
}

make_model2_neg <- function(data_model_gold_gold_to_study2) {
  prep <- data_model_gold_gold_to_study2 |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(confounds == "True" | stringr::str_detect(type, "UKB")) |>
    dplyr::filter(
      stringr::str_detect(replacement, "False", TRUE) |
        stringr::str_detect(sim, "gold")
    ) |>
    dplyr::filter(stringr::str_detect(type, "UKB_SMALL", TRUE)) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type),
      n_sub = dplyr::case_when(
        sim == "gold" & stringr::str_detect(measure, "^f.") ~ "ukb",
        sim == "gold" & stringr::str_detect(measure, "^f.", TRUE) ~ "hcp",
        TRUE ~ as.character(n_sub)
      ) |>
        factor(
          levels = c(
            "20",
            "40",
            "60",
            "80",
            "100",
            "hcp",
            "1000",
            "10000",
            "ukb"
          ),
          labels = c("20", "40", "60", "80", "100", "hcp", "1k", "10k", "ukb"),
          ordered = TRUE
        ),
      v = (a * b) / ((a + b)^2 * (a + b + 1))
    )

  avgs <- prep |>
    dplyr::summarise(
      avg = mean(avg),
      v = mean(v),
      .by = c(task, measure, n_sub)
    ) |>
    dplyr::mutate(
      dataset = dplyr::if_else(
        stringr::str_detect(measure, "^f."),
        "ukb",
        "hcpya"
      )
    )

  ss <- avgs |>
    dplyr::filter(n_sub %in% c("ukb", "hcp")) |>
    dplyr::filter(avg > 0) |>
    dplyr::distinct(measure, task)

  a <- prep |>
    dplyr::anti_join(ss) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = avg)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = type),
      outliers = FALSE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(group = measure, color = dataset),
      data = dplyr::anti_join(avgs, ss),
      alpha = 0.2
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::ylab(
      "Rate of Significant Rank Correlation"
    ) +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside")
    )

  b <- prep |>
    dplyr::filter(!n_sub %in% c("hcp", "ukb")) |>
    dplyr::anti_join(ss) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = v)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = type),
      outliers = FALSE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(group = measure, color = dataset),
      data = dplyr::anti_join(avgs, ss) |>
        dplyr::filter(!n_sub %in% c("hcp", "ukb")),
      alpha = 0.2
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::ylab(
      "Rate of Significant Rank Correlation"
    ) +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside")
    )

  a +
    b +
    patchwork::plot_layout(nrow = 2) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        legend.margin = ggplot2::margin(0, 0, 0, 0), # turned off for alignment
        legend.justification.top = "left",
        legend.justification.left = "bottom",
        legend.justification.bottom = "right",
        legend.justification.inside = c(1, 0),
        legend.box = "horizontal"
      )
}

make_model3 <- function(
  data_model_gold_gold_to_study3,
  data_model_study_to_study3
) {
  prep <- data_model_gold_gold_to_study3 |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(confounds == "True" | stringr::str_detect(type, "UKB")) |>
    dplyr::filter(
      stringr::str_detect(replacement, "False", TRUE)
    ) |>
    dplyr::filter(stringr::str_detect(type, "UKB_SMALL", TRUE)) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type),
      n_sub = factor(
        n_sub,
        levels = c(
          "20",
          "40",
          "60",
          "80",
          "100",
          "1000",
          "10000"
        ),
        labels = c("20", "40", "60", "80", "100", "1k", "10k"),
        ordered = TRUE
      )
    )

  avgs <- prep |>
    dplyr::summarise(
      .estimate = mean(.estimate),
      .by = c(task, measure, n_sub)
    ) |>
    dplyr::mutate(
      dataset = dplyr::if_else(
        stringr::str_detect(measure, "^f."),
        "ukb",
        "hcpya"
      )
    )

  a <- prep |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = .estimate)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = type),
      outliers = FALSE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(color = dataset, group = measure),
      alpha = 0.1,
      data = avgs
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::ylab(
      "Product-Moment Correlation of Coefficients\n(Samples to Gold)"
    ) +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside")
    )

  prep_b <- data_model_study_to_study3 |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(stringr::str_detect(replacement, "False", TRUE)) |>
    dplyr::filter(
      stringr::str_detect(type, "UKB_SMALL", TRUE),
      confounds == "True" | stringr::str_detect(type, "UKB")
    ) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type),
      n_sub = factor(
        n_sub,
        levels = c(
          "20",
          "40",
          "60",
          "80",
          "100",
          "1000",
          "10000"
        ),
        labels = c("20", "40", "60", "80", "100", "1k", "10k"),
        ordered = TRUE
      )
    )

  avgs_b <- prep_b |>
    dplyr::summarise(
      .estimate = mean(.estimate),
      .by = c(n_sub, task, measure)
    ) |>
    dplyr::mutate(
      dataset = dplyr::if_else(
        stringr::str_detect(measure, "^f."),
        "ukb",
        "hcpya"
      )
    )

  b <- prep_b |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = .estimate),
      alpha = 0.5
    ) +
    ggplot2::facet_wrap(~task, nrow = 2, scales = "free_x") +
    ggplot2::geom_boxplot(ggplot2::aes(fill = type), outliers = FALSE) +
    ggplot2::geom_line(
      ggplot2::aes(group = measure, color = dataset),
      alpha = 0.1,
      data = avgs_b
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::ylab("ICC(C,1) of Predictions") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside")
    )

  a +
    b +
    patchwork::plot_layout(ncol = 1) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        legend.margin = ggplot2::margin(0, 0, 0, 0), # turned off for alignment
        legend.justification.top = "left",
        legend.justification.left = "bottom",
        legend.justification.bottom = "right",
        legend.justification.inside = c(1, 0),
        legend.location = "plot",
        legend.box = "horizontal"
      )
}


make_model_model <- function(
  data_model_gold_gold_to_study,
  data_model_gold_gold_to_study2,
  data_model_gold_gold_to_study3,
  data_model_study_to_study,
  data_model_study_to_study3
) {
  a <- data_model_gold_gold_to_study |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(confounds == "False") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0"),
      sim == "simulation"
    ) |>
    ggplot2::ggplot(aes(x = n_sub, y = avg, color = type)) +
    ggplot2::facet_wrap(~model, scales = "free_x") +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      linewidth = 0.5,
      width = 0,
      alpha = 0.5
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = avg - 2 * sem, ymax = avg + 2 * sem),
      linewidth = 3,
      width = 0
    ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(x = n_sub, y = avg, fill = type),
      data = dplyr::filter(
        data_model_gold_gold_to_study,
        task == "EMOTION",
        measure %in% c("PMAT24_A_CR", "f.20016.2.0"),
        sim == "gold",
        confounds == "False"
      ),
      pch = 21,
      color = "gold"
    ) +
    ggplot2::scale_x_log10("N Sub") +
    ggplot2::ylab("Average Rank Correlation (CI)\nPrediction-Truth (gF)") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme(legend.position = "bottom")

  b <- data_model_study_to_study |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(confounds == "False", method == "consistency") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0")
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = icc, color = type),
      alpha = 0.5
    ) +
    ggplot2::facet_wrap(~model) +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper)) +
    ggplot2::xlab("N Sub") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::ylab("ICC(C,1) of Predictions")

  a2 <- data_model_gold_gold_to_study2 |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(confounds == "False") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0"),
      sim == "simulation"
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = avg, color = type)) +
    ggplot2::facet_wrap(~model, scales = "free_x") +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower, ymax = upper),
      width = 0
    ) +
    ggplot2::geom_jitter(
      mapping = ggplot2::aes(x = n_sub, y = avg, fill = type),
      data = dplyr::filter(
        data_model_gold_gold_to_study2,
        measure %in% c("PMAT24_A_CR", "f.20016.2.0"),
        sim == "gold",
        confounds == "False"
      ),
      pch = 21,
      color = "gold",
      alpha = 0.3,
      height = 0.05,
      width = 0
    ) +
    ggplot2::scale_x_log10("N Sub") +
    ggplot2::ylab(
      "Rate of Significant Rank Correlation\nFor Fluid Intelligence Prediction"
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme(legend.position = "bottom")

  a3 <- data_model_gold_gold_to_study3 |>
    dplyr::filter(confounds == "False") |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0")
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = .estimate, color = type)) +
    ggplot2::facet_wrap(~model, scales = "free_x") +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .lower, ymax = .upper)) +
    ggplot2::scale_x_log10("N Sub") +
    ggplot2::ylab(
      "Product-Moment Correlation of Coefficients\n(Samples to Gold)"
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme(legend.position = "bottom")

  b3 <- data_model_study_to_study3 |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(confounds == "False") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0")
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = .estimate, color = type),
      alpha = 0.5
    ) +
    ggplot2::facet_wrap(~model) +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = .lower, ymax = .upper)) +
    ggplot2::xlab("N Sub") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::ylab("ICC(C,1) of Coefficients")

  a +
    b +
    a2 +
    a3 +
    b3 +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}


make_all_cog <- function(
  data_model_gold_gold_to_study,
  data_model_study_to_study
) {
  .data_model_gold_gold_to_study <- data_model_gold_gold_to_study |>
    dplyr::filter(
      model == "RIDGE_CV",
      stringr::str_detect(type, "UKB_SMALL", TRUE),
      stringr::str_detect(replacement, "False", TRUE) | sim == "gold"
    ) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type),
      n_sub = dplyr::case_when(
        sim == "gold" & stringr::str_detect(measure, "^f.") ~ "ukb",
        sim == "gold" & stringr::str_detect(measure, "^f.", TRUE) ~ "hcp",
        TRUE ~ as.character(n_sub)
      ),
      n_sub = factor(
        n_sub,
        levels = c(
          "20",
          "40",
          "60",
          "80",
          "100",
          "hcp",
          "1000",
          "10000",
          "ukb"
        ),
        labels = c("20", "40", "60", "80", "100", "hcp", "1k", "10k", "ukb"),
        ordered = TRUE
      )
    )

  avgs <- .data_model_gold_gold_to_study |>
    dplyr::summarise(avg = mean(avg), .by = c(task, measure, n_sub)) |>
    dplyr::mutate(
      dataset = dplyr::if_else(
        stringr::str_detect(measure, "^f."),
        "ukb",
        "hcpya"
      )
    )

  a <- .data_model_gold_gold_to_study |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = avg)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = type), outliers = FALSE) +
    ggplot2::geom_line(
      ggplot2::aes(group = measure, color = dataset),
      alpha = 0.1,
      data = avgs
    ) +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::ylab(
      "Average Rank Correlation"
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside"),
    )

  avgs_b <- data_model_study_to_study |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(stringr::str_detect(replacement, "False", TRUE)) |>
    dplyr::filter(
      stringr::str_detect(type, "UKB_SMALL", TRUE),
      stringr::str_detect(method, "consis"),
      confounds == "True" | stringr::str_detect(type, "UKB")
    ) |>
    dplyr::summarise(icc = mean(icc), .by = c(n_sub, task, measure)) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      n_sub = factor(
        n_sub,
        levels = c("20", "40", "60", "80", "100", "1000", "10000"),
        labels = c("20", "40", "60", "80", "100", "1k", "10k"),
        ordered = TRUE
      ),
      dataset = dplyr::if_else(
        stringr::str_detect(measure, "^f."),
        "ukb",
        "hcpya"
      )
    )

  b <- data_model_study_to_study |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(stringr::str_detect(replacement, "False", TRUE)) |>
    dplyr::filter(stringr::str_detect(type, "UKB_SMALL", TRUE)) |>
    dplyr::filter(
      confounds == "True" | stringr::str_detect(type, "UKB"),
      method == "consistency"
    ) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type),
      n_sub = factor(
        n_sub,
        levels = c("20", "40", "60", "80", "100", "1000", "10000"),
        labels = c("20", "40", "60", "80", "100", "1k", "10k"),
        ordered = TRUE
      )
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = icc),
      alpha = 0.5
    ) +
    ggplot2::facet_wrap(~task, nrow = 2, scales = "free_x") +
    ggplot2::geom_boxplot(ggplot2::aes(fill = type), outliers = FALSE) +
    ggplot2::geom_line(
      ggplot2::aes(group = measure, color = dataset),
      alpha = 0.1,
      data = avgs_b
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::ylab("ICC(C,1) of Predictions") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside")
    )

  a +
    b +
    patchwork::plot_layout(nrow = 2) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_minimal(base_size = 8) &
    ggplot2::theme(
      legend.margin = ggplot2::margin(0, 0, 0, 0), # turned off for alignment
      legend.justification.top = "left",
      legend.justification.left = "bottom",
      legend.justification.bottom = "right",
      legend.justification.inside = c(1, 0),
      legend.location = "plot",
      legend.box = "horizontal"
    )
}

make_model_all_icc <- function(data_model_study_to_study, type) {
  data_model_study_to_study |>
    dplyr::filter(
      stringr::str_detect(task, "EMOTION", TRUE),
      confounds == "True"
    ) |>
    dplyr::filter(type == .env$type) |>
    dplyr::mutate(max_avg = max(icc), .by = c(measure)) |>
    dplyr::mutate(
      measure = stringr::str_replace_all(measure, "_", "\\\\_"),
      measure = factor(measure),
      measure = forcats::fct_reorder(measure, max_avg)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = measure)) +
    ggplot2::facet_wrap(~task) +
    ggplot2::geom_raster(ggplot2::aes(fill = icc)) +
    ggplot2::scale_fill_viridis_c(option = "turbo", name = "ICC(C,1)") +
    ggplot2::xlab("N Sub") +
    ggplot2::ylab("Instrument") +
    ggplot2::theme_minimal(base_size = 7) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = unit(8, "pt")
    )
}


make_tikz <- function(p, file, width, height) {
  ggplot2::ggsave(
    file,
    p,
    tikzDevice::tikz,
    width = width,
    height = height,
    standAlone = TRUE
  )
  file
}

plot_sigmas <- function(.data) {
  .data |>
    tidyr::pivot_longer(c(sigma2_b, sigma2_w, var.data)) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = value, color = type)
    ) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(model ~ name) +
    ggplot2::scale_color_viridis_d(option = "turbo")
}

make_model_sigmas <- function(data_model_study_to_study) {
  data_model_study_to_study |>
    dplyr::filter(
      method == "consistency",
      task == "EMOTION",
      confounds == "False"
    ) |>
    dplyr::select(n_sub:var.data) |>
    plot_sigmas()
}

make_model_sigmas3 <- function(data_model_study_to_study3) {
  data_model_study_to_study3 |>
    dplyr::filter(task == "EMOTION", confounds == "False") |>
    dplyr::select(type:n_sub, sigma2_b:var.data) |>
    na.omit() |>
    plot_sigmas()
}

make_model_model_ukb <- function(
  data_model_gold_gold_to_study,
  data_model_gold_gold_to_study2,
  data_model_gold_gold_to_study3,
  data_model_study_to_study,
  data_model_study_to_study3
) {
  a <- data_model_gold_gold_to_study |>
    dplyr::filter(stringr::str_detect(type, "UKB")) |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(confounds == "False") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0"),
      sim == "simulation"
    ) |>
    ggplot2::ggplot(aes(x = n_sub, y = avg, color = type)) +
    ggplot2::facet_wrap(~model, scales = "free_x") +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(
      ymin = avg - 2 * sem,
      ymax = avg + 2 * sem
    )) +
    ggplot2::ylab("Average Rank Correlation (CI)\nPrediction-Truth (gF)") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme(legend.position = "bottom")

  b <- data_model_study_to_study |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(confounds == "False", method == "consistency") |>
    dplyr::filter(stringr::str_detect(type, "UKB")) |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0")
    ) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = icc, color = type),
      alpha = 0.5
    ) +
    ggplot2::facet_wrap(~model) +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper)) +
    ggplot2::xlab("N Sub") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::ylab("ICC(C,1) of Predictions")

  a2 <- data_model_gold_gold_to_study2 |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(confounds == "False") |>
    dplyr::filter(stringr::str_detect(type, "UKB")) |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0"),
      sim == "simulation"
    ) |>
    ggplot2::ggplot(aes(x = n_sub, y = avg, color = type)) +
    ggplot2::facet_wrap(~model, scales = "free_x") +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
    ggplot2::ylab(
      "Rate of Significant Rank Correlation\nFor Fluid Intelligence Prediction"
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme(legend.position = "bottom")

  a3 <- data_model_gold_gold_to_study3 |>
    dplyr::filter(stringr::str_detect(type, "UKB")) |>
    dplyr::filter(confounds == "False") |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0")
    ) |>
    ggplot2::ggplot(aes(x = n_sub, y = .estimate, color = type)) +
    ggplot2::facet_wrap(~model, scales = "free_x") +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(aes(ymin = .lower, ymax = .upper)) +
    ggplot2::ylab(
      "Product-Moment Correlation of Coefficients\n(Samples to Gold)"
    ) +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::theme(legend.position = "bottom")

  b3 <- data_model_study_to_study3 |>
    dplyr::filter(stringr::str_detect(type, "UKB")) |>
    dplyr::filter(task == "EMOTION") |>
    dplyr::filter(confounds == "False") |>
    dplyr::filter(
      measure %in% c("PMAT24_A_CR", "f.20016.2.0")
    ) |>
    ggplot2::ggplot(aes(x = n_sub, y = .estimate, color = type), alpha = 0.5) +
    ggplot2::facet_wrap(~model) +
    ggplot2::geom_line() +
    ggplot2::geom_errorbar(aes(ymin = .lower, ymax = .upper)) +
    ggplot2::xlab("N Sub") +
    ggplot2::scale_color_viridis_d(option = "turbo") +
    ggplot2::ylab("ICC(C,1) of Coefficients")

  a +
    b +
    a2 +
    a3 +
    b3 +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

make_model_sigmas_ukb <- function(data_model_study_to_study) {
  data_model_study_to_study |>
    dplyr::filter(
      method == "consistency",
      task == "EMOTION",
      confounds == "False",
      stringr::str_detect(type, "UKB")
    ) |>
    dplyr::select(n_sub:var.data) |>
    plot_sigmas() +
    ggplot2::ggtitle("Variances for Predictions")
}

make_model_sigmas3_ukb <- function(data_model_study_to_study3) {
  data_model_study_to_study3 |>
    dplyr::filter(
      task == "EMOTION",
      confounds == "False",
      stringr::str_detect(type, "UKB")
    ) |>
    plot_sigmas() +
    ggplot2::ggtitle("Variances for Coefficients")
}

make_table_of_studies_without_peaks <- function(study_peaks, dst) {
  study_peaks |>
    dplyr::distinct(type, Task, iter, n_sub, threshold) |>
    dplyr::count(type, Task, n_sub, threshold) |>
    dplyr::filter(n < 100) |>
    dplyr::mutate(
      `Proportion With Peak` = n / 100,
      Task = stringr::str_to_lower(Task)
    ) |>
    dplyr::select(-n) |>
    gt::gt() |>
    gt::as_latex() |>
    as.character() |>
    readr::write_lines(dst)
  dst
}

make_modelroi <- function(data_modelroi_gold_gold_to_study) {
  .data_model_gold_gold_to_study <- data_modelroi_gold_gold_to_study |>
    dplyr::filter(
      model == "RIDGE_CV",
      stringr::str_detect(type, "UKB_SMALL", TRUE),
      stringr::str_detect(replacement, "False", TRUE) | sim == "gold",
      stringr::str_detect(type, "UKB")
    ) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type),
      n_sub = dplyr::case_when(
        sim == "gold" & stringr::str_detect(measure, "^f.") ~ "ukb",
        sim == "gold" & stringr::str_detect(measure, "^f.", TRUE) ~ "hcp",
        TRUE ~ as.character(n_sub)
      ),
      n_sub = factor(
        n_sub,
        levels = c(
          "20",
          "40",
          "60",
          "80",
          "100",
          "hcp",
          "1000",
          "10000",
          "ukb"
        ),
        labels = c("20", "40", "60", "80", "100", "hcp", "1k", "10k", "ukb"),
        ordered = TRUE
      )
    )

  avgs <- .data_model_gold_gold_to_study |>
    dplyr::summarise(avg = mean(avg), .by = c(task, measure, n_sub)) |>
    dplyr::mutate(
      dataset = dplyr::if_else(
        stringr::str_detect(measure, "^f."),
        "ukb",
        "hcpya"
      )
    )

  .data_model_gold_gold_to_study |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = avg)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = type), outliers = FALSE) +
    ggplot2::geom_line(
      ggplot2::aes(group = measure, color = dataset),
      alpha = 0.2,
      data = avgs
    ) +
    ggplot2::scale_color_manual(values = c(viridisLite::turbo(4)[3])) +
    ggplot2::scale_fill_manual(values = c(viridisLite::turbo(4)[3])) +
    ggplot2::ylab(
      "Average Rank Correlation"
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::theme_minimal()
}

make_model_r2 <- function(data_model_gold_gold_to_study_r2) {
  prep <- data_model_gold_gold_to_study_r2 |>
    dplyr::filter(model == "RIDGE_CV") |>
    dplyr::filter(confounds == "True" | stringr::str_detect(type, "UKB")) |>
    dplyr::filter(
      stringr::str_detect(replacement, "False", TRUE) |
        stringr::str_detect(sim, "gold")
    ) |>
    dplyr::filter(stringr::str_detect(type, "UKB_SMALL", TRUE)) |>
    dplyr::mutate(
      task = stringr::str_to_lower(task),
      type = stringr::str_to_lower(type),
      n_sub = dplyr::case_when(
        sim == "gold" & stringr::str_detect(measure, "^f.") ~ "ukb",
        sim == "gold" & stringr::str_detect(measure, "^f.", TRUE) ~ "hcp",
        TRUE ~ as.character(n_sub)
      ),
      n_sub = factor(
        n_sub,
        levels = c(
          "20",
          "40",
          "60",
          "80",
          "100",
          "hcp",
          "1000",
          "10000",
          "ukb"
        ),
        labels = c(
          "20",
          "40",
          "60",
          "80",
          "100",
          "hcp",
          "1k",
          "10k",
          "ukb"
        ),
        ordered = TRUE
      ),
      avg = dplyr::if_else(avg < 0, 0, avg)
    )

  avgs <- prep |>
    dplyr::summarise(avg = mean(avg), .by = c(task, measure, n_sub)) |>
    dplyr::mutate(
      dataset = dplyr::if_else(
        stringr::str_detect(measure, "^f."),
        "ukb",
        "hcpya"
      )
    )

  prep |>
    ggplot2::ggplot(ggplot2::aes(x = n_sub, y = avg)) +
    ggplot2::facet_wrap(~task, scales = "free_x", nrow = 2) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = type), outliers = FALSE) +
    ggplot2::geom_line(
      ggplot2::aes(color = dataset, group = measure),
      data = avgs,
      alpha = 0.2
    ) +
    ggplot2::xlab("N Sub") +
    ggplot2::ylab("Average Coefficient of Determination") +
    ggplot2::coord_cartesian(ylim = c(0, NA)) +
    ggplot2::scale_color_manual(values = c("blue", viridisLite::turbo(4)[3])) +
    ggplot2::scale_fill_viridis_d(option = "turbo") +
    ggplot2::guides(
      colour = ggplot2::guide_legend(position = "inside"),
      fill = ggplot2::guide_legend(position = "inside"),
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      legend.margin = ggplot2::margin(0, 0, 0, 0), # turned off for alignment
      legend.justification.top = "left",
      legend.justification.left = "bottom",
      legend.justification.bottom = "right",
      legend.justification.inside = c(1, 0),
      legend.location = "plot",
      legend.box = "horizontal"
    )
}

write_regions <- function(
  rois_pop,
  rois_pop_ukb,
  dst = "analyses/tables/top_ten_regions.tsv"
) {
  dplyr::bind_rows(rois_pop, rois_pop_ukb) |>
    dplyr::filter(n_parcels == 400) |>
    dplyr::mutate(
      r = dplyr::row_number(dplyr::desc(abs(estimate))),
      .by = c(Task, n_parcels, type)
    ) |>
    dplyr::filter(r < 11) |>
    dplyr::mutate(estimate = statistic / sqrt(parameter)) |>
    dplyr::select(type, Task, label, rank = r, estimate) |>
    dplyr::arrange(type, Task, rank) |>
    readr::write_tsv(dst)
  dst
}

write_model_performance <- function(
  data_model_gold_gold_to_study,
  data_model_gold_gold_to_study_r2,
  dst = "analyses/tables/measures.tsv"
) {
  cors <- data_model_gold_gold_to_study |>
    dplyr::filter(
      replacement == "True",
      confounds == "True" | type == "UKB",
      model == "RIDGE_CV"
    ) |>
    dplyr::select(type, task, n_sub, measure, avg)

  r2 <- data_model_gold_gold_to_study_r2 |>
    dplyr::filter(
      replacement == "True",
      confounds == "True" | type == "UKB",
      model == "RIDGE_CV"
    ) |>
    dplyr::select(type, task, n_sub, measure, avg)

  dplyr::bind_rows(list(cor = cors, r2 = r2), .id = "m") |>
    tidyr::pivot_wider(names_from = m, values_from = avg) |>
    dplyr::mutate(task = stringr::str_to_lower(task)) |>
    dplyr::arrange(type, measure, task, n_sub) |>
    readr::write_tsv(dst)
  dst
}

write_peaks <- function(
  gold_tested,
  n_peaks = 10,
  dst = "analyses/tables/top_ten_peaks.tsv"
) {
  gold_tested |>
    dplyr::filter(n_parcels == 400) |>
    dplyr::mutate(d = abs(statistic) / sqrt(parameter)) |>
    dplyr::select(Task, n_parcels, label, type, d) |>
    dplyr::slice_max(
      order_by = d,
      by = c(type, Task, label),
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::slice_max(
      order_by = d,
      by = c(type, Task),
      n = n_peaks,
      with_ties = FALSE
    ) |>
    dplyr::mutate(
      rank = rank(abs(d) * -1, ties.method = "first"),
      .by = c(type, Task)
    ) |>
    dplyr::arrange(type, Task, rank, d) |>
    readr::write_tsv(dst)
  dst
}

make_peak_avg_bysize <- function(study_to_gold_distances, glm_pop2) {
  peaks <- study_to_gold_distances |>
    dplyr::filter(type == "VOL") |>
    dplyr::summarise(
      avg_d = mean(d),
      .by = c(type, Task, n_sub, x, y, z, threshold)
    )

  mapping <- to_tbl(MNITemplate::getMNIPath("Brain", res = "2mm")) |>
    mask() |>
    mask_gray() |>
    mask_atlas() |>
    dplyr::mutate(index = 1:dplyr::n()) |>
    dplyr::select(-value)

  glm_pop2 |>
    dplyr::filter(type == "VOL") |>
    dplyr::mutate(
      glm = stringr::str_remove(glm, "/dcl01/smart/data/psadil/meta/"),
      data = purrr::map(
        glm,
        ~ duckplyr::read_parquet_duckdb(.x) |> dplyr::collect()
      )
    ) |>
    dplyr::select(-glm) |>
    tidyr::unnest(data) |>
    dplyr::mutate(g = abs(pe) / sigma * correct_d(n_sub)) |>
    dplyr::select(-z, -pe, -sigma, -n_sub) |>
    dplyr::left_join(mapping, by = dplyr::join_by(index)) |>
    dplyr::full_join(peaks, by = dplyr::join_by(type, Task, x, y, z)) |>
    na.omit() |>
    dplyr::mutate(g = cut(g, breaks = c(0, 0.1, .3, .5, Inf))) |>
    dplyr::summarise(avg_d = mean(avg_d), .by = c(g, threshold))
}


make_peak_avg_bynetwork <- function(study_to_gold_distances, glm_pop2, at) {
  peaks <- study_to_gold_distances |>
    dplyr::filter(type == "VOL") |>
    dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
    dplyr::mutate(
      `Network Name` = dplyr::if_else(
        is.na(`Network Name`) & !is.na(label),
        "subcortical",
        `Network Name`
      )
    ) |>
    dplyr::filter(!is.na(`Network Name`)) |>
    dplyr::summarise(
      avg_d = mean(d),
      .by = c(type, Task, n_sub, `Network Name`, threshold, iter)
    )

  mapping <- to_tbl(MNITemplate::getMNIPath("Brain", res = "2mm")) |>
    mask() |>
    mask_gray() |>
    mask_atlas() |>
    dplyr::mutate(index = 1:dplyr::n()) |>
    dplyr::select(-value)

  eff_size <- glm_pop2 |>
    dplyr::filter(type == "VOL") |>
    dplyr::mutate(
      glm = stringr::str_remove(glm, "/dcl01/smart/data/psadil/meta/"),
      data = purrr::map(
        glm,
        ~ duckplyr::read_parquet_duckdb(.x) |> dplyr::collect()
      )
    ) |>
    dplyr::select(-glm) |>
    tidyr::unnest(data) |>
    dplyr::mutate(g = abs(pe) / sigma * correct_d(n_sub)) |>
    dplyr::select(-z, -pe, -sigma, -n_sub) |>
    dplyr::left_join(mapping, by = dplyr::join_by(index)) |>
    dplyr::left_join(at, by = dplyr::join_by(x, y, z)) |>
    dplyr::mutate(
      `Network Name` = dplyr::if_else(
        is.na(`Network Name`) & !is.na(label),
        "subcortical",
        `Network Name`
      )
    ) |>
    dplyr::filter(!is.na(`Network Name`)) |>
    dplyr::summarise(
      g = mean(g),
      .by = c(type, Task, `Network Name`)
    )

  peaks |>
    dplyr::left_join(
      eff_size,
      by = dplyr::join_by(type, Task, `Network Name`)
    ) |>
    dplyr::summarise(
      avg_d = mean(avg_d),
      g = unique(g),
      .by = c(type, Task, n_sub, `Network Name`, threshold)
    ) |>
    dplyr::filter(`Network Name` == "somatomotor", threshold > 0)
}
