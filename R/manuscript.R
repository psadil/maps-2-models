make_roi <- function(
    data_roi_study_to_gold,
    data_roi_study_to_study,
    data_roi_sub_to_sub,
    file) {
  a <- data_roi_study_to_gold |>
    dplyr::filter(n_parcels == 400) |>
    ggplot2::ggplot(
      ggplot2::aes(x = n_sub, y = prop, group = label, color = abs(d))
    ) +
    ggplot2::geom_point(alpha = 0.2) +
    ggplot2::geom_line(alpha = 0.2) +
    facet_wrap(~Task) +
    scale_y_continuous(
      "Proportion Simulations w/\nActivity in Most Active ROI",
      limits = c(0, 1),
      labels = c(0, 0.5, 1),
      breaks = c(0, 0.5, 1)
    ) +
    scale_x_continuous(
      "N Sub",
      breaks = c(40, 80),
      labels = c(40, 80)
    ) +
    scale_color_viridis_c(
      "Abs. Effect Size",
      option = "turbo",
      limits = c(0, NA),
      n.breaks = 3
    )

  b <- data_roi_study_to_study |>
    dplyr::filter(
      forcats::fct_match(n_parcels, "N Parcels: 400")
    ) |>
    ggplot(aes(y = n_sub, x = phi, group = n_sub)) +
    facet_wrap(~Task) +
    ggdist::stat_dots(quantiles = 50) +
    ylab("N Sub") +
    scale_x_continuous(
      "Phi Coefficient\n(Study-Study)",
      limits = c(-0.25, 1),
      breaks = c(-0.25, 0.5),
      labels = c(-0.25, 0.5)
    )

  cc <- data_roi_sub_to_sub |>
    dplyr::mutate(
      Task = factor(Task),
      Task = forcats::fct_rev(Task)
    ) |>
    dplyr::filter(n_parcels == 400) |>
    ggplot(aes(y = Task, x = rho)) +
    ggdist::stat_dotsinterval(
      quantiles = 100
    ) +
    scale_x_continuous(
      "Product-Moment Correlation\n(Sub-Sub)",
      limits = c(-0.75, 1),
      breaks = c(-0.75, 0, 0.75),
      labels = c(-0.75, 0, 0.75)
    )

  a + b + cc +
    patchwork::plot_layout(ncol = 1, heights = c(1, 1.5, 1)) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    ggplot2::theme_gray(base_size = 8) +
      ggplot2::theme(
        legend.position = "bottom",
        legend.key.size = unit(8, "pt")
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
    active_threshold = 0.02) {
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
      n_sub, iter, label, `Label Name`,
      `Full component name`, n_parcels
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
    ggplot(aes(x = n_sub, y = prop, group = l)) +
    geom_ribbon(
      aes(
        ymin = qbeta(0.05 / 2, 5, 100 - 5 + 1),
        ymax = qbeta(1 - 0.05 / 2, 5, 100 - 5)
      ),
      alpha = 0.01,
      fill = "lightblue",
      linetype = "dashed",
      color = "black"
    ) +
    geom_point(show.legend = FALSE, alpha = 0.2) +
    geom_line(show.legend = FALSE, alpha = 0.2) +
    facet_wrap(~n_parcels) +
    scale_y_continuous(
      "Proportion Simulations w/\nActivity in Most Active ROI",
      limits = c(0, 0.12)
    ) +
    xlab("N Sub") +
    theme_gray(base_size = 12)
}

make_prop_active_most_active_roi_ptfce <- function(data_roi_study_to_gold) {
  data_roi_study_to_gold |>
    ggplot(aes(x = n_sub, y = prop, group = label, color = abs(d))) +
    geom_point(alpha = 0.2) +
    geom_line(alpha = 0.2) +
    facet_grid(n_parcels ~ Task) +
    scale_y_continuous(
      "Proportion Simulations w/\nActivity in Most Active ROI",
      limits = c(0, 1),
      labels = c(0, 0.5, 1),
      breaks = c(0, 0.5, 1)
    ) +
    scale_x_continuous(
      "N Sub",
      breaks = c(40, 80),
      labels = c(40, 80)
    ) +
    scale_color_viridis_c(
      "Abs. Effect Size",
      option = "turbo",
      limits = c(0, NA),
      n.breaks = 3
    ) +
    theme(legend.position = "bottom")
}

make_peaks <- function(
    data_peak_study_to_gold,
    data_peak_study_to_study,
    data_peak_sub_to_sub) {
  a <- data_peak_study_to_gold |>
    ggplot(aes(x = within, group = peak, y = n_simulations, color = Value)) +
    geom_point(alpha = 0.2) +
    geom_line(alpha = 0.2) +
    facet_grid(n_sub ~ Task) +
    scale_y_continuous(
      "Proportion Simulations w/\nPeak in Radius",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c(0, 0.5, 1)
    ) +
    scale_x_continuous(
      "Radius (mm)",
      limits = c(0, 20),
      breaks = c(0, 10, 20),
      labels = c(0, 10, 20)
    ) +
    scale_color_viridis_c(
      "Peak Height",
      option = "turbo",
      limits = c(0, NA),
      n.breaks = 3
    )

  b <- data_peak_study_to_study |>
    ggplot(aes(y = n_sub, x = d)) +
    facet_wrap(~Task) +
    ggdist::stat_dots(quantiles = 100) +
    ylab("N Sub") +
    scale_x_continuous(
      "Distance Between\nHighest Peaks\n(Study-Study)"
    )

  cc <- data_peak_sub_to_sub |>
    dplyr::mutate(
      Task = factor(Task),
      Task = forcats::fct_rev(Task)
    ) |>
    ggplot(aes(y = Task, x = d)) +
    ggdist::stat_dots(quantiles = 100) +
    ylab(NULL) +
    scale_x_continuous(
      "Distance Between\nHighest Peaks\n(Sub-Sub)"
    )

  a + b + cc +
    patchwork::plot_layout(nrow = 1, widths = c(4, 2, 1)) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    theme_gray(base_size = 6) +
      theme(
        legend.position = "bottom",
        legend.key.size = unit(6, "pt")
      )
}

make_peak_bysize <- function(space, data_topo_gold, gold_peaks, at) {
  gold_peaks_ <- gold_peaks |>
    dplyr::select(Task, m) |>
    tidyr::unnest(m) |>
    dplyr::left_join(at) |>
    dplyr::group_by(Task, label) |>
    dplyr::slice_max(
      order_by = Value,
      n = 1,
      with_ties = FALSE
    ) |> # grab highest peak from each label
    dplyr::group_by(Task) |>
    dplyr::slice_max(
      order_by = Value,
      n = 10,
      with_ties = FALSE
    ) |> # grab highest 10 peaks (distinct labels)
    dplyr::ungroup() |>
    dplyr::distinct(Task, x, y, z, Value)

  space |>
    dplyr::mutate(Value = Value / sqrt(n_pop)) |>
    dplyr::filter(Value > 0.2) |>
    dplyr::group_by(Value, n_sub, corrp_thresh, Task) |>
    dplyr::summarise(d = mean(d, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(
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
    ) |>
    ggplot(aes(y = d, x = Value)) +
    geom_point(alpha = 0.1, shape = 20) +
    facet_grid(`N Sub` ~ Task) +
    scale_x_continuous(
      "Gold Standard Peak Cohen's d",
      breaks = c(0, 1),
      labels = c(0, 1)
    ) +
    scale_y_log10("avg dist(Gold Standard Peak, Study Peak) (mm)") +
    theme_gray(base_size = 8)
}

make_peak_bynetwork <- function(space, data_topo_gold, gold_peaks, at) {
  gold_peaks_ <- gold_peaks |>
    dplyr::select(Task, m) |>
    tidyr::unnest(m) |>
    dplyr::left_join(at) |>
    dplyr::group_by(Task, label) |>
    dplyr::slice_max(
      order_by = Value,
      n = 1,
      with_ties = FALSE
    ) |> # grab highest peak from each label
    dplyr::group_by(Task) |>
    dplyr::slice_max(
      order_by = Value,
      n = 10,
      with_ties = FALSE
    ) |> # grab highest 10 peaks (distinct labels)
    dplyr::ungroup() |>
    dplyr::distinct(Task, x, y, z, Value)

  space |>
    dplyr::mutate(Value = Value / sqrt(n_pop)) |>
    dplyr::filter(Value > 0.2, !is.na(`Network Name`)) |>
    dplyr::group_by(n_sub, corrp_thresh, Task, `Network Name`, iter) |>
    dplyr::summarise(
      d = mean(d, na.rm = TRUE),
      Value = mean(Value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
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
    ) |>
    dplyr::filter(!is.na(d)) |>
    ggplot(aes(y = `Network Name`, x = d, color = Value)) +
    geom_boxplot(outlier.shape = NA) +
    scattermore::geom_scattermore(
      pointsize = 5,
      position = position_jitter(width = 0),
      alpha = 0.5
    ) +
    facet_grid(`N Sub` ~ Task) +
    scale_color_viridis_c(
      option = "turbo",
      guide = guide_colorbar("Cohen's d")
    ) +
    ylab("Network") +
    scale_x_continuous("avg dist(Gold Standard Peak, Study Peak) (mm)") +
    theme_gray(base_size = 8)
}

make_topo <- function(
    data_topo_gold,
    data_topo_gold_to_study,
    data_topo_study_to_study,
    data_topo_sub_to_sub) {
  a <- data_topo_gold |>
    ggplot() +
    facet_wrap(~Task) +
    scattermore::geom_scattermore(aes(x = cope, y = sigma), alpha = 0.1, pointsize = 1) +
    geom_function(fun = function(x) sqrt(20) * x / qt(0.001, 19), xlim = c(-200, 0)) +
    geom_function(fun = function(x) sqrt(100) * x / qt(0.001, 99), xlim = c(-200, 0), color = "gray50") +
    geom_function(fun = function(x) sqrt(100) * x / qt(0.999, 99), xlim = c(0, 200), color = "gray50") +
    geom_function(fun = function(x) sqrt(20) * x / qt(0.999, 19), xlim = c(0, 200)) +
    scale_y_continuous(
      limits = c(0, 150)
    ) +
    scale_x_continuous(
      limits = c(-200, 200),
      labels = c(-150, 0, 150),
      breaks = c(-150, 0, 150)
    ) +
    xlab(expression(beta ~ mean)) +
    ylab(expression(beta ~ SD))

  b <- data_topo_gold_to_study |>
    dplyr::mutate(`N Sub` = factor(n_sub)) |>
    ggplot(aes(x = rho, y = `N Sub`)) +
    facet_wrap(~Task) +
    geom_boxplot(outlier.shape = NA, size = .1) +
    geom_point(
      position = position_jitter(width = 0),
      shape = 20,
      size = 0.1,
      alpha = 0.1
    ) +
    xlab("Rank Correlation\n(Gold to Study)")

  cc <- data_topo_study_to_study |>
    ggplot(aes(x = rr, y = `N Sub`, color = Task)) +
    geom_line(aes(group = Task)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper)) +
    xlab("Pairwise Rank Correlation\n(Study to Study)") +
    theme(
      legend.position = "bottom"
    )

  d <- data_topo_sub_to_sub |>
    dplyr::filter(!is.na(rho), stringr::str_detect(task, "EMOTION", TRUE)) |>
    ggplot(aes(x = rho, y = task)) +
    ggdist::stat_dots(
      quantiles = 100
    ) +
    ylab("Task") +
    xlab("Pairwise Product-Moment Correlation\n(Sub to Sub)")

  a + b + cc + d +
    patchwork::plot_layout(ncol = 1) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    theme_gray(base_size = 8) +
      theme(
        legend.position = "bottom",
        legend.key.size = unit(8, "pt")
      )
}

make_prop_effect_size <- function(data_topo_gold) {
  data_topo_gold |>
    dplyr::filter(!is.na(d)) |>
    dplyr::count(Task, d, name = "N") |>
    dplyr::group_by(Task) |>
    dplyr::mutate(Proportion = N / sum(N)) |>
    ggplot(aes(x = Task, fill = d, y = Proportion)) +
    geom_col(position = "dodge") +
    guides(fill = guide_legend("Cohen's d")) +
    xlab(NULL) +
    theme_gray(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}

make_topo_bynetwork <- function(pop_cor_region) {
  pop_cor_region |>
    dplyr::mutate(
      f = atanh(rho),
      `Network Name` =
        dplyr::if_else(is.na(`Network Name`) & !is.na(label), "subcortical", `Network Name`)
    ) |>
    dplyr::filter(!is.na(rho) & is.finite(f)) |>
    dplyr::group_by(Task, n_sub, ContrastName, method, `Network Name`, iter) |>
    dplyr::summarise(
      f = mean(f),
      N = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      rr = tanh(f),
      `N Sub` = factor(n_sub)
    ) |>
    ggplot(aes(y = rr, x = `N Sub`, color = `Network Name`)) +
    facet_wrap(~Task) +
    geom_boxplot(outlier.alpha = 0.25) +
    scale_color_viridis_d(option = "turbo") +
    ylab("Rank Correlation with Reference") +
    theme(legend.position = "bottom")
}

make_model <- function(
    data_model_gold_gold_to_study,
    data_model_study_to_study,
    data_model_sub_to_sub) {
  a <- data_model_gold_gold_to_study |>
    dplyr::filter(confounds == "True") |>
    dplyr::filter(
      measure == "PMAT24_A_CR",
      type == "simulation",
      stringr::str_detect(task, "EMOTION", TRUE)
    ) |>
    ggplot(aes(x = n_sub, y = avg)) +
    facet_wrap(~task) +
    geom_line() +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    geom_point(
      mapping = aes(x = n_sub, y = statistic_rep),
      color = "gold",
      data = dplyr::filter(
        data_model_gold_gold_to_study,
        measure == "PMAT24_A_CR",
        type == "gold",
        stringr::str_detect(task, "EMOTION", TRUE),
        confounds == "True"
      )
    ) +
    scale_x_log10("N Sub") +
    ylab("Average Rank Correlation (CI)\nPrediction-Truth (gF)") +
    theme(legend.position = "bottom")

  b <- data_model_study_to_study |>
    dplyr::filter(confounds == "True") |>
    dplyr::filter(
      measure == "PMAT24_A_CR",
      stringr::str_detect(task, "EMOTION", TRUE)
    ) |>
    ggplot(aes(x = n_sub, y = icc, color = type), alpha = 0.5) +
    facet_wrap(~task) +
    geom_line() +
    geom_errorbar(aes(ymin = lower, ymax = upper)) +
    xlab("N Sub") +
    ylab("ICC")

  cc <- data_model_sub_to_sub |>
    dplyr::filter(confounds) |>
    dplyr::filter(stringr::str_detect(task, "EMOTION", TRUE)) |>
    ggplot(aes(x = r, y = task)) +
    ggdist::stat_dots(quantiles = 100) +
    ylab("Task") +
    xlab("Pairwise Rank Correlation of Features\n(Sub to Sub)")

  a + b + cc +
    patchwork::plot_layout(ncol = 1) +
    patchwork::plot_annotation(tag_levels = "a", tag_suffix = ")") &
    theme_gray(base_size = 8) +
      theme(
        legend.position = "bottom",
        legend.key.size = unit(8, "pt")
      )
}

make_all_cog <- function(data_model_gold_gold_to_study) {
  data_model_gold_gold_to_study |>
    dplyr::filter(
      type == "simulation",
      stringr::str_detect(task, "EMOTION", TRUE)
    ) |>
    dplyr::mutate(max_avg = max(avg), .by = c(measure)) |>
    dplyr::mutate(
      measure = stringr::str_replace_all(measure, "_", "\\\\_"),
      measure = factor(measure),
      measure = forcats::fct_reorder(measure, max_avg)
    ) |>
    ggplot(aes(x = n_sub, y = measure)) +
    facet_wrap(~task) +
    geom_raster(aes(fill = avg)) +
    scale_fill_viridis_c(option = "turbo", name = "Rank\nCorrelation") +
    xlab("N Sub") +
    ylab("Instrument") +
    theme_gray(base_size = 7) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(8, "pt")
    )
}

make_model_all_icc <- function(data_model_study_to_study, type) {
  data_model_study_to_study |>
    dplyr::filter(stringr::str_detect(task, "EMOTION", TRUE)) |>
    dplyr::filter(type == .env$type) |>
    dplyr::mutate(max_avg = max(icc), .by = c(measure)) |>
    dplyr::mutate(
      measure = stringr::str_replace_all(measure, "_", "\\\\_"),
      measure = factor(measure),
      measure = forcats::fct_reorder(measure, max_avg)
    ) |>
    ggplot(aes(x = n_sub, y = measure)) +
    facet_wrap(~task) +
    geom_raster(aes(fill = icc)) +
    scale_fill_viridis_c(option = "turbo", name = "ICC(C,1)") +
    xlab("N Sub") +
    ylab("Instrument") +
    theme_gray(base_size = 7) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(8, "pt")
    )
}


make_tikz <- function(p, file, width, height) {
  ggsave(
    file,
    p,
    tikzDevice::tikz,
    width = width,
    height = height,
    standAlone = TRUE
  )
  file
}
