make_data_peak_study_to_gold <- function(
  study_to_gold_distances,
  at,
  gold_tested,
  n_peaks = 10,
  dscalar_file = "data-raw/Schaefer2018_400Parcels_7Networks_order.dscalar.nii"
) {
  dscalar <- ciftiTools::read_cifti(dscalar_file)
  # curvature <- read_cifti_with_mwall("data-raw/HCP_S1200_GroupAvg_v1/S1200.All.curvature_MSMAll.32k_fs_LR.dscalar.nii")

  cii <- ciftiTools::read_cifti("data-raw/palm/extrema.dscalar.nii")

  subcortext_at <- tibble::tibble(
    index = seq_along(cii$meta$subcort$labels),
    label.cii = as.character(cii$meta$subcort$labels),
    structure = "subcort"
  )

  labels <- read_cifti_labels()

  # curvatures <- as_tibble(cbind(rowMeans(curvature$data$cortex_left), rowMeans(curvature$data$cortex_right))) |>
  #   rename(cortex_left = V1, cortex_right = V2) |>
  #   mutate(index = row_number()) |>
  #   pivot_longer(starts_with("cortex"), names_to = "structure", values_to = "curvature")

  gold_most <- gold_tested |>
    dplyr::mutate(d = statistic / sqrt(parameter)) |>
    dplyr::select(Task, n_parcels, label, type, d)

  ats <- tibble::as_tibble(cbind(
    dscalar$data$cortex_left,
    dscalar$data$cortex_right
  )) |>
    dplyr::rename(cortex_left = V1, cortex_right = V2) |>
    dplyr::mutate(index = dplyr::row_number()) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("cortex"),
      names_to = "structure",
      values_to = "label.cii"
    ) |>
    dplyr::left_join(labels, by = dplyr::join_by(label.cii == index)) |>
    dplyr::select(-label.cii) |>
    dplyr::rename(label.cii = label) |>
    # left_join(curvatures, by = join_by(index, structure)) |>
    dplyr::bind_rows(subcortext_at)

  gold_ranks <- study_to_gold_distances |>
    dplyr::distinct(type, Task, index, x, y, z, Value, structure) |>
    dplyr::left_join(dplyr::select(at, x, y, z, label)) |>
    dplyr::left_join(ats) |>
    dplyr::mutate(label = dplyr::if_else(is.na(label), label.cii, label)) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::slice_max(
      order_by = Value,
      by = c(type, Task, label),
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::slice_max(
      order_by = Value,
      by = c(type, Task),
      n = n_peaks,
      with_ties = FALSE
    ) |>
    dplyr::mutate(
      rank = rank(abs(Value) * -1, ties.method = "first"),
      .by = c(type, Task)
    ) |>
    dplyr::select(type, Task, index, x, y, z, rank, label)

  study_to_gold_distances |>
    dplyr::right_join(gold_ranks) |>
    dplyr::summarize(
      within_2 = any(d < 2),
      within_4 = any(d < 4),
      within_6 = any(d < 6),
      within_8 = any(d < 8),
      within_10 = any(d < 10),
      within_20 = any(d < 20),
      .by = c(Task, n_sub, rank, type, iter, index, x, y, z, label, threshold)
    ) |>
    dplyr::summarize(
      dplyr::across(
        tidyselect::starts_with("within"),
        sum
      ),
      .by = c(Task, n_sub, rank, type, index, x, y, z, label, threshold)
    ) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("within"),
      values_to = "n_simulations",
      names_to = "within",
      names_pattern = "within_([[:digit:]]+)",
      names_transform = as.integer
    ) |>
    dplyr::mutate(
      n_simulations = n_simulations / 100,
      n_sub = glue::glue("N Sub: {n_sub}"),
      n_sub = factor(
        n_sub,
        levels = c(
          "N Sub: 20",
          "N Sub: 40",
          "N Sub: 60",
          "N Sub: 80",
          "N Sub: 100"
        )
      )
    ) |>
    dplyr::left_join(gold_most, by = dplyr::join_by(Task, type, label))
}


make_data_peak_study_to_study <- function(
  study_to_gold_distances,
  at,
  n_peaks = 10,
  dscalar_file = "data-raw/Schaefer2018_400Parcels_7Networks_order.dscalar.nii"
) {
  dscalar <- ciftiTools::read_cifti(dscalar_file)
  # curvature <- read_cifti_with_mwall("data-raw/HCP_S1200_GroupAvg_v1/S1200.All.curvature_MSMAll.32k_fs_LR.dscalar.nii")

  cii <- ciftiTools::read_cifti("data-raw/palm/extrema.dscalar.nii")

  subcortext_at <- tibble::tibble(
    index = seq_along(cii$meta$subcort$labels),
    label.cii = cii$meta$subcort$labels,
    structure = "subcort"
  )
  n_parcels <- unique(at$n_parcels)
  n_networks <- unique(at$n_networks)

  labels <- readr::read_delim(
    here::here(
      "data-raw",
      "Parcellations",
      "MNI",
      "fsleyes_lut",
      glue::glue(
        "Schaefer2018_{n_parcels}Parcels_{n_networks}Networks_order.lut"
      )
    ),
    col_names = c("dscalar", "R", "G", "B", "label.cii"),
    show_col_types = FALSE
  ) |>
    dplyr::select(dscalar, label.cii)

  ats <- tibble::as_tibble(
    cbind(
      dscalar$data$cortex_left,
      dscalar$data$cortex_right
    )
  ) |>
    dplyr::rename(cortex_left = V1, cortex_right = V2) |>
    dplyr::mutate(index = dplyr::row_number()) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("cortex"),
      names_to = "structure",
      values_to = "dscalar"
    ) |>
    dplyr::left_join(labels, by = dplyr::join_by(dscalar)) |>
    dplyr::bind_rows(subcortext_at)

  gold_ranks <- study_to_gold_distances |>
    dplyr::distinct(type, Task, index, x, y, z, Value, structure) |>
    dplyr::left_join(dplyr::select(at, x, y, z, label)) |>
    dplyr::left_join(ats) |>
    dplyr::mutate(label = dplyr::if_else(is.na(label), label.cii, label)) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::slice_max(
      order_by = Value,
      by = c(type, Task, label),
      n = 1,
      with_ties = FALSE
    ) |>
    dplyr::slice_max(
      order_by = Value,
      by = c(type, Task),
      n = n_peaks,
      with_ties = FALSE
    ) |>
    dplyr::mutate(
      rank = rank(abs(Value) * -1, ties.method = "first"),
      .by = c(type, Task)
    ) |>
    dplyr::select(type, Task, index, x, y, z, rank, label, structure)

  vol_ds <- study_to_gold_distances |>
    dplyr::right_join(
      gold_ranks |> dplyr::filter(type %in% c("VOL", "UKB")),
      by = dplyr::join_by(type, Task, index, x, y, z)
    ) |>
    dplyr::select(
      n_sub,
      x = x.study,
      y = y.study,
      z = z.study,
      iter,
      Task,
      label,
      type,
      threshold,
      rank
    ) |>
    dplyr::group_nest(Task, n_sub, label, type, threshold, rank) |>
    dplyr::mutate(
      d = purrr::map(
        data,
        ~ .x |>
          dplyr::arrange(iter) |>
          dplyr::select(-iter) |>
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

  surf_ds <- study_to_gold_distances |>
    dplyr::select(-x, -y, -z, -`Cluster Index`) |>
    dplyr::right_join(
      gold_ranks |> dplyr::filter(!(type %in% c("VOL", "UKB"))),
      by = dplyr::join_by(type, Task, index, structure)
    ) |>
    dplyr::select(
      n_sub,
      iter,
      Task,
      label,
      index = index.study,
      type,
      x = x.study,
      y = y.study,
      z = z.study,
      threshold,
      rank
    ) |>
    dplyr::group_nest(Task, n_sub, label, type, threshold, rank) |>
    dplyr::mutate(
      surf = dplyr::if_else(
        stringr::str_detect(label, "left"),
        "_L_",
        "_R_"
      ),
      surface = purrr::pmap_chr(
        list(task = Task, type = type, surf = surf),
        function(task, type, surf) {
          fs::dir_ls(
            Sys.getenv("PALMDIR"),
            glob = glue::glue("*iter-0*{task}*{type}*{surf}*mid*")
          )
        }
      ),
      corrected_areas = stringr::str_replace(surface, "midthickness", "area"),
      distances = purrr::pmap(
        list(.d = data, surface = surface, corrected_areas = corrected_areas),
        .get_surf_dist_pairs
      )
    ) |>
    dplyr::select(-surface, -corrected_areas, -surf, -data) |>
    tidyr::unnest(distances)

  dplyr::bind_rows(vol_ds, surf_ds) |>
    dplyr::mutate(n_sub = factor(n_sub))
}

.get_surf_dist_pairs <- function(.d, surface, corrected_areas) {
  distance_list <- list()
  for (iter in seq_len(nrow(.d))) {
    # -1 because wb_command indexing starts at 0
    # and all stored values are starting at 1
    index <- .d$index[[iter]] - 1
    if (is.na(index)) {
      out <- dplyr::select(.d, x, y, z) |>
        dist() |>
        as.matrix() |>
        corrr::as_cordf() |>
        corrr::shave() |>
        corrr::stretch(na.rm = TRUE) |>
        dplyr::rename(d = r)
      return(out)
    } else {
      distances <- tempfile(fileext = ".shape.gii")
      check <- ciftiTools::run_wb_cmd(
        glue::glue(
          "-surface-geodesic-distance {surface} {index} {distances} -corrected-areas {corrected_areas}"
        ),
        intern = FALSE
      )
      if (!check) {
        stop("Failed to calculate geodesic distance")
      }
      gii <- gifti::read_gifti(distances)
      to_append <- tibble::tibble(
        x = iter,
        y = seq_len(nrow(.d)),
        d = gii$data[[1]][.d$index]
      ) |>
        list()
      fs::file_delete(distances)
    }
    distance_list <- append(distance_list, to_append)
  }

  dplyr::bind_rows(distance_list) |>
    dplyr::filter(x < y) |>
    dplyr::mutate(dplyr::across(c(x, y), as.character))
}

make_peaks_gold_table <- function(data_peak_study_to_gold, dst) {
  data_peak_study_to_gold |>
    dplyr::filter(n_parcels == 400) |>
    dplyr::distinct(type, Task, label, rank) |>
    dplyr::arrange(type, Task, rank) |>
    dplyr::mutate(label = stringr::str_remove(label, "7Networks_")) |>
    dplyr::mutate(Task = stringr::str_to_lower(Task)) |>
    gt::gt() |>
    gt::sub_missing() |>
    gt::as_latex() |>
    as.character() |>
    readr::write_lines(dst)
}
