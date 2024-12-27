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

make_data_peak_study_to_study <- function(at, space, gold_peaks) {
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
  
  space |>
    dplyr::filter(!is.na(study_ind)) |> # can happen in no voxels pass threshold
    dplyr::inner_join(gold_peaks_) |>
    dplyr::select(n_sub, x = x.study, y = y.study, z = z.study, iter, Task, peak_i, corrp_thresh) |>
    dplyr::group_nest(Task, n_sub, peak_i, corrp_thresh) |>
    dplyr::mutate(
      n_sub = factor(n_sub),
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
}

bootstrap_r <- function(.data, times=2000){
  
  boots <- boot::boot(
    data = .data$r,
    statistic = function(x,i) mean(x[i]), 
    R = 2000)
  
  ci <- boot::boot.ci(boots, type="perc")
  
  tibble::tibble(
    rr = boots$t0,
    lower = ci$percent[4],
    upper = ci$percent[5]
  )
}

make_data_topo_study_to_study <- function(pairwise, times=2000) {
  cors <- pairwise |>
    dplyr::filter(y > x) |>
    dplyr::group_nest(Task, n_sub, ContrastName, method) |>
    dplyr::mutate(
      `N Sub` = factor(n_sub),
      boots = purrr::map(
        data, bootstrap_r
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(boots)
  
}




make_data_peak_study_to_gold <- function(at, gold_peaks, space, data_topo_gold) {
  at_ <- at |>
    dplyr::mutate(dplyr::across(c(x, y, z), as.integer))
    
  gold_peaks_ <- gold_peaks |>
    dplyr::select(Task, m) |>
    tidyr::unnest(m) |>
    dplyr::left_join(at_, by = dplyr::join_by(x,y,z)) |>
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
    dplyr::distinct(Task, x, y, z) |>
    dplyr::left_join(data_topo_gold, by = dplyr::join_by(x,y,z,Task)) |>
    dplyr::mutate(
      Value = cope / sigma * correct_d(n_sub)
    ) |>
    dplyr::select(Task, x, y, z, Value)
  
  space |>
    dplyr::filter(
      !is.na(study_ind),
      corrp_thresh == 0.95
    ) |>
    dplyr::inner_join(gold_peaks_, by = dplyr::join_by(x,y,z,Task)) |>
    dplyr::group_by(Task, n_sub, x, y, z, iter) |>
    dplyr::summarize(
      within_2 = any(d < 2),
      within_4 = any(d < 4),
      within_6 = any(d < 6),
      within_8 = any(d < 8),
      within_10 = any(d < 10),
      within_20 = any(d < 20),
      .groups = "drop"
    ) |>
    dplyr::group_by(Task, n_sub, x, y, z) |>
    dplyr::summarize(
      dplyr::across(
        tidyselect::starts_with("within"),
        sum
      ),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("within"),
      values_to = "n_simulations",
      names_to = "within",
      names_pattern = "within_([[:digit:]]+)",
      names_transform = as.integer
    ) |>
    dplyr::left_join(gold_peaks_, by = dplyr::join_by(Task, x, y, z)) |>
    tidyr::unite(col = "peak", x, y, z) |>
    dplyr::mutate(
      n_simulations = n_simulations / 100,
      n_sub = glue::glue("N Sub: {n_sub}"),
      n_sub = factor(
        n_sub,
        levels = c("N Sub: 20", "N Sub: 40", "N Sub: 60", "N Sub: 80", "N Sub: 100")
      )
    ) 
}

make_data_topo_gold_to_study <- function(tfce, tfce_pop, storage_dir, method = "spearman") {
  # tfce <- tfce |>
  #   dplyr::filter(stringr::str_detect(tfce_corrp_tstat, glue::glue("flags-{ContrastName}_tfce")))
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
    dplyr::mutate(test = cope / sigma * correct_d(tfce_pop$n_sub)) |>
    dplyr::select(x, y, z, test)
  
  tfce |>
    dplyr::select(Task, CopeNumber, ContrastName, iter, n_sub) |>
    dplyr::bind_cols(
      dplyr::left_join(study, test, by = c("x", "y", "z")) |>
        dplyr::summarise(rho = cor(study, test, method = .env$method, use = "complete.obs"))
    ) |>
    dplyr::mutate(method = .env$method)
}
