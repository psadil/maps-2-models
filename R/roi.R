
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
        ~stats::t.test(value ~ 1, data = .x) |>
          broom::tidy()),
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
      iter = 0) 
}


roi_from_cifti <- function(file, n_parcels){
  parc <- load_parc(n_parcels=n_parcels)
  parc <- ciftiTools::parc_add_subcortex(parc)
  xii <- ciftiTools::read_xifti(file)
  ciftiTools::apply_parc(xii, parc, FUN=mean, na.rm=TRUE) |>
    tibble::as_tibble(rownames = "label") |>
    dplyr::rename(value = V1) |>
    dplyr::filter(!label == "???") |>
    na.omit()
}


roi_from_nifti <- function(file, n_parcels){
  
  at <- make_atlas_full(n_parcels=n_parcels)
  
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
  if (nrow(vol) > 0){
    vol <- vol |>
      dplyr::mutate(VOL = purrr::map(VOL, roi_from_nifti, n_parcels))  |>
      tidyr::unnest(VOL)
  }else{
    vol <- dplyr::select(vol, -VOL)
  }
  
  surface <- test |>
    dplyr::select(-VOL, -MSMALL) |>
    na.omit() |>
    dplyr::mutate(SURFACE = purrr::map(SURFACE, roi_from_cifti, n_parcels)) |>
    tidyr::unnest(SURFACE)
  
  dplyr::bind_rows(
    list(MSMALL=msmall, VOL=vol, SURFACE=surface), 
    .id = "type") |>
    dplyr::select(
      type, Task, CopeNumber, ContrastName, sub, label, value
    ) |>
    dplyr::mutate(
      n_parcels = .env$n_parcels
    )
}


make_data_roi_study_to_gold <- function(
    gold_tested, 
    rois_tested, 
    data_topo_gold, 
    at_list) {
  
  gold_most <- gold_tested |>
    dplyr::mutate(
      r = dplyr::row_number(dplyr::desc(abs(estimate))),
      .by = c(Task, n_parcels, type)
    ) |>
    dplyr::filter(r < 11) |>
    dplyr::select(Task, n_parcels, label, type)
  
  gold_roi <- data_topo_gold |>
    dplyr::left_join(
      at_list, 
      by = dplyr::join_by(x, y, z), 
      relationship = "many-to-many") |>
    dplyr::semi_join(gold_most, by = dplyr::join_by(Task, n_parcels, label)) |>
    dplyr::filter(!is.na(label)) |>
    dplyr::summarise(
      d = mean(cope / sigma * correct_d(n_sub)),
      .by = c(Task, label, n_parcels)
    )
  
  
  rois_tested |>
    dplyr::semi_join(gold_most, by = dplyr::join_by(Task, n_parcels, label, type)) |>
    dplyr::summarise(
      prop = mean(active),
      .by = c(Task, n_parcels, label, n_sub, type)
    ) |>
    dplyr::left_join(gold_roi)
}

make_data_roi_study_to_study <- function(rois_tested) {
  rois_tested |>
    dplyr::select(Task, n_parcels, n_sub, iter, active, label) |>
    dplyr::collect() |>
    dplyr::group_nest(Task, n_parcels, n_sub) |>
    dplyr::mutate(
      data = purrr::map(
        data,
        ~ .x |>
          dplyr::select(label, active, iter) |>
          tidyr::pivot_wider(names_from = label, values_from = active) |>
          dplyr::arrange(iter)
      ),
      phi = purrr::map(
        data,
        ~ dplyr::select(.x, -iter) |>
          as.matrix() |>
          sim.phi() |>
          corrr::as_cordf() |>
          corrr::shave() |>
          corrr::stretch(na.rm = TRUE) |>
          dplyr::rename(phi = r)
      )
    ) |>
    dplyr::mutate(
      n_sub = factor(
        n_sub,
        levels = c(20, 40, 60, 80, 100),
        labels = c(20, 40, 60, 80, 100)
      ),
      n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
      n_parcels = forcats::fct_relabel(
        n_parcels,
        .fun = ~ glue::glue("N Parcels: {.x}") |>
          as.character()
      )
    ) |>
    dplyr::select(-data) |>
    tidyr::unnest(phi)
}

make_data_roi_sub_to_sub <- function(rois_pop) {
  rois_pop |>
    dplyr::select(Task, n_parcels, Z, label, sub) |>
    dplyr::collect() |>
    dplyr::group_nest(Task, n_parcels) |>
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