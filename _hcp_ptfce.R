Sys.setenv(TAR_PROJECT = "hcp_ptfce")

library(targets)
library(tarchetypes)
library(rlang)

source(here::here("R", "ale.R"))
source(here::here("R", "spatial.R"))
source(here::here("R", "tfce.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "loading.R"))
source(here::here("R", "poster.R"))
source(here::here("R", "hcp.R"))
source(here::here("R", "ptfce.R"))
source(here::here("R", "figures.R"))
source(here::here("R", "pairwise.R"))

Sys.setenv(
  NIIDIR = here::here("data-raw","hcp-niis-ptfce")
) # explicitly avoiding tracking this

targets::tar_option_set(
  format = "qs", 
  storage="worker", 
  packages = c("oro.nifti"))
library(future.callr)
plan(callr)


list(
  tar_group_by(
    contrasts,
    readr::read_csv(here::here("data-raw", "hcp", "contrasts.csv")) |>
      dplyr::filter(
        !Task=="EMOTION",
        (Task=="WM" & CopeNumber %in%  c(11)) |
          (Task=="GAMBLING" & CopeNumber %in%  c(6)) |
          (Task=="MOTOR" & CopeNumber %in% c(7)) |
          (Task=="LANGUAGE" & CopeNumber %in%  c(3)) |
          (Task=="SOCIAL" & CopeNumber %in%  c(3)) |
          (Task=="RELATIONAL" & CopeNumber %in%  c(3))
      ),
    Task, CopeNumber
  ),
  tar_target(
    test,
    contrasts |>
      dplyr::mutate(
        avail = purrr::map2(
          Task, CopeNumber,
          ~fs::dir_ls(
            "data-raw/hcp", 
            recurse = TRUE, 
            glob = glue::glue("*{.x}*cope{.y}.feat*cope1*"))),
        avail = purrr::map(
          avail,
          ~.x[
            stringr::str_detect(
              .x, 
              stringr::str_c(not_avail(), collapse = "|"), 
              negate = TRUE)])), 
    pattern = map(contrasts)
  ),
  tar_target(n_sub, c(20, 40, 60, 80, 100)),
  tar_target(iter, seq_len(100)),
  tar_target(
    tfce,
    test |>
      dplyr::mutate(
        tmp = purrr::map2(
          avail, ContrastName,
          ~do_ptfce(
            .x, 
            n_sub=n_sub, 
            iter=iter, 
            storage_dir=Sys.getenv("NIIDIR"), 
            flags=.y,
            resample=TRUE))) |>
      tidyr::unnest(tmp),
    pattern = cross(n_sub, iter, map(test))
  ),
  tar_target(
    tfce_pop,
    test |>
      dplyr::mutate(
        tmp = purrr::map2(
          avail, ContrastName,
          ~do_ptfce(
            .x, 
            n_sub=length(.x), 
            iter=0, 
            storage_dir=Sys.getenv("NIIDIR"), 
            flags=.y,
            resample = FALSE,
            enhance = FALSE))) |>
      tidyr::unnest(tmp),
    pattern = map(test)),
  tar_target(
    active0, 
    tfce |> dplyr::mutate(tmp = purrr::map(ptfce, get_active_ptfce)), 
    pattern = map(tfce)),
  tar_target(
    active, 
    tidyr::unnest(active0, tmp) |> dplyr::select(-ptfce, -avail), 
    format = format_arrow_table()),
  tar_target(at, make_atlas_full()),
  tar_target(n_parcels, c(200, 400, 600, 800, 1000)),
  tar_target(
    at_list, 
    make_atlas_full(n_parcels=n_parcels), 
    pattern = map(n_parcels)),
  tar_target(
    gold_peaks,
    dplyr::mutate(
      tfce_pop,
      m = purrr::map(
        .data$ptfce,
        ~get_ptfce_maxes_pop(q=.x, cluster_thresh=0.001, minextent=0)))
  ),
  tar_target(corrp_thresh, c(0.95)),
  tar_target(
    maxes,
    dplyr::mutate(
      tfce,
      m = purrr::map(
        .data$ptfce,
        ~get_ptfce_maxes(q=.x, corrp_thresh=corrp_thresh, minextent=0)),
      corrp_thresh = corrp_thresh),
    pattern = cross(map(tfce), corrp_thresh)
  ),  
  tar_target(
    space0,
    dplyr::left_join(
      maxes, 
      gold_peaks, 
      by=c("Task", "CopeNumber", "ContrastName", "tar_group"), 
      suffix=c(".study", ".ref")) |>
      dplyr::mutate(
        augmented = purrr::map2(
          m.study, m.ref,
          ~augment_distance(
            study = .x,
            reference = .y,
            vox_mm = 2
          )
        )
      ) |>
      dplyr::select(n_sub=n_sub.study, augmented, iter=iter.study, corrp_thresh, Task, CopeNumber, n_pop=n_sub.ref) |>
      tidyr::unnest(augmented)
  ),
  tar_target(space, add_labels(space=space0, at=at)),
  tar_target(
    tfce_null,
    do_ptfce(
      fs::dir_ls("data-raw/Fake2B"),
      n_sub=n_sub, 
      iter=iter, 
      storage_dir=here::here("data-raw/Fake2B-ptfce"), 
      flags="2BK-0BK",
      resample=TRUE),
    pattern = cross(n_sub, iter)
  ),
  tar_target(
    tfce_pop_null,
    do_ptfce(
      fs::dir_ls("data-raw/Fake2B"),
      n_sub=length(fs::dir_ls("data-raw/Fake2B")), 
      iter=0, 
      storage_dir=here::here("data-raw/Fake2B-ptfce"), 
      flags="2BK-0BK",
      resample=FALSE,
      enhance=FALSE)),
  tar_target(
    active0_null, 
    tfce_null |> 
      dplyr::mutate(tmp = purrr::map(ptfce, get_active_ptfce)), 
    pattern = map(tfce_null)),
  tar_target(
    active_null, 
    active0_null |>
      dplyr::select(-copes) |>
      tidyr::unnest(tmp) |> 
      dplyr::select(-ptfce), 
    format = format_arrow_table()),
  tar_target(
    rois,
    test |>
      dplyr::mutate(
        tmp = purrr::map(
          avail, 
          ~do_roi(
            .x, 
            n_sub = n_sub, 
            iter = iter, 
            at = at_list,
            resample=TRUE))) |>
      tidyr::unnest(tmp) |>
      dplyr::select(-avail),
    pattern = cross(n_sub, iter, map(test), at_list),
    format = "parquet"
  ),
  tar_target(
    rois_pop,
    test |>
      dplyr::mutate(
        tmp = purrr::map(
          avail,
          ~do_roi(
            .x, 
            iter = 0, 
            at = at_list))) |>
      tidyr::unnest(tmp) |>
      dplyr::select(-avail),
    pattern = cross(map(test), at_list),
    format = "parquet"
  ),
  tar_target(
    gold_tested,
    test_roi(rois_pop, Task, n_parcels),
    pattern = map(rois_pop),
    format = "parquet"
  ),
  tar_target(
    rois_tested,
    test_roi(rois, Task, n_parcels, iter),
    pattern = map(rois),
    format = "parquet"
  ),
  tar_target(
    rois_tested2,
    rois_tested, 
    format = format_arrow_table()
  ),
  tar_target(
    rois_pop2,
    rois_pop, 
    format = format_arrow_table()
  ),
  tar_target(
    data_topo_gold_to_study,
    make_data_topo_gold_to_study(tfce, tfce_pop, method="spearman"),
    cross(map(tfce))),
  tar_target(ContrastNames, contrasts$ContrastName),
  tar_target(
    pairwise, 
    cor_pairwise_ptfce(tfce, ContrastNames, n_sub, method="spearman"), 
    cross(n_sub, ContrastNames)),
  tar_target(
    data_topo_study_to_study, 
    make_data_topo_study_to_study(pairwise=pairwise), 
    format = "parquet"),
  tar_target(dataset, here::here("data-raw/hcp-parquet")),
  tarchetypes::tar_group_by(
    subtask, 
    get_subs_tasks(dataset), 
    sub,
    format = "parquet"),
  tar_target(
    rhos, 
    do_loo_cor(dataset, dplyr::select(subtask, sub, task)), 
    map(subtask),  
    format = "parquet"),
  tar_target(
    maxes_no_thresh,
    dplyr::mutate(
      tfce,
      m = purrr::map(
        .data$ptfce,
        ~get_ptfce_maxes_pop(q=.x, minextent=0, cluster_thresh = 0.0001)),
      corrp_thresh = 0),
    pattern = cross(map(tfce))
  ),  
  tar_target(
    space0_no_thresh,
    dplyr::left_join(
      maxes_no_thresh, 
      gold_peaks, 
      by=c("Task", "CopeNumber", "ContrastName", "tar_group"), 
      suffix=c(".study", ".ref")) |>
      dplyr::mutate(
        augmented = purrr::map2(
          m.study, m.ref,
          ~augment_distance(
            study = .x,
            reference = .y,
            vox_mm = 2
          )
        )
      ) |>
      dplyr::select(
        n_sub=n_sub.study, 
        augmented, 
        iter=iter.study, 
        corrp_thresh, 
        Task, 
        CopeNumber, 
        n_pop=n_sub.ref) |>
      tidyr::unnest(augmented)
  ),
  tar_target(space_no_thresh, add_labels(space=space0_no_thresh, at=at)),
  tar_target(
    data_peak_study_to_study,
    make_data_peak_study_to_study(
      at=at, 
      space=dplyr::bind_rows(space_no_thresh, space), 
      gold_peaks=gold_peaks), 
    format = "parquet"),
  tar_target(
    data_roi_study_to_gold,
    make_data_roi_study_to_gold(
      gold_tested=gold_tested,
      rois_tested=rois_tested2), 
    format = "parquet"),
  tar_target(
    data_roi_study_to_study,
    make_data_roi_study_to_study(rois_tested=rois_tested2), 
    format = "parquet"),
  tar_target(
    data_roi_sub_to_sub,
    make_data_roi_sub_to_sub(rois_pop=rois_pop2), 
    format = "parquet"),
  tar_target(
    data_peak_study_to_gold,
    make_data_peak_study_to_gold(
      at=at, 
      gold_peaks=gold_peaks,
      space=space), 
    format = "parquet"),
  tar_target(
    zstats, 
    test |>
      tidyr::unnest(avail) |> 
      dplyr::mutate(zstat=stringr::str_replace(avail, "/cope1.nii", "/zstat1.nii")) |>
      dplyr::select(Task, CopeNumber, ContrastName, zstat),
    format = "parquet"),
  tar_target(
    sub_peaks,
    dplyr::mutate(
      zstats,
      m = purrr::map(
        .data$zstat,
        ~get_ptfce_maxes_sub(zstat=.x, cluster_thresh=0.001, minextent=0))),
    pattern = map(zstats)
  ),
  tar_target(
    space_sub,
  dplyr::left_join(
    sub_peaks, 
    gold_peaks, 
    by=c("Task", "CopeNumber", "ContrastName"), 
    suffix=c(".study", ".ref")) |>
    dplyr::mutate(
      augmented = purrr::map2(
        m.study, m.ref,
        ~augment_distance(
          study = .x,
          reference = .y,
          vox_mm = 2
        )
      )
    ) |>
      dplyr::select(augmented, Task, CopeNumber) |>
      tidyr::unnest(augmented)
  ),
  tar_target(
    data_peak_sub_to_sub,
    make_data_peak_sub_to_sub(
      at=at, 
      space_sub=space_sub, 
      gold_peaks=gold_peaks), 
    format = "parquet")
)
