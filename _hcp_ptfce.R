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
            resample=FALSE,
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
    tidyr::unnest(active0_null, tmp) |> dplyr::select(-ptfce), 
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
    pop_cor2,
    cor_w_pop2(tfce, tfce_pop, method="pearson"),
    cross(map(tfce))),
)
