Sys.setenv(TAR_PROJECT = "hcp")

library(targets)
library(tarchetypes)
library(rlang)

source(here::here("R", "ale.R"))
source(here::here("R", "spatial.R"))
source(here::here("R", "tfce.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "loading.R"))
source(here::here("R", "poster.R"))

Sys.setenv(
  NIIDIR = here::here("data-raw","hcp-niis")
) # explicitly avoiding tracking this

targets::tar_option_set(format = "qs", storage="worker", retrieval="worker", error="continue")
library(future.callr)
plan(callr)


list(
  tar_group_by(
    contrasts,
    readr::read_csv(here::here("data-raw", "hcp", "contrasts.csv")) |>
      dplyr::filter(
        !Task=="EMOTION",
          (Task=="WM" & CopeNumber %in%  c(9, 10, 11)) |
          (Task=="GAMBLING" & CopeNumber %in%  c(1, 2, 3)) |
          (Task=="MOTOR" & CopeNumber %in% c(1, 7, 8)) |
          (Task=="LANGUAGE" & CopeNumber %in%  c(1, 2, 3)) |
          (Task=="SOCIAL" & CopeNumber %in%  c(1, 2, 3)) |
          (Task=="RELATIONAL" & CopeNumber %in%  c(1, 2, 3))
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
            glob = glue::glue("*{.x}*cope{.y}.feat*cope1*")))), 
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
          ~do_tfce2(.x, n_sub=n_sub, iter=iter, n=1000, storage_dir=Sys.getenv("NIIDIR"), flags=.y))) |>
      tidyr::unnest(tmp),
    pattern = cross(n_sub, iter, map(test))
  ),
  tar_target(
    tfce_pop,
    test |>
      dplyr::mutate(
        tmp = purrr::map2(
          avail, ContrastName,
          ~do_tfce_pop2(.x, n_sub=length(.x), iter=0, storage_dir=Sys.getenv("NIIDIR"), flags=.y))) |>
      tidyr::unnest(tmp),
    pattern = map(test)),
  tar_target(corrp_thresh, c(0.01, 0.95)),
  tar_target(
    maxes,
    dplyr::mutate(
      tfce,
      m = purrr::map2(
        .data$tfce_corrp_tstat, .data$tfce_tstat,
        ~get_tfce_maxes(corrp=.x, tstat=.y, corrp_thresh=corrp_thresh, minextent=0)),
      corrp_thresh = corrp_thresh),
    pattern = cross(map(tfce), corrp_thresh)
  ),
  tar_target(
    gold_peaks,
    dplyr::mutate(
      tfce_pop,
      m = purrr::map(
        .data$tstat,
        ~get_tfce_maxes_pop(tstat=.x, cluster_thresh=0.001, minextent=0)))
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
            reference = .y
          )
        )
      ) |>
      dplyr::select(n_sub=n_sub.study, augmented, iter=iter.study, corrp_thresh, Task, CopeNumber, n_pop=n_sub.ref) |>
      tidyr::unnest(augmented)
  ),
  tar_target(at, make_atlas_full()),
  tar_target(space, add_labels(space=space0, at=at))
)
