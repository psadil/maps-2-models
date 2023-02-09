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
source(here::here("R", "hcp.R"))

Sys.setenv(
  NIIDIR = here::here("data-raw","hcp-niis")
) # explicitly avoiding tracking this

targets::tar_option_set(format = "qs", storage="worker")
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
        # (Task=="WM" & CopeNumber %in%  c(9, 10, 11)) |
        # (Task=="GAMBLING" & CopeNumber %in%  c(1, 2, 3)) |
        # (Task=="MOTOR" & CopeNumber %in% c(1, 7, 21)) |
        # (Task=="LANGUAGE" & CopeNumber %in%  c(1, 2, 3)) |
        # (Task=="SOCIAL" & CopeNumber %in%  c(1, 2, 3)) |
        # (Task=="RELATIONAL" & CopeNumber %in%  c(1, 2, 3))
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
  tar_target(corrp_thresh, c(0.95)),
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
            reference = .y,
            vox_mm = 2
          )
        )
      ) |>
      dplyr::select(n_sub=n_sub.study, augmented, iter=iter.study, corrp_thresh, Task, CopeNumber, n_pop=n_sub.ref) |>
      tidyr::unnest(augmented)
  ),
  tar_target(at, make_atlas_full()),
  tar_target(space, add_labels(space=space0, at=at)),
  tar_target(
    prop, 
    get_active(tfce), 
    pattern = map(tfce)),
  tar_target(propall, prop, format = format_arrow_table()),
  tar_target(
    prop0, 
    get_active0(tfce), 
    pattern = map(tfce)),
  tar_target(prop0all, prop0, format = format_arrow_table()),
  tar_target(pop_d, get_pop_d(tfce_pop)),
  tar_target(method, c("spearman", "pearson")),
  tar_target(pop_cor, cor_w_pop(tfce, tfce_pop, method=method), cross(map(tfce), method)),
  tar_target(ContrastNames, contrasts$ContrastName), 
  tar_target(pairwise, cor_pairwise(tfce, ContrastNames, n_sub, method=method), cross(n_sub, ContrastNames, method)),
  tar_target(gray0, get_all_gray(tfce), pattern = map(tfce)),
  tar_target(gray, gray0, format = format_arrow_table()),
  tar_target(pop_cor_region, cor_w_pop_by_region(tfce=tfce, tfce_pop=tfce_pop, at=at, method="spearman"), map(tfce)),
  tar_target(pop_cor2, cor_w_pop2(tfce, tfce_pop, method="spearman"), map(tfce)),
  tar_target(eff_by_roi, get_eff_by_roi(test, at=at), format = format_arrow_table()),
  tar_target(n_parcels, c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)),
  tar_target(at_list, make_atlas_full(n_parcels=n_parcels), pattern = map(n_parcels)),
  tar_target(space_, add_labels(space=space0, at=at_list), pattern = map(at_list)),
  tar_target(space2, space_, format = format_arrow_table())
)
