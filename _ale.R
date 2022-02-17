# Sys.setenv(TAR_PROJECT = "ale")

library(targets)
library(tarchetypes)
library(rlang)
library(magrittr)

source(here::here("R", "ale.R"))
source(here::here("R", "spatial.R"))

library(future)
# library(future.callr)
# plan(callr)
library(future.batchtools)
plan(batchtools_sge, template = "tools/sge.tmpl")


list(
  tar_target(
    avail,
    here::here("data-raw", "cope5.txt"),
    format = "file"),
  tar_target(
    train,
    readr::read_lines(avail, num_threads=1, n_max=10000),
    format = "qs"),
  tar_target(
    test,
    readr::read_lines(avail, num_threads=1, skip=10000)),
  tar_target(n_sub, c(10, 20, 50, 100, 200, 500, 1000)),
  tar_target(n_study, c(1)),
  tar_target(iter, seq_len(50)),
  tar_target(
    cope_files,
    prep_cope_tbl(train, n_sub=n_sub, iter=iter, n_study=n_study),
    pattern = cross(n_sub, iter, n_study),
    format = "fst_tbl"),
  tar_target(
    small,
    dplyr::filter(cope_files, .data$n_sub %in% .env$n_sub),
    format = "fst_tbl",
    pattern = slice(n_sub, c(1,2,3))),
  tar_group_by(
    small_group,
    small,
    iter, n_sub, n_study, study,
    format = "fst_tbl"),
  tar_target(
    med,
    dplyr::filter(cope_files, .data$n_sub %in% .env$n_sub),
    format = "fst_tbl",
    pattern = slice(n_sub, c(4,5,6))),
  tar_group_by(
    med_group,
    med,
    iter, n_sub, n_study, study,
    format = "fst_tbl"),
  tar_target(
    large,
    dplyr::filter(cope_files, .data$n_sub %in% .env$n_sub),
    format = "fst_tbl",
    pattern = slice(n_sub, c(7))),
  tar_group_by(
    large_group,
    large,
    iter, n_sub, n_study, study,
    format = "fst_tbl"),
  tar_target(
    z_img_small,
    do_z(small_group),
    pattern = map(small_group),
    format = "fst_tbl",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge,
          template = "tools/sge.tmpl",
          resources = list(mem_free = "3G"))))),
  tar_target(
    z_img_med,
    do_z(med_group),
    pattern = map(med_group),
    format = "fst_tbl",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge,
          template = "tools/sge.tmpl",
          resources = list(mem_free = "20G"))))),
  tar_target(
    z_img_large,
    do_z(large_group),
    pattern = map(large_group),
    format = "fst_tbl",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge,
          template = "tools/sge.tmpl",
          resources = list(mem_free = "50G"))))),
  tar_target(
    z_img,
    dplyr::bind_rows(z_img_small, z_img_med, z_img_large),
    format = "fst_tbl"
    ),
  # tar_target(
  #   t_img,
  #   do_t(cope_files2),
  #   pattern = map(cope_files2),
  #   format = "fst_tbl",
  #   resources = tar_resources(
  #     future = tar_resources_future(
  #       plan = tweak(
  #         batchtools_sge,
  #         template = "tools/sge.tmpl",
  #         resources = list(mem_free = "15G"))))),
  # tar_target(
  #   ale_py_script,
  #   fs::path("python", "ale.py"),
  #   format = "file"),
  #  tar_target(
  #    t_img2,
  #    t_img |>
  #      dplyr::group_by(iter, n_study, n_sub) |>
  #      tar_group(),
  #    iteration = "group",
  #    format = "fst_tbl"),
  #  tar_target(
  #    ale,
  #    do_ale_py(t_img2, python_source = here::here(ale_py_script), condaenv = "meta"),
  #    cue = tar_cue(depend = FALSE),
  #    pattern = map(t_img2),
  #    format = "fst_tbl",
  #    resources = tar_resources(
  #      future = tar_resources_future(
  #        plan = tweak(
  #          batchtools_sge,
  #          template = "tools/sge.tmpl",
  #          resources = list(mem_free = "5G"))))),
  tar_target(
    z_pop,
    calc_z(test),
    format = "qs",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge,
          template = "tools/sge.tmpl",
          resources = list(mem_free = "200G"))))),
  tar_target(
    at,
    fst::read_fst(here::here("data-raw", "atlas.fst")),
    format = "fst_tbl"),
  tar_target(threshold, c(3.1)),
  tar_target(
      maxes,
      z_img |>
        dplyr::mutate(
          m = purrr::map(.data$z, ~get_maxes(.x, threshold=.env$threshold, minextent=0, at = at))),
      pattern = cross(threshold, map(z_img)),
      format = "qs",
      resources = tar_resources(
        future = tar_resources_future(
          plan = tweak(
            batchtools_sge,
            template = "tools/sge.tmpl",
            resources = list(mem_free = "1G"))))),
    tar_target(
      gold_peaks,
      get_maxes(neurobase::tempimg(z_pop), threshold=3.1, minextent=0, at=at),
      format = "fst_tbl"),
    tar_target(
     space,
     maxes |>
       dplyr::mutate(
         augmented = purrr::map(
           .data$m, 
           augment_distance, 
           reference=gold_peaks |> dplyr::filter(!is.na(label)))) |>
       dplyr::select(n_sub, n_study, augmented, iter) |>
       tidyr::unnest(augmented),
     format = "fst_tbl")
   # tar_target(
   #   comparison,
   #    tidy_ale(ale),
   #    pattern = map(ale),
   #    format = "fst_tbl"
   #  ),
  #  tar_target(
  #    ibma_py_script,
  #    fs::path("python", "ibma.py"),
  #    format = "file"),
  #  tar_target(
  #    ibma,
  #    do_ibma_py(ale, python_source = here::here(ibma_py_script), condaenv = "meta"),
  #    pattern = map(ale),
  #    format = "fst_tbl"),
  #  tar_target(
  #    ibma_summary,
  #    tidy_ibma(ibma),
  #    pattern = map(ibma),
  #    format = "fst_tbl")

)
