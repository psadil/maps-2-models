# Sys.setenv(TAR_PROJECT = "ale")

library(targets)
library(rlang)
library(magrittr)

source(here::here("R", "ale.R"))

library(future)
# library(future.callr)
# plan(callr)
library(future.batchtools)
plan(batchtools_sge, template = "tools/sge.tmpl")



list(
  tar_target(
    avail,
    here::here("data-raw", "avail.txt"),
    format = "file"),
  tar_target(
    feat_dirs,
    readr::read_lines(avail, num_threads=1),
    format = "qs"),
  tar_target(n_sub, c(10, 20, 30, 50, 100, 250, 500, 1000)),
  tar_target(n_study, c(1)),
  tar_target(iter, seq_len(50)),
  tar_target(cope5_index, seq_len(10000)),
  tar_target(
    cope5,
    apply_reg_cope(feat_dirs[cope5_index], tar_path()),
    format = "file", error="continue",
    pattern = map(cope5_index)),
  tar_target(
    cope_files,
    prep_cope_tbl(cope5, n_sub=n_sub, iter=iter, n_study=n_study),
    pattern = cross(n_sub, iter, n_study),
    format = "fst_tbl"),
  tar_target(
    cope_files2,
    cope_files %>%
      dplyr::group_by(iter, n_study, n_sub, study) %>%
      tar_group(), 
    iteration = "group",
    format = "fst_tbl"),
  tar_target(
    t_img,
    do_t(cope_files2),
    pattern = map(cope_files2),
    format = "fst_tbl",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge,
          template = "tools/sge.tmpl",
          resources = list(mem_free = "15G"))))),
  tar_target(
    ale_py_script,
    fs::path("python", "ale.py"),
    format = "file"),
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
    calc_z(cope5),
    format = "qs",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge, 
          template = "tools/sge.tmpl", 
          resources = list(mem_free = "30G")))))
#  tar_target(
#    comparison,
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
