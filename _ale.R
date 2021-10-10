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
  tar_target(n_sub, c(5, 10, 20)),
  tar_target(n_study, c(5, 10, 20, 30)),
  tar_target(iter, seq_len(5)),
  tar_target(cope5_index, seq_len(600)),
  tar_target(
    cope5,
    apply_reg_cope(feat_dirs[cope5_index], tar_path()),
    format = "file",
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
    clusters,
    calc_clusters(cope_files2, pthresh = 0.01),
    cue = tar_cue(depend = FALSE),
    pattern = map(cope_files2),
    format = "fst_tbl",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge, 
          template = "tools/sge.tmpl", 
          resources = list(mem_free = "10G"))))),
  tar_target(
    clusters2,
    clusters %>% 
      dplyr::group_by(index, study, n_sub, iter, n_study) %>%
      dplyr::filter(Value == max(Value)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(n_sub, iter, n_study) %>%
      tar_group(),
    iteration = "group",
    format = "fst_tbl"),
  tar_target(
    ale,
    do_ale(clusters2, p=0.05, perm=1000, clust=0.05),
    pattern = map(clusters2),
    cue = tar_cue(depend = FALSE),
    format = "file",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge, 
          template = "tools/sge.tmpl", 
          resources = list(mem_free = "10G"))))),
  tar_target(
    z_pop,
    calc_z(cope5),
    format = "qs",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge, 
          template = "tools/sge.tmpl", 
          resources = list(mem_free = "20G"))))),
  tar_target(
    comparison,
    avg_by_clust(ale, z_pop),
    pattern = map(ale),
    format = "fst_tbl"
  )
)
