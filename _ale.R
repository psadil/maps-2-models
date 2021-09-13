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
    # get_stats_dirs(here::here("data-raw","task_fMRI")),
    # get_stats_dirs(fs::path("/dcl01", "smart","data", "UKBiobank", "task_fMRI")),
    readr::read_lines(avail, num_threads=1),
    format = "qs"),
  tar_target(n_sub, c(10, 20, 50, 100)),
  tar_target(study, seq_len(5)),
  tar_target(iter, seq_len(1)),
  tar_target(cope5_index, seq_len(100)),
  tar_target(
    cope5,
    apply_reg_cope(feat_dirs[cope5_index], tar_path()),
    format = "file",
    pattern = map(cope5_index)),
  tar_target(
    cope_files,
    prep_cope_tbl(cope5, n = n_sub, study=study, iter=iter), 
    format = "fst_tbl",
    pattern = cross(n_sub, study, iter)),
  tar_target(
    clusters,
    calc_clusters(cope_files, pthresh = 0.001),
    format = "fst_tbl",
    pattern = map(cope_files),
    resources = tar_resources(
      future = tar_resources_future(
        resources = list(mem_free = 10)))),
  tar_target(
    clusters_grouped,
    clusters %>%
      dplyr::group_by(index, study, n_sub, iter) %>%
      dplyr::filter(Value == max(Value)) %>%
      dplyr::group_by(n_sub, iter) %>%
      tar_group(),
    iteration = "group",
    format = "qs"),
  tar_target(
    ale,
    do_ale(clusters_grouped, p=0.001),
    pattern = map(clusters_grouped),
    format = "file",
    iteration = "vector",
    resources = tar_resources(
      future = tar_resources_future(
        resources = list(mem_free = 10)))),
  tar_target(
    clusters_grouped_index,
    clusters %>%
      dplyr::group_by(n_sub, iter) %>%
      dplyr::mutate(
        index.bak = index,
        index = kmeans(cbind(x,y,z), dplyr::n_distinct(index), nstart = 10)$cluster) %>%
      dplyr::group_by(n_sub, iter, index) %>%
      tar_group(),
    iteration = "group",
    format = "qs"),
  tar_target(
    ale_index,
    do_ale(clusters_grouped_index, p=0.05, perm=1000, clust=0.05),
    pattern = map(clusters_grouped_index),
    format = "file",
    iteration = "vector",
    resources = tar_resources(
      future = tar_resources_future(
        resources = list(mem_free = 10)))),
  tar_target(
    z_pop,
    calc_z(cope5), 
    format = "qs")
)
