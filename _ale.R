# Sys.setenv(TAR_PROJECT = "ale")

library(targets)
library(rlang)
library(magrittr)

source(here::here("R", "ale.R"))

list(
  tar_target(
    stats_dirs,
    # get_stats_dirs(here::here("data-raw","task_fMRI")),
    get_stats_dirs(fs::path("/dcl01", "smart","data", "UKBiobank", "task_fMRI")),
    format = "file"),
  tar_target(
    n_sub,
    c(4, 6)),
  tar_target(
    study,
    !!1:2
  ),
  # tar_target(i, seq_len(10)),
  tar_target(cope5_index, !!1:6),
  tar_target(
    cope5,
    apply_reg_cope(stats_dirs[cope5_index], tar_path()),
    format = "file",
    pattern = map(cope5_index)
  ),
  tar_target(
    cope_files,
    prep_cope_tbl(cope5, n = n_sub, study=study),
    format = "fst_tbl",
    pattern = cross(n_sub, study)
  ),
  tar_target(
    clusters,
    calc_clusters(cope_files, pthresh = 0.5),
    format = "fst_tbl",
    pattern = map(cope_files)),
  tar_target(
    clusters_grouped,
    clusters |>
      dplyr::group_by(n_sub) |>
      tar_group(),
    iteration = "group"),
  tar_target(
    ale,
    do_ale(clusters_grouped, foci_stem=tar_path(), p=0.001),
    pattern = map(clusters_grouped),
    format = "file",
    iteration = "vector"
  )
)
