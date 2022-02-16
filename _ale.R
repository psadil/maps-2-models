# Sys.setenv(TAR_PROJECT = "ale")

library(targets)
library(tarchetypes)
library(rlang)
library(magrittr)

source(here::here("R", "ale.R"))

library(future)
# library(future.callr)
# plan(callr)
library(future.batchtools)
plan(batchtools_sge, template = "tools/sge.tmpl")

get_sizes <- function(cls){
  sizes <- neurobase::img_indices(cls$osize, add_values = TRUE) |>
    tibble::as_tibble() |>
    dplyr::rename(cluster_size = value)
  neurobase::img_indices(cls$oindex, add_values = TRUE) |>
    tibble::as_tibble() |>
    dplyr::filter(value > 0) |>
    dplyr::rename(`Cluster Index` = value) |>
    dplyr::left_join(sizes, by = c("x","y","z")) |>
    dplyr::distinct(`Cluster Index`, cluster_size)
  
}

get_maxes <- function(
  niifile, 
  threshold=2, 
  mask=fslr::mni_fname(mm="2", brain = TRUE, mask = TRUE), 
  at=atlas,
  solitary = TRUE,
  minextent = 10){
  m <- neurobase::readnii(mask)
  z <- neurobase::readnii(niifile) |>
    neurobase::mask_img(mask=m) 
  
  cls1 <- fslr::fslcluster(
    z, 
    threshold=threshold,
    opts = glue::glue("--volume={sum(m)} --minextent={minextent}"))
  
  cls2 <- fslr::fslcluster(
    z * -1, 
    threshold=threshold,
    opts = glue::glue("--volume={sum(m)} --minextent={minextent}"))
  
  pos <- fslr::read_cluster_table(cls1$olmax) |>
    tibble::as_tibble() |>
    dplyr::mutate(sign = "positive") |>
    dplyr::left_join(get_sizes(cls1), by = "Cluster Index")
  neg <- fslr::read_cluster_table(cls2$olmax) |>
    tibble::as_tibble() |>
    dplyr::mutate(sign = "negative") |>
    dplyr::left_join(get_sizes(cls2), by = "Cluster Index")
  
  out <- dplyr::bind_rows(pos, neg) |>
    dplyr::left_join(at, by = c("x","y","z"))
  if (solitary){
    out <- out |>
      dplyr::group_by(label) |>
      dplyr::filter(Value == max(Value)) |>
      dplyr::ungroup() 
  }
  return(out)
}


augment_distance <- function(study, reference){
  neg <- fuzzyjoin::distance_left_join(
    study |> dplyr::filter(sign == "negative"),
    reference |> dplyr::filter(sign == "negative") |> dplyr::select(-`Cluster Index`, -Value, -cluster_size, -n_voxels, -volume, -sign),
    distance_col = "d",
    max_dist = 100,
    by = c("x", "y", "z")) 
  pos <- fuzzyjoin::distance_left_join(
    study |> dplyr::filter(sign == "positive"),
    reference |> dplyr::filter(sign == "positive") |> dplyr::select(-`Cluster Index`, -Value, -cluster_size, -n_voxels, -volume, -sign),
    distance_col = "d",
    max_dist = 100,
    by = c("x", "y", "z")) 
  
  out <- dplyr::bind_rows(neg, pos) |>
    dplyr::select(-x.y, -y.y, -z.y, x=x.x, y=y.x, z=z.x, label=label.x, ref_label=label.y)
  
  return(out)
}

list(
  tar_target(
    avail,
    here::here("data-raw", "avail.txt"),
    format = "file"),
  tar_target(
    feat_dirs,
    readr::read_lines(avail, num_threads=1),
    format = "qs"),
  tar_target(n_sub, c(10, 20, 50, 100, 250, 500, 1000)),
  tar_target(n_study, c(1)),
  tar_target(iter, seq_len(50)),
  tar_target(
    cope5,
    apply_reg_cope(feat_dirs, tar_path()),
    format = "file", 
    pattern = map(feat_dirs),
    storage = "worker",
    retrieval = "worker", error="continue",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge,
          template = "tools/sge.tmpl",
          resources = list(mem_free = "512M"))))),
  tar_target(
    cope_files,
    prep_cope_tbl(cope5, n_sub=n_sub, iter=iter, n_study=n_study),
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
    calc_z(cope5),
    format = "qs",
    resources = tar_resources(
      future = tar_resources_future(
        plan = tweak(
          batchtools_sge, 
          template = "tools/sge.tmpl", 
          resources = list(mem_free = "150G")))))
  #tar_target(
  # at,
  # fst::read_fst(here::here("data-raw", "atlas.fst")),
  # format = "fst_tbl")
  # tar_target(
  #   maxes,
  #   z_img |>
  #     dplyr::mutate(
  #       maxes = purrr::map(z_stat, ~get_maxes(.x, solitary = FALSE, threshold=2.3, minextent=0, at = at))),
  #   pattern = map(z_img),
  #   format = "qs")
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
