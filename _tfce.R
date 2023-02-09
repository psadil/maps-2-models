Sys.setenv(TAR_PROJECT = "tfce")

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
  NIIDIR = here::here("data-raw","niis"),
  AVAIL = here::here("data-raw", "copes")) # explicitly avoiding tracking this

#options(clustermq.scheduler = "multiprocess")

targets::tar_option_set(format = "qs", storage="worker", retrieval="worker")
library(future.callr)
plan(callr)


list(
  tar_target(
    train,
    readr::read_lines(Sys.getenv("AVAIL"), n_max=10000)
  ),
  tar_target(
    test,
    readr::read_lines(Sys.getenv("AVAIL"), skip=10000)),
  tar_target(n_sub, c(10, 20, 50, 100, 200, 500)),
  tar_target(iter, seq_len(100)),
  tar_target(
    tfce,
    do_tfce(train, n_sub=n_sub, iter=iter, n=1000, storage_dir=Sys.getenv("NIIDIR")),
    pattern = cross(n_sub, iter)
  ),
  tar_target(
    tfce_pop,
    do_tfce_pop(test, n_sub=length(test), iter=0, storage_dir=Sys.getenv("NIIDIR"))),
  tar_target(tfce_cor, do_cor(tfce, tfce_pop)),
  tar_target(corrp_thresh, c(0.01, 0.95)),
  tar_target(
    maxes,
    dplyr::mutate(
      tfce,
      m = purrr::map2(
        .data$tfce_corrp_tstat, stringr::str_replace(.data$tfce_tstat,"/_tstat", "_tstat"),
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
    maxes |>
      dplyr::mutate(
        augmented = purrr::map(
          .data$m,
          ~augment_distance(
            study = .x,
            reference = gold_peaks$m[[1]],
            vox_mm = 2.4
          )
        )
      ) |>
      dplyr::select(n_sub, augmented, iter, corrp_thresh) |>
      tidyr::unnest(augmented)
  ),
  tar_target(at, make_atlas_full()),
  tar_target(space, add_labels(space=space0, at=at)),
  tar_target(cor2, get_cor2(tfce_cor, tfce_pop)),
  tar_target(center, get_center(tfce_cor, tfce_pop)),
  tar_target(sigmacope_bias, make_sigmacope_bias_tbl(tfce_cor, tfce_pop)),
  tar_target(pop_d, make_pop_d(tfce_pop)),
  tar_target(iters, make_iters(tfce_cor, pop_d)),
  tar_target(big, make_big(tfce_cor, pop_d, prop=1), format = "parquet"),
  tar_target(
    key_file,
    here::here("data-raw", "Data_Dictionary_Showcase.csv"),
    format = "file"),
  tar_target(key, load_key(key_file)),
  tar_target(lit_file, here::here("data-raw","lit-review.xlsx"), format = "file"),
  tar_target(lit, load_lit(lit_file)),
  tar_target(fname, here::here("data-raw", "ukb.parquet"), format = "file"),
  tar_target(ukb, load_ukb_parquet_for_corsplit(fname)),
  tar_target(n_for_cor, c(10, 20, 50, 100, 200, 500)),
  tar_target(i_for_cor, seq_len(1000)),
  tar_target(
    splitcors, 
    split_and_correlate(ukb, n=n_for_cor, i=i_for_cor, y=fluid), 
    pattern = cross(n_for_cor, i_for_cor),
    deployment = "main"),
  tar_target(brainvolfluid, include_keys(splitcors, key)),
  tar_target(sample_study, tfce_cor$d[[310]] |> dplyr::left_join(pop_d)),
  tar_target(ukb_cor, get_gold_bvf_cor(ukb, test)),
  tar_target(sampling_distr, get_cor_sampling_distr(10, 500, rho=ukb_cor)),
  tarchetypes::tar_render(poster, "analyses/ohbm2022/poster3.rmd")
)
