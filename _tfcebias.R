Sys.setenv(TAR_PROJECT = "bias")

library(targets)
library(tarchetypes)
library(rlang)

source(here::here("R", "ale.R"))
source(here::here("R", "spatial.R"))
source(here::here("R", "tfce.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "loading.R"))

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
  tar_target(n_sub, c(10)),
  tar_target(iter, seq_len(2)),
  tar_target(train_files, list(i=iter, f=sample(train, n_sub, replace=TRUE)), pattern = cross(n_sub, iter)),
  tar_target(
    tfce,
    do_tfcebias(train_files$f, n_sub=length(train_files$f), iter=train_files$i, n=1000, storage_dir=Sys.getenv("NIIDIR")),
    pattern = map(train_files)
  ),
  tar_target(
    tfce_pop,
    do_tfce_pop(test, n_sub=length(test), iter=0, storage_dir=Sys.getenv("NIIDIR"))),
  tar_target(tfce_cor, do_cor(tfce, tfce_pop))
)
