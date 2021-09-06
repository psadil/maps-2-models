# Sys.setenv(TAR_PROJECT = "ale")

library(targets)

source(here::here("R", "loading.R"))
source(here::here("R", "utils.R"))

# End this file with a list of target objects.
list(
  tar_target(
    key_file,
    here::here("data-raw","Data_Dictionary_Showcase.csv"),
    format = "file"),
  tar_target(
    key,
    load_key(key_file),
    format = "fst_tbl"),
  tar_target(
    lit_file,
    here::here("data-raw","lit-review.xlsx"),
    format = "file"),
  tar_target(
    lit,
    load_lit(lit_file),
    format = "fst_tbl"),
  tar_target(
    fname,
    here::here("data-raw", "ukb37917.tab"),
    format = "file"),
  tar_target(
    ukb, 
    load_tab(fname, key, n_max=2e5),
    format = "fst_tbl"),
  tar_target(
    n,
    c(10, 20, 40, 100, 1000, 3000)),
  tar_target(i, seq_len(1000)),
  tar_target(
    ukb_split,
    sample_ukb(ukb, n=n, i=i),
    pattern = cross(n, i),
    format = "fst_tbl"),
  tar_target(
    cors,
    correlate(ukb_split, `20016`),
    format = "fst_tbl"),
  tar_target(
    out,
    include_keys(cors, key),
    format = "fst_tbl"),
  tar_target(
    rough,
    plot_rough2(out, lit),
    format = "qs"
  )
)

