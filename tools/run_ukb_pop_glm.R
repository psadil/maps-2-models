
source(here::here("R", "spatial.R"))
source(here::here("R", "tfce.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "hcp.R"))
source(here::here("R", "ptfce.R"))
source(here::here("R", "figures.R"))
source(here::here("R", "model.R"))
source(here::here("R", "manuscript.R"))
source(here::here("R", "roi.R"))
source(here::here("R", "topo.R"))
source(here::here("R", "ukb.R"))
source(here::here("R", "cifti.R"))


Sys.setenv(DUCKPLYR_TEMP_DIR = "/fastscratch/myscratch/pssadil")

ex <- to_tbl("data-raw/ukb_mni/1000043.nii.gz") |>
  mask() |>
  mask_gray() |>
  mask_atlas() |>
  dplyr::select(-value) |>
  dplyr::mutate(index = dplyr::row_number()) |>
  duckplyr::as_duckdb_tibble()

# files <- fs::dir_ls("_targets/objects", glob = "*ukb_gray_*")
# arrow::open_dataset(files) |>
#   dplyr::right_join(ex, by = dplyr::join_by(x, y, z)) |>
#   dplyr::summarise(
#     pe = mean(value), 
#     sigma = sd(value), 
#     n_sub = dplyr::n_distinct(sub),
#     .by = c(index)) |> 
#   dplyr::collect() |>
#   dplyr::mutate(
#     z = (sqrt(n_sub) * pe / sigma) |> # t_stat -> z
#       pt(n_sub - 1, log.p = TRUE, lower.tail = FALSE) |>
#       qnorm(log.p = TRUE, lower.tail = FALSE)
#   ) |>
#   arrow::write_parquet("data-raw/ukb_glm_pop.parquet")

duckplyr::db_exec("PRAGMA max_temp_directory_size = '512GiB'")


files <- fs::dir_ls("_targets/objects", glob = "*ukb_gray_*")
duckplyr::read_parquet_duckdb(files) |>
  dplyr::summarise(
    pe = mean(value), 
    sigma = sd(value), 
    n_sub = dplyr::n_distinct(sub),
    .by = c(x, y, z)) |> 
  dplyr::right_join(ex, by = dplyr::join_by(x, y, z)) |>
  dplyr::select(-x, -y, -z) |>
  dplyr::collect() |>
  dplyr::mutate(
    z = (sqrt(n_sub) * pe / sigma) |> # t_stat -> z
      pt(n_sub - 1, log.p = TRUE, lower.tail = FALSE) |>
      qnorm(log.p = TRUE, lower.tail = FALSE)
  ) |>
  arrow::write_parquet("data-raw/ukb_glm_pop.parquet")



