library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)


# What proportion of simulated studies have active voxels within the ROIs within 
# that are most active in the gold standard?
Sys.setenv(TAR_PROJECT = "hcp_ptfce")

gold_tested <- targets::tar_read(gold_tested) |>
  mutate(
    d = statistic / sqrt(n_sub),
    active = abs(d) > 0.2) |>
  mutate(r = row_number(desc(abs(estimate))), .by = c(Task, n_parcels))


gold_most <- gold_tested |>
  filter(r < 11) |>
  select(Task, n_parcels, label) |>
  group_nest(n_parcels) |>
  mutate(
    at = map(
      n_parcels, 
      ~make_atlas_full2(n_parcels=.x) |>
        distinct(index, label, `Full component name`)),
    data = map2(data, at, left_join)
  ) |>
  select(-at) |>
  unnest(data)

readr::write_csv(gold_most, "analyses/figures/gold_most.csv")


gold_most |>
  filter(n_parcels==400) |>
  select(-index, -n_parcels, -`Full component name`) |>
  print(n=100) 
