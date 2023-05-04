library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)
library(pTFCE)
library(oro.nifti)

source(here::here("R", "ale.R"))
source(here::here("R", "spatial.R"))
source(here::here("R", "tfce.R"))
source(here::here("R", "utils.R"))
source(here::here("R", "loading.R"))
source(here::here("R", "poster.R"))
source(here::here("R", "hcp.R"))
source(here::here("R", "ptfce.R"))


# What proportion of simulated studies have active voxels within the ROIs within 
# that are most active in the gold standard?

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

at <- targets::tar_read(at) |>
  mutate(across(c(x,y,z), as.integer))

at_list <- targets::tar_read(at_list) |>
  mutate(across(c(x,y,z), as.integer))

out <- targets::tar_read(active) |>
  dplyr::left_join(at_list) |>
  distinct(
    n_sub, iter, Task, label, `Label Name`,
    `Full component name`, n_parcels) |>
  collect() |>
  count(
    n_sub, Task, label, `Label Name`,
    `Full component name`, n_parcels) |>
  mutate(prop = n / 100)

get_max <- function(q){
  x <- qs::qread(q)
  to_tbl0(x$Z, measure = "Z") |> mask()
}

# these are the regions whose activation is the "strongest"
# so, assumption is that researcher picks the task to activate
# this specific region
gold <- targets::tar_read(tfce_pop) |>
  select(-avail) |>
  dplyr::mutate(tmp = purrr::map(ptfce, get_max)) |>
  select(-ptfce) |>
  tidyr::unnest(tmp) |>
  left_join(at_list, multiple = "all") |>
  group_by(
    Task, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  summarise(
    Z = mean(Z),
    .groups = "drop"
  ) |>
  group_by(Task, n_parcels) |>
  slice_max(order_by=Z, n=20) |>
  ungroup() |>
  filter(!is.na(label)) |> 
  group_by(Task) |>
  mutate(
    l=label |> 
      factor() |> 
      as.numeric() |>
      factor()) |>
  ungroup()


# regions with at least one voxel active
out |>
  semi_join(distinct(gold, Task, l, label, n_parcels)) |>
  right_join(
    distinct(gold, Task, l, label, n_parcels) |> 
      crossing(distinct(out, n_sub))) |>
  filter(n_parcels %in% c(200, 400, 600, 800, 1000)) |>
  mutate(
    propr = if_else(is.na(prop), 0, prop),
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}"))
  ) |>
  ggplot(aes(x=n_sub, y=prop, group=l)) +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_grid(n_parcels~Task) +
  scale_y_continuous(
    "Proportion Simulations w/ Activity in Most Active ROI",
    limits = c(0, 1)
  ) +
  xlab("N Sub") +
  theme_gray(base_size = 12)

ggsave("analyses/figures/prop-active-most-active-roi-ptfce.png", width = 10, height = 6.5)

## NULL

# regions with at least one voxel active
out_null <- targets::tar_read(active_null) |>
  dplyr::left_join(at_list) |>
  distinct(
    n_sub, iter, label, `Label Name`,
    `Full component name`, n_parcels) |>
  collect() |>
  count(
    n_sub, label, `Label Name`,
    `Full component name`, n_parcels) |>
  mutate(prop = n / 100)

gold_null <- targets::tar_read(tfce_pop_null) |>
  dplyr::mutate(tmp = purrr::map(ptfce, get_max)) |>
  select(-ptfce) |>
  tidyr::unnest(tmp) |>
  left_join(at_list, multiple = "all") |>
  group_by(
    label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  summarise(
    Z = mean(Z),
    .groups = "drop"
  ) |>
  group_by(n_parcels) |>
  slice_max(order_by=Z, n=20) |>
  ungroup() |>
  filter(!is.na(label)) |> 
  mutate(
    l=label |> 
      factor() |> 
      as.numeric() |>
      factor()) |>
  ungroup()

out_null |>
  semi_join(distinct(gold_null, l, label, n_parcels)) |>
  right_join(
    distinct(gold_null, label, n_parcels, l) |> 
      crossing(distinct(out_null, n_sub))) |>
  filter(n_parcels %in% c(200, 400, 600, 800, 1000)) |>
  mutate(
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}")),
    propr = if_else(is.na(prop), 0, prop)
  ) |>
  mutate(prop = dplyr::if_else(is.na(prop), 0, prop)) |>
  ggplot(aes(x=n_sub, y=prop, group=l)) +
  geom_ribbon(
    aes(ymin= 0.0073, ymax = 0.0926), 
    alpha=0.01, 
    fill="lightblue",
    linetype = "dashed",
    color = "black") +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_wrap(~n_parcels) +
  scale_y_continuous(
    "Proportion Simulations w/ Activity in Most Active ROI",
    limits = c(0, 0.1)
  ) +
  xlab("N Sub") +
  theme_gray(base_size = 12)

ggsave("analyses/figures/prop-active-most-active-roi_ptfce-null.png", width = 10, height = 6.5)
