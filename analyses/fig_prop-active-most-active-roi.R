library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

# What proportion of simulated studies have active voxels within the ROIs within 
# that are most active in the gold standard?

Sys.setenv(TAR_PROJECT = "hcp")

at <- targets::tar_read(at) |>
  mutate(across(c(x,y,z), as.integer))

at_list <- targets::tar_read(at_list) |>
  mutate(across(c(x,y,z), as.integer))

out <- targets::tar_read(prop0all) |>
  dplyr::left_join(at_list) |>
  distinct(
    n_sub, iter, Task, ContrastName, CopeNumber, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  collect() |>
  count(
    n_sub, Task, ContrastName, CopeNumber, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  mutate(prop = n / 100)

# these are the regions whose activation is the "strongest"
# so, assumption is that researcher picks the task to activate
# this specific region
gold <- targets::tar_read(pop_d) |>
  left_join(at_list) |>
  group_by(
    Task, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  summarise(
    d = mean(cope/sigma),
    .groups = "drop"
  ) |>
  group_by(Task, n_parcels) |>
  slice_max(order_by=d, n=10) |>
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
  right_join(distinct(gold, Task, l, label, n_parcels)) |>
  filter(n_parcels %in% c(200, 400, 600, 800, 1000)) |>
  mutate(
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}"))
  ) |>
  ggplot(aes(x=n_sub, y=prop, group=l)) +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_grid(n_parcels~Task) +
  scale_y_continuous(
    "Proportion Simulations w/ Activity in Most Active ROI",
    limits = c(0.25, 1)
  ) +
  xlab("N Sub") +
  theme_gray(base_size = 12)

ggsave("analyses/figures/prop-active-most-active-roi.png", width = 10, height = 6.5)
