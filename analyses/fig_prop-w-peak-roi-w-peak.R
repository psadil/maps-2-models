#' What proportion of simulated studies have active voxels within the ROIs 
#' whose average effect size in the gold standard is the highest?

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp")

targets::tar_load(space)

at <- targets::tar_read(at) |>
  mutate(across(c(x,y,z), as.integer))

at_list <- targets::tar_read(at_list) |>
  mutate(across(c(x,y,z), as.integer))


gold2 <- targets::tar_read(gold_peaks) |>
  select(Task, m) |>
  unnest(m) |>
  left_join(at_list) |>
  filter(!is.na(label)) |> 
  group_by(Task, label, n_parcels, n_networks) |>
  slice_max(order_by = Value, n = 1) |>
  group_by(Task, n_parcels, n_networks) |>
  slice_max(order_by = Value, n = 10, with_ties = FALSE) |>
  ungroup() |>
  group_by(Task) |>
  mutate(
    l=label |> 
      factor() |> 
      as.numeric() |>
      factor()) |>
  ungroup()


space |>
  filter(corrp_thresh==0.95) |>
  select(Task, iter, n_sub, x=x.study, y=y.study, z=z.study) |>
  left_join(select(at_list, -`Label Name`, -hemi, -n_voxels, -volume)) |>
  semi_join(distinct(gold2, Task, label, n_parcels)) |>
  distinct(Task, iter, label, n_sub, n_parcels) |>
  right_join(distinct(gold2, Task, label, n_parcels) |> crossing(distinct(space, n_sub))) |>
  group_by(Task, label, n_sub, n_parcels) |>
  summarise(
    prop = sum(!is.na(iter)) / 100,
    .groups = "drop"
  ) |>
  filter(n_parcels %in% c(200, 400, 600, 800, 1000)) |>
  mutate(
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}"))
  ) |>
  ggplot(aes(x=n_sub, y=prop, group=label)) +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_grid(n_parcels~Task) +
  scale_y_continuous(
    "Proportion Simulations with Peak in ROI with Peak",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1)
  ) +
  xlab("N Sub") +
  theme_gray(base_size = 12)

ggsave("analyses/figures/prop-w-peak-roi-w-peak.png", width = 10, height = 6.5)
