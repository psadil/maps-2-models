library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(c(data_roi_study_to_gold, data_roi_study_to_study, data_roi_sub_to_sub))

a <- data_roi_study_to_gold |>
  filter(n_parcels==400) |>
  ggplot(aes(x=n_sub, y=prop, group=label, color=abs(d))) +
  geom_point(alpha=0.2) +
  geom_line(alpha=0.2) +
  facet_wrap(~Task) +
  scale_y_continuous(
    "Proportion Simulations w/\nActivity in Most Active ROI",
    limits = c(0, 1),
    labels = c(0, 0.5, 1),
    breaks = c(0, 0.5, 1)
  ) +
  scale_x_continuous(
    "N Sub",
    breaks = c(40, 80),
    labels = c(40, 80)
  ) +
  scale_color_viridis_c(
    "Abs. Effect Size",
    option = "turbo",
    limits = c(0, NA),
    n.breaks = 3) 

b <- data_roi_study_to_study |>
  filter(fct_match(n_parcels, "N Parcels: 400")) |>
  ggplot(aes(y=n_sub, x=phi, group=n_sub)) +
  facet_wrap(~Task) +
  ggdist::stat_dots(quantiles = 50) +
  ylab("N Sub") +
  scale_x_continuous(
    "Phi Coefficient\n(Study-Study)",
    limits = c(-0.25, 1),
    breaks = c(-0.25, 0.5),
    labels = c(-0.25, 0.5))

cc <- data_roi_sub_to_sub |>
  mutate(
    Task = factor(Task),
    Task = forcats::fct_rev(Task)) |>
  filter(fct_match(n_parcels, 400)) |>
  ggplot(aes(y=Task, x=rho)) +
  ggdist::stat_dotsinterval(
    quantiles = 100) +
  scale_x_continuous(
    "Product-Moment Correlation\n(Sub-Sub)",
    limits = c(-0.75, 1),
    breaks = c(-0.75, 0, 0.75),
    labels = c(-0.75, 0, 0.75))

p <- a + b + cc +
  plot_layout(ncol = 1, heights = c(1, 1.5, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme_gray(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(8, "pt"))

ggsave(
  "analyses/figures/roi.png", 
  p,
  width = 3, 
  height = 6,
  device = ragg::agg_png)
