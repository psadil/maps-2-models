library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(c(data_peak_study_to_gold, data_peak_study_to_study, data_peak_sub_to_sub))

a <- data_peak_study_to_gold |>
  ggplot(aes(x=within, group=peak, y=n_simulations, color=Value)) +
  geom_point(alpha=0.2) +
  geom_line(alpha=0.2) +
  facet_grid(n_sub~Task) +
  scale_y_continuous(
    "Proportion Simulations w/\nPeak in Radius",
    limits = c(0,1),
    breaks = c(0,0.5,1),
    labels = c(0,0.5,1)
  ) +
  scale_x_continuous(
    "Radius (mm)",
    limits = c(0, 20),
    breaks = c(0, 10, 20),
    labels = c(0, 10, 20)
  ) +
  scale_color_viridis_c(
    "Peak Height",
    option = "turbo",
    limits = c(0, NA),
    n.breaks = 3) 

b <- data_peak_study_to_study |>
  ggplot(aes(y=n_sub, x=d)) +
  facet_wrap(~Task) +
  ggdist::stat_dots(quantiles = 100) +
  ylab("N Sub") +
  scale_x_continuous(
    "Distance Between\nHighest Peaks\n(Study-Study)") 

cc <- data_peak_sub_to_sub |>
  mutate(
    Task = factor(Task),
    Task = forcats::fct_rev(Task)) |>
  ggplot(aes(y=Task, x=d)) +
  ggdist::stat_dots(quantiles = 100) +
  ylab(NULL) +
  scale_x_continuous(
    "Distance Between\nHighest Peaks\n(Sub-Sub)") 


# p <- a + b + cc +
#   plot_layout(design = c(
#     area(1, 1, 6, 3),
#     area(1, 3, 5, 5),
#     area(5, 3, 6, 5)
#   )) &
#   theme_gray(base_size = 8) +
#   theme(
#     legend.position = "bottom",
#     legend.key.size = unit(8, "pt"))

p <- a + b + cc +
  plot_layout(nrow = 1, widths = c(4, 2, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme_gray(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(8, "pt"))

ggsave(
  "analyses/figures/peaks.png", 
  p,
  width = 9, 
  height = 4,
  device = ragg::agg_png)
