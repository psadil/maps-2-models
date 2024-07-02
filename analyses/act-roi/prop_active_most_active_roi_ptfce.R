library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(data_roi_study_to_gold)

p <- data_roi_study_to_gold |>
  ggplot(aes(x = n_sub, y = prop, group = label, color = abs(d))) +
  geom_point(alpha = 0.2) +
  geom_line(alpha = 0.2) +
  facet_grid(n_parcels~Task) +
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
    n.breaks = 3
  ) +
  theme(legend.position = "bottom")

tikzDevice::tikz(
  "analyses/figures/prop-active-most-active-roi-ptfce.tex",
  width=6,
  height=6.5)
p
dev.off()

