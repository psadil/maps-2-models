library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "pairwise")

targets::tar_load(rhos)

rhos |>
  filter(!is.na(rho), stringr::str_detect(task, "EMOTION", TRUE)) |>
  summarise(
    rho = mean(rho),
    .by = c(probe, task)
  ) |>
  ggplot(aes(x=rho, y=task)) +
  ggdist::stat_dots(quantiles = 50) +
  ylab("Task")  +
  xlab("Pairwise (Sub-Sub) Avgerage Correlation")

ggsave(
  "analyses/figures/topo_sub-to-sub_rho.png", 
  device = ragg::agg_png,
  width = 4,
  height = 4)


rhos |>
  filter(!is.na(rho), stringr::str_detect(task, "EMOTION", TRUE)) |>
  ggplot(aes(x=rho, y=task)) +
  ggdist::stat_dots(quantiles = 1000) +
  ylab("Task")  +
  xlab("Pairwise (Sub-Sub) Correlation")

ggsave(
  "analyses/figures/topo_sub-to-sub_rho-noavg.png", 
  device = ragg::agg_png,
  width = 4,
  height = 4)
