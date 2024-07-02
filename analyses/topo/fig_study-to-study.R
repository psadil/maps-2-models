library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

source(here::here("R", "utils.R"))
source(here::here("R", "tfce.R"))
source(here::here("R", "updates.R"))
source(here::here("R", "poster.R"))
source(here::here("R", "hcp.R"))

Sys.setenv(TAR_PROJECT = "hcp")


targets::tar_load(data_topo_study_to_study)

data_topo_study_to_study |>
  ggplot(aes(x = rr, y = `N Sub`, color = Task)) +
  geom_point() +
  geom_line(aes(group = Task)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper)) +
  xlab("Pairwise Product-Moment Correlation") +
  theme_gray(base_size = 9) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  "analyses/figures/topo_study-to-study_rho.png",
  device = ragg::agg_png,
  width = 4,
  height = 4
)
