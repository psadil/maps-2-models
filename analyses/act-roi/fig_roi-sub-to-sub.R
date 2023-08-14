library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

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


out |>
  ggplot(aes(y=n_parcels, x=rho)) +
  facet_wrap(~Task) +
  ggdist::stat_dots(quantiles = 100) +
  ylab("N Parcels") +
  scale_x_continuous(
    "Product-Moment Correlation",
    limits = c(-0.75, 1),
    breaks = c(-0.75, 0, 0.75),
    labels = c(-0.75, 0, 0.75)) +
  theme_gray(base_size = 9)

ggsave(
  "analyses/figures/prop-activation-roi-sub-to-sub_rho.png", 
  width = 3, 
  height = 3)

