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

targets::tar_load(data_roi_study_to_study)

data_roi_study_to_study |>
  ggplot(aes(y=n_sub, x=phi, group=n_sub)) +
  facet_grid(n_parcels~Task) +
  ggdist::stat_dots(quantiles = 50) +
  ylab("N Sub") +
  scale_x_continuous(
    "Phi Coefficient",
    limits = c(-0.5, 1),
    breaks = c(-0.5, 0, 0.5, 1),
    labels = c(-0.5, 0, 0.5, 1))

ggsave("analyses/figures/prop-activation-roi-study-to-study_phi.png", width = 10, height = 6.5)

