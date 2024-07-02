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

targets::tar_load(pop_cor2)
pop_cor2 |>
  mutate(`N Sub` = factor(n_sub)) |>
  ggplot(aes(x = rho, y = `N Sub`)) +
  facet_wrap(~Task) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(alpha = 0.05) +
  xlab("Product-Moment Correlation with Reference") +
  theme_gray(base_size = 9)

ggsave(
  "analyses/figures/topo_gold-to-study_rho.png",
  device = ragg::agg_png,
  width = 4,
  height = 4
)
