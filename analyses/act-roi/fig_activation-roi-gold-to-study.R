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

targets::tar_load(data_roi_study_to_gold)

data_roi_study_to_gold |>
  ggplot(aes(x=n_sub, y=prop, group=label, color=abs(d))) +
  geom_point(alpha=0.2) +
  geom_line(alpha=0.2) +
  facet_grid(n_parcels~Task) +
  scale_y_continuous(
    "Proportion Simulations w/ Activity in Most Active ROI",
    limits = c(0, 1)
  ) +
  scale_color_viridis_c(
    "Abs. Effect Size",
    option = "turbo") +
  xlab("N Sub") +
  theme_gray(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("analyses/figures/prop-activation-roi-gold-to-study.png", width = 10, height = 6.5)



data_roi_study_to_gold |>
  filter(n_parcels==400) |>
  mutate(expected = purrr::map2_dbl(n_sub, d, ~power.t.test(n=.x, delta=.y)$power)) |>
  ggplot(aes(x=expected, y=prop, group=label, color=abs(d))) +
  geom_point(alpha=0.2) +
  geom_line(alpha=0.2) +
  facet_wrap(~Task) +
  scale_y_continuous(
    "Proportion Simulations w/ Activity in Most Active ROI",
    limits = c(0, 1)
  ) +
  scale_color_viridis_c(
    "Abs. Effect Size",
    option = "turbo") +
  xlab("Expected Power") +
  geom_abline() +
  theme_gray(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("analyses/figures/prop-activation-roi-gold-to-study_expected-power.png", width = 5, height = 5)

# gold_tested <- targets::tar_read(gold_tested) |>
#   mutate(
#     d = statistic / sqrt(n_sub),
#     active = abs(d) > 0.2)


allact <- rois_tested |>
  left_join(
    select(gold_tested, Task, n_parcels, label, active),
    by = join_by(Task, n_parcels, label)) |>
  summarise(
    phi = phi(active.x, active.y),
    F1 = dice(active.x, active.y),
    .by = c(Task, n_parcels, iter, n_sub)
  ) |>
  mutate(
    n_sub = factor(
      n_sub,
      levels = c(20, 40, 60, 80, 100),
      labels = c(20, 40, 60, 80, 100)),
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}"))
  ) 


allact |>
  ggplot(aes(y=n_sub, x=phi, group=n_sub)) +
  facet_grid(n_parcels~Task) +
  ggdist::stat_dots(quantiles = 20) +
  # ggdist::stat_halfeye(size=1) +
  ylab("N Sub") +
  scale_x_continuous(
    limits = c(-0.2, 1),
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1))

allact |>
  ggplot(aes(y=n_sub, x=F1, group=n_sub)) +
  facet_grid(n_parcels~Task) +
  ggdist::stat_dots(quantiles = 20) +
  # ggdist::stat_halfeye(size=1) +
  ylab("N Sub") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1))

ggsave("analyses/figures/prop-activation-roi-gold-to-study-F1.png", width = 10, height = 6.5)


allact |>
  filter(fct_match(n_parcels, "N Parcels: 400")) |>
  ggplot(aes(y=n_sub, x=F1, group=n_sub)) +
  facet_wrap(~Task) +
  ggdist::stat_dots(quantiles = 20) +
  ylab("N Sub") +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1))

ggsave(
  "analyses/figures/prop-activation-roi-gold-to-study-F1400.png", 
  width = 3, 
  height = 3)
