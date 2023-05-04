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

phi <- function(x, y){
  tp <- mean(x & y)
  fp <- mean(!x & y)
  tn <- mean(!x & !y)
  fn <- mean(x & !y)
  if(((tp&fp) == 0) | ((tp & fn) == 0) | ((fn & tn) == 0) | ((fp & tn) == 0)){
    return(0)
  }
  cor(x, y)
}

sim.phi <- function(X, ...){
  n <- nrow(X)
  D <- array(0, dim=c(n, n))
  for (i in 1:(n - 1)) {
    for (j in i:n) {
      D[i,j] <- phi(X[i,], X[j,])
    }
  }
  D <- D + t(D)
  return(D)
}

# What proportion of simulated studies have active voxels within the ROIs within 
# that are most active in the gold standard?
Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(rois_tested2)

out <- rois_tested2 |>
  select(Task, n_parcels, n_sub, iter, active, label) |>
  collect() |>
  group_nest(Task, n_parcels, n_sub) |>
  mutate(
    data = map(
      data, 
      ~.x |> 
        select(label, active, iter) |> 
        pivot_wider(names_from = label, values_from=active) |> 
        arrange(iter)),
    phi = map(
      data,
      ~select(.x, -iter) |> 
        as.matrix() |> 
        sim.phi() |>
        corrr::as_cordf() |> 
        corrr::shave() |> 
        corrr::stretch(na.rm = TRUE) |>
        rename(phi = r))
  ) |>
  mutate(
    n_sub = factor(
      n_sub,
      levels = c(20, 40, 60, 80, 100),
      labels = c(20, 40, 60, 80, 100)),
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}"))
  ) |>
  select(-data) |>
  unnest(phi)


out |>
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


out |>
  filter(fct_match(n_parcels, "N Parcels: 400")) |>
  ggplot(aes(y=n_sub, x=phi, group=n_sub)) +
  facet_wrap(~Task) +
  ggdist::stat_dots(quantiles = 50) +
  ylab("N Sub") +
  scale_x_continuous(
    "Phi Coefficient",
    limits = c(-0.5, 1),
    breaks = c(-0.5, 0, 0.5, 1),
    labels = c(-0.5, 0, 0.5, 1))

ggsave("analyses/figures/prop-activation-roi-study-to-study_phi400.png", width = 3, height = 3)
