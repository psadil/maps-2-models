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

sim.rho <- function(X, ...){
  n <- nrow(X)
  D <- array(0, dim=c(n, n))
  for (i in 1:(n - 1)) {
    for (j in i:n) {
      D[i,j] <- cor(X[i,], X[j,])
    }
  }
  D <- D + t(D)
  return(D)
}

rois_pop <- targets::tar_read(rois_pop2) |>
  mutate(d = Z / sd) |>
  select(Task, n_parcels, d, label, sub) |>
  collect()

out <- rois_pop |>
  collect() |>
  group_nest(Task, n_parcels) |>
  mutate(
    data = map(
      data, 
      ~.x |> 
        select(label, d, sub) |> 
        pivot_wider(names_from = label, values_from=d) |> 
        arrange(sub)),
    rho = map(
      data,
      ~select(.x, -sub) |> 
        as.matrix() |> 
        sim.rho() |>
        corrr::as_cordf() |> 
        corrr::shave() |> 
        corrr::stretch(na.rm = TRUE) |>
        rename(rho = r))
  ) |>
  mutate(n_parcels = factor(n_parcels)) |>
  select(-data) |>
  unnest(rho)

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

