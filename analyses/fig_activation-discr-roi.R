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

jaccard <- function(x, y) {
  denom <- sum(c(sum(x), sum(y), -sum(x & y)))
  if (denom == 0){
    return(-Inf)
  }
  sum(x & y) / denom
}

dice <- function(x, y){
  2 * sum(x & y) / (sum(x) + sum(y))
}

dist.jaccard <- function(X, ...){
  n <- nrow(X)
  D <- array(0, dim=c(n, n))
  for (i in 1:(n - 1)) {
    for (j in i:n) {
      D[i,j] <- 1 - jaccard(X[i,], X[j,])
    }
  }
  D <- D + t(D)
  return(D)
}

sim.dice <- function(X, ...){
  n <- nrow(X)
  D <- array(0, dim=c(n, n))
  for (i in 1:(n - 1)) {
    for (j in i:n) {
      D[i,j] <- dice(X[i,], X[j,])
    }
  }
  D <- D + t(D)
  return(D)
}

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

at <- targets::tar_read(at) |>
  mutate(across(c(x,y,z), as.integer))

at_list <- targets::tar_read(at_list) |>
  mutate(across(c(x,y,z), as.integer))


active <- targets::tar_read(active) |>
  dplyr::left_join(at_list) |>
  filter(!is.na(label)) |>
  distinct(
    n_sub, iter, Task, label, `Label Name`,
    `Full component name`, n_parcels) |>
  collect()

notactive <- at_list |>
  distinct(
    label, `Label Name`, 
    `Full component name`, n_parcels) |>
  crossing(n_sub = unique(active$n_sub), iter = unique(active$iter), Task = unique(active$Task)) |>
  anti_join(
    active, 
    by = join_by(label, `Label Name`, `Full component name`, n_parcels, n_sub, iter, Task))


allact <- bind_rows(list(`TRUE`=active, `FALSE`=notactive), .id = "above") |>
  mutate(above = as.logical(above)) |>
  filter((n_parcels / 100) %% 2 == 0) |>
  group_nest(n_sub, n_parcels, Task) |>
  mutate(
    data = map(
      data, 
      ~.x |> 
        select(label, above, iter) |> 
        pivot_wider(names_from = label, values_from=above) |> 
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
  ) 


allact |>
  mutate(
    n_sub = factor(
      n_sub,
      levels = c(20, 40, 60, 80, 100),
      labels = c(20, 40, 60, 80, 100)),
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}"))
  ) |>
  select(-data) |>
  unnest(phi) |>
  ggplot(aes(y=n_sub, x=phi, group=n_sub)) +
  facet_grid(n_parcels~Task) +
  ggdist::stat_dots(quantiles = 50) +
  # ggdist::stat_halfeye(size=1) +
  ylab("N Sub") +
  scale_x_continuous(
    limits = c(-0.2, 1),
    breaks = c(0, 0.5, 1),
    labels = c(0, 0.5, 1))




# regions with at least one voxel active
out |>
  semi_join(distinct(gold, Task, l, label, n_parcels)) |>
  right_join(
    distinct(gold, Task, l, label, n_parcels) |> 
      crossing(distinct(out, n_sub))) |>
  filter(n_parcels %in% c(200, 400, 600, 800, 1000)) |>
  mutate(
    propr = if_else(is.na(prop), 0, prop),
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}"))
  ) |>
  ggplot(aes(x=n_sub, y=prop, group=l)) +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_grid(n_parcels~Task) +
  scale_y_continuous(
    "Proportion Simulations w/ Activity in Most Active ROI",
    limits = c(0, 1)
  ) +
  xlab("N Sub") +
  theme_gray(base_size = 12)

ggsave("analyses/figures/prop-active-most-active-roi-ptfce.png", width = 10, height = 6.5)

## NULL

# regions with at least one voxel active
out_null <- targets::tar_read(active_null) |>
  dplyr::left_join(at_list) |>
  distinct(
    n_sub, iter, label, `Label Name`,
    `Full component name`, n_parcels) |>
  collect() |>
  count(
    n_sub, label, `Label Name`,
    `Full component name`, n_parcels) |>
  mutate(prop = n / 100)

gold_null <- targets::tar_read(tfce_pop_null) |>
  dplyr::mutate(tmp = purrr::map(ptfce, get_max)) |>
  select(-ptfce) |>
  tidyr::unnest(tmp) |>
  left_join(at_list, multiple = "all") |>
  group_by(
    label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  summarise(
    Z = mean(Z),
    .groups = "drop"
  ) |>
  group_by(n_parcels) |>
  slice_max(order_by=Z, n=20) |>
  ungroup() |>
  filter(!is.na(label)) |> 
  mutate(
    l=label |> 
      factor() |> 
      as.numeric() |>
      factor()) |>
  ungroup()

out_null |>
  semi_join(distinct(gold_null, l, label, n_parcels)) |>
  right_join(
    distinct(gold_null, label, n_parcels, l) |> 
      crossing(distinct(out_null, n_sub))) |>
  filter(n_parcels %in% c(200, 400, 600, 800, 1000)) |>
  mutate(
    n_parcels = factor(n_parcels, levels = unique(n_parcels) |> sort()),
    n_parcels = fct_relabel(n_parcels, .fun = ~glue::glue("N Parcels: {.x}")),
    propr = if_else(is.na(prop), 0, prop)
  ) |>
  mutate(prop = dplyr::if_else(is.na(prop), 0, prop)) |>
  ggplot(aes(x=n_sub, y=prop, group=l)) +
  geom_ribbon(
    aes(ymin= 0.0073, ymax = 0.0926), 
    alpha=0.01, 
    fill="lightblue",
    linetype = "dashed",
    color = "black") +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_wrap(~n_parcels) +
  scale_y_continuous(
    "Proportion Simulations w/ Activity in Most Active ROI",
    limits = c(0, 0.1)
  ) +
  xlab("N Sub") +
  theme_gray(base_size = 12)

ggsave("analyses/figures/prop-active-most-active-roi_ptfce-null.png", width = 10, height = 6.5)
