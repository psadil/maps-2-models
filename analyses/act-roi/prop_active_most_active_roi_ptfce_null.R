library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

source("R/tfce.R")

get_max <- function(q){
  x <- qs::qread(q)
  to_tbl0(x$Z, measure = "Z") |> mask()
}

# What proportion of simulated studies have active voxels within the ROIs within 
# that are most active in the gold standard?

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(c(at_list, active_null, tfce_pop_null, iter, gold_tested))

n_sims <- n_distinct(iter)

# regions with at least one voxel active
out_null <- active_null |>
  collect() |>
  dplyr::left_join(at_list, by = join_by(x,y,z), relationship = "many-to-many") |>
  distinct(
    n_sub, iter, label, `Label Name`,
    `Full component name`, n_parcels) |>
  count(
    n_sub, label, `Label Name`,
    `Full component name`, n_parcels) |>
  mutate(prop = n / n_sims)

# gold_null <- tfce_pop_null |>
#   dplyr::mutate(tmp = purrr::map(ptfce, get_max)) |>
#   select(-ptfce) |>
#   tidyr::unnest(tmp) |>
#   left_join(at_list, multiple = "all") |>
#   group_by(
#     label, `Label Name`,
#     `Network Name`, `Full component name`, n_parcels) |>
#   summarise(
#     Z = mean(Z),
#     .groups = "drop"
#   ) |>
#   group_by(n_parcels) |>
#   slice_max(order_by=Z, n=10) |>
#   ungroup() |>
#   filter(!is.na(label)) |> 
#   mutate(
#     l=label |> 
#       factor() |> 
#       as.numeric() |>
#       factor()) |>
#   ungroup()

active_threshold <- 0.02
gold_null <- gold_tested |>
  filter(Task=="WM") |>
  dplyr::mutate(
    d = statistic / sqrt(n_sub),
    active = abs(d) > active_threshold
  ) |>
  dplyr::mutate(
    r = dplyr::row_number(dplyr::desc(abs(estimate))),
    .by = c(Task, n_parcels)
  ) |>
  dplyr::filter(r < 11) |>
  dplyr::select(Task, n_parcels, label, d) |>
  mutate(l=label |>
           factor() |>
           as.numeric() |>
           factor())

p <- out_null |>
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
    aes(
      ymin=qbeta(0.05/2, 5, 100-5+1), 
      ymax=qbeta(1-0.05/2, 5, 100-5)
      ), 
    alpha=0.01, 
    fill="lightblue",
    linetype = "dashed",
    color = "black") +
  geom_point(show.legend = FALSE, alpha=0.2) +
  geom_line(show.legend = FALSE, alpha=0.2) +
  facet_wrap(~n_parcels) +
  scale_y_continuous(
    "Proportion Simulations w/\nActivity in Most Active ROI",
    limits = c(0, 0.12)
  ) +
  xlab("N Sub") +
  theme_gray(base_size = 12)


tikzDevice::tikz(
  "analyses/figures/prop-active-most-active-roi-ptfce-null.tex",
  width=5,
  height=3)
p
dev.off()

