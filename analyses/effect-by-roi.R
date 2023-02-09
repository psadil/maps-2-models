library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)
source(here::here("R","utils.R"))
source(here::here("R","tfce.R"))
source(here::here("R","updates.R"))
source(here::here("R","poster.R"))



eff_by_roi <- targets::tar_read(eff_by_roi, store=here::here("_hcp")) |>
  group_by(
    label, `Label Name`, `Network Name`, `Full component name`, hemi, n_voxels, volume, 
    Task, CopeNumber, ContrastName, avail) |>
  summarise(cope = mean(cope), .groups = "drop") |>
  collect()

# at <- targets::tar_read(at, store=here::here("_hcp")) |>
#   mutate(across(.cols = c(x, y, z), as.integer))
# out <- targets::tar_read(gray, store=here::here("_hcp")) |>
#   left_join(at, by = c("x", "y", "z")) |>
#   group_by(
#     label, `Label Name`, `Network Name`, `Full component name`, hemi, n_voxels, volume, 
#     Task, CopeNumber, ContrastName, n_sub, iter) |>
#   summarise(cope = mean(cope), sigma=mean(sigma), .groups = "drop") |>
#   collect()

eff_by_roi |>
  filter(Task == "MOTOR") |>
  mutate(
    label = str_remove(label, "7Networks|7Networks_"),
    label = str_extract(label, "[[:digit:]]+"),
    label = factor(as.numeric(label))) |>
  ggplot(aes(y=label, x=cope, color=hemi)) +
  facet_wrap(~`Network Name`, scales = "free") +
  ggridges::geom_density_ridges(
    jittered_points = TRUE,
    position = ggridges::position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 1, point_alpha = 0.1, alpha = 0.7
  )


eff_by_roi |>
  filter(stringr::str_detect(Task, "MOTOR")) |>
  mutate(
    `Network Name` = if_else(is.na(`Network Name`), "subcortical", `Network Name`), 
    label = str_remove(label, "7Networks|7Networks_"),
    label = str_extract(label, "(?<=_)[[:digit:]]+$|[[:alpha:]]+$"),
    label = factor(label)) |>
  ggplot(aes(x=label, y=cope, color=hemi)) +
  facet_wrap(~`Network Name`, scales = "free") +
  geom_violin(alpha=0.1) +
  coord_flip()
