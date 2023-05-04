library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp")

# what proportion of studies have a peak that is within x mm of the gold standard?

targets::tar_load(space)

at <- targets::tar_read(at) |>
  mutate(across(c(x,y,z), as.integer))

at_list <- targets::tar_read(at_list) |>
  mutate(across(c(x,y,z), as.integer))

gold_peaks <- targets::tar_read(gold_peaks) |>
  select(Task, m) |>
  unnest(m) |>
  left_join(at) |>
  group_by(Task, label) |>
  slice_max(order_by = Value, n=1, with_ties = FALSE) |> # grab highest peak from each label
  group_by(Task) |>
  slice_max(order_by = Value, n=10, with_ties = FALSE) |> # grab highest 10 peaks (distinct labels)
  ungroup() |>
  distinct(Task, x, y, z) 

space |>
  filter(
    !is.na(study_ind), # not sure how this happens...
    corrp_thresh==0.95) |> 
  semi_join(gold_peaks) |>
  group_by(Task, n_sub, x, y, z, iter) |>
  summarize(
    within_2 = any(d < 2),
    within_4 = any(d < 4),
    within_6 = any(d < 6),
    within_8 = any(d < 8),
    within_10 = any(d < 10),
    within_20 = any(d < 20),
    .groups = "drop"
  ) |>
  group_by(Task, n_sub, x, y, z) |>
  summarize(across(starts_with("within"), sum), .groups="drop") |>
  pivot_longer(
    starts_with("within"),
    values_to = "n_simulations",
    names_to = "within", 
    names_pattern = "within_([[:digit:]]+)", 
    names_transform = as.integer) |>
  unite(col="peak", x, y, z) |>
  mutate(
    n_simulations = n_simulations / 100,
    n_sub = glue::glue("N Sub: {n_sub}"),
    n_sub = factor(
      n_sub, 
      levels = c("N Sub: 20", "N Sub: 40", "N Sub: 60", "N Sub: 80", "N Sub: 100"))) |>
  ggplot(aes(x=within, group=peak, y=n_simulations)) +
  geom_point(alpha=0.2) +
  geom_line(alpha=0.2) +
  facet_grid(n_sub~Task) +
  scale_y_continuous(
    "Proportion Simulations with Peak in Radius",
    limits = c(0,1),
    breaks = c(0,0.5,1),
    labels = c(0,0.5,1)
  ) +
  scale_x_continuous(
    "Radius (mm)",
    limits = c(0, 20)
  ) +
  theme_gray(base_size = 14)

ggsave("analyses/figures/prop-studies-within-mm.png", width = 9, height = 6)

