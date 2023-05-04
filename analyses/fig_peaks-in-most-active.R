
#' How often are the highest peaks in the gold standard located in ROIs that 
#' have the largest average effect size?

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

Sys.setenv(TAR_PROJECT = "hcp")


# most active
gold <- targets::tar_read(pop_d) |>
  left_join(at_list) |>
  group_by(
    Task, label, `Label Name`,
    `Network Name`, `Full component name`, n_parcels) |>
  summarise(
    d = mean(cope/sigma),
    .groups = "drop"
  ) |>
  group_by(Task, n_parcels) |>
  slice_max(order_by=d, n=10) |>
  ungroup() |>
  filter(!is.na(label)) |>
  mutate(most_active = TRUE)

# all peaks
gold2 <- targets::tar_read(gold_peaks) |>
  select(Task, m) |>
  unnest(m) |>
  left_join(at_list) |>
  filter(!is.na(label)) |> 
  group_by(Task, label, n_parcels, n_networks) |>
  slice_max(order_by = Value, n = 1, with_ties = FALSE) |>
  ungroup() |>
  distinct(Task, n_parcels, n_networks, label) |>
  mutate(local_peak = TRUE)


c_on_act <- left_join(gold, gold2) |>
  mutate(local_peak = if_else(is.na(local_peak), FALSE, TRUE)) |>
  group_by(Task, n_parcels) |>
  summarise(
    local_peak = mean(local_peak),    
    .groups = "drop"
  ) 

c_on_peak <- left_join(gold2, gold) |>
  mutate(most_active = if_else(is.na(most_active), FALSE, TRUE)) |>
  group_by(Task, n_parcels) |>
  summarise(
    most_active = mean(most_active),
    .groups = "drop"
  )


left_join(c_on_act, c_on_peak) |>
  pivot_longer(
    c(local_peak, most_active), 
    names_to = "Condition",
    values_to = "Proportion") |>
  mutate(
    Condition = case_when(
      Condition %in% "local_peak" ~ "Local Peak",
      Condition %in% "most_active" ~ "Most Active"
    )
  ) |>
  ggplot(aes(x=n_parcels, y=Proportion, color=Condition)) +
  geom_line() +
  facet_wrap(~Task) +
  scale_y_continuous(
    "Proportion",
    limits = c(0,1)) +
  scale_x_continuous(
    "N Parcels",
    n.breaks = 6) +
  theme_gray(base_size = 16)

ggsave("analyses/figures/peaks-in-most-active.png", width=10, height = 6)

