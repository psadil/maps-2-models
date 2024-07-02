library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")


dataset <- "/home/ubuntu/mnt/meta/act_preds/data/out"

d <- arrow::open_dataset(dataset) |>
  dplyr::select(sub, fold, n_sub, task, y_hat) |>
  dplyr::collect() |>
  group_nest(n_sub, task) |>
  mutate(
    rex = map(
      data,
      ~ ReX::lme_ICC_2wayR(data = .x$y_hat, subID = .x$sub, session = .x$fold) |>
        as_tibble()
    )
  ) |>
  select(-data) |>
  unnest(rex)


d |>
  select(n_sub, task, starts_with("ICC")) |>
  pivot_longer(starts_with("ICC"), names_to = "Model", values_to = "ICC") |>
  mutate(
    Measure = case_match(
      Model,
      "ICC.a" ~ "Agreement",
      "ICC.c" ~ "Consistency",
      "ICCk.a" ~ "Agreement",
      "ICCk.c" ~ "Consistency"
    ),
    Model = case_match(
      Model,
      "ICC.a" ~ "ICC(2,1)",
      "ICC.c" ~ "ICC(2,1)",
      "ICCk.a" ~ "ICC(2,k)",
      "ICCk.c" ~ "ICC(2,k)"
    )
  ) |>
  filter(Model == "ICC(2,k)") |>
  ggplot(aes(x = n_sub, y = ICC)) +
  geom_line(aes(linetype = Measure)) +
  facet_wrap(~task) +
  xlab("N Sub") +
  ylab("ICC(2,k)")
