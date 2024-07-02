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
  mutate(across(c(x, y, z), as.integer))

at_list <- targets::tar_read(at_list) |>
  mutate(across(c(x, y, z), as.integer))

gold_peaks <- targets::tar_read(gold_peaks) |>
  select(Task, m) |>
  unnest(m) |>
  left_join(at) |>
  group_by(Task, label) |>
  slice_max(order_by = Value, n = 1, with_ties = FALSE) |> # grab highest peak from each label
  group_by(Task) |>
  slice_max(order_by = Value, n = 10, with_ties = FALSE) |> # grab highest 10 peaks (distinct labels)
  ungroup() |>
  distinct(Task, x, y, z) |>
  mutate(peak_i = 1:n(), .by = Task)

tmp <- space |>
  filter(
    !is.na(study_ind), # not sure how this happens... maybe graymatter mask vs cereb?
    corrp_thresh == 0.01
  ) |>
  left_join(gold_peaks) |>
  select(n_sub, x = x.study, y = y.study, z = z.study, iter, Task, peak_i) |>
  group_nest(Task, n_sub, peak_i) |>
  mutate(
    d = map(
      data,
      ~ .x |>
        arrange(iter) |>
        select(-iter) |>
        dist() |>
        as.matrix() |>
        corrr::as_cordf() |>
        corrr::shave() |>
        corrr::stretch(na.rm = TRUE) |>
        rename(d = r)
    )
  ) |>
  select(-data) |>
  mutate(n_sub = factor(n_sub)) |>
  unnest(d)


tmp |>
  filter(!is.na(n_sub)) |>
  ggplot(aes(y = n_sub, x = d)) +
  facet_wrap(~Task) +
  ggdist::stat_dots(quantiles = 50) +
  ylab("N Sub") +
  scale_x_continuous(
    "Distance Between Highest Peaks"
  ) +
  theme_gray(base_size = 9)

ggsave(
  "analyses/figures/prop-peak-vox-study-to-study_euclid.png",
  width = 3,
  height = 3
)
