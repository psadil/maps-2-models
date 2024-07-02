library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(c(space, data_topo_gold, gold_peaks, at))

gold_peaks_ <- gold_peaks |>
  dplyr::select(Task, m) |>
  tidyr::unnest(m) |>
  dplyr::left_join(at) |>
  dplyr::group_by(Task, label) |>
  dplyr::slice_max(
    order_by = Value,
    n = 1,
    with_ties = FALSE
  ) |> # grab highest peak from each label
  dplyr::group_by(Task) |>
  dplyr::slice_max(
    order_by = Value,
    n = 10,
    with_ties = FALSE
  ) |> # grab highest 10 peaks (distinct labels)
  dplyr::ungroup() |>
  dplyr::distinct(Task, x, y, z, Value)


p <- space |>
  mutate(Value = Value / sqrt(n_pop)) |>
  filter(Value > 0.2) |>
  group_by(Value, n_sub, corrp_thresh, Task) |>
  summarise(d = mean(d, na.rm = TRUE), .groups = "drop") |>
  mutate(
    `N Sub` = glue::glue("N Sub: {n_sub}"),
    `N Sub` = factor(
      `N Sub`,
      levels = c(
        "N Sub: 20",
        "N Sub: 40",
        "N Sub: 60",
        "N Sub: 80",
        "N Sub: 100"
      )
    )
  ) |>
  ggplot(aes(y = d, x = Value)) +
  geom_point(alpha = 0.1, shape = 20) +
  facet_grid(`N Sub` ~ Task) +
  scale_x_continuous(
    "Gold Standard Peak Cohen's d",
    breaks = c(0, 1),
    labels = c(0, 1)
  ) +
  scale_y_log10("avg dist(Gold Standard Peak, Study Peak) (mm)") +
  theme_gray(base_size = 8)

tikzDevice::tikz(
  "analyses/figures/peak-bysize.tex",
  width=5,
  height=4)
p
dev.off()


