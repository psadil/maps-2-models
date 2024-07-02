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
  filter(Value > 0.2, !is.na(`Network Name`)) |>
  group_by(n_sub, corrp_thresh, Task, `Network Name`, iter) |>
  summarise(
    d = mean(d, na.rm = TRUE),
    Value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  ) |>
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
  filter(!is.na(d)) |>
  ggplot(aes(y = `Network Name`, x = d, color = Value)) +
  geom_boxplot(outlier.shape = NA) +
  scattermore::geom_scattermore(
    pointsize = 5, 
    position = position_jitter(width = 0),
    alpha = 0.5) +
  facet_grid(`N Sub` ~ Task) +
  scale_color_viridis_c(
    option = "turbo",
    guide = guide_colorbar("Cohen's d")
  ) +
  ylab("Network") +
  scale_x_continuous("avg dist(Gold Standard Peak, Study Peak) (mm)") +
  theme_gray(base_size = 8)

ggsave(
  "analyses/figures/peak-bynetwork.png",
  p,
  width = 6,
  height = 4,
  device = ragg::agg_png
)

tikzDevice::tikz(
  "analyses/figures/peak-bynetwork.tex",
  width=6,
  height=4)
p
dev.off()

