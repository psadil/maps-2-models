library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(data_model_gold_gold_to_study)

p <- data_model_gold_gold_to_study |>
  filter(type == "simulation", stringr::str_detect(task, "EMOTION", TRUE)) |>
  mutate(max_avg = max(avg), .by = c(measure)) |>
  mutate(
    measure = stringr::str_replace_all(measure, "_", "\\\\_"),
    measure = factor(measure),
    measure = forcats::fct_reorder(measure, max_avg)
  ) |>
  ggplot(aes(x = n_sub, y = measure)) +
  facet_wrap(~task) +
  geom_raster(aes(fill = avg)) +
  scale_fill_viridis_c(option = "turbo", name = "Rank\nCorrelation") +
  xlab("N Sub") +
  ylab("Instrument") +
  theme_gray(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(8, "pt")
  )

ggsave(
  "analyses/figures/all_cog.png",
  p,
  width = 6,
  height = 8,
  device = ragg::agg_png
)

tikzDevice::tikz(
  "analyses/figures/all_cog.tex",
  width=6,
  height=8)
p
dev.off()
