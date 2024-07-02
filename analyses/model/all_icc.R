library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(data_model_study_to_study)

a <- data_model_study_to_study |>
  filter(stringr::str_detect(task, "EMOTION", TRUE)) |>
  filter(type == "agreement") |>
  mutate(max_avg = max(icc), .by = c(measure)) |>
  mutate(
    measure = stringr::str_replace_all(measure, "_", "\\\\_"),
    measure = factor(measure),
    measure = forcats::fct_reorder(measure, max_avg),
  ) |>
  ggplot(aes(x = n_sub, y = measure)) +
  facet_wrap(~task) +
  geom_raster(aes(fill = icc)) +
  scale_fill_viridis_c(option = "turbo", name = "ICC(A,1)") +
  xlab("N Sub") +
  ylab("Instrument") +
  theme_gray(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(8, "pt")
  )


b <- data_model_study_to_study |>
  filter(stringr::str_detect(task, "EMOTION", TRUE)) |>
  filter(type == "consistency") |>
  mutate(max_avg = max(icc), .by = c(measure)) |>
  mutate(
    measure = stringr::str_replace_all(measure, "_", "\\\\_"),
    measure = factor(measure),
    measure = forcats::fct_reorder(measure, max_avg)
  ) |>
  ggplot(aes(x = n_sub, y = measure)) +
  facet_wrap(~task) +
  geom_raster(aes(fill = icc)) +
  scale_fill_viridis_c(option = "turbo", name = "ICC(C,1)") +
  xlab("N Sub") +
  ylab("Instrument") +
  theme_gray(base_size = 7) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(8, "pt")
  )



ggsave(
  "analyses/figures/model_all_agreement.png",
  a,
  width = 6,
  height = 8,
  device = ragg::agg_png
)


ggsave(
  "analyses/figures/model_all_consistency.png",
  b,
  width = 6,
  height = 8,
  device = ragg::agg_png
)

tikzDevice::tikz(
  "analyses/figures/model_all_agreement.tex",
  width=6,
  height=8)
a
dev.off()

tikzDevice::tikz(
  "analyses/figures/model_all_consistency.tex",
  width=6,
  height=8)
b
dev.off()
