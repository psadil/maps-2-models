library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(
  c(
    data_model_gold_gold_to_study,
    data_model_study_to_study,
    data_model_sub_to_sub
  )
)

a <- data_model_gold_gold_to_study |>
  filter(confounds=="True") |>
  filter(measure == "PMAT24_A_CR", type == "simulation", stringr::str_detect(task, "EMOTION", TRUE)) |>
  ggplot(aes(x = n_sub, y = avg)) +
  facet_wrap(~task) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_point(
    mapping = aes(x = n_sub, y = statistic_rep),
    color = "gold",
    data = filter(
      data_model_gold_gold_to_study,
      measure == "PMAT24_A_CR",
      type == "gold",
      stringr::str_detect(task, "EMOTION", TRUE),
      confounds=="True"
    )
  ) +
  scale_x_log10("N Sub") +
  ylab("Average Rank Correlation (CI)\nPrediction-Truth (gF)") +
  theme(legend.position = "bottom")

b <- data_model_study_to_study |>
  filter(confounds=="True") |>
  filter(measure == "PMAT24_A_CR", stringr::str_detect(task, "EMOTION", TRUE)) |>
  ggplot(aes(x = n_sub, y = icc, color = type), alpha = 0.5) +
  facet_wrap(~task) +
  geom_line() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  xlab("N Sub") +
  ylab("ICC")

# b <- data_model_study_to_study2 |>
#   ggplot(aes(x=s, y=n_sub)) +
#   facet_wrap(~task) +
#   ggdist::stat_dots(quantiles=100) +
#   scale_x_continuous(
#     "Prediction Standard Deviation\n(Study to Study)",
#     limits = c(0, NA)) +
#   ylab("N Sub")

cc <- data_model_sub_to_sub |>
  filter(confounds=="True") |>
  filter(stringr::str_detect(task, "EMOTION", TRUE)) |>
  ggplot(aes(x = r, y = task)) +
  ggdist::stat_dots(quantiles = 100) +
  ylab("Task") +
  xlab("Pairwise Rank Correlation of Features\n(Sub to Sub)")

p <- a + b + cc +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme_gray(base_size = 8) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(8, "pt")
    )

ggsave(
  "analyses/figures/model.png",
  p,
  width = 3.25,
  height = 6,
  device = ragg::agg_png
)

tikzDevice::tikz(
  "analyses/figures/model.tex",
  width=3.25,
  height=6)
p
dev.off()

# "Predictions for fluid intelligence were significant with all of the considered tasks..."
data_model_gold_gold_to_study |> 
  filter(measure=="PMAT24_A_CR", n_sub>200)

