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
    data_model_sub_to_sub))

a <- data_model_gold_gold_to_study |> 
  ggplot(aes(x=n_sub, y=avg)) +
  facet_wrap(~task) +
  geom_line(aes(group=type)) +
  geom_point() +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd)) +
  scale_x_log10("N Sub") +
  ylab("Average Rank Correlation (+-SD)\nPrediction-Truth (gF)") +
  theme(legend.position = "bottom")

b <- data_model_study_to_study |>
  filter(Model == "ICC(2,1)") |>
  ggplot(aes(x=n_sub, y=ICC)) +
  geom_line(aes(linetype=Measure)) +
  facet_wrap(~task) +
  xlab("N Sub") +
  ylab("ICC(2,1)")

# b <- data_model_study_to_study2 |>
#   ggplot(aes(x=s, y=n_sub)) +
#   facet_wrap(~task) +
#   ggdist::stat_dots(quantiles=100) +
#   scale_x_continuous(
#     "Prediction Standard Deviation\n(Study to Study)",
#     limits = c(0, NA)) +
#   ylab("N Sub")

cc <- data_model_sub_to_sub |>
  ggplot(aes(x=r, y=task)) +
  ggdist::stat_dots(quantiles = 100) +
  ylab("Task") +
  xlab("Pairwise Rank Correlation of Features\n(Sub to Sub)")

p <- a + b + cc +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme_gray(base_size = 8) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(8, "pt"))

ggsave(
  "analyses/figures/model.png", 
  p,
  width = 3.25, 
  height = 6,
  device = ragg::agg_png)
