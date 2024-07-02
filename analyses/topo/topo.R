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
    data_topo_gold,
    data_topo_gold_to_study,
    data_topo_study_to_study,
    data_topo_sub_to_sub
  )
)

a <- data_topo_gold |>
  ggplot() +
  facet_wrap(~Task) +
  scattermore::geom_scattermore(aes(x = cope, y = sigma), alpha = 0.1, pointsize = 1) +
  geom_function(fun = function(x) sqrt(20) * x / qt(0.001, 19), xlim = c(-200, 0)) +
  geom_function(fun = function(x) sqrt(100) * x / qt(0.001, 99), xlim = c(-200, 0), color = "gray50") +
  geom_function(fun = function(x) sqrt(100) * x / qt(0.999, 99), xlim = c(0, 200), color = "gray50") +
  geom_function(fun = function(x) sqrt(20) * x / qt(0.999, 19), xlim = c(0, 200)) +
  scale_y_continuous(
    limits = c(0, 150)
  ) +
  scale_x_continuous(
    limits = c(-200, 200),
    labels = c(-150, 0, 150),
    breaks = c(-150, 0, 150)
  ) +
  xlab(expression(beta ~ mean)) +
  ylab(expression(beta ~ SD))

b <- data_topo_gold_to_study |>
  mutate(`N Sub` = factor(n_sub)) |>
  ggplot(aes(x = rho, y = `N Sub`)) +
  facet_wrap(~Task) +
  geom_boxplot(outlier.shape = NA, size = .1) +
  geom_point(
    position = position_jitter(width = 0),
    shape = 20,
    size = 0.1,
    alpha = 0.1
  ) +
  xlab("Rank Correlation\n(Gold to Study)")

cc <- data_topo_study_to_study |>
  ggplot(aes(x = rr, y = `N Sub`, color = Task)) +
  geom_line(aes(group = Task)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper)) +
  xlab("Pairwise Rank Correlation\n(Study to Study)") +
  theme(
    legend.position = "bottom"
  )

d <- data_topo_sub_to_sub |>
  filter(!is.na(rho), stringr::str_detect(task, "EMOTION", TRUE)) |>
  ggplot(aes(x = rho, y = task)) +
  ggdist::stat_dots(
    quantiles = 100
  ) +
  ylab("Task") +
  xlab("Pairwise Product-Moment Correlation\n(Sub to Sub)")

p <- a + b + cc + d +
  plot_layout(ncol = 1) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme_gray(base_size = 8) +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(8, "pt")
    )

ggsave(
  "analyses/figures/topo.png",
  p,
  width = 3.25,
  height = 7,
  device = ragg::agg_png
)

tikzDevice::tikz(
  "analyses/figures/topo.tex",
  width=3.25,
  height=7)
p
dev.off()


data_topo_study_to_study |>
  group_by(task) |>
  summarise(
    q10 = quantile(rho, 0.1, na.rm = TRUE),
    q90 = quantile(rho, 0.9, na.rm = TRUE)
  ) |>
  mutate(di = q90 - q10) |>
  arrange(di)
