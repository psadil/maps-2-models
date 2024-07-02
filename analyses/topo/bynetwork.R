library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)

Sys.setenv(TAR_PROJECT = "hcp_ptfce")

targets::tar_load(pop_cor_region)

p <- pop_cor_region |>
  mutate(
    f = atanh(rho),
    `Network Name` = if_else(is.na(`Network Name`) & !is.na(label), "subcortical", `Network Name`)
  ) |>
  filter(!is.na(rho) & is.finite(f)) |>
  group_by(Task, n_sub, ContrastName, method, `Network Name`, iter) |>
  summarise(
    f = mean(f),
    N = n(),
    .groups = "drop"
  ) |>
  mutate(
    rr = tanh(f),
    `N Sub` = factor(n_sub)
  ) |>
  ggplot(aes(y = rr, x = `N Sub`, color = `Network Name`)) +
  facet_wrap(~Task) +
  geom_boxplot(outlier.alpha = 0.25) +
  scale_color_viridis_d(option = "turbo") +
  ylab("Rank Correlation with Reference") +
  theme(legend.position = "bottom")

ggsave(
  "analyses/figures/topo-bynetwork.png",
  p,
  width = 6,
  height = 4,
  device = ragg::agg_png
)


tikzDevice::tikz(
  "analyses/figures/topo-bynetwork.tex",
  width=6,
  height=4)
p
dev.off()

