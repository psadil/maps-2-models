library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(patchwork)
library(forcats)
source(here::here("R","utils.R"))
source(here::here("R","tfce.R"))
source(here::here("R","updates.R"))
source(here::here("R","poster.R"))
source(here::here("R","hcp.R"))

Sys.setenv(TAR_PROJECT = "hcp")

targets::tar_load(pop_cor_region)

pop_cor_region |>
  mutate(
    f = atanh(rho),
    `Network Name` = if_else(is.na(`Network Name`) & !is.na(label), "subcortical", `Network Name`)) |>
  filter(!is.na(rho) & is.finite(f)) |>
  group_by(Task, n_sub, ContrastName, method, `Network Name`) |>
  summarise(
    f = mean(f),
    N = n(),
    .groups = "drop"
  ) |>
  mutate(
    `N Sub` = factor(n_sub),
    SE = 1/sqrt(N-3),
    rr = tanh(f),
    lower = tanh(f - 1.96*SE),
    upper = tanh(f + 1.96*SE)) |>
  ggplot(aes(y=rr, x=`N Sub`, color = `Network Name`)) +
  facet_wrap(~Task) +
  geom_point(alpha=0.25) +
  geom_errorbar(aes(ymin=lower, ymax=upper)) +
  ylab("Correlation with Reference")


p <- pop_cor_region |>
  mutate(
    f = atanh(rho),
    `Network Name` = if_else(is.na(`Network Name`) & !is.na(label), "subcortical", `Network Name`)) |>
  filter(!is.na(rho) & is.finite(f)) |>
  group_by(Task, n_sub, ContrastName, method, `Network Name`, iter) |>
  summarise(
    f = mean(f),
    N = n(),
    .groups = "drop"
  ) |>
  mutate(
    rr = tanh(f),
    `N Sub` = factor(n_sub)) |>
  ggplot(aes(y=rr, x=`N Sub`, color = `Network Name`)) +
  facet_wrap(~Task) +
  geom_boxplot(outlier.alpha = 0.25) +
  scale_color_viridis_d(option="turbo") +
  ylab("Rank Correlation with Reference") +
  theme(legend.position = "bottom")

ggsave("analyses/cor-by-network.png", p)

targets::tar_load(at)


tmp <- pop_cor_region |>
  mutate(
    f = atanh(rho),
    `Network Name` = if_else(is.na(`Network Name`) & !is.na(label), "subcortical", `Network Name`)) |>
  filter(!is.na(rho) & is.finite(f)) |>
  group_by(Task, n_sub, ContrastName, method, `Label Name`) |>
  summarise(
    f = mean(f),
    N = n(),
    .groups = "drop"
  ) |>
  mutate(
    rr = tanh(f),
    `N Sub` = factor(n_sub)) |>
  right_join(filter(at, (z %% 5) == 0)) |>
  filter(between(z, 30, 70))

plot_cor_region_task <- function(tmp, task) {
  tmp |>
    filter(.data$Task==.env$task) |>
    ggplot(aes(x=x, y=y)) +
    coord_fixed() +
    facet_grid(`N Sub` ~ z, labeller = label_both) +
    geom_raster(
      aes(fill=value),
      data = 
        to_tbl(MNITemplate::getMNIPath(what="Brain", res = "2mm")) |>
        filter(z %in% unique(tmp$z)), show.legend = FALSE) +
    scale_fill_distiller(
      type = "seq",
      direction = -1,
      palette = "Greys") +
    ggnewscale::new_scale_fill() +
    geom_raster(aes(fill=rr), alpha = 0.5) +
    scale_fill_viridis_c(
      option="turbo",
      limits = c(0, 1)) +
    labs(fill="correlation") +
    theme_void() +
    theme(
      legend.position = "bottom"
    )  +
    ggtitle(task)
  
}

ps <- map(unique(tmp$Task), ~plot_cor_region_task(tmp, .x))

for (x in ps){
  ggsave(
    glue::glue("analyses/cor_region/{x$labels$title}.png"),
    x,
    device = ragg::agg_png,
    dpi = 600)
}


targets::tar_load(pop_cor2)

